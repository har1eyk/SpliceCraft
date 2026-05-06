#!/usr/bin/env python3
"""commercialsaas_inspect — reverse-engineering aid for CommercialSaaS .dna files.

The CommercialSaaS binary format is undocumented. BioPython parses 5 packet
types out of ~15+. To round-trip with CommercialSaaS we have to discover the
rest. This script walks any .dna file's packet stream and prints what
we find so we can identify unknown packets, decode their payloads, and
build up a documented catalog.

Usage
-----
    python3 scripts/commercialsaas_inspect.py path/to/file.dna
    python3 scripts/commercialsaas_inspect.py path/to/file.dna --dump 0x07
    python3 scripts/commercialsaas_inspect.py --catalog

Modes
-----
- default: enumerate every packet — type byte, length, content sniff
  (XML root tag, ASCII vs binary, hex preview).
- ``--dump <type_byte>``: write the raw payload of every packet matching
  that type to ``./commercialsaas_dump_<type>_<n>.bin`` for offline inspection.
- ``--catalog``: print the catalog of packet types we've documented so
  far. Edit this file's ``KNOWN_PACKETS`` dict as we learn more.

Format primitives
-----------------
Each packet is a TLV: 1 byte type + 4 bytes big-endian length + N bytes
payload. The cookie packet (0x09) MUST come first; everything after it
is in implementation-defined order.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from struct import unpack


# Packet types we've documented (so far). Update as we learn more.
# Sources:
#  - BioPython Bio/SeqIO/CommercialSaaSIO.py (cookie, dna, primers, notes, features)
#  - Reverse-engineering against real .dna files (in progress)
KNOWN_PACKETS: dict[int, dict[str, str]] = {
    0x00: {
        "name":   "DNA",
        "format": "1-byte flags + N-byte ASCII sequence. "
                  "Flag bit 0x01 = circular. Other bits unknown "
                  "(methylation? damage? read-only?).",
    },
    0x02: {
        "name":   "EnzymeSet",
        "format": "Binary table of 12-byte records — appears to be "
                  "the per-file enzyme cut catalog (which enzymes "
                  "cut the sequence + how many times). Encoding "
                  "TBD — needs cross-reference against the .dna's "
                  "actual cut counts.",
    },
    0x03: {
        "name":   "EnzymeDefs",
        "format": "Binary with embedded ASCII recognition strings "
                  "separated by commas (e.g., 'GACNNNNNNGTC,"
                  "CCAGGCCTGG'). Likely the IUPAC site list for the "
                  "current enzyme set. Encoding TBD.",
    },
    0x05: {
        "name":   "Primers",
        "format": "UTF-8 XML, root <Primers>. Each <Primer> has "
                  "<BindingSite> children with location, strand, "
                  "annealedBases, meltingTemperature.",
    },
    0x06: {
        "name":   "Notes",
        "format": "UTF-8 XML, root <Notes>. Children: <Type>, "
                  "<LastModified> (YYYY.MM.DD), <AccessionNumber>, "
                  "<Comments>, others.",
    },
    0x07: {
        "name":   "History",
        "format": "**xz-compressed** UTF-8 XML, root <HistoryTree>. "
                  "Decompress with `lzma.decompress(payload)`. "
                  "Schema: one top-level <Node> (the current plasmid) "
                  "with recursive child <Node> elements representing "
                  "parent fragments. Each <Node> has attrs name, "
                  "type='DNA', seqLen, strandedness, ID, circular, "
                  "operation ('insertFragment', 'replace', 'insert', "
                  "...), resurrectable. Children: <RegeneratedSite>, "
                  "<HistoryColors>, <InputSummary>, <Features>, and "
                  "nested <Node> for parents.",
    },
    0x08: {
        "name":   "AddProps",
        "format": "UTF-8 XML, root <AdditionalSequenceProperties>. "
                  "Children: <UpstreamStickiness>, "
                  "<DownstreamStickiness>, <UpstreamModification> "
                  "(e.g., FivePrimePhosphorylation), "
                  "<DownstreamModification>.",
    },
    0x09: {
        "name":   "Cookie",
        "format": "8-byte ASCII 'CommercialSaaS' + 3 unsigned shorts: "
                  "seqType, exportVersion, importVersion. Always first.",
    },
    0x0A: {
        "name":   "Features",
        "format": "UTF-8 XML, root <Features>. Each <Feature> has "
                  "type, directionality (1=fwd, 2=rev), name, "
                  "<Segment> children with range='start-end' "
                  "(1-based), and <Q name=...> qualifier elements.",
    },
    0x0B: {
        "name":   "Alignment",
        "format": "Binary alignment record. Multiple of these in one "
                  "file — one per aligned sequence. Encoding TBD; "
                  "starts with a 4-byte ID followed by per-position "
                  "data. Cross-reference with <AlignableSequences> "
                  "(packet 0x11) for the human-readable index.",
    },
    0x0D: {
        "name":   "EnzymeVisSet",
        "format": "Binary, contains an ASCII enzyme-set name string "
                  "(e.g., 'Unsaved Enzyme Set'). Encoding TBD.",
    },
    0x0E: {
        "name":   "CustomEnzymes",
        "format": "UTF-8 XML, root <CustomEnzymeSets>. Each "
                  "<CustomEnzymeSet> has type, name, enzymeNames "
                  "(space-separated list).",
    },
    0x10: {
        "name":   "AlignmentData",
        "format": "Binary. Header bytes 'ae 5a 54 52 0d 0a 1a 0a' "
                  "(non-PNG; CommercialSaaS-internal magic). Payload is "
                  "indexed by a leading 4-byte ID — matches the "
                  "<AlignableSequences> entries one-to-one.",
    },
    0x11: {
        "name":   "Alignable",
        "format": "UTF-8 XML, root <AlignableSequences>. Children: "
                  "<Sequence> with name, ID, sortOrder, trimmedRange. "
                  "Index for the binary alignment packets (0x0B, 0x10).",
    },
    0x1B: {
        "name":   "GZipped",
        "format": "**gzip-compressed** payload (magic 1f 8b 08). "
                  "Content TBD — possibly a thumbnail / preview map / "
                  "rendered image.",
    },
    0x1C: {
        "name":   "EnzymeVis",
        "format": "UTF-8 XML, root <EnzymeVisibilities>. Per-enzyme "
                  "visibility flags for the current view.",
    },
}


def iter_packets(data: bytes):
    """Yield ``(packet_idx, type_byte, length, payload_bytes, file_offset)``
    tuples for every packet in the .dna byte stream. Stops on EOF or on
    a malformed length header (so partial files don't raise)."""
    offset = 0
    idx = 0
    n = len(data)
    while offset < n:
        if offset + 5 > n:
            print(f"  ! truncated header at offset {offset} "
                  f"({n - offset} bytes remaining)",
                  file=sys.stderr)
            return
        type_byte = data[offset]
        length = unpack(">I", data[offset + 1:offset + 5])[0]
        payload_start = offset + 5
        payload_end = payload_start + length
        if payload_end > n:
            print(f"  ! packet {idx} type=0x{type_byte:02X} length={length} "
                  f"runs past EOF (only {n - payload_start} bytes available)",
                  file=sys.stderr)
            return
        yield idx, type_byte, length, data[payload_start:payload_end], offset
        offset = payload_end
        idx += 1


def sniff_payload(payload: bytes) -> tuple[str, str]:
    """Return (kind, preview) where kind ∈ {'xml', 'ascii', 'binary'}.
    Best-effort detection — we want a quick eyeball for what's inside."""
    if not payload:
        return ("empty", "")
    # XML sniff: starts with optional BOM/whitespace, then '<'.
    head = payload[:32]
    stripped = head.lstrip(b" \t\r\n\xef\xbb\xbf")
    if stripped[:1] == b"<":
        try:
            txt = payload.decode("utf-8")
            m = re.match(r"\s*<\s*(\w+)", txt)
            root = m.group(1) if m else "?"
            preview = txt[:200].replace("\n", " ")
            if len(txt) > 200:
                preview += "…"
            return ("xml", f"<{root}> root  ·  {preview}")
        except UnicodeDecodeError:
            pass
    # ASCII text (e.g., the DNA sequence).
    try:
        txt = payload.decode("ascii")
        if all(c.isprintable() or c in "\r\n\t" for c in txt):
            preview = txt[:64]
            if len(txt) > 64:
                preview += "…"
            return ("ascii", f"{len(txt)} chars  ·  {preview!r}")
    except UnicodeDecodeError:
        pass
    # Binary — hex preview.
    hex_preview = payload[:32].hex(" ")
    if len(payload) > 32:
        hex_preview += " …"
    return ("binary", hex_preview)


def cmd_enumerate(path: Path, *, max_payload_preview: int = 200) -> None:
    """Walk every packet in `path` and print a one-line summary per
    packet. The first packet should be the cookie; any other type at
    position 0 is a red flag."""
    raw = path.read_bytes()
    print(f"# {path}  ({len(raw):,} bytes)")
    print(f"# {'idx':>3}  {'off':>8}  {'type':>4}  "
          f"{'len':>9}  {'name':<10}  kind     preview")
    for idx, type_byte, length, payload, file_off in iter_packets(raw):
        info = KNOWN_PACKETS.get(type_byte, {})
        name = info.get("name", "(unknown)")
        kind, preview = sniff_payload(payload)
        if len(preview) > max_payload_preview:
            preview = preview[:max_payload_preview] + "…"
        print(f"  {idx:3d}  {file_off:>8d}  0x{type_byte:02X}  "
              f"{length:>9,d}  {name:<10}  {kind:<7}  {preview}")


def cmd_dump(path: Path, target_type: int, *, out_dir: Path) -> None:
    """Write every payload matching `target_type` to a file in
    `out_dir`. Useful when you want to feed an unknown XML into a
    pretty-printer or save a hex dump to study byte-by-byte."""
    out_dir.mkdir(parents=True, exist_ok=True)
    raw = path.read_bytes()
    n = 0
    stem = path.stem
    for idx, type_byte, length, payload, _ in iter_packets(raw):
        if type_byte != target_type:
            continue
        out = out_dir / f"{stem}_{idx:02d}_0x{type_byte:02X}.bin"
        out.write_bytes(payload)
        n += 1
        kind, preview = sniff_payload(payload)
        print(f"wrote {out} ({length:,} bytes, {kind}: {preview[:80]})")
    if n == 0:
        print(f"no packets of type 0x{target_type:02X} found in {path}",
              file=sys.stderr)


def cmd_catalog() -> None:
    print("# Documented CommercialSaaS .dna packet types")
    for byte, info in sorted(KNOWN_PACKETS.items()):
        print(f"\n0x{byte:02X}  {info['name']}")
        print(f"      {info['format']}")
    print("\n# Unknown packet types (encountered in real .dna files but "
          "not documented):")
    print("# Run the enumerator on a sample file to find them; add "
          "entries to KNOWN_PACKETS as their structure is decoded.")


def _parse_hex_byte(s: str) -> int:
    s = s.strip().lower()
    if s.startswith("0x"):
        s = s[2:]
    return int(s, 16)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Enumerate / dump CommercialSaaS "
                                              ".dna packets.")
    p.add_argument("file", nargs="?", type=Path,
                    help="Path to the .dna file to inspect.")
    p.add_argument("--dump", metavar="HEX",
                    help="Write every packet matching this type byte "
                         "(e.g. 0x07) to ./commercialsaas_dump/. Use with file.")
    p.add_argument("--out", type=Path, default=Path("commercialsaas_dump"),
                    help="Output directory for --dump (default: "
                         "./commercialsaas_dump).")
    p.add_argument("--catalog", action="store_true",
                    help="Print the documented packet-type catalog and "
                         "exit.")
    args = p.parse_args(argv)
    if args.catalog:
        cmd_catalog()
        return 0
    if not args.file:
        p.print_help(sys.stderr)
        return 2
    if not args.file.exists():
        print(f"file not found: {args.file}", file=sys.stderr)
        return 1
    if args.dump:
        cmd_dump(args.file, _parse_hex_byte(args.dump), out_dir=args.out)
    else:
        cmd_enumerate(args.file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
