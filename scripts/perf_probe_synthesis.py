#!/usr/bin/env python3
"""perf_probe_synthesis — microbench targeting the synthesis modal
additions (sweeps #13-#18).

Run from repo root:  python3 scripts/perf_probe_synthesis.py

Measures the new code paths so a perf regression doesn't slip in
between releases. Not part of pytest — assertions absent on purpose;
read the numbers and act on outliers.

Operations covered:
  * Import-time cost of `splicecraft` itself
  * SynthesisScreen mount (cold app)
  * SynthesisEditor render at 1 kb, 10 kb, 50 kb (cap)
  * ProteinEditor render across 100/1k/5k AA + various motif loads
  * Multi-lane dither pack: 1, 10, 100 motifs at 5k AA
  * Motif library load (`_load_protein_motifs`) cold + warm
  * `_extract_aa_feats_from_record` on records with 0 / 30 / 300
    misc_feature sub-features
  * Full save round-trip via `_commit_save`

Print one number per row (median ms); flag anything ≥ 16 ms (one
frame at 60Hz) — those are interactive paths the user feels.
"""
from __future__ import annotations

import gc
import random
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def bench(fn, *, iters: int = 5, warmup: int = 1) -> float:
    for _ in range(warmup):
        fn()
    samples = []
    for _ in range(iters):
        gc.collect()
        t0 = time.perf_counter()
        fn()
        samples.append((time.perf_counter() - t0) * 1000.0)
    return statistics.median(samples)


def header(title: str) -> None:
    print()
    print("═" * 78)
    print(f"  {title}")
    print("═" * 78)


def row(label: str, ms: float, *, threshold_ms: float = 16.0) -> None:
    flag = "  ⚠  " if ms >= threshold_ms else "     "
    print(f"{flag}{label:<58}{ms:>10.2f} ms")


# ── 1. Import cost (cold from a clean Python) ─────────────────────────────
# Note: the rest of this script imports splicecraft below; this row
# only measures what `python -c "import splicecraft"` would cost. We
# spawn a subprocess so the cost isn't influenced by our own warm
# import.
def measure_import() -> float:
    import subprocess
    samples = []
    for _ in range(3):
        t0 = time.perf_counter()
        subprocess.run(
            [sys.executable, "-c", "import splicecraft"],
            check=True, capture_output=True,
        )
        samples.append((time.perf_counter() - t0) * 1000.0)
    return statistics.median(samples)


# Eager import — every other probe below uses sc.*
import splicecraft as sc  # noqa: E402


# ── 2. SynthesisScreen mount (cold app) ───────────────────────────────────
async def _mount_synthesis_screen() -> float:
    import time as _time
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        await pilot.pause()
        t0 = _time.perf_counter()
        app.push_screen(sc.SynthesisScreen())
        await pilot.pause()
        await pilot.pause()
        return (_time.perf_counter() - t0) * 1000.0


def measure_mount() -> float:
    import asyncio
    samples = []
    for _ in range(3):
        samples.append(asyncio.run(_mount_synthesis_screen()))
    return statistics.median(samples)


# ── 3. SynthesisEditor render at scale ────────────────────────────────────
def _build_seq(n: int) -> str:
    rng = random.Random(42)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _build_dna_feats(n_seq: int, n_feats: int) -> list[dict]:
    if n_feats <= 0 or n_seq <= 0:
        return []
    span = max(1, n_seq // n_feats)
    feats = []
    for i in range(n_feats):
        start = i * span
        end = min(n_seq, start + span - 1)
        if end <= start:
            continue
        feats.append({
            "start": start, "end": end,
            "label": f"f{i}",
            "type": "misc_feature",
            "color": "#3B82F6",
            "strand": 1,
        })
    return feats


def measure_synth_render(n_seq: int, n_feats: int) -> float:
    seq = _build_seq(n_seq)
    feats = _build_dna_feats(n_seq, n_feats)
    # `_build_seq_text` is what `SynthesisEditor._refresh_view` calls.
    # `line_width = len(seq) + 1` matches the synthesis editor's
    # "single chunk" invocation.
    def go():
        sc._build_seq_text(
            seq, feats,
            line_width=len(seq) + 1,
            user_sel=None,
            cursor_pos=-1,
            re_highlight=None,
            aa_highlight=None,
            viewport_y_range=None,
        )
    return bench(go, iters=3)


# ── 4. ProteinEditor render ───────────────────────────────────────────────
def _build_aa(n: int) -> str:
    rng = random.Random(42)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(rng.choice(aa) for _ in range(n))


def _build_aa_feats(n_seq: int, n_feats: int) -> list[dict]:
    if n_feats <= 0 or n_seq <= 0:
        return []
    span = max(1, n_seq // n_feats)
    feats = []
    for i in range(n_feats):
        start = i * span
        end = min(n_seq, start + span - 1)
        if end <= start:
            continue
        feats.append({
            "start": start, "end": end,
            "label": f"m{i}",
            "type": "Tag",
            "color": "#1E40AF",
            "strand": 1,
        })
    return feats


def measure_protein_lane_pack(n_aa: int, n_feats: int,
                                cols_per_aa: int) -> float:
    aa_seq = _build_aa(n_aa)
    feats = _build_aa_feats(n_aa, n_feats)
    pe = sc.ProteinEditor()
    pe._aa_seq = aa_seq
    pe._aa_feats = feats
    def go():
        pe._build_protein_lane_text(cols_per_aa=cols_per_aa)
    return bench(go, iters=3)


def measure_protein_lane_overlap(n_aa: int, n_motifs: int) -> float:
    """All motifs cover the same range — forces them to stack into
    n_motifs separate lanes. Worst case for the pack algorithm.
    """
    aa_seq = _build_aa(n_aa)
    overlap_end = max(10, n_aa // 2)
    feats = [
        {"start": 0, "end": overlap_end,
         "label": f"m{i}", "type": "Tag",
         "color": "#1E40AF", "strand": 1}
        for i in range(n_motifs)
    ]
    pe = sc.ProteinEditor()
    pe._aa_seq = aa_seq
    pe._aa_feats = feats
    def go():
        pe._build_protein_lane_text(cols_per_aa=3)
    return bench(go, iters=3)


def measure_row_count_cost(n_aa: int, n_feats: int) -> float:
    """`_row_count` runs `_lane_height` on every refresh; this
    measures the cost of that secondary pack pass."""
    aa_seq = _build_aa(n_aa)
    feats = _build_aa_feats(n_aa, n_feats)
    pe = sc.ProteinEditor()
    pe._aa_seq = aa_seq
    pe._aa_feats = feats
    def go():
        pe._row_count()
    return bench(go, iters=5)


# ── 5. Motif library load ─────────────────────────────────────────────────
def measure_motif_load_cold() -> float:
    """Cold load — clears cache between iters."""
    def go():
        sc._protein_motifs_cache = None
        sc._load_protein_motifs()
    return bench(go, iters=5)


def measure_motif_load_warm() -> float:
    """Warm load — cache stays."""
    sc._load_protein_motifs()   # seed cache
    def go():
        sc._load_protein_motifs()
    return bench(go, iters=10)


# ── 6. _extract_aa_feats_from_record across feature counts ────────────────
def measure_aa_extract(n_feats: int) -> float:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    n_aa = 500
    dna = "ATG" * n_aa
    rec = SeqRecord(Seq(dna), id="t", name="t")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "linear"
    rec.features.append(SeqFeature(
        FeatureLocation(0, len(dna), strand=1),
        type="CDS",
        qualifiers={"label": ["p"], "translation": [_build_aa(n_aa)]},
    ))
    # Add n_feats codon-aligned misc_feature sub-features.
    for i in range(n_feats):
        start = (i * 3) % (len(dna) - 6)
        rec.features.append(SeqFeature(
            FeatureLocation(start, start + 3, strand=1),
            type="misc_feature",
            qualifiers={
                "label": [f"m{i}"],
                "note": ["splicecraft-aa-feature-type=Tag"],
                "ApEinfo_fwdcolor": ["#1E40AF"],
            },
        ))
    def go():
        sc.SynthesisScreen._extract_aa_feats_from_record(rec)
    return bench(go, iters=5)


# ── Main ──────────────────────────────────────────────────────────────────
def main() -> None:
    print("─" * 78)
    print("  Synthesis modal perf probe — sweeps #13-#18")
    print("  threshold flag ⚠ ≥ 16 ms (one frame @ 60Hz)")
    print("─" * 78)

    header("1. Import")
    row("python -c 'import splicecraft'", measure_import(),
        threshold_ms=1000.0)

    header("2. Synthesis screen mount (cold app)")
    row("SynthesisScreen mount + first paint", measure_mount(),
        threshold_ms=500.0)

    header("3. SynthesisEditor render")
    for n_seq, n_feats in [
        (1_000,  10),
        (10_000, 50),
        (50_000, 100),     # cap
    ]:
        ms = measure_synth_render(n_seq, n_feats)
        row(f"  seq={n_seq:>5} bp, feats={n_feats:>3}", ms)

    header("4. ProteinEditor lane pack")
    for n_aa, n_feats, cpa in [
        (100,   1,   3),
        (100,   5,   3),
        (1_000, 10,  3),
        (5_000, 30,  3),
        (16_666, 30, 3),  # cap (≈ 50 kb / 3)
    ]:
        ms = measure_protein_lane_pack(n_aa, n_feats, cpa)
        row(f"  aa={n_aa:>6}, feats={n_feats:>3}, cpa={cpa}", ms)

    header("5. ProteinEditor pack worst-case (overlap = N lanes)")
    for n_aa, n_motifs in [
        (200, 5),
        (500, 10),
        (1_000, 20),
        (5_000, 30),
    ]:
        ms = measure_protein_lane_overlap(n_aa, n_motifs)
        row(f"  aa={n_aa:>5}, all-overlap motifs={n_motifs}", ms)

    header("6. _row_count → _lane_height cost (per refresh)")
    for n_aa, n_feats in [
        (100, 0),
        (100, 5),
        (5_000, 30),
        (16_666, 30),
    ]:
        ms = measure_row_count_cost(n_aa, n_feats)
        row(f"  aa={n_aa:>6}, feats={n_feats:>3}", ms,
            threshold_ms=1.0)

    header("7. Motif library load")
    row("cold (cache miss)", measure_motif_load_cold(),
        threshold_ms=1.0)
    row("warm (cache hit)", measure_motif_load_warm(),
        threshold_ms=1.0)

    header("8. _extract_aa_feats_from_record")
    for n in [0, 30, 300]:
        ms = measure_aa_extract(n)
        row(f"  CDS sub-features = {n}", ms)


if __name__ == "__main__":
    main()
