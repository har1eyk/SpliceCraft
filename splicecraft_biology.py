"""splicecraft_biology — pure biology primitives extracted from
splicecraft.py as a controlled test of the single-file seam.

The single-file rule (entire app in splicecraft.py) is intentional;
see CLAUDE.md + docs/architecture.md. This module is the first
deliberate extraction, scoped to entities that:

1. Have NO `PlasmidApp` coupling (no `self.notify`, no Textual,
   no reactive attrs).
2. Are pure functions or top-level constants.
3. Are imported back into `splicecraft.py` so external callers
   keep `splicecraft._rc(...)` etc. unchanged.

If this stays clean (no cross-imports back into splicecraft.py, no
test churn), it's the precedent for future extractions. See
CONTRIBUTING.md's "three-test rule" for the criteria.

Sacred invariants this module owns:
  #3  — `_rc` handles full IUPAC via `_IUPAC_COMP` (not just ACGT).
  #4  — `_iupac_pattern` is bounded-LRU cached in `_PATTERN_CACHE`.
  #8  — `_feat_len` returns `(total - start) + end` when `end < start`
       so wrap features have the right length.

These are the same invariants documented in `CLAUDE.md`; the
extraction did not change them, only their physical home.
"""
from __future__ import annotations

import functools
import re
from collections import OrderedDict


# ── IUPAC + reverse complement ────────────────────────────────────────────


_IUPAC_RE: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "W": "[AT]", "S": "[CG]",
    "M": "[AC]", "K": "[GT]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


# Pattern cache (sacred invariant #4). Bounded LRU so a long-lived
# process scanning many recognition sites can't grow the cache
# indefinitely; the catalog is ~120 enzymes (palindromic + RC variants
# < 256), so 256 is comfortably above steady state. Public dict so
# tests can `.clear()` and inspect membership.
_PATTERN_CACHE: "OrderedDict[str, re.Pattern[str]]" = OrderedDict()
_PATTERN_CACHE_MAX = 256


def _iupac_pattern(site: str) -> "re.Pattern[str]":
    # Sweep #22: case-fold cache key. Patterns ARE always built from
    # `site.upper()` internally, so two calls with `"gaattc"` and
    # `"GAATTC"` produce identical regex objects but pre-fix occupied
    # two separate slots — wasting one cap unit per mixed-case
    # variant. Normalize to uppercase before lookup AND store.
    #
    # 2026-05-27 (audit-5 restriction M3): reject unknown characters
    # instead of silently letting them into the regex. Pre-fix any
    # char not in `_IUPAC_RE` fell through to `c` itself — a custom
    # enzyme site like ``"GAATU"`` (RNA U typo) compiled to a literal
    # ``U`` pattern that never matched DNA, and a stray regex
    # metacharacter (``*``, ``(``, ``?``) became part of the pattern.
    # User-defined enzyme sites are an attack surface; validate here.
    key = site.upper()
    pat = _PATTERN_CACHE.get(key)
    if pat is not None:
        _PATTERN_CACHE.move_to_end(key)
        return pat
    bad = [c for c in key if c not in _IUPAC_RE]
    if bad:
        raise ValueError(
            f"recognition site {site!r} contains non-IUPAC "
            f"character(s) {', '.join(repr(c) for c in bad[:6])}"
            f"{' (truncated)' if len(bad) > 6 else ''}"
        )
    pat = re.compile("".join(_IUPAC_RE[c] for c in key))
    _PATTERN_CACHE[key] = pat
    if len(_PATTERN_CACHE) > _PATTERN_CACHE_MAX:
        _PATTERN_CACHE.popitem(last=False)
    return pat


_IUPAC_COMP = str.maketrans(
    "ACGTRYWSMKBDHVN",
    "TGCAYRWSKMVHDBN",
)

# Case-preserving ACGT complement used by the sequence-panel renderer.
_DNA_COMP_PRESERVE_CASE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


# Tiny LRU on `_rc` — the cached value of a 200 kb plasmid's RC is
# itself ~200 kb, so a 4-entry cap is enough to cover the working set
# (current sequence + a couple of recent rotations / undo snapshots)
# without ballooning RAM. Benches show 35–113× speedup on cache hit
# for cosmid-size sequences (scripts/perf_probe.py).
@functools.lru_cache(maxsize=4)
def _rc(seq: str) -> str:
    return seq.upper().translate(_IUPAC_COMP)[::-1]


# ── Wrap-aware coordinate helpers ─────────────────────────────────────────


def _feat_len(start: int, end: int, total: int) -> int:
    """Circular-aware feature length. A wrap feature (end < start) is
    (total - start) + end bp long; a linear feature is end - start."""
    return (total - start) + end if end < start else end - start


def _seq_len(record) -> int:
    """Length of ``record.seq`` in bp, or 0 if the record has no
    sequence attached. BioPython's ``SeqRecord.seq`` is typed as
    ``Seq | MutableSeq | None`` because the dataclass allows records
    without sequences (rare — e.g. annotation-only GenBank views).
    SpliceCraft always loads records with sequence content, but
    routing length lookups through this helper sidesteps the
    ``"None" is not assignable to "Sized"`` pyright noise at every
    ``len(rec.seq)`` call site without an inline None guard."""
    seq = getattr(record, "seq", None)
    return len(seq) if seq is not None else 0


def _slice_circular(seq: str, start: int, end: int) -> str:
    """Circular-aware slice. If end > start this is a normal slice; if
    end < start the slice wraps the origin and returns seq[start:] + seq[:end].
    end == start is treated as empty (not "wrap whole plasmid") — callers
    that want the latter should pass explicit boundaries. Used by the
    primer-design helpers so a region straddling the origin can be
    primer-designed without special casing at every call site.
    """
    if end >= start:
        return seq[start:end]
    return seq[start:] + seq[:end]


# What is intentionally NOT extracted (yet):
#
#   _scan_restriction_sites + _scan_restriction_sites_impl — depend on
#   `_NEB_ENZYMES` (giant dict), `_RESTR_SCAN_CACHE`, and several module-
#   level worker-side helpers. Not pure.
#
#   _translate_cds — depends on `_GENETIC_CODE` and Biopython's translate.
#   Could move but increases the import surface meaningfully.
#
#   _feat_bounds — touches Biopython's `SeqFeature` / `CompoundLocation`
#   types; extracting would force splicecraft_biology to import Biopython
#   eagerly at module load, costing ~250 ms of startup time on every
#   `splicecraft-cli` call (which only imports for type hints today).
#
#   _bp_in — a METHOD on PlasmidMap, not a module-level function.
#
# Future extractions are welcome but must pass the three-test rule
# in CONTRIBUTING.md: no PlasmidApp coupling, reduces complexity at
# the call site, every existing test passes unchanged.
