#!/usr/bin/env python3
"""perf_probe_render — keystroke-to-paint micro-bench for the seq-panel.

Measures `_build_seq_text` cost (the bulk of the per-keystroke work)
across realistic scenarios:
  - Steady-state cache hit (rotation to same origin twice)
  - Cursor move
  - Selection change
  - Per-char edit (worst case: changes everything)
  - Rotation simulation (different `id(seq)` but same content)

Run:  python3 scripts/perf_probe_render.py
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
import splicecraft as sc  # noqa: E402


def bench(fn, *, iters: int = 20, warmup: int = 2) -> float:
    """Return median per-call ms across `iters` runs after `warmup`."""
    for _ in range(warmup):
        fn()
    samples = []
    for _ in range(iters):
        gc.collect()
        t0 = time.perf_counter()
        fn()
        samples.append((time.perf_counter() - t0) * 1000.0)
    return statistics.median(samples)


def synth_feats(seq_len: int, n: int = 20, *, seed: int = 0) -> list[dict]:
    rng = random.Random(seed)
    feats = []
    palette = ["#FF8855", "#55AAFF", "#33CC66", "#FFCC00", "#CC55FF"]
    for i in range(n):
        start = rng.randint(0, max(0, seq_len - 500))
        end = min(seq_len, start + rng.randint(100, 500))
        feats.append({
            "type":   rng.choice(["CDS", "promoter", "terminator", "misc_feature"]),
            "start":  start,
            "end":    end,
            "strand": rng.choice([1, -1]),
            "label":  f"f{i}",
            "color":  rng.choice(palette),
        })
    return feats


def synth_seq(bp: int, *, seed: int = 0) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(bp))


def clear_render_caches() -> None:
    sc._BUILD_SEQ_CACHE.clear()
    sc._CHUNK_LAYOUT_CACHE.clear()
    sc._CHUNK_STATIC_CACHE.clear()
    sc._CHUNK_OVERLAY_CACHE.clear()


def cache_state() -> str:
    return (
        f"build={len(sc._BUILD_SEQ_CACHE)}/{sc._BUILD_SEQ_CACHE_MAX}, "
        f"layout={len(sc._CHUNK_LAYOUT_CACHE)}/{sc._CHUNK_LAYOUT_CACHE_MAX}, "
        f"static={len(sc._CHUNK_STATIC_CACHE)}/{sc._CHUNK_STATIC_CACHE_MAX}, "
        f"overlay={len(sc._CHUNK_OVERLAY_CACHE)}/{sc._CHUNK_OVERLAY_CACHE_MAX}"
    )


def header(title: str) -> None:
    print()
    print("═" * 78)
    print(f"  {title}")
    print("═" * 78)


def probe_steady_state() -> None:
    header("Steady-state: same args, repeated — measures cache-hit fast path")
    for bp in [2_686, 10_000, 50_000]:
        seq = synth_seq(bp)
        feats = synth_feats(bp)
        clear_render_caches()
        # warm
        sc._build_seq_text(seq, feats, line_width=60, cursor_pos=100)
        t = bench(
            lambda: sc._build_seq_text(seq, feats, line_width=60, cursor_pos=100),
            iters=20, warmup=2,
        )
        print(f"  seq={bp:>6} bp, 20 feats : {t:>8.2f} ms   ({cache_state()})")


def probe_cursor_move() -> None:
    header("Cursor move: same seq + feats, cursor offset varies — overlay cost")
    for bp in [2_686, 10_000, 50_000]:
        seq = synth_seq(bp)
        feats = synth_feats(bp)
        clear_render_caches()
        # warm with cursor at 0
        sc._build_seq_text(seq, feats, line_width=60, cursor_pos=0)
        counter = {"i": 0}

        def move_cursor():
            counter["i"] = (counter["i"] + 1) % bp
            sc._build_seq_text(seq, feats, line_width=60, cursor_pos=counter["i"])

        t = bench(move_cursor, iters=20, warmup=2)
        print(f"  seq={bp:>6} bp, 20 feats : {t:>8.2f} ms   ({cache_state()})")


def probe_selection_change() -> None:
    header("Selection change: same seq + feats, sel_range varies — overlay invalidation")
    for bp in [2_686, 10_000, 50_000]:
        seq = synth_seq(bp)
        feats = synth_feats(bp)
        clear_render_caches()
        sc._build_seq_text(seq, feats, line_width=60, sel_range=(0, 100))
        counter = {"i": 0}

        def change_sel():
            counter["i"] = (counter["i"] + 1) % (bp - 100)
            sc._build_seq_text(seq, feats, line_width=60,
                                sel_range=(counter["i"], counter["i"] + 100))

        t = bench(change_sel, iters=20, warmup=2)
        print(f"  seq={bp:>6} bp, 20 feats : {t:>8.2f} ms   ({cache_state()})")


def probe_per_char_edit() -> None:
    header("Per-char edit: every call has a DIFFERENT seq — worst case")
    for bp in [2_686, 10_000, 50_000]:
        feats = synth_feats(bp)
        clear_render_caches()
        # warm
        sc._build_seq_text(synth_seq(bp, seed=999), feats, line_width=60)
        counter = {"i": 0}
        base_seq = synth_seq(bp, seed=0)

        def edit_one_char():
            counter["i"] = (counter["i"] + 1) % bp
            mutated = (base_seq[:counter["i"]]
                       + ("A" if base_seq[counter["i"]] != "A" else "T")
                       + base_seq[counter["i"] + 1:])
            sc._build_seq_text(mutated, feats, line_width=60)

        t = bench(edit_one_char, iters=10, warmup=2)
        print(f"  seq={bp:>6} bp, 20 feats : {t:>8.2f} ms   ({cache_state()})")


def probe_rotation_simulation() -> None:
    """A rotation creates a fresh-allocated `disp_seq = seq[origin:] + seq[:origin]`
    and a fresh-allocated `disp_feats` list. The content is structurally
    similar to the un-rotated version (same bases, shifted features) but
    `id()` differs. This probe simulates that to expose the id-keyed cache
    miss the outer view-cache currently pays on every rotation step."""
    header("Rotation: same content, fresh allocations — exposes id() vs hash() cache cost")
    for bp in [2_686, 10_000, 50_000]:
        seq = synth_seq(bp)
        feats = synth_feats(bp)
        clear_render_caches()
        # warm with origin 0
        sc._build_seq_text(seq, feats, line_width=60)
        counter = {"origin": 0}

        def rotate_by_1():
            counter["origin"] = (counter["origin"] + 1) % bp
            origin = counter["origin"]
            disp_seq = seq[origin:] + seq[:origin]
            disp_feats = [
                {**f,
                 "start": (f["start"] - origin) % bp,
                 "end":   (f["end"]   - origin) % bp}
                for f in feats
            ]
            sc._build_seq_text(disp_seq, disp_feats, line_width=60)

        t = bench(rotate_by_1, iters=10, warmup=2)
        print(f"  seq={bp:>6} bp, 20 feats : {t:>8.2f} ms   ({cache_state()})")


def probe_lru_workload() -> None:
    """A user opens 5 plasmids and rotates among them. With cache size 4
    the 5th evicts the 1st. Measures cost of cache-misses-as-eviction."""
    header("LRU workload: 5 plasmids cycled — exposes cache-size limits")
    for bp in [10_000, 50_000]:
        plasmids = [(synth_seq(bp, seed=i), synth_feats(bp, seed=i)) for i in range(5)]
        clear_render_caches()
        # warm each once
        for seq, feats in plasmids:
            sc._build_seq_text(seq, feats, line_width=60)
        counter = {"i": 0}

        def hop():
            counter["i"] = (counter["i"] + 1) % 5
            seq, feats = plasmids[counter["i"]]
            sc._build_seq_text(seq, feats, line_width=60, cursor_pos=10)

        t = bench(hop, iters=20, warmup=2)
        print(f"  seq={bp:>6} bp × 5 plasmids : {t:>8.2f} ms   ({cache_state()})")


def main() -> int:
    print(f"SpliceCraft perf_probe_render — {sc.__version__}")
    print(f"Python {sys.version.split()[0]} on {sys.platform}")

    probe_steady_state()
    probe_cursor_move()
    probe_selection_change()
    probe_per_char_edit()
    probe_rotation_simulation()
    probe_lru_workload()

    print()
    print("Done. Each row: median of 10-20 runs after warmup. Lower is better.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
