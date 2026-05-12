"""
test_alignment_overlay — linear-map alignment overlay helpers.

Covers `_alignment_to_target_segments` and `_alignment_to_target_letters`
— the pure functions that classify each target column of a pairwise
alignment as match / mismatch / gap. They drive the linear-view
alignment lanes (blue / red / gray bars + letters) that overlay
sequencing-read pile-ups on the plasmid map.

Sacred behaviours under test:

  * Target-resolution coordinates — target gaps (insertions in the
    query) consume no target column and don't subdivide the surrounding
    state.
  * Three-state classification — match / mismatch / gap match the
    user-spec 3-color scheme.
  * Case-insensitive matching — the bp comparison ignores case.
  * `t_start` offset — for local alignments with non-zero target
    offsets, segments shift accordingly.
"""
from __future__ import annotations

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# _alignment_to_target_segments
# ═══════════════════════════════════════════════════════════════════════════════

class TestAlignmentToTargetSegments:
    def test_all_match(self):
        assert sc._alignment_to_target_segments("ATGC", "ATGC") == [
            (0, 4, "match"),
        ]

    def test_all_mismatch(self):
        assert sc._alignment_to_target_segments("TTTT", "AAAA") == [
            (0, 4, "mismatch"),
        ]

    def test_match_mismatch_match(self):
        # col 2: G≠C mismatch; flanks all match
        assert sc._alignment_to_target_segments("ATGC", "ATCC") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
            (3, 4, "match"),
        ]

    def test_query_deletion_makes_gap_segment(self):
        # query has 2-bp deletion against target
        assert sc._alignment_to_target_segments("AT--GC", "ATCCGC") == [
            (0, 2, "match"),
            (2, 4, "gap"),
            (4, 6, "match"),
        ]

    def test_target_gap_invisible_at_target_resolution(self):
        # target has a 2-col gap (insertion in query) — consumes zero
        # target columns; surrounding state continues unbroken
        assert sc._alignment_to_target_segments("ATXXGC", "AT--GC") == [
            (0, 4, "match"),
        ]

    def test_target_gap_then_state_change(self):
        # query insertion immediately followed by a mismatch — the
        # mismatch starts at the target column right after the
        # insertion (insertion contributes no target position)
        assert sc._alignment_to_target_segments("ATXG", "AT-C") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
        ]

    def test_t_start_offset(self):
        assert sc._alignment_to_target_segments("ATGC", "ATGC", t_start=100) == [
            (100, 104, "match"),
        ]

    def test_case_insensitive(self):
        assert sc._alignment_to_target_segments("atgc", "ATGC") == [
            (0, 4, "match"),
        ]

    def test_empty(self):
        assert sc._alignment_to_target_segments("", "") == []

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="differ in length"):
            sc._alignment_to_target_segments("ATGC", "ATG")

    def test_complex_mixed_states(self):
        # M M MM M gap M  →  segments at positions 0..2, 2..3, 3..4, 4..5, 5..6
        assert sc._alignment_to_target_segments("ATGC-G", "ATCCAG") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
            (3, 4, "match"),
            (4, 5, "gap"),
            (5, 6, "match"),
        ]

    def test_consecutive_runs_coalesce(self):
        # 3 consecutive matches → one segment, not three
        result = sc._alignment_to_target_segments("AAAA", "AAAA")
        assert len(result) == 1
        assert result[0] == (0, 4, "match")

    def test_all_gap(self):
        assert sc._alignment_to_target_segments("----", "ATGC") == [
            (0, 4, "gap"),
        ]


# ═══════════════════════════════════════════════════════════════════════════════
# _alignment_to_target_letters
# ═══════════════════════════════════════════════════════════════════════════════

class TestAlignmentToTargetLetters:
    def test_all_match(self):
        assert sc._alignment_to_target_letters("ATGC", "ATGC") == {
            0: ("A", "match"),
            1: ("T", "match"),
            2: ("G", "match"),
            3: ("C", "match"),
        }

    def test_mismatch_letter_is_query_base(self):
        # target ATGC, query ATGT — col 3 query says T, target says C
        letters = sc._alignment_to_target_letters("ATGT", "ATGC")
        assert letters[3] == ("T", "mismatch")

    def test_gap_letter_is_dash(self):
        letters = sc._alignment_to_target_letters("AT-G", "ATCG")
        assert letters[2] == ("-", "gap")

    def test_target_gap_skipped(self):
        # target column 1 is a gap — query base at that column never
        # makes it into the per-target dict
        letters = sc._alignment_to_target_letters("ATXG", "A-TG")
        assert letters == {
            0: ("A", "match"),
            1: ("X", "mismatch"),
            2: ("G", "match"),
        }

    def test_t_start_offset(self):
        assert sc._alignment_to_target_letters("AT", "AT", t_start=50) == {
            50: ("A", "match"),
            51: ("T", "match"),
        }

    def test_case_insensitive(self):
        # lowercase query, uppercase target — match classification holds
        letters = sc._alignment_to_target_letters("atgc", "ATGC")
        for pos in range(4):
            _, state = letters[pos]
            assert state == "match"

    def test_empty(self):
        assert sc._alignment_to_target_letters("", "") == {}

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="differ in length"):
            sc._alignment_to_target_letters("AT", "ATG")


# ═══════════════════════════════════════════════════════════════════════════════
# Cross-check: segments and letters agree on state
# ═══════════════════════════════════════════════════════════════════════════════

class TestSegmenterLetterConsistency:
    """The two helpers walk the same gapped strings with the same
    classification — every target column in `letters` must fall inside
    exactly one segment of the matching state.
    """

    @pytest.mark.parametrize("aq,at", [
        ("ATGC",       "ATGC"),
        ("ATGT",       "ATGC"),
        ("AT--GC",     "ATCCGC"),
        ("ATXG",       "AT-C"),
        ("ATGC-G",     "ATCCAG"),
        ("----",       "ATGC"),
        ("ATCGATCG",   "ATCGTTCG"),  # one mismatch in the middle
    ])
    def test_consistency(self, aq, at):
        segs = sc._alignment_to_target_segments(aq, at)
        letters = sc._alignment_to_target_letters(aq, at)
        for t_pos, (_letter, state) in letters.items():
            matching = [
                s for s in segs if s[0] <= t_pos < s[1] and s[2] == state
            ]
            assert len(matching) == 1, (
                f"t_pos={t_pos} state={state!r} not covered by any "
                f"matching segment in {segs!r}"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Registration + lifecycle hardening
# ═══════════════════════════════════════════════════════════════════════════════
# Generation-counter race guards prevent in-flight workers from
# resurrecting cleared alignments. The two assertions below back the
# `_alignments_generation` contract:
#   * `_clear_alignments` ALWAYS bumps the counter (even when the band
#     is already empty) so a worker that hadn't registered yet still
#     gets poisoned by a "preemptive" clear.
#   * `_register_alignment` refuses degenerate input (empty aligned
#     strings) — those would paint nothing and surface as a phantom
#     row.

TERMINAL_SIZE = (160, 48)


class TestAlignmentLifecycle:
    """Pilot-driven tests for the register/clear contract."""

    async def test_clear_bumps_generation_when_non_empty(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Seed one alignment so clear has work to do.
            app._alignments = [{
                "name": "fake", "query_label": "q", "target_label": "t",
                "target_record": tiny_record,
                "result": {"aligned_q": "A", "aligned_t": "A"},
                "aligned_q": "A", "aligned_t": "A",
                "t_start": 0, "segments": [(0, 1, "match")],
                "t_lo": 0, "t_hi": 1, "letters": None,
            }]
            gen_before = app._alignments_generation
            app._clear_alignments()
            assert app._alignments == []
            assert app._alignments_generation == gen_before + 1

    async def test_clear_bumps_generation_when_already_empty(
            self, tiny_record, isolated_library):
        """Empty-band clear still bumps the counter — workers that
        started before the clear must still see the bump and refuse
        to register, even if there was nothing visible to clear."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            gen_before = app._alignments_generation
            assert app._alignments == []
            app._clear_alignments()
            assert app._alignments_generation == gen_before + 1

    async def test_register_rejects_empty_aligned_strings(
            self, tiny_record, isolated_library):
        """Degenerate `_pairwise_align` results (empty aligned_q /
        aligned_t) MUST NOT register — they'd surface as a phantom
        zero-width row."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            assert app._alignments == []
            # All three empty-string permutations should be refused.
            for aq, at in (("", ""), ("ATGC", ""), ("", "ATGC")):
                app._register_alignment(
                    name="empty",
                    query_label="q",
                    target_label="t",
                    target_record=tiny_record,
                    result={"aligned_q": aq, "aligned_t": at},
                )
            assert app._alignments == []

    async def test_register_succeeds_with_valid_result(
            self, tiny_record, isolated_library):
        """Sanity: a valid result lands in the band."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._register_alignment(
                name="ok",
                query_label="q",
                target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            assert len(app._alignments) == 1
            entry = app._alignments[0]
            assert entry["segments"] == [(0, 4, "match")]
            assert entry["t_lo"] == 0 and entry["t_hi"] == 4
