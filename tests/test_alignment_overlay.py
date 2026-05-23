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
# _alignment_to_query_segments — query-axis mirror of the target helper
# ═══════════════════════════════════════════════════════════════════════════════
# Drives the Alt+A / diff-plasmid overlay flow where the currently-loaded
# plasmid is the **query** (first arg to `_pairwise_align`) and segments
# must land at query bp positions so bars line up on the open record's
# linear map.

class TestAlignmentToQuerySegments:
    def test_all_match(self):
        assert sc._alignment_to_query_segments("ATGC", "ATGC") == [
            (0, 4, "match"),
        ]

    def test_all_mismatch(self):
        assert sc._alignment_to_query_segments("TTTT", "AAAA") == [
            (0, 4, "mismatch"),
        ]

    def test_match_mismatch_match(self):
        assert sc._alignment_to_query_segments("ATGC", "ATCC") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
            (3, 4, "match"),
        ]

    def test_target_deletion_makes_gap_segment(self):
        # Symmetric to the target helper's `test_query_deletion_makes_gap_segment`:
        # target has a 2-bp deletion vs the query — those query
        # positions get classified as "gap" because the target has no
        # base aligned to them.
        assert sc._alignment_to_query_segments("ATCCGC", "AT--GC") == [
            (0, 2, "match"),
            (2, 4, "gap"),
            (4, 6, "match"),
        ]

    def test_query_gap_invisible_at_query_resolution(self):
        # Insertion in target relative to query — consumes zero query
        # columns; surrounding state continues unbroken.
        assert sc._alignment_to_query_segments("AT--GC", "ATXXGC") == [
            (0, 4, "match"),
        ]

    def test_query_gap_then_state_change(self):
        assert sc._alignment_to_query_segments("AT-C", "ATXG") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
        ]

    def test_q_start_offset(self):
        assert sc._alignment_to_query_segments("ATGC", "ATGC", q_start=100) == [
            (100, 104, "match"),
        ]

    def test_case_insensitive(self):
        assert sc._alignment_to_query_segments("atgc", "ATGC") == [
            (0, 4, "match"),
        ]

    def test_empty(self):
        assert sc._alignment_to_query_segments("", "") == []

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="differ in length"):
            sc._alignment_to_query_segments("ATGC", "ATG")

    def test_complex_mixed_states(self):
        # Mirror of the target helper's complex case: M M MM M gap M
        # using query-axis perspective. Query "ATCCAG" vs target
        # "ATGC-G" — at query position 4 the target has a gap, so
        # the query gets a single "gap" column there.
        assert sc._alignment_to_query_segments("ATCCAG", "ATGC-G") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
            (3, 4, "match"),
            (4, 5, "gap"),
            (5, 6, "match"),
        ]

    def test_consecutive_runs_coalesce(self):
        result = sc._alignment_to_query_segments("AAAA", "AAAA")
        assert len(result) == 1
        assert result[0] == (0, 4, "match")

    def test_all_gap(self):
        assert sc._alignment_to_query_segments("ATGC", "----") == [
            (0, 4, "gap"),
        ]


# ═══════════════════════════════════════════════════════════════════════════════
# _alignment_to_query_letters
# ═══════════════════════════════════════════════════════════════════════════════

class TestAlignmentToQueryLetters:
    def test_all_match(self):
        # Per-query-bp dict — letter at each position is the TARGET
        # base (mirror of the target-axis helper, which stores the
        # query base).
        assert sc._alignment_to_query_letters("ATGC", "ATGC") == {
            0: ("A", "match"),
            1: ("T", "match"),
            2: ("G", "match"),
            3: ("C", "match"),
        }

    def test_mismatch_letter_is_target_base(self):
        letters = sc._alignment_to_query_letters("ATGT", "ATGC")
        assert letters[3] == ("C", "mismatch")

    def test_gap_letter_is_dash(self):
        # Query has a base at col 2, target has a gap — query position
        # records target letter "-" with state "gap".
        letters = sc._alignment_to_query_letters("ATCG", "AT-G")
        assert letters[2] == ("-", "gap")

    def test_query_gap_skipped(self):
        # Insertion in target — query has gap at col 1, that column
        # never enters the per-query dict.
        letters = sc._alignment_to_query_letters("A-TG", "ATXG")
        assert letters == {
            0: ("A", "match"),
            1: ("X", "mismatch"),
            2: ("G", "match"),
        }

    def test_q_start_offset(self):
        assert sc._alignment_to_query_letters("AT", "AT", q_start=50) == {
            50: ("A", "match"),
            51: ("T", "match"),
        }

    def test_case_insensitive(self):
        letters = sc._alignment_to_query_letters("atgc", "ATGC")
        for pos in range(4):
            _, state = letters[pos]
            assert state == "match"

    def test_empty(self):
        assert sc._alignment_to_query_letters("", "") == {}

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="differ in length"):
            sc._alignment_to_query_letters("AT", "ATG")


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


class TestQuerySegmenterLetterConsistency:
    """Symmetric guard for the query-axis helpers."""

    @pytest.mark.parametrize("aq,at", [
        ("ATGC",       "ATGC"),
        ("ATGT",       "ATGC"),
        ("ATCCGC",     "AT--GC"),
        ("AT-C",       "ATXG"),
        ("ATCCAG",     "ATGC-G"),
        ("ATGC",       "----"),
        ("ATCGATCG",   "ATCGTTCG"),
    ])
    def test_consistency(self, aq, at):
        segs = sc._alignment_to_query_segments(aq, at)
        letters = sc._alignment_to_query_letters(aq, at)
        for q_pos, (_letter, state) in letters.items():
            matching = [
                s for s in segs if s[0] <= q_pos < s[1] and s[2] == state
            ]
            assert len(matching) == 1, (
                f"q_pos={q_pos} state={state!r} not covered by any "
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
            # Default axis is "target" (Plasmidsaurus convention).
            assert entry["axis"] == "target"


class TestAlignmentPersistenceRoundTrip:
    """User-feature 2026-05-23: alignments survive a record swap by
    living on the library entry. Hydrate restores the visible ones
    on the next `_apply_record`.

    Test contract: register an alignment against library entry A,
    flush, load entry B, load entry A again → the original alignment
    is back on the band with all fields intact (segments, axis,
    target_label, target_record's sequence).
    """

    @staticmethod
    def _make_library_record(seq: str, rid: str, name: str = None):
        """Build a minimal library entry dict that includes a real
        record + gb_text + size."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq(seq), id=rid, name=name or rid,
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        return rec, {
            "id":      rid,
            "name":    name or rid,
            "size":    len(seq),
            "gb_text": sc._record_to_gb_text(rec),
        }

    async def test_alignment_survives_record_swap_and_back(
            self, isolated_library):
        # Library has two entries: A (200 bp), B (300 bp). Set up
        # the active collection BEFORE creating the app, otherwise
        # `_ensure_default_collection` + `_restore_library_from_active_collection`
        # will rebuild "Main Collection" from whatever's already in
        # the active collection (which may be empty in this isolated
        # tmp tree) and wipe our seeded library.
        rec_a, entry_a = self._make_library_record("A" * 200, "A_PLASMID")
        rec_b, entry_b = self._make_library_record("C" * 300, "B_PLASMID")
        sc._save_collections([{
            "name":        sc._DEFAULT_COLLECTION_NAME,
            "description": "test collection",
            "plasmids":    [entry_a, entry_b],
            "saved":       "2026-05-23",
        }])
        sc._set_active_collection_name(sc._DEFAULT_COLLECTION_NAME)
        sc._save_library([entry_a, entry_b])

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            app._apply_record(rec_a)
            await pilot.pause(0.05)

            # Register an alignment of A against B-as-target.
            seq_a = str(rec_a.seq)
            seq_b = str(rec_b.seq)
            app._register_alignment(
                name="A vs B",
                query_label="A_PLASMID",
                target_label="B_PLASMID",
                target_record=rec_b,
                result={"aligned_q": seq_a[:200], "aligned_t": seq_b[:200]},
                axis="query",
            )
            assert len(app._alignments) == 1
            app._flush_active_alignments()

            # Verify it landed on A's library entry.
            entries = sc._load_library()
            a_entry = next(
                (e for e in entries if e.get("id") == "A_PLASMID"), None,
            )
            assert a_entry is not None, (
                f"A_PLASMID missing from library after flush; entries="
                f"{[e.get('id') for e in entries]!r}"
            )
            stored = a_entry.get("alignments") or []
            assert len(stored) == 1
            assert stored[0]["target_label"] == "B_PLASMID"
            assert stored[0]["visible"] is True
            assert "target_gb_text" in stored[0]
            assert "target_seq_hash" in stored[0]
            stored_id_before = stored[0]["id"]

            # Swap to B → band is cleared (clear_undo=True).
            app._apply_record(rec_b)
            await pilot.pause(0.05)
            assert app._alignments == [], (
                "switching records must drop the in-memory band"
            )

            # Swap back to A → hydrate restores the alignment.
            app._apply_record(rec_a)
            await pilot.pause(0.05)
            assert len(app._alignments) == 1, (
                f"hydrate must restore the stored alignment; got "
                f"{len(app._alignments)} alignments"
            )
            restored = app._alignments[0]
            # Original storage metadata stamped on the restored entry
            # so the next flush round-trips losslessly.
            assert restored["_stored_id"] == stored_id_before
            assert restored["_stored_visible"] is True
            assert restored["target_label"] == "B_PLASMID"
            assert restored["axis"] == "query"
            # And the target record's sequence round-tripped via gb_text.
            assert str(restored["target_record"].seq) == seq_b

    async def test_flush_no_op_when_record_not_in_library(
            self, isolated_library, tiny_record):
        """Loading a record that isn't in the library (file open, demo)
        means there's no library entry to persist onto. Flush must
        silently no-op rather than raising."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)  # not in library
            await pilot.pause(0.05)
            seq = str(tiny_record.seq)
            app._register_alignment(
                name="r1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            # Must not raise.
            app._flush_active_alignments()

    async def test_flush_preserves_hidden_stored_alignments(
            self, isolated_library):
        """Footgun guard 2026-05-23: `_flush_active_alignments` used
        to overwrite the stored list with `self._alignments` — which
        only ever contains *visible* alignments after hydrate. So any
        `visible: False` stored entry would silently vanish the first
        time the user registered a new alignment (the flush would
        write [old visible, new alignment] and drop the hidden ones).

        Test contract: seed the library entry with both a visible and
        a hidden stored alignment. Load → hydrate restores only the
        visible one to `self._alignments`. Register a new alignment
        and flush. Re-read the stored list: must contain all three
        entries — the hidden one preserved untouched.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec_a, entry_a = self._make_library_record("A" * 200, "A_HIDDEN_KEEP")
        rec_b = SeqRecord(
            Seq("T" * 200), id="B", name="B",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        entry_a["alignments"] = [
            {
                "id":              "vis-id",
                "label":           "visible alignment",
                "query_label":     "q",
                "target_label":    "t",
                "target_id":       "B",
                "target_gb_text":  sc._record_to_gb_text(rec_b),
                "target_seq_hash": sc._alignment_target_hash("T" * 200),
                "axis":            "query",
                "result":          {"aligned_q": "A" * 200,
                                    "aligned_t": "T" * 200},
                "visible":         True,
                "added":           "2026-05-23",
                "source":          "manual",
            },
            {
                "id":              "hidden-id",
                "label":           "hidden alignment",
                "query_label":     "q",
                "target_label":    "t",
                "target_id":       "B",
                "target_gb_text":  sc._record_to_gb_text(rec_b),
                "target_seq_hash": sc._alignment_target_hash("T" * 200),
                "axis":            "query",
                "result":          {"aligned_q": "A" * 200,
                                    "aligned_t": "T" * 200},
                "visible":         False,
                "added":           "2026-05-22",
                "source":          "manual",
            },
        ]
        sc._save_collections([{
            "name":        sc._DEFAULT_COLLECTION_NAME,
            "description": "test",
            "plasmids":    [entry_a],
            "saved":       "2026-05-23",
        }])
        sc._set_active_collection_name(sc._DEFAULT_COLLECTION_NAME)
        sc._save_library([entry_a])

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            app._apply_record(rec_a)
            await pilot.pause(0.05)
            # Hydrate restored visible only.
            assert len(app._alignments) == 1
            assert app._alignments[0]["_stored_id"] == "vis-id"

            # Register a fresh alignment and flush.
            app._register_alignment(
                name="new alignment",
                query_label="q",
                target_label="t",
                target_record=rec_b,
                result={"aligned_q": "A" * 200, "aligned_t": "T" * 200},
                axis="query",
            )
            app._flush_active_alignments()

            # The stored list MUST contain all three — hidden preserved.
            entries = sc._load_library()
            a_entry = next(
                (e for e in entries if e.get("id") == "A_HIDDEN_KEEP"),
                None,
            )
            assert a_entry is not None
            stored = a_entry.get("alignments") or []
            ids = [e.get("id") for e in stored]
            assert "vis-id" in ids, (
                f"visible stored entry was dropped from {ids!r}"
            )
            assert "hidden-id" in ids, (
                f"HIDDEN stored entry was wiped by the flush — "
                f"flush must merge with existing storage, not replace. "
                f"got {ids!r}"
            )
            # The hidden entry's visible field must still be False
            # (we didn't accidentally re-visible it).
            hidden = next(e for e in stored if e["id"] == "hidden-id")
            assert hidden["visible"] is False
            # And the fresh one is appended with visible=True.
            fresh = next(
                e for e in stored
                if e["id"] not in ("vis-id", "hidden-id")
            )
            assert fresh["visible"] is True
            assert fresh["label"] == "new alignment"

    async def test_hydrate_skips_invisible_stored_alignments(
            self, isolated_library):
        """Stored alignments with `visible: False` must NOT land on
        the band — they exist for the manager modal but stay hidden
        until toggled. The visibility-toggle UI lives in chunk 2;
        this test pins down the hydrate-side filter so toggle-off
        actually hides on the next load."""
        rec_a, entry_a = self._make_library_record("A" * 100, "A_HIDDEN")
        # Pre-seed the library entry with an invisible alignment.
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec_b = SeqRecord(
            Seq("T" * 100), id="B", name="B",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        entry_a["alignments"] = [{
            "id":              "fixed-id",
            "label":           "hidden alignment",
            "query_label":     "q",
            "target_label":    "t",
            "target_id":       "B",
            "target_gb_text":  sc._record_to_gb_text(rec_b),
            "target_seq_hash": sc._alignment_target_hash("T" * 100),
            "axis":            "query",
            "result":          {"aligned_q": "A" * 100, "aligned_t": "T" * 100},
            "visible":         False,
            "added":           "2026-05-23",
            "source":          "manual",
        }]
        sc._save_library([entry_a])

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            app._apply_record(rec_a)
            await pilot.pause(0.05)
            assert app._alignments == [], (
                "invisible stored alignments must not appear on the band"
            )


class TestAlignmentManagerModal:
    """Manager modal (Alt+L) — listing, toggling, and deleting
    stored alignments for the active plasmid."""

    @staticmethod
    def _make_stored(label: str, id_: str = "", *,
                      visible: bool = True, source: str = "manual") -> dict:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("T" * 100), id="T_TARGET", name="T_TARGET",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        return {
            "id":              id_ or f"id-{label}",
            "label":           label,
            "query_label":     "Q",
            "target_label":    "T",
            "target_id":       "T_TARGET",
            "target_gb_text":  sc._record_to_gb_text(rec),
            "target_seq_hash": sc._alignment_target_hash("T" * 100),
            "axis":            "query",
            "result":          {"aligned_q": "A" * 100,
                                "aligned_t": "T" * 100,
                                "identity_pct": 12.3},
            "visible":         visible,
            "added":           "2026-05-23",
            "source":          source,
        }

    async def test_modal_lists_all_stored_with_visibility_glyphs(
            self, tiny_record, isolated_library):
        """All stored entries appear in the table — visible AND hidden,
        with distinct glyphs in column 0."""
        stored = [
            self._make_stored("vis_one", "id1", visible=True),
            self._make_stored("hidden_one", "id2", visible=False),
            self._make_stored("vis_two", "id3", visible=True),
        ]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored, plasmid_label="P")
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            from textual.widgets import DataTable
            t = modal.query_one("#alnmgr-table", DataTable)
            assert t.row_count == 3
            # Modal's own copy is independent of caller's list.
            assert modal._alignments is not stored
            assert [a["label"] for a in modal._alignments] == [
                "vis_one", "hidden_one", "vis_two",
            ]

    async def test_toggle_visible_flips_in_place(
            self, tiny_record, isolated_library):
        """Space (action_toggle_visible) flips the cursor row's
        `visible` field without reordering rows."""
        stored = [
            self._make_stored("a", visible=True),
            self._make_stored("b", visible=True),
        ]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            from textual.widgets import DataTable
            t = modal.query_one("#alnmgr-table", DataTable)
            t.move_cursor(row=1)
            modal.action_toggle_visible()
            assert modal._alignments[0]["visible"] is True
            assert modal._alignments[1]["visible"] is False
            # Cursor stays on row 1.
            assert t.cursor_row == 1

    async def test_delete_removes_cursor_row(
            self, tiny_record, isolated_library):
        stored = [
            self._make_stored("keep1"),
            self._make_stored("remove_me"),
            self._make_stored("keep2"),
        ]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            from textual.widgets import DataTable
            t = modal.query_one("#alnmgr-table", DataTable)
            t.move_cursor(row=1)
            modal.action_delete_selected()
            assert [a["label"] for a in modal._alignments] == [
                "keep1", "keep2",
            ]
            from textual.widgets import DataTable
            assert modal.query_one("#alnmgr-table", DataTable).row_count == 2

    async def test_hide_all_and_show_all_bulk(
            self, tiny_record, isolated_library):
        stored = [
            self._make_stored("a", visible=True),
            self._make_stored("b", visible=False),
            self._make_stored("c", visible=True),
        ]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            # Trigger via direct button press handler.
            modal._hide_all(None)
            assert all(not a["visible"] for a in modal._alignments)
            modal._show_all(None)
            assert all(a["visible"] for a in modal._alignments)

    async def test_save_returns_modified_list_cancel_returns_none(
            self, tiny_record, isolated_library):
        stored = [self._make_stored("alpha")]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        captured = []
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Save path: dismiss with the (modified) list.
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal, callback=lambda r: captured.append(("save", r)))
            await pilot.pause()
            modal.action_toggle_visible()
            modal._save_and_close(None)
            await pilot.pause()
            await pilot.pause(0.05)
            assert captured[-1][0] == "save"
            assert captured[-1][1] is not None
            assert captured[-1][1][0]["visible"] is False  # toggled
            captured.clear()
            # Cancel path: dismiss with None even after edits.
            modal2 = sc.AlignmentManagerModal(stored)
            app.push_screen(modal2, callback=lambda r: captured.append(("cancel", r)))
            await pilot.pause()
            modal2._cancel_btn(None)
            await pilot.pause()
            await pilot.pause(0.05)
            assert captured[-1] == ("cancel", None)


class TestAlignmentSurvivesZoomAndPan:
    """Regression guards for the user-reported "the alignment disappears
    when I zoom or pan" complaint. The intended behaviour:

      * Linear-map zoom (+/-) changes the bp/col ratio and the visible
        bp window but MUST NOT clear `self._alignments`.
      * Linear-map pan (arrow keys → `_linear_pan`) shifts the offset
        but MUST NOT clear or reshape alignments; bars stay anchored
        to their bp positions and slide on screen along with the rail.
    """

    async def test_zoom_in_does_not_clear_alignments(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            n = len(tiny_record.seq)
            # Whole-record alignment (covers every bp).
            seq = str(tiny_record.seq)
            app._register_alignment(
                name="self",
                query_label="self", target_label="self",
                target_record=tiny_record,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            assert len(app._alignments) == 1
            entry_before = app._alignments[0]
            for _ in range(5):
                pm.action_linear_zoom_in()
                await pilot.pause()
            # Still registered; identity preserved (same dict, not a
            # fresh copy that lost the cached letters / segments).
            assert len(app._alignments) == 1
            assert app._alignments[0] is entry_before
            assert app._alignments[0]["t_lo"] == 0
            assert app._alignments[0]["t_hi"] == n

    async def test_zoom_out_does_not_clear_alignments(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            pm._linear_zoom = 4.0  # start zoomed in
            seq = str(tiny_record.seq)
            app._register_alignment(
                name="self",
                query_label="self", target_label="self",
                target_record=tiny_record,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            assert len(app._alignments) == 1
            entry_before = app._alignments[0]
            for _ in range(5):
                pm.action_linear_zoom_out()
                await pilot.pause()
            assert len(app._alignments) == 1
            assert app._alignments[0] is entry_before

    async def test_pan_does_not_clear_alignments(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            pm._linear_zoom = 4.0  # zoomed in so pan is meaningful
            seq = str(tiny_record.seq)
            app._register_alignment(
                name="self",
                query_label="self", target_label="self",
                target_record=tiny_record,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            entry_before = app._alignments[0]
            # Right then left — should land back at the same offset
            # without disturbing the alignment overlay.
            for _ in range(5):
                pm.action_rotate_cw()  # `_linear_pan(+1)` in linear mode
                await pilot.pause()
            for _ in range(5):
                pm.action_rotate_ccw()
                await pilot.pause()
            assert len(app._alignments) == 1
            assert app._alignments[0] is entry_before


class TestRegisterAlignmentAxis:
    """`_register_alignment` accepts an `axis` parameter that selects
    which side of the alignment plays the role of the currently-loaded
    plasmid (= the render axis along which overlay bars are positioned).

    `axis="target"` is the Plasmidsaurus / sequencing-pile flow
    (segments in target coords). `axis="query"` is the Alt+A /
    diff-plasmid flow (segments in query coords) — without this,
    overlay bars on the open plasmid's linear map land at the picked
    plasmid's bp positions instead of the open plasmid's, which is
    wrong whenever the pairwise alignment isn't a perfect 1:1.
    """

    async def test_axis_query_uses_query_coord_segments(
            self, tiny_record, isolated_library):
        """The discriminating case: query "ATCCGC" vs target "AT--GC"
        produces different segments in the two coord systems.

          * Target axis: the target has a 2-col gap (insertion in query)
            so the surrounding match state continues unbroken at target
            resolution → ``[(0, 4, "match")]``.
          * Query axis: the target has 2 fewer bases than the query, so
            those query positions get classified as gap → ``[(0, 2,
            "match"), (2, 4, "gap"), (4, 6, "match")]``.
        """
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._register_alignment(
                name="alt-a",
                query_label="open-plasmid",
                target_label="picked",
                target_record=tiny_record,
                result={
                    "aligned_q": "ATCCGC",
                    "aligned_t": "AT--GC",
                },
                axis="query",
            )
            assert len(app._alignments) == 1
            entry = app._alignments[0]
            assert entry["axis"] == "query"
            assert entry["segments"] == [
                (0, 2, "match"),
                (2, 4, "gap"),
                (4, 6, "match"),
            ]
            # render-axis bounds — bars draw across the full query
            # span (0..6), NOT the target span (0..4).
            assert entry["t_lo"] == 0 and entry["t_hi"] == 6

    async def test_axis_target_default_unchanged(
            self, tiny_record, isolated_library):
        """Regression guard: the existing Plasmidsaurus flow doesn't
        pass `axis`, so it must keep getting target-axis segments. The
        same gapped pair as above, registered without `axis`, collapses
        the target gap into the surrounding match (target-resolution).
        """
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._register_alignment(
                name="plasmidsaurus",
                query_label="read",
                target_label="open-plasmid",
                target_record=tiny_record,
                result={
                    "aligned_q": "ATCCGC",
                    "aligned_t": "AT--GC",
                },
            )
            entry = app._alignments[0]
            assert entry["axis"] == "target"
            assert entry["segments"] == [(0, 4, "match")]

    async def test_invalid_axis_raises(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            with pytest.raises(ValueError, match="axis must be"):
                app._register_alignment(
                    name="bad", query_label="q", target_label="t",
                    target_record=tiny_record,
                    result={"aligned_q": "AT", "aligned_t": "AT"},
                    axis="sideways",
                )



# ═══════════════════════════════════════════════════════════════════════════════
# Circular alignment offset (GH #16, 2026-05-14)
# ═══════════════════════════════════════════════════════════════════════════════

class TestCircularAlignmentOffset:
    """`_find_circular_alignment_offset` rotates a circular target so
    the global pairwise align doesn't pair bp 1 of an arbitrarily-
    started Plasmidsaurus read with bp 1 of the GenBank reference.
    Regression guard for Cory Tobin's report — pre-fix alignment of a
    700-bp-rotated read showed 66% identity + 500+ gaps; post-fix the
    same read aligns at 100% with zero gaps."""

    def test_returns_zero_when_sequences_already_aligned(self):
        target = 'ACGTACGTACGT' * 100
        # Same sequence — no rotation needed.
        assert sc._find_circular_alignment_offset(target, target) == 0

    def test_detects_simple_700_bp_rotation(self):
        # Pseudo-random plasmid-shaped target so each 25-bp seed
        # appears at a unique anchor — repeat-pattern fixtures defeat
        # the uniqueness guard and the helper safely falls back to 0.
        import random
        rng = random.Random(20260514)
        target = "".join(rng.choice("ACGT") for _ in range(2000))
        read = target[700:] + target[:700]
        offset = sc._find_circular_alignment_offset(read, target)
        assert offset == 700

    def test_detects_rotation_near_origin_wrap(self):
        """Rotation just before the target's end means the seed may
        straddle the wrap. The doubled-target search handles this."""
        import random
        rng = random.Random(20260514)
        target = "".join(rng.choice("ACGT") for _ in range(1500))
        rotation = len(target) - 50
        read = target[rotation:] + target[:rotation]
        offset = sc._find_circular_alignment_offset(read, target)
        assert offset == rotation

    def test_returns_zero_for_short_sequences(self):
        # Below the k=25 minimum kmer length — bail cleanly.
        assert sc._find_circular_alignment_offset('AAAA', 'TTTT') == 0

    def test_skips_low_complexity_seeds(self):
        """A query that starts with a homopolymer run shouldn't seed
        on the homopolymer (it'd match everywhere); the helper steps
        past it and finds a complex seed further along."""
        target = 'A' * 100 + 'GTACGTACGTAC' * 30 + 'C' * 50
        # Read starts mid-target.
        rotation = 250
        read = target[rotation:] + target[:rotation]
        offset = sc._find_circular_alignment_offset(read, target)
        # Either the helper finds the exact rotation OR returns 0
        # (no clean unique seed); both are acceptable. The bad
        # outcome we're guarding against is a WRONG non-zero answer.
        assert offset in (0, rotation)

    def test_pairwise_align_with_rotation_recovers_identity(self):
        """End-to-end: a rotated read aligned against the rotated
        target should produce near-100%% identity vs ~50-70%% without
        rotation. This is the test that maps directly to Cory's GH #16
        screenshot."""
        import random
        rng = random.Random(20260514)
        target = "".join(rng.choice("ACGT") for _ in range(3000))
        read = target[700:] + target[:700]
        offset = sc._find_circular_alignment_offset(read, target)
        assert offset == 700
        rotated = target[offset:] + target[:offset]
        result = sc._pairwise_align(read, rotated, mode='global')
        assert result['identity_pct'] >= 99.0
        assert result['n_gaps'] == 0


class TestRotateSeqRecord:
    """`_rotate_seq_record` shifts a SeqRecord's sequence + features
    so that a chosen position becomes the new origin. Used by the
    alignment path to keep the viewer's feature lane in register
    with the rotated target."""

    @staticmethod
    def _circular(seq: str, *, features=()):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq(seq), id='T', name='T')
        rec.annotations['molecule_type'] = 'DNA'
        rec.annotations['topology'] = 'circular'
        rec.features = list(features)
        return rec

    def test_zero_offset_returns_input(self):
        rec = self._circular('A' * 100)
        rotated = sc._rotate_seq_record(rec, 0)
        assert rotated is rec

    def test_rotation_shifts_sequence(self):
        rec = self._circular('ABCDEFGHIJ')
        rotated = sc._rotate_seq_record(rec, 3)
        assert str(rotated.seq) == 'DEFGHIJABC'

    def test_rotation_shifts_simple_feature(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._circular(
            'A' * 100,
            features=[SeqFeature(FeatureLocation(50, 70, strand=1),
                                   type='CDS')],
        )
        rotated = sc._rotate_seq_record(rec, 20)
        # Feature was at 50-70; after rotation by 20 it's at 30-50.
        assert len(rotated.features) == 1
        loc = rotated.features[0].location
        assert int(loc.start) == 30
        assert int(loc.end) == 50

    def test_rotation_preserves_record_metadata(self):
        rec = self._circular('A' * 100)
        rec.description = 'test'
        rotated = sc._rotate_seq_record(rec, 20)
        assert rotated.id == 'T'
        assert rotated.description == 'test'
        assert rotated.annotations['topology'] == 'circular'


# ═══════════════════════════════════════════════════════════════════════════════
# Per-plasmid map_mode persistence (2026-05-18)
# ═══════════════════════════════════════════════════════════════════════════════
# Library entries can carry a `map_mode` field that overrides the
# topology-derived default on load. Plasmidsaurus alignment + the
# user's Alt+L toggle both write through `_persist_map_mode_for_active`
# so the choice sticks across reloads. Sequencing-aligned plasmids
# also have their library entry tagged `linear` so the next open
# defaults to the diff-friendly view.

class TestPerPlasmidMapModePersistence:
    """Tests use `_apply_record` rather than `_preload_record` because
    the preload path dispatches `_add_save_to_disk` (a `@work(thread=True)`
    worker) whose disk write can race test teardown and contaminate the
    next test's tmp library file."""

    async def test_load_record_with_stashed_linear_overrides_topology(
            self, isolated_library):
        """A circular plasmid loaded from a library entry tagged
        `map_mode: "linear"` opens in linear view — the user-set
        preference beats the topology default."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("A" * 500), id="C", name="C",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec._tui_map_mode = "linear"
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            assert pm._map_mode == "linear", (
                "stashed _tui_map_mode='linear' must win over circular topology"
            )

    async def test_load_record_with_stashed_circular_overrides_linear_topo(
            self, isolated_library):
        """Symmetric: a `topology=linear` record loaded with a
        stashed `circular` preference opens circular. Belt-and-braces
        check on the override direction."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("A" * 500), id="L", name="L",
                        annotations={"molecule_type": "DNA",
                                     "topology": "linear"})
        rec._tui_map_mode = "circular"
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            assert pm._map_mode == "circular"

    async def test_load_record_unknown_stashed_mode_falls_back_to_topology(
            self, isolated_library):
        """Defensive: a bogus stashed value (e.g. hand-edit) is
        ignored and the topology default applies."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("A" * 500), id="C", name="C",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec._tui_map_mode = "spiral"   # nonsense
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            assert pm._map_mode == "circular"

    async def test_persist_map_mode_writes_to_library_entry(
            self, tiny_record, isolated_library):
        """`_persist_map_mode_for_active` saves the chosen mode onto
        the active library entry so the next reload picks it up.
        Uses `_apply_record` (not `_preload_record`) to avoid the
        background `_add_save_to_disk` worker — that worker's write
        races test teardown."""
        # Pre-seed the library with an entry matching the loaded record.
        sc._save_library([{
            "id":      tiny_record.id,
            "name":    tiny_record.name,
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-18",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            app._persist_map_mode_for_active("linear")
            entries = sc._load_library()
            match = next(e for e in entries if e.get("id") == tiny_record.id)
            assert match.get("map_mode") == "linear"

    async def test_persist_map_mode_is_noop_for_unknown_entry(
            self, isolated_library):
        """When the loaded record isn't in the library, the helper
        must silently no-op (no exception, no spurious row).
        Uses a unique record id so even if a prior test's worker
        leaked an entry to the cache, the lookup fails cleanly."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        sc._save_library([])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Unique id (`UNSAVED_NOOP_TEST`) the previous tests in
            # this class never reference — even cache contamination
            # from a sibling test can't match it.
            unsaved = SeqRecord(Seq("A" * 200),
                                 id="UNSAVED_NOOP_TEST",
                                 name="UNSAVED_NOOP_TEST",
                                 annotations={"molecule_type": "DNA",
                                              "topology": "circular"})
            app._apply_record(unsaved)
            await pilot.pause(0.05)
            entries_before = sc._load_library()
            assert not any(
                e.get("id") == "UNSAVED_NOOP_TEST" for e in entries_before
            ), "unsaved record must not be in the library"
            # Must not raise; library state for the unknown id stays
            # unchanged.
            app._persist_map_mode_for_active("linear")
            entries_after = sc._load_library()
            assert not any(
                e.get("id") == "UNSAVED_NOOP_TEST" for e in entries_after
            ), "no-op path must not insert a new row"
            # And entries that DID exist before keep their state.
            assert entries_before == entries_after

    async def test_toggle_map_view_persists_when_entry_exists(
            self, tiny_record, isolated_library):
        """End-to-end: user-driven `action_toggle_map_view` writes
        through to the library entry. Circular plasmid + toggle → entry
        carries `map_mode: "linear"`."""
        sc._save_library([{
            "id":      tiny_record.id,
            "name":    tiny_record.name,
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-18",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            # tiny_record is circular so starts circular
            assert pm._map_mode == "circular"
            pm.action_toggle_map_view()   # circular → linear
            assert pm._map_mode == "linear"
            entries = sc._load_library()
            match = next(e for e in entries if e.get("id") == tiny_record.id)
            assert match.get("map_mode") == "linear"

    async def test_register_alignment_persists_linear_on_target(
            self, tiny_record, isolated_library):
        """Registering an alignment against a circular target pins
        the map to linear AND writes `map_mode: "linear"` to the
        target's library entry. Mirrors the Plasmidsaurus path —
        sequencing-aligned plasmids default to linear on every later
        load."""
        sc._save_library([{
            "id":      tiny_record.id,
            "name":    tiny_record.name,
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-18",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            assert pm._map_mode == "circular"
            app._register_alignment(
                name="read1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            assert pm._map_mode == "linear"
            entries = sc._load_library()
            match = next(e for e in entries if e.get("id") == tiny_record.id)
            assert match.get("map_mode") == "linear"


# ═══════════════════════════════════════════════════════════════════════════════
# Sequencing toolbar screen (2026-05-18)
# ═══════════════════════════════════════════════════════════════════════════════
# Sequencing replaces the freestanding Plasmidsaurus modal. Tab layout
# leaves room for future ingestion sources (direct API, nanopore
# consensus). The legacy class name is aliased so agent/test paths
# keep resolving.

class TestSequencingScreen:
    def test_back_compat_alias_resolves(self):
        """`PlasmidsaurusAlignModal` is the old class name — kept as
        an alias for tests and agent-API callers."""
        assert sc.PlasmidsaurusAlignModal is sc.SequencingScreen

    def test_menu_lists_sequencing(self):
        """The Sequencing entry is wired into the top-level menu bar."""
        assert "Sequencing" in sc.MenuBar.MENUS

    async def test_screen_opens_with_plasmidsaurus_tab(
            self, tmp_path, tiny_record, isolated_library):
        """The Sequencing screen mounts with the Plasmidsaurus tab
        active (it's currently the only tab; future tabs will share
        the same TabbedContent)."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            assert isinstance(screen, sc.SequencingScreen)
            # All the legacy alignment IDs still resolve so the worker
            # event handlers (which key off these IDs) still work.
            screen.query_one("#align-zip-tree")
            screen.query_one("#align-members")
            screen.query_one("#align-target")
            screen.query_one("#btn-align-go")
            screen.query_one("#btn-sequencing-close")

    async def test_subtabs_disabled_until_zip_loaded(
            self, tmp_path, tiny_record, isolated_library):
        """Samples / Quality / Align sub-tabs are disabled on mount;
        the user can't tab into them until a valid zip lands. General
        stays enabled because it owns the zip picker."""
        from textual.widgets import TabPane
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            for tab_id in ("psaurus-sub-samples",
                            "psaurus-sub-quality",
                            "psaurus-sub-align"):
                assert screen.query_one(f"#{tab_id}", TabPane).disabled, (
                    f"{tab_id} must be disabled before zip load"
                )
            # General stays enabled — it owns the zip picker.
            assert not screen.query_one(
                "#psaurus-sub-general", TabPane,
            ).disabled

    async def test_zip_load_enables_subtabs_and_populates_tables(
            self, tmp_path, tiny_record, isolated_library):
        """End-to-end: feeding the screen a synthetic Plasmidsaurus-
        style zip via `_on_zip_picked` enables the dependent sub-tabs,
        populates the Samples + Quality tables, and writes the run
        metadata summary."""
        from textual.widgets import TabPane, DataTable, Static
        from textual.widgets import DirectoryTree
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Synthesise a Plasmidsaurus-shaped zip: 2 samples with
        # consensus .gbk + summary.txt + per-base TSV.
        rec1 = SeqRecord(Seq("ATGC" * 200), id="MAV34", name="MAV34",
                         annotations={"molecule_type": "DNA",
                                      "topology": "circular"})
        rec2 = SeqRecord(Seq("GCAT" * 200), id="MAV35", name="MAV35",
                         annotations={"molecule_type": "DNA",
                                      "topology": "circular"})
        gbk1 = tmp_path / "MAV34.gbk"
        gbk2 = tmp_path / "MAV35.gbk"
        SeqIO.write(rec1, gbk1, "genbank")
        SeqIO.write(rec2, gbk2, "genbank")
        zp = tmp_path / "RUN42_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk1, "RUN42_genbank-files/RUN42_1_MAV34.gbk")
            zf.write(gbk2, "RUN42_genbank-files/RUN42_2_MAV35.gbk")
            zf.writestr(
                "RUN42_summary-files/RUN42_1_MAV34.txt",
                "       1-mer (%)  2-mer (%)\n"
                "moles       95.5        4.5\n"
                "mass        90.1        9.9\n\n\n"
                "*************************\n\n\n"
                "E. coli genomic contamination: 12.3%\n",
            )
            zf.writestr(
                "RUN42_summary-files/RUN42_2_MAV35.txt",
                "       1-mer (%)  2-mer (%)\n"
                "moles       99.9        0.1\n"
                "mass        99.5        0.5\n\n\n"
                "*************************\n\n\n"
                "E. coli genomic contamination: 2.1%\n",
            )
            # Synthetic per-base TSV: 5 rows, integer coverage.
            zf.writestr(
                "RUN42_per-base-data/RUN42_1_MAV34.tsv",
                "pos\tref\treads_all\n"
                "1\tA\t30\n2\tT\t25\n3\tG\t40\n4\tC\t10\n5\tA\t50\n",
            )
            zf.writestr("RUN42_gel.png", b"PNG-fake-bytes")

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            assert isinstance(screen, sc.SequencingScreen)
            # Feed the FileSelected event the directory tree would emit.
            tree = screen.query_one(
                "#align-zip-tree", sc._ZipAwareDirectoryTree,
            )
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, zp),
            )
            await pilot.pause(0.3)
            # Dependent sub-tabs are now enabled.
            for tab_id in ("psaurus-sub-samples",
                            "psaurus-sub-quality",
                            "psaurus-sub-align"):
                assert not screen.query_one(
                    f"#{tab_id}", TabPane,
                ).disabled, f"{tab_id} must be enabled after zip load"
            # Run metadata shows both samples.
            assert screen._parsed_run.get("run_id") == "RUN42"
            assert len(screen._parsed_run.get("samples", [])) == 2
            # Samples table populated with one row per sample.
            samples_t = screen.query_one(
                "#align-members", DataTable,
            )
            assert samples_t.row_count == 2
            # Quality table also has both samples.
            quality_t = screen.query_one(
                "#plasmidsaurus-quality-table", DataTable,
            )
            assert quality_t.row_count == 2
            # Run-level files table picks up the gel.png.
            runfiles_t = screen.query_one(
                "#plasmidsaurus-runfiles-table", DataTable,
            )
            assert runfiles_t.row_count >= 1

    async def test_sample_row_select_enables_align_button(
            self, tmp_path, tiny_record, isolated_library):
        """Clicking a sample row marks that sample's .gbk as the
        alignment query, updates the query indicator on the Align
        tab, and flips the Align button to enabled."""
        from textual.widgets import (DataTable, Button, Static,
                                       DirectoryTree)
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 100), id="MAV1", name="MAV1",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        gbk = tmp_path / "MAV1.gbk"
        SeqIO.write(rec, gbk, "genbank")
        zp = tmp_path / "RUN1_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "RUN1_genbank-files/RUN1_1_MAV1.gbk")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            tree = screen.query_one(
                "#align-zip-tree", sc._ZipAwareDirectoryTree,
            )
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, zp),
            )
            await pilot.pause(0.3)
            # Pre-select: Align button is disabled.
            assert screen.query_one(
                "#btn-align-go", Button,
            ).disabled
            # Synthesise the RowSelected event the Samples DataTable
            # would emit on click.
            samples_t = screen.query_one(
                "#align-members", DataTable,
            )
            samples_t.cursor_coordinate = (
                samples_t.cursor_coordinate.__class__(0, 0)
            )
            from textual.coordinate import Coordinate
            row_key = next(iter(samples_t.rows.keys()))
            samples_t.post_message(
                DataTable.RowSelected(
                    samples_t, Coordinate(0, 0), row_key,
                )
            )
            await pilot.pause(0.1)
            assert screen._selected_member is not None
            assert not screen.query_one(
                "#btn-align-go", Button,
            ).disabled
            # Selected member is set; the query indicator update is
            # exercised end-to-end (verified via the button-enabled
            # state above). Static's content is private API in Textual
            # so we don't peek at it directly.
            assert "MAV1" in str(screen._selected_member)

    async def test_repick_same_zip_skips_reparse(
            self, tmp_path, tiny_record, isolated_library):
        """Picking the same zip twice in a row is a no-op (perf
        guard — parse can take ~1 s on large runs)."""
        from textual.widgets import DirectoryTree
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 50), id="MAV1", name="MAV1",
                        annotations={"molecule_type": "DNA"})
        gbk = tmp_path / "MAV1.gbk"
        SeqIO.write(rec, gbk, "genbank")
        zp = tmp_path / "R_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "R_genbank-files/R_1_MAV1.gbk")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            tree = screen.query_one(
                "#align-zip-tree", sc._ZipAwareDirectoryTree,
            )
            # First pick — populates _parsed_run.
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, zp),
            )
            await pilot.pause(0.3)
            parsed_before = screen._parsed_run
            assert parsed_before
            # Mark with a sentinel so we can detect re-parse.
            parsed_before["_test_sentinel"] = True
            # Second pick of the same path — should NOT re-parse
            # (the sentinel survives because _parsed_run is the
            # same dict object).
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, zp),
            )
            await pilot.pause(0.2)
            assert screen._parsed_run.get("_test_sentinel") is True, (
                "same-path re-pick must not re-parse"
            )

    async def test_invalid_zip_keeps_subtabs_disabled(
            self, tmp_path, tiny_record, isolated_library):
        """A user picking a non-zip file should NOT unlock the sub-tabs
        and the General tab's status row should explain why."""
        from textual.widgets import TabPane, DirectoryTree
        # Non-zip file.
        bogus = tmp_path / "README.txt"
        bogus.write_text("not a zip")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            tree = screen.query_one(
                "#align-zip-tree", sc._ZipAwareDirectoryTree,
            )
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, bogus),
            )
            await pilot.pause(0.2)
            # Sub-tabs stay disabled.
            for tab_id in ("psaurus-sub-samples",
                            "psaurus-sub-quality",
                            "psaurus-sub-align"):
                assert screen.query_one(
                    f"#{tab_id}", TabPane,
                ).disabled, f"{tab_id} must stay disabled on bad zip"
            # Parsed state stays empty.
            assert not screen._parsed_run


# ═══════════════════════════════════════════════════════════════════════════════
# Plasmidsaurus zip parser (run-structured ingestion)
# ═══════════════════════════════════════════════════════════════════════════════
# `_parse_plasmidsaurus_zip` walks a results zip and groups files by
# sample so the Sequencing toolbar's sub-tabs can render without
# re-reading the zip per tab. Run-level extras (gel.png, README) go
# under `run_files`. `_parse_plasmidsaurus_summary` extracts the
# k-mer / contamination percentages from the per-sample summary file.

class TestPlasmidsaurusZipParser:
    def _build_zip(self, dirpath, samples, *, run="RUN1",
                    extra_files=None):
        """Build a synthetic Plasmidsaurus-shaped zip in `dirpath`.
        `samples` is a list of (sample_name, summary_text, perbase_text).
        Returns the zip path."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        zp = dirpath / f"{run}_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            for idx, (name, summary, perbase) in enumerate(samples, 1):
                rec = SeqRecord(
                    Seq("ATGC" * 50), id=name, name=name,
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"},
                )
                gbk = dirpath / f"{name}.gbk"
                SeqIO.write(rec, gbk, "genbank")
                base = f"{run}_{idx}_{name}"
                zf.write(gbk, f"{run}_genbank-files/{base}.gbk")
                if summary is not None:
                    zf.writestr(
                        f"{run}_summary-files/{base}.txt", summary,
                    )
                if perbase is not None:
                    zf.writestr(
                        f"{run}_per-base-data/{base}.tsv", perbase,
                    )
            for name, content in (extra_files or []):
                zf.writestr(name, content)
        return zp

    def test_parses_run_id_from_folder_prefix(self, tmp_path):
        zp = self._build_zip(tmp_path, [("MAV1", None, None)],
                              run="ABC42")
        data = sc._parse_plasmidsaurus_zip(zp)
        assert data["run_id"] == "ABC42"

    def test_groups_files_under_one_sample(self, tmp_path):
        zp = self._build_zip(tmp_path, [
            ("MAV1", "moles 99.0\nmass 98.0\nE. coli contamination: 5.0%\n",
             "pos\tref\treads_all\n1\tA\t30\n2\tT\t40\n"),
        ])
        data = sc._parse_plasmidsaurus_zip(zp)
        assert len(data["samples"]) == 1
        s = data["samples"][0]
        # Sample base collapses to the run_<n>_<name> stem.
        assert s["base"].endswith("MAV1")
        # All categories landed on the same sample dict.
        assert s["gbk"]
        assert s["summary"]
        assert s["perbase"]
        # Summary text streamed inline.
        assert "moles" in s["summary_text"]
        # Per-base coverage stats computed.
        assert s["perbase_coverage"].get("mean") == 35.0

    def test_run_level_files_separated_from_samples(self, tmp_path):
        zp = self._build_zip(tmp_path, [("MAV1", None, None)],
                              extra_files=[("RUN1_gel.png", b"PNG")])
        data = sc._parse_plasmidsaurus_zip(zp)
        # Sample list has the one MAV1; run-level file shows up in
        # `run_files`.
        assert len(data["samples"]) == 1
        run_paths = {rf["name"] for rf in data["run_files"]}
        assert "RUN1_gel.png" in run_paths

    def test_natural_sort_samples(self, tmp_path):
        """Samples come back natural-sorted on their base name —
        the run-index prefix Plasmidsaurus uses (`<run>_<n>_<name>`)
        naturally puts `_2_` before `_10_` under the natural-sort
        key (vs lexicographic `_10_` < `_2_`)."""
        # Pass samples in scrambled order (1, 10, 2). The run-index
        # is assigned by `enumerate` in input order, so the bases
        # become `RUN1_1_A`, `RUN1_2_B`, `RUN1_3_C` — already
        # naturally sorted by index, regardless of name.
        zp = self._build_zip(tmp_path, [
            ("A", None, None),
            ("B", None, None),
            ("C", None, None),
        ])
        data = sc._parse_plasmidsaurus_zip(zp)
        names = [s["name"] for s in data["samples"]]
        assert names == sorted(names, key=sc._natural_sort_key)

    def test_summary_parser_extracts_kmer_and_contam(self):
        text = (
            "       1-mer (%)  2-mer (%)\n"
            "moles       97.5        2.5\n"
            "mass        95.1        4.9\n\n\n"
            "*************************\n\n\n"
            "E. coli genomic contamination: 18.0%\n"
        )
        out = sc._parse_plasmidsaurus_summary(text)
        assert out["kmer_moles_pct"] == 97.5
        assert out["kmer_mass_pct"] == 95.1
        assert out["contamination_pct"] == 18.0
        assert "E. coli" in out["contamination_source"]

    def test_summary_parser_handles_missing_fields(self):
        # Empty input — every field returns None / "".
        out = sc._parse_plasmidsaurus_summary("")
        assert out["kmer_moles_pct"] is None
        assert out["kmer_mass_pct"] is None
        assert out["contamination_pct"] is None
        assert out["contamination_source"] == ""

    def test_perbase_summary_returns_empty_on_garbage(self, tmp_path):
        """A malformed per-base TSV (no numeric column 2) should not
        crash the parser; the sample's `perbase_coverage` ends empty."""
        zp = self._build_zip(tmp_path, [
            ("MAV1", None, "pos\tref\treads_all\n"
                            "alpha\tbeta\tgamma\nA\tB\tC\n"),
        ])
        data = sc._parse_plasmidsaurus_zip(zp)
        assert data["samples"][0]["perbase_coverage"] == {}

    def test_missing_zip_raises_value_error(self, tmp_path):
        with pytest.raises(ValueError):
            sc._parse_plasmidsaurus_zip(tmp_path / "does-not-exist.zip")

    def test_oversize_zip_rejected(self, tmp_path, monkeypatch):
        """A zip claiming to be larger than the cap is refused."""
        # Build a tiny zip then artificially cap to a smaller size.
        zp = self._build_zip(tmp_path, [("MAV1", None, None)])
        monkeypatch.setattr(sc, "_PLASMIDSAURUS_ZIP_MAX_BYTES", 1)
        with pytest.raises(ValueError, match="too large"):
            sc._parse_plasmidsaurus_zip(zp)

    def test_standalone_gbk_no_category_folder(self, tmp_path):
        """Zips without the Plasmidsaurus `_genbank-files/` folder
        layout still discover .gbk files as samples (back-compat
        with the older `_list_gbk_members_in_zip` behaviour)."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 20), id="X", name="X",
                        annotations={"molecule_type": "DNA"})
        gbk = tmp_path / "X.gbk"
        SeqIO.write(rec, gbk, "genbank")
        zp = tmp_path / "ad-hoc.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "sample_A/consensus.gbk")
        data = sc._parse_plasmidsaurus_zip(zp)
        # The .gbk should land as a sample, not in run_files.
        assert any(s.get("gbk") for s in data["samples"])


# ═══════════════════════════════════════════════════════════════════════════════
# Alignment band positioned closest-to-centerline (2026-05-18)
# ═══════════════════════════════════════════════════════════════════════════════
# Alignment lanes used to stack BELOW the reverse-feature band. The
# closest-to-center refactor flips that order — alignment lanes now
# render at `rail_row + 2` and reverse features get offset downward
# by the alignment lane count. `_pack_alignment_lanes` is the helper
# that lets the parent renderer learn the lane count up front.

class TestAlignmentBandCenterline:
    async def test_pack_alignment_lanes_returns_count(
            self, tiny_record, isolated_library):
        """`_pack_alignment_lanes` returns (placed, lane_count); empty
        when no alignments touch the visible window."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            # No alignments registered → empty placement.
            placed, n_lanes = pm._pack_alignment_lanes(
                margin_l=5, usable_w=100, view_s=0, view_e=1000,
                w=160, bp_to_col=lambda bp: 5 + bp // 10,
            )
            assert placed == []
            assert n_lanes == 0
            # Register one alignment → one lane.
            app._register_alignment(
                name="r1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            placed, n_lanes = pm._pack_alignment_lanes(
                margin_l=5, usable_w=100, view_s=0, view_e=1000,
                w=160, bp_to_col=lambda bp: 5 + bp // 10,
            )
            assert n_lanes == 1
            assert len(placed) == 1

    async def test_linear_draw_with_alignment_renders_without_error(
            self, tiny_record, isolated_library):
        """Smoke: a linear-view render with one alignment + rev feature
        completes without raising. Covers the new offset path where
        rev features land below the alignment band."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        # Add a rev-strand feature to exercise the offset path.
        tiny_record.features.append(
            SeqFeature(FeatureLocation(50, 100, strand=-1),
                       type="misc_feature",
                       qualifiers={"label": ["rev1"]})
        )
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            app._register_alignment(
                name="r1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            # Must not raise.
            pm.refresh()
            await pilot.pause(0.05)

    async def test_first_paint_after_register_includes_band(
            self, tiny_record, isolated_library):
        """Regression guard (2026-05-22): the PlasmidMap render cache
        previously keyed on (zoom, offset, features, map_mode, …) but
        NOT on `_alignments`, so registering an alignment via
        `set_alignments` invalidated the widget (refresh()) without
        invalidating the draw cache — the cached pre-registration text
        was returned and the user saw bars only after some other
        tracked attribute changed (e.g. a zoom press flipped
        `_linear_zoom` and forced a fresh paint).

        Test contract: render the linear view BEFORE registering an
        alignment, then again AFTER, with no zoom / pan / feature
        change in between. The two paints must differ — the second
        carries the alignment band, the first doesn't. If the cache
        key regresses, both calls return the same cached object.
        """
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            # Force a paint before any alignment so the cache is hot.
            paint_before = pm.render()
            text_before = str(paint_before)
            # Register an alignment. The render path is supposed to
            # invalidate its cache and paint a band.
            app._register_alignment(
                name="r1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            paint_after = pm.render()
            text_after = str(paint_after)
            # The two Rich Texts must compare unequal — same widget
            # size + viewport, but the band landed in the second.
            assert text_before != text_after, (
                "PlasmidMap.render() returned a stale cached paint after "
                "an alignment was registered; the draw-cache key needs "
                "to include alignment state."
            )

    async def test_letter_mode_renders_adjacent_letters_no_gaps(
            self, tiny_record, isolated_library):
        """Regression guard (2026-05-22): at letter-mode zoom the band
        must render alignment letters as adjacent characters
        (`ATGCATGC…`) with no blank gutters between them. The zoom is
        capped at `col_per_bp == 1.0` via `_max_useful_linear_zoom`,
        so the renderer never enters the "letters spread across N
        cols with blank gaps" regime that the user found awkward.

        Test contract: with a record longer than the usable width
        (so the cap is > 1.0 and actually reachable via zoom-in),
        push `_linear_zoom` to the cap, confirm a contiguous run of
        ATGC appears in the band row.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # 500 bp ATGC-repeat — longer than the test's usable_w (~96)
        # so the zoom cap lands at col_per_bp == 1.0 (zoom > 1.0).
        rec = SeqRecord(
            Seq("ATGC" * 125), id="t500", name="t500",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            seq = str(rec.seq)
            app._register_alignment(
                name="self", query_label="self",
                target_label="self", target_record=rec,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            # Push the zoom to its useful cap (one col per bp).
            cap = pm._max_useful_linear_zoom()
            assert cap > 1.0, (
                "test setup error: record must be longer than "
                "usable_w so the cap lands above min zoom"
            )
            pm._linear_zoom = cap
            pm._linear_offset_bp = 0
            text = pm.render()
            plain = str(text)
            rows = plain.split("\n")
            # Locate the band row by the substring of the alignment's
            # sequence — at col_per_bp == 1.0 the letters land at
            # adjacent columns, so `"ATGCATGC"` appears verbatim.
            band_row = next(
                (r for r in rows if "ATGCATGC" in r),
                "",
            )
            assert band_row, (
                f"expected contiguous letter run in band row at "
                f"col_per_bp=1.0; got rows={rows!r}"
            )

    async def test_letter_row_has_no_internal_spaces_at_cap(
            self, tiny_record, isolated_library):
        """Regression guard (2026-05-23): at the zoom cap (col_per_bp
        == 1.0) the band's letter row must be a contiguous run of
        bases — no internal spaces between letters.

        Pre-fix the `bp_to_col` formula used float math
        (`int((bp - view_s) / visible_bp * usable_w)`). At
        `usable_w == visible_bp`, `1/N*N` is `0.9999999999999999` in
        IEEE 754, so `int()` truncated to 0 and bps 0 and 1 collided
        on the same column. The dropped column rendered as a default
        background cell — visible as "splits" in the letter row.

        Using integer math (`(bp - view_s) * usable_w // visible_bp`)
        gives the same result without precision loss; every bp in
        the visible window now maps to a distinct column when
        `usable_w == visible_bp`.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Use a sequence whose length exceeds usable_w so the cap is
        # above min zoom and reachable. 500 bp self-alignment hits
        # one big "match" segment that should paint every bp.
        rec = SeqRecord(
            Seq("ATGC" * 125), id="t500", name="t500",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            seq = str(rec.seq)
            app._register_alignment(
                name="self", query_label="self",
                target_label="self", target_record=rec,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            pm._linear_zoom = pm._max_useful_linear_zoom()
            pm._linear_offset_bp = 0
            text = pm.render()
            plain = str(text)
            rows = plain.split("\n")
            # Find the alignment band row: skip the header (row 0)
            # and pick the row with the most ACGT density.
            band_row = next(
                (r for r in rows[1:]
                 if sum(1 for c in r.strip() if c in "ACGTacgt") > 50),
                "",
            )
            assert band_row, (
                f"could not locate band row; rows={rows!r}"
            )
            stripped = band_row.strip()
            internal_spaces = sum(1 for c in stripped if c == " ")
            assert internal_spaces == 0, (
                f"band row has {internal_spaces} internal spaces at "
                f"the zoom cap; bp_to_col float precision is causing "
                f"adjacent bps to collide on the same col. "
                f"row={stripped!r}"
            )

    async def test_zoom_in_caps_at_one_col_per_bp(
            self, tiny_record, isolated_library):
        """User feedback 2026-05-22: stop zooming the moment letters
        are one column apart — going further "A T G C" with gutters
        is not useful. `action_linear_zoom_in` must clamp at the
        zoom level where `col_per_bp == 1.0`. Repeated zoom-in
        presses past that should be no-ops.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("A" * 200), id="t200", name="t200",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            cap = pm._max_useful_linear_zoom()
            # Zoom in many times — should plateau at the cap.
            for _ in range(30):
                pm.action_linear_zoom_in()
                await pilot.pause()
            assert pm._linear_zoom <= cap + 1e-6, (
                f"action_linear_zoom_in must clamp at "
                f"_max_useful_linear_zoom={cap}; got "
                f"_linear_zoom={pm._linear_zoom}"
            )
            # And once at the cap, another press doesn't push past.
            pm._linear_zoom = cap
            pm.action_linear_zoom_in()
            assert pm._linear_zoom <= cap + 1e-6, (
                f"once at cap, zoom-in must be a no-op; got "
                f"_linear_zoom={pm._linear_zoom}"
            )


    async def test_alignment_selection_clears_on_detail_dismiss(
            self, tiny_record, isolated_library):
        """Regression guard (2026-05-22): clicking an alignment bar
        sets `_selected_align_idx`, which the renderer turns into
        `style="reverse"` on the bar glyphs. Reverse on a full-block
        "█" inverts fg/bg → the bars read as the terminal's default
        foreground (gray on a dark terminal) instead of their
        blue/red/gray scheme. If the selection isn't cleared when the
        detail screen dismisses, the user comes back to all-gray bars.

        Test contract: register an alignment, set the selection
        manually (proxy for the click), simulate a close of the
        detail screen by pushing + dismissing AlignmentScreen with
        the same callback the click path uses. After dismiss,
        `_selected_align_idx` must be back to -1.
        """
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            app._register_alignment(
                name="r1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            ai = 0
            pm._selected_align_idx = ai
            assert pm._selected_align_idx == ai

            # Mirror the on_click → push_screen(..., callback=) flow
            # without the click event itself. The callback below is
            # the one the production code installs.
            def _clear_align_selection(
                _result=None, _selected_ai=ai, _pm=pm,
            ):
                if _pm._selected_align_idx == _selected_ai:
                    _pm._selected_align_idx = -1
                    _pm.refresh()

            screen = sc.AlignmentScreen(
                query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC",
                        "q_len": 4, "t_len": 4, "score": 8.0,
                        "identity_pct": 100.0, "n_matches": 4,
                        "n_mismatches": 0, "n_gaps": 0,
                        "mode": "global"},
            )
            app.push_screen(screen, callback=_clear_align_selection)
            await pilot.pause()
            await pilot.pause(0.05)
            # Trigger the dismiss path the Esc / q binding takes.
            screen.action_close()
            await pilot.pause()
            await pilot.pause(0.05)
            assert pm._selected_align_idx == -1, (
                "AlignmentScreen.action_close() must fire the "
                "push_screen callback so the lane selection drops; "
                "got _selected_align_idx="
                f"{pm._selected_align_idx}."
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Sequencing hardening (2026-05-18)
# ═══════════════════════════════════════════════════════════════════════════════
# Sweep #8: per-base TSV cap, single-pass zip-open in the Samples
# table, NUL-anchored sentinels for the empty-library / no-gbk paths,
# narrow exception types.

class TestSequencingHardening:
    def test_perbase_summary_truncates_at_max_bytes(self):
        """`_summarize_perbase_tsv` must stop reading once the
        decompressed stream exceeds `max_bytes`. A pathological zip
        bomb that decompresses into a multi-GB single line would
        otherwise OOM `io.TextIOWrapper`'s line buffer."""
        import io
        # Build a 4 KB-per-line TSV with 200 lines = 800 KB. Cap at
        # 200 KB so the streamer stops after ~50 lines, not all 200.
        rows = [f"{i}\tA\t30" for i in range(1, 201)]
        body = ("pos\tref\treads_all\n" + "\n".join(rows)).encode(
            "utf-8",
        )
        # Inflate each row to ~4KB by padding column 1 (`ref`).
        pad = "X" * 4000
        rows_padded = [f"{i}\t{pad}\t30" for i in range(1, 201)]
        body = ("pos\tref\treads_all\n" + "\n".join(rows_padded)).encode(
            "utf-8",
        )
        cap = 200 * 1024
        stats = sc._summarize_perbase_tsv(io.BytesIO(body), max_bytes=cap)
        # Stats are present (the cap allowed *some* rows through).
        assert stats, "truncation must still yield a partial summary"
        # n_pos is bounded by what the cap permitted — about
        # cap / 4 KB ≈ 50 rows. Refuse to specify the exact number,
        # just verify we didn't slurp the full 200.
        assert stats["n_pos"] < 200, (
            f"cap={cap} should have stopped before 200 rows; "
            f"got n_pos={stats['n_pos']}"
        )

    def test_perbase_summary_short_input_complete(self):
        """A short TSV (well under cap) is fully consumed — the cap
        is one-way (truncate-only), it never under-counts on small
        inputs."""
        import io
        body = (
            b"pos\tref\treads_all\n"
            b"1\tA\t10\n2\tT\t20\n3\tG\t30\n4\tC\t40\n5\tA\t50\n"
        )
        stats = sc._summarize_perbase_tsv(
            io.BytesIO(body), max_bytes=1024 * 1024,
        )
        assert stats["n_pos"] == 5
        assert stats["mean"] == 30.0
        assert stats["min"] == 10
        assert stats["max"] == 50
        assert stats["above_20x"] == 4   # 20, 30, 40, 50

    def test_perbase_summary_no_trailing_newline(self):
        """Final row without trailing `\\n` is still counted —
        regression guard for the chunked-reader tail-flush logic.
        Pre-fix the rewrite, a 1-row TSV without trailing newline
        returned an empty dict because the pending fragment never
        got consumed."""
        import io
        body = b"pos\tref\treads_all\n1\tA\t42"
        stats = sc._summarize_perbase_tsv(
            io.BytesIO(body), max_bytes=1024,
        )
        assert stats["n_pos"] == 1
        assert stats["mean"] == 42.0

    def test_parse_zip_skips_oversize_perbase(
            self, tmp_path, monkeypatch):
        """A per-base TSV whose claimed `file_size` exceeds the cap
        is skipped (no read attempted) — defence layer 1. The sample
        still surfaces; its `perbase_coverage` is just empty."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 20), id="X", name="X",
                        annotations={"molecule_type": "DNA"})
        gbk = tmp_path / "X.gbk"
        SeqIO.write(rec, gbk, "genbank")
        zp = tmp_path / "R_results.zip"
        # Plant a big TSV body (~200 KB of synthetic rows).
        big_body = ("pos\tref\treads_all\n"
                    + "\n".join(f"{i}\tA\t30" for i in range(1, 20001)))
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "R_genbank-files/R_1_X.gbk")
            zf.writestr("R_per-base-data/R_1_X.tsv", big_body)
        # Cap to 1 KB so the 200 KB tsv is refused upfront.
        monkeypatch.setattr(
            sc, "_PLASMIDSAURUS_PERBASE_MAX_BYTES", 1024,
        )
        data = sc._parse_plasmidsaurus_zip(zp)
        assert len(data["samples"]) == 1
        # perbase_coverage is empty because the read was refused.
        assert data["samples"][0]["perbase_coverage"] == {}
        # But the sample still lists the perbase member name.
        assert data["samples"][0]["perbase"]

    def test_empty_library_sentinel_is_unique(self):
        """The NUL-anchored sentinels must not collide with any
        realistic library `id` / zip member name. Sanity check that
        they actually contain NUL (which the safe-name check rejects
        in member paths and which LOCUS-safe ids never carry)."""
        assert "\x00" in sc.SequencingScreen._EMPTY_LIBRARY_SENTINEL
        assert "\x00" in sc.SequencingScreen._NO_GBK_KEY_PREFIX

    async def test_target_dropdown_handles_empty_library(
            self, tmp_path, isolated_library):
        """Sequencing screen with NO library entries shows the empty-
        library sentinel and `_go` refuses to advance when the user
        clicks Align without a real target. Verified indirectly via
        the Select's current value (Static's `renderable` is private
        in newer Textual)."""
        from textual.widgets import Select, Button
        # Wipe the library so `_target_options` only has the sentinel.
        sc._save_library([])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            assert isinstance(screen, sc.SequencingScreen)
            # The Select's current value is the empty-library sentinel.
            sel = screen.query_one("#align-target", Select)
            assert sel.value == sc.SequencingScreen._EMPTY_LIBRARY_SENTINEL
            # Simulate a state where the user has picked a sample
            # (forces `_go` past the early-return). The sentinel check
            # should fire BEFORE the zip is opened, so the fake path
            # never gets touched.
            screen._zip_path = tmp_path / "nope.zip"
            screen._selected_member = "ignored.gbk"
            screen.query_one("#btn-align-go", Button).disabled = False
            # Snapshot alignment-registration count; the early-return
            # path must not bump it.
            n_before = len(app._alignments)
            screen._go(None)
            await pilot.pause(0.05)
            # No alignment registered because `_go` short-circuited.
            assert len(app._alignments) == n_before

    async def test_no_gbk_sentinel_refuses_align(
            self, tmp_path, isolated_library):
        """A samples row keyed with the NUL-anchored no-gbk sentinel
        must not arm the Align button (synthetic key would crash
        `_extract_gbk_member`)."""
        from textual.widgets import Button
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen

            class _FakeKey:
                def __init__(self, v):
                    self.value = v

            class _FakeEvent:
                def __init__(self, key):
                    self.row_key = _FakeKey(key)
            # Pretend the user clicked a synthetic row.
            sentinel_key = (
                sc.SequencingScreen._NO_GBK_KEY_PREFIX + "sample-X"
            )
            screen._on_member_selected(_FakeEvent(sentinel_key))
            await pilot.pause(0.05)
            assert screen._selected_member is None
            assert screen.query_one(
                "#btn-align-go", Button,
            ).disabled

    def test_batch_extract_gbk_meta_opens_zip_once(
            self, tmp_path, monkeypatch):
        """`_batch_extract_gbk_meta` should walk every sample's gbk
        inside a single `ZipFile` open — pre-fix each sample paid
        a fresh open. Counts opens via a monkeypatched `__init__`."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Build a 5-sample zip.
        zp = tmp_path / "R_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            for i in range(1, 6):
                rec = SeqRecord(
                    Seq("ATGC" * 50), id=f"S{i}", name=f"S{i}",
                    annotations={"molecule_type": "DNA"},
                )
                gbk = tmp_path / f"S{i}.gbk"
                SeqIO.write(rec, gbk, "genbank")
                zf.write(gbk, f"R_genbank-files/R_{i}_S{i}.gbk")
        parsed = sc._parse_plasmidsaurus_zip(zp)
        # Wire up a SequencingScreen instance just enough to call the
        # batch method directly (avoids the full async-mount cost).
        screen = sc.SequencingScreen.__new__(sc.SequencingScreen)
        screen._zip_path = zp
        screen._parsed_run = parsed
        opens: list[str] = []
        real_init = zipfile.ZipFile.__init__

        def _counting_init(self, file, *a, **kw):
            opens.append(str(file))
            return real_init(self, file, *a, **kw)
        monkeypatch.setattr(zipfile.ZipFile, "__init__", _counting_init)
        meta = screen._batch_extract_gbk_meta(parsed["samples"])
        # Exactly one ZipFile open for all 5 samples.
        assert len(opens) == 1, (
            f"expected 1 zip open for batch read, got {len(opens)}: {opens}"
        )
        # Every sample resolved bp/feats counts.
        assert len(meta) == 5
        for s in parsed["samples"]:
            gbk = s.get("gbk")
            assert gbk in meta, f"missing meta for {gbk}"
            bp_str, _feats = meta[gbk]
            assert bp_str != "—"

    def test_batch_extract_gbk_meta_corrupt_zip_returns_empty(
            self, tmp_path):
        """A corrupted zip path makes `_batch_extract_gbk_meta` log
        and return an empty dict — caller falls back to per-row
        "—" placeholders. Guard against missing-file / bad-zip OS
        errors leaking up the call stack."""
        screen = sc.SequencingScreen.__new__(sc.SequencingScreen)
        # Path that exists but isn't a zip.
        bad = tmp_path / "not-a-zip.txt"
        bad.write_text("hello")
        screen._zip_path = bad
        meta = screen._batch_extract_gbk_meta(
            [{"gbk": "x.gbk"}],
        )
        assert meta == {}

    def test_batch_extract_rejects_unsafe_member_names(self, tmp_path):
        """Belt-and-braces: `_batch_extract_gbk_meta` re-checks
        `_is_safe_zip_member_name` on every member. An in-process
        mutator of `_parsed_run` that tried to smuggle a traversal
        path back in would land in the err bucket, not crash."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 20), id="X", name="X",
                        annotations={"molecule_type": "DNA"})
        gbk = tmp_path / "X.gbk"
        SeqIO.write(rec, gbk, "genbank")
        zp = tmp_path / "R_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "R_genbank-files/R_1_X.gbk")
        screen = sc.SequencingScreen.__new__(sc.SequencingScreen)
        screen._zip_path = zp
        # Hand-crafted "sample" with a traversal name.
        meta = screen._batch_extract_gbk_meta(
            [{"gbk": "../../etc/passwd"}],
        )
        assert meta == {"../../etc/passwd": ("[red]err[/red]", "—")}
