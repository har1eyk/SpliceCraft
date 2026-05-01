"""
test_new_plasmid — NewPlasmidModal + `_annotate_seq_from_feature_library`.

Covers:
  - The pure helper `_annotate_seq_from_feature_library` finds substring
    matches on both strands, respects the `min_overlap` skip, and handles
    circular wrap.
  - Modal compose mounts cleanly at the baseline terminal size.
  - "Create" button validates input (empty / non-IUPAC bases reject)
    and dismisses with a SeqRecord on success.
  - "Annotate from library" populates SeqFeatures on the resulting record.
  - "Annotate via BLAST" stub surfaces the Phase-2 pending message
    without crashing.

The boundary regression for these modals is in `test_modal_boundaries.py`;
this file owns the **functional** contract.
"""
from __future__ import annotations

import splicecraft as sc


TERMINAL_SIZE = (160, 48)


# ═══════════════════════════════════════════════════════════════════════════════
# Pure helper: _annotate_seq_from_feature_library
# ═══════════════════════════════════════════════════════════════════════════════

class TestAnnotateFromFeatureLibrary:

    def _seed_library(self, entries: list[dict]) -> None:
        """Seed the feature library through the same atomic save the
        rest of the app uses, so the next `_load_features` reads it."""
        sc._save_features(entries)

    def test_empty_sequence_returns_empty(self):
        self._seed_library([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATGAAA", "strand": 1, "color": "",
        }])
        assert sc._annotate_seq_from_feature_library("") == []

    def test_no_library_returns_empty(self):
        # Empty library → nothing to match.
        self._seed_library([])
        assert sc._annotate_seq_from_feature_library("ATGAAATCG") == []

    def test_finds_forward_strand_hit(self):
        # 30 bp library entry, embedded in a longer template. Library
        # entries < 12 bp are skipped, so we use 30 bp here.
        self._seed_library([{
            "name":         "lacZ_30",
            "feature_type": "CDS",
            "sequence":     "A" * 30,
            "strand":       1,
            "color":        "",
        }])
        template = "GGGG" + "A" * 30 + "TTTT"
        hits = sc._annotate_seq_from_feature_library(template)
        assert len(hits) == 1
        h = hits[0]
        assert h["start"] == 4
        assert h["end"]   == 34
        assert h["strand"] == 1
        assert h["name"]   == "lacZ_30"

    def test_finds_reverse_strand_hit(self):
        # Library entry "AAAATTTGGGCCCAAAATCGAACGT" (25 bp, non-palindromic).
        # Template embeds its reverse-complement, so a -1 strand hit
        # should land at the embed position.
        lib_seq = "AAAATTTGGGCCCAAAATCGAACGT"
        self._seed_library([{
            "name":         "thing",
            "feature_type": "misc_feature",
            "sequence":     lib_seq,
            "strand":       1,
            "color":        "",
        }])
        rc = sc._rc(lib_seq)
        template = "GGG" + rc + "TTT"
        hits = sc._annotate_seq_from_feature_library(template)
        assert any(h["strand"] == -1 and h["start"] == 3 for h in hits)

    def test_palindromic_entry_only_counted_once(self):
        # A perfect palindrome on its own RC. Should appear only once,
        # otherwise every restriction-site-shaped library entry would
        # double-count.
        # 12 bp palindrome (min_overlap default).
        pal = "AAAAGGCCTTTT"
        assert sc._rc(pal) == pal
        self._seed_library([{
            "name":         "pal_test",
            "feature_type": "misc_feature",
            "sequence":     pal,
            "strand":       1,
            "color":        "",
        }])
        template = "GGG" + pal + "TTT"
        hits = sc._annotate_seq_from_feature_library(template)
        assert len(hits) == 1

    def test_short_entries_skipped(self):
        # 10 bp < default min_overlap=12 → skipped. Otherwise a 6 bp
        # restriction site library entry would explode the hit count
        # on any pasted sequence.
        self._seed_library([{
            "name":         "short_one",
            "feature_type": "misc_feature",
            "sequence":     "AAAAAGGGGG",   # 10 bp
            "strand":       1,
            "color":        "",
        }])
        template = "TTT" + "AAAAAGGGGG" + "TTT"
        assert sc._annotate_seq_from_feature_library(template) == []

    def test_circular_wrap_hit_records_end_past_total(self):
        # A 14 bp library entry that, on a circular template of length
        # 20, straddles the origin: bp 15..29 → bp 15..20 + 0..9.
        # The helper should append seq[:longest-1] internally and emit
        # a single hit with end > n.
        lib_seq = "GGGGAAAATTTTAA"           # 14 bp
        self._seed_library([{
            "name":         "wrap_one",
            "feature_type": "misc_feature",
            "sequence":     lib_seq,
            "strand":       1,
            "color":        "",
        }])
        # Construct a 20 bp circular template with the lib_seq starting
        # at bp 12 and wrapping (12..20 + 0..6).
        # First 12 chars are filler, last 8 chars = lib_seq[:8],
        # head 6 chars = lib_seq[8:14].
        template = lib_seq[8:14] + "XXXXXX" + lib_seq[:8]
        # template is 20 bp; lib_seq lives at [12,) with 6 bp wrapping
        # back to [0, 6). Because 'X' isn't a valid IUPAC base in
        # _sanitize_bases, swap to 'C'.
        template = lib_seq[8:14] + "CCCCCC" + lib_seq[:8]
        assert len(template) == 20
        hits = sc._annotate_seq_from_feature_library(
            template, circular=True,
        )
        # Expect at least one hit at start=12 with end past 20.
        wrap_hits = [h for h in hits if h["start"] == 12 and h["end"] > 20]
        assert wrap_hits, f"expected wrap hit, got: {hits}"


# ═══════════════════════════════════════════════════════════════════════════════
# Modal compose + interaction
# ═══════════════════════════════════════════════════════════════════════════════

class TestNewPlasmidModalFlow:

    async def test_modal_mounts(self, tiny_record, isolated_library):
        """Smoke: pushing the modal doesn't raise and the key widgets
        are present."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.NewPlasmidModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            for wid in ("#newplas-name", "#newplas-seq",
                        "#btn-newplas-create",
                        "#btn-newplas-annot-lib",
                        "#btn-newplas-annot-blast",
                        "#btn-newplas-cancel"):
                modal.query_one(wid)   # raises NoMatches if missing

    async def test_create_with_valid_seq_returns_record(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        result_holder: dict = {}
        def _capture(res):
            result_holder["got"] = res
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.NewPlasmidModal(), callback=_capture)
            await pilot.pause()
            modal = app.screen
            modal.query_one("#newplas-name", sc.Input).value = "myplas"
            modal.query_one("#newplas-seq",  sc.TextArea).text = "ATGCATGCATGC"
            await pilot.pause()
            modal.query_one("#btn-newplas-create", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.05)
        result = result_holder.get("got")
        assert result is not None, "modal didn't dismiss on Create"
        rec = result["record"]
        assert str(rec.seq) == "ATGCATGCATGC"
        assert rec.annotations.get("topology") == "circular"
        # No annotated features should exist on a plain Create call.
        assert all(f.type == "source" for f in rec.features) or \
               len(rec.features) == 0

    async def test_create_rejects_empty_seq(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        result_holder: dict = {}
        def _capture(res):
            result_holder["got"] = res
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.NewPlasmidModal(), callback=_capture)
            await pilot.pause()
            modal = app.screen
            # Don't fill the TextArea — Create should reject and the
            # status label should mention "at least one base".
            modal.query_one("#btn-newplas-create", sc.Button).press()
            await pilot.pause()
            status = modal.query_one("#newplas-status", sc.Static)
            # Status renderable is a Rich Text-ish; render it for the assertion.
            txt = str(status.render())
            assert "base" in txt.lower() or "paste" in txt.lower()
        # Modal still up, no dismissal happened.
        assert "got" not in result_holder

    async def test_create_rejects_non_iupac(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NewPlasmidModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#newplas-seq", sc.TextArea).text = "ATGZZZ"
            await pilot.pause()
            modal.query_one("#btn-newplas-create", sc.Button).press()
            await pilot.pause()
            status = modal.query_one("#newplas-status", sc.Static)
            txt = str(status.render())
            # Sanitizer message includes "non-IUPAC" — close enough.
            assert "iupac" in txt.lower() or "z" in txt.lower()

    async def test_annotate_from_library_creates_features(
            self, tiny_record, isolated_library):
        # Seed library with a 30 bp entry and expect one feature on the
        # output record.
        sc._save_features([{
            "name":         "marker30",
            "feature_type": "CDS",
            "sequence":     "ATG" * 10,
            "strand":       1,
            "color":        "",
        }])
        template = "C" * 6 + ("ATG" * 10) + "C" * 6
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        result_holder: dict = {}
        def _capture(res):
            result_holder["got"] = res
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NewPlasmidModal(), callback=_capture)
            await pilot.pause()
            modal = app.screen
            modal.query_one("#newplas-seq", sc.TextArea).text = template
            await pilot.pause()
            modal.query_one("#btn-newplas-annot-lib", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.05)
        result = result_holder.get("got")
        assert result is not None
        rec = result["record"]
        # At least one non-source feature, and it has type CDS.
        feats = [f for f in rec.features if f.type != "source"]
        assert len(feats) >= 1
        assert feats[0].type == "CDS"

    async def test_blast_button_rejects_empty_seq(
            self, tiny_record, isolated_library):
        # Engine wired: empty paste → "Paste at least one base." reject.
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NewPlasmidModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#btn-newplas-annot-blast", sc.Button).press()
            await pilot.pause()
            status = modal.query_one("#newplas-status", sc.Static)
            txt = str(status.render())
            assert "base" in txt.lower() or "paste" in txt.lower()


class TestBlastModalFlow:
    """BlastModal is now engine-wired (Phase 2 done). Verify its
    scaffolding: program/source selects mount, Run with empty query
    reports an error, Run with a query renders engine output, and
    Build Database reports a real summary."""

    async def test_modal_mounts(self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            for wid in ("#blast-query", "#blast-program", "#blast-source",
                        "#btn-blast-run", "#btn-blast-build",
                        "#btn-blast-close"):
                modal.query_one(wid)

    async def test_run_with_empty_query_complains(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#btn-blast-run", sc.Button).press()
            await pilot.pause()
            status = modal.query_one("#blast-status", sc.Static)
            txt = str(status.render())
            assert "query" in txt.lower() or "paste" in txt.lower()

    async def test_run_executes_engine_against_empty_db(
            self, tiny_record, isolated_library):
        # No collections seeded → engine reports an empty database
        # rather than crashing or returning the old Phase-2 stub.
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#blast-query", sc.TextArea).text = (
                "ATGAAATCGGAATTCATGAAA"
            )
            await pilot.pause()
            modal.query_one("#btn-blast-run", sc.Button).press()
            await pilot.pause()
            results = modal.query_one("#blast-results", sc.Static)
            txt = str(results.render())
            # Engine wired: results are either "no hits" or an empty-DB
            # notice, but never the old "Phase 2 pending" string.
            assert "pending" not in txt.lower()

    async def test_build_database_reports_summary(
            self, tiny_record, isolated_library):
        # Engine builds the (possibly empty) DB and reports its summary
        # rather than the Phase-2 stub.
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#btn-blast-build", sc.Button).press()
            await pilot.pause()
            status = modal.query_one("#blast-status", sc.Static)
            txt = str(status.render())
            assert "pending" not in txt.lower()
            # Summary mentions BLASTN + 'subjects' even for an empty DB.
            assert "blastn" in txt.lower()
            assert "subject" in txt.lower()
