"""
test_intron_render — the intron special-feature render (the ``╱╲`` zigzag).

The intron is a special feature type: instead of the ``▒`` block bar every
other feature uses, it paints an alternating box-drawing zigzag so
exon→intron→exon reads at a glance (introns aren't transcribed). The
translation side (splicing introns out via ``_exons``) is covered by
``test_circular_math``; these tests pin the *render* contract, which was
otherwise untested:

  * alternation is keyed on ABSOLUTE column parity, so a chunk-spanning
    intron stays seamless across a line-wrap boundary (no phase shift);
  * the type match is case-insensitive (CommercialSaaS writes "Intron");
  * a zero-length / off-chunk intron paints nothing (no degenerate art).
"""
from __future__ import annotations

import splicecraft as sc


def _blank(width: int) -> list:
    return [(" ", "") for _ in range(width)]


def _chars(arr, lo: int, hi: int) -> list:
    return [arr[c][0] for c in range(lo, hi)]


class TestIntronZigzag:
    def test_zigzag_parity_keyed_on_absolute_column(self):
        # intron [5, 11) in chunk [0, 20). Parity: even abs col → ╲, odd → ╱.
        arr = _blank(20)
        feat = {"start": 5, "end": 11, "strand": 1, "type": "intron",
                "color": "grey"}
        sc._paint_feature_bar(arr, feat, 0, 20)
        assert _chars(arr, 5, 11) == ["╱", "╲", "╱", "╲", "╱", "╲"]
        # outside the intron span is untouched
        assert arr[4][0] == " " and arr[11][0] == " "

    def test_zigzag_seamless_across_chunk_boundary(self):
        # The SAME intron [8, 14) rendered in two adjacent chunks must keep
        # one continuous zigzag — parity follows the absolute bp, not the
        # chunk-local column, so col 9→10 doesn't repeat a glyph.
        feat = {"start": 8, "end": 14, "strand": 1, "type": "intron",
                "color": "grey"}
        a0 = _blank(10)
        sc._paint_feature_bar(a0, feat, 0, 10)     # abs cols 8, 9
        a1 = _blank(10)
        sc._paint_feature_bar(a1, feat, 10, 20)    # abs cols 10..13
        # abs 8(even)=╲ 9(odd)=╱ | 10(even)=╲ 11(odd)=╱ 12=╲ 13=╱
        assert a0[8][0] == "╲" and a0[9][0] == "╱"
        assert [a1[c][0] for c in range(4)] == ["╲", "╱", "╲", "╱"]

    def test_type_match_is_case_insensitive(self):
        arr = _blank(12)
        feat = {"start": 2, "end": 6, "strand": 1, "type": "Intron",
                "color": "grey"}
        sc._paint_feature_bar(arr, feat, 0, 12)
        assert all(arr[c][0] in ("╱", "╲") for c in range(2, 6))

    def test_zero_length_intron_paints_nothing(self):
        arr = _blank(12)
        feat = {"start": 6, "end": 6, "strand": 1, "type": "intron",
                "color": "grey"}
        sc._paint_feature_bar(arr, feat, 0, 12)
        assert all(ch == " " for ch, _ in arr)

    def test_intron_clipped_to_chunk_window(self):
        # intron [3, 30) but chunk only spans [0, 10): paints cols 3..9 only.
        arr = _blank(10)
        feat = {"start": 3, "end": 30, "strand": 1, "type": "intron",
                "color": "grey"}
        sc._paint_feature_bar(arr, feat, 0, 10)
        assert all(arr[c][0] in ("╱", "╲") for c in range(3, 10))
        assert all(arr[c][0] == " " for c in range(0, 3))

    def test_intron_uses_feature_color(self):
        arr = _blank(10)
        feat = {"start": 2, "end": 6, "strand": 1, "type": "intron",
                "color": "magenta"}
        sc._paint_feature_bar(arr, feat, 0, 10)
        assert all(arr[c][1] == "magenta" for c in range(2, 6))
