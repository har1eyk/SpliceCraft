# Rotation Ranker Synthetic-Pair Sweep

- Total pairs: **20**
- Plain (no rot/RC): 1
- With rotation: 18
- With RC: 12

## Results

| Ranker | Correct | Pct |
|---|---|---|
| Current `(n_matches, ungapped_identity_pct)` | 18 | 90.0% |
| Alternative `(score, n_matches)` | 17 | 85.0% |

- Disagreements: 1
- Δ: **-1**

## Verdict

TIE / WITHIN NOISE — current ranker is fine; revisit only if real-world data exposes a specific regression.
