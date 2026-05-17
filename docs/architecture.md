# Architecture

How SpliceCraft is organized, why, and how to navigate it.

## The single-file rule

SpliceCraft is a single-file Python app (`splicecraft.py`,
~61,000 lines) on Textual + Biopython, with a stdlib-only sidecar
`splicecraft_cli.py` for agent-API clients. The single-file layout is
**intentional** — no import puzzles, everything is greppable from one
place, one totally-ordered source of truth for behaviour.

`grep -n "^class \|^def " splicecraft.py` gives an authoritative live
map.

Test files are 1:1 named after the subsystem they cover.

## File layout (top-to-bottom)

`splicecraft.py` is laid out by concern, top-down:

1. Imports + persistence helpers + path resolution
2. Logging primitives (`_log_event`, `_action_log`, `_timed`)
3. Enzyme catalog + IUPAC + scanner + 2D feature packer + seq-panel
   renderer
4. GenBank I/O
5. `_Canvas` / `_BrailleCanvas` / `PlasmidMap` / `FeatureSidebar`
6. `LibraryPanel`
7. `SequencePanel`
8. Core modals
9. Grammars + settings
10. Codon registry + Kazusa + mutagenesis
11. Feature-library workbench
12. Parts bin
13. Domesticator + Constructor
14. Mutagenize modal
15. Primer design
16. Small modals
17. `PlasmidApp` (controller, keybindings, undo stashes, autosave,
    `@work` threads)
18. `main()`

## The agent handover document

`CLAUDE.md` at the repo root is the contributor-and-agent handover
document. It contains:

- **41 numbered sacred invariants** — biology correctness, persistence
  contracts, concurrency rules, UI conventions, lock-file
  hardening, the `.dna` writer's expected packet inventory. Touching
  invariant code without updating its regression test will trip the
  test in under two seconds.
- **Known pitfalls** — wrap features, `id()` cache keys, Textual
  reactive auto-invalidation, the `_source_path` survival rule across
  in-place edits, Primer3's linear-only constraint, etc.
- **Persistent user preferences** — the conventions for adding a new
  `settings.json` toggle (4 mechanical steps).
- **Pairwise alignment + Plasmidsaurus ingestion** — the two-stage
  pipeline (size caps + alignment).
- **Architecture pointers** + grep recipes.

**Read it before touching the rendering layer, record pipeline,
primer design, or any persisted-data save path.**

## Why one file

The constraint started as a personal preference and has held up under
scrutiny:

- **Greppability.** Every callsite, every state mutation, every error
  path is reachable with one grep. No "find usages" gymnastics
  across packages.
- **No import puzzles.** New contributors don't need to learn the
  module layout. New subsystems can land without any package-graph
  thinking.
- **Editor responsiveness.** Modern editors handle 60k LoC files
  fine; the trade-off is that the LSP / type-checker re-checks the
  whole file on edit, which is acceptable when the suite runs in
  ~5 minutes anyway.
- **Refactor cost is real.** The cost is internalised, not externalised.
  See [V1_GATE.md](
  https://github.com/Binomica-Labs/SpliceCraft/blob/master/V1_GATE.md)
  soft gate S3 / S6 for the long-term plan.

When the constraint will be reconsidered: when the file passes
~100k LoC, or when a subsystem with no `PlasmidApp` coupling
appears that benefits clearly from extraction. See `CONTRIBUTING.md`
for the three-test rule on extractions.

## Test pyramid

| Suite                                  | Use                                              | Runtime |
|----------------------------------------|--------------------------------------------------|---------|
| `tests/test_dna_sanity.py`             | Inner loop while iterating on biology code       | < 2 s   |
| `tests/test_commercialsaas_io.py`      | When touching `.dna` reader / writer             | ~30 s   |
| `tests/test_agent_api.py`              | When touching `_h_*` endpoints                   | ~45 s   |
| `tests/test_smoke.py`                  | End-to-end + update / restore flows              | ~90 s   |
| `tests/test_perf_regression.py`        | Best-of-N regression gates (perf-baseline.json)  | ~3 s    |
| `tests/test_cli_client.py`             | splicecraft-cli sidecar                          | ~2 s    |
| Full suite (`pytest -n auto -q`)       | Before commit, before release                    | ~5 min  |

All tests run offline against synthetic `SeqRecord`s and monkeypatched
data paths; an autouse fixture in `tests/conftest.py` guarantees no
test can write to real user files.

## Concurrency model

- **`@work(thread=True)` workers** for everything heavier than the
  16 ms frame budget: restriction scan, BLAST, HMMscan, Primer3,
  pairwise align, Gibson sim, GB cycle, classifier digest.
- **Stale-record guard**: workers capture
  `_record_load_counter` at entry; the post-work callback bails out
  if the canvas has moved on.
- **Worker exclusivity**: `@work(exclusive=True, group=...)` for
  user-driven modal flows so a click-spam can't pile up superseded
  work.
- **RLock-protected saves**: `_cache_lock` serialises every `_save_*`
  + cache reassignment so two concurrent saves can't land
  `os.replace A→B` while caches land `B→A`. RLock because save
  chains nest.

## Observability

Every save / load / migration / network / lock / shutdown emits a
structured event via `_log_event(event, **fields)`:

```
app.<area>.<verb>      # user actions
op.<area>.<verb>       # operations
<noun>.<verb>          # state (save.ok, record.loaded, migration.step)
```

Sequence content is **never** logged — `_repr_for_log` truncates and
tags anything DNA-shaped.

Design target: **user pastes log → AI parses → patch shipped same
loop**. The agent-friendly event taxonomy is documented in `CLAUDE.md`
invariant #42.
