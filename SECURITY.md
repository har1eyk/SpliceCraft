# Security Policy

SpliceCraft is a desktop / terminal tool that reads scientific file
formats (`.gb` / `.gbk` / `.genbank` / `.dna` / `.fasta` / `.zip`)
and contacts a small number of public scientific endpoints (NCBI E-utils,
PyPI, Kazusa) on the user's behalf. Most of the threat model is
"the file or the response I just loaded is hostile."

## Supported versions

Only the latest `0.x.y` release on PyPI is supported. Older `0.x` lines
are not patched — upgrade with `splicecraft update` (which atomically
snapshots your data first; see invariant #39 in `CLAUDE.md`).

| Version | Supported            |
|---------|----------------------|
| Latest  | :white_check_mark:   |
| Earlier | :x:                  |

This will change at v1.0.0 — see `V1_GATE.md` for the gate, and the
maintenance policy will be revised in the same release.

## Reporting a vulnerability

**Please do not file public issues for security bugs.**

Email `scocioba@gmail.com` with `[splicecraft-security]` in the
subject line. Include:

- SpliceCraft version (`splicecraft --version`)
- Python version + OS
- A minimal reproduction (input file, command, expected vs. observed)
- Whether the bug requires a user action (open file, follow link) or
  fires on launch / on autosave / on a background timer
- An asset hint: code execution, file disclosure outside `DATA_DIR`,
  data corruption, denial-of-service, network egress to a third party,
  etc.

You should expect an acknowledgement within **5 business days**. If
a fix is straightforward, expect a PyPI release on the next batched
release (cadence has historically been a few times per week — see
`CHANGELOG.md`). For coordinated disclosure, propose a date in the
initial email; the default ceiling is **90 days** from
acknowledgement.

Credit will be offered in the changelog entry unless you ask
otherwise.

## Scope

### In scope

The defensive posture described in `CLAUDE.md` is what we maintain.
Bugs that **bypass any of the following** are in scope:

| Surface                | Defence                                                                         | Sacred invariant |
|------------------------|---------------------------------------------------------------------------------|------------------|
| NCBI XML responses     | `_safe_xml_parse` rejects DOCTYPE / ENTITY before `ET.fromstring`               | #8, #19          |
| `.dna` history packets | Same `_safe_xml_parse` path                                                     | #19              |
| `.dna` LZMA history    | Streaming decompress with `max_length` cap (no decompression bomb expansion)    | #21              |
| `.dna` sidecar paths   | `_dna_sidecar_path` strips `..`, dot-only, NUL; appends 8-char SHA-1 prefix     | #22, #41         |
| JSON library files     | `_safe_load_json` size-capped at 1 GB                                           | #23              |
| Agent-API uploads      | `_h_load_file` 50 MB cap, `force=true` override; per-handler caps documented    | #24              |
| PyPI / NCBI / Kazusa   | `resp.read(MAX + 1)` + bail-if-exceeded                                          | #20              |
| Agent-API write paths  | `_check_agent_write_path` walks full ancestor chain via `resolve()` divergence   | #41              |
| Agent-API read paths   | `_check_agent_read_dir` via `lstat` + `S_ISDIR`                                  | #41              |
| Atomic writes          | `_atomic_write_text` / `_atomic_write_bytes` + `_fsync_parent_dir` after rename | #31, #41         |
| Export extension       | `_check_export_extension` whitelists `.gb` / `.gbk` / `.gff3` / `.fasta` etc.   | #41              |
| Path sanitisation      | `_sanitize_path` refuses `~user` (user-enumeration oracle)                       | #41              |
| Concurrent saves       | `_cache_lock` (RLock) wraps every `_save_*` + cache reassignment                 | #41              |

Examples of in-scope vulnerabilities:

- XXE or billion-laughs through any XML ingest path
- Path traversal through `.dna` sidecar naming, GenBank record IDs,
  agent-API endpoints, or export targets
- Symlink escape from `DATA_DIR` (e.g. a symlink in the data dir
  redirecting a `_save_*` write to `~/.ssh/authorized_keys`)
- Decompression bombs through `.dna` or Plasmidsaurus `.zip` ingest
- Lock-file races leading to silent data loss
- Pre-update snapshot bypass (allowing a bad release to wipe data
  without recovery)
- TOCTOU between `_safe_save_json`'s ancestor check and `os.replace`
- Cache-poisoning across the deepcopy-on-read+save contract
  (invariant #17)
- Code execution via the agent-API (`splicecraft --agent-api`) from
  a non-localhost origin

### Out of scope

- **Local attacker with arbitrary disk write.** If they can write to
  `~/.local/share/splicecraft/`, they can already replace the
  library; SpliceCraft does not defend against post-compromise
  persistence in the user's own home directory.
- **The agent-API listening on a non-localhost interface.** It binds
  to `127.0.0.1` by default and trust is single-tenant by design.
  Exposing it on a LAN is not supported.
- **Bugs in upstream dependencies** (Biopython, Textual, Primer3,
  pyhmmer, platformdirs). Report those to their respective projects.
  We're happy to receive a heads-up so we can pin a workaround.
- **Resource exhaustion on intentionally pathological input** larger
  than the documented cap (e.g. a 200 MB GenBank file — we cap at
  64 MB and surface an error).
- **NCBI / PyPI / Kazusa serving compromised responses.** We defend
  against malformed responses but not against an upstream
  compromise. We do verify size caps and parse defensively.

## Defence-in-depth surfaces you may also exercise

`scripts/perf_probe.py` and `scripts/perf_probe_render.py` can be
used to gauge cost ceilings without writing exploit code. The CLI
client (`splicecraft_cli.py`) has its own 50 MB response cap and
1 KB token-file cap — those are intentional and tested.

## What we will not do

- Issue CVEs ourselves. If a vulnerability merits one, please open
  the CVE through MITRE; we will mirror the ID in the changelog.
- Ship a security advisory through GHSA without your acknowledgement.
- Email-bomb your inbox. One acknowledgement, one heads-up on the
  release that contains the fix.

## Hall of fame

Reporters will be listed here as fixes ship. Be the first.
