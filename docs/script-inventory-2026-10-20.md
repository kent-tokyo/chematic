# Maintenance script inventory (2026-10-20)

This is a path and reference audit of `origin/main` at `e842ae41`, followed by
the first layout cleanup. It is not a claim that every historical script is
still executed. `scripts/README.md` describes the retention policy and main
entry points.

| Measure | Before | After this cleanup |
|---|---:|---:|
| Tracked files under `scripts/` | 597 | 593 |
| Direct children of `scripts/` | 516 | 512 |
| Files under `scripts/tests/` | 81 | 81 |
| Files with a literal `scripts/<path>` reference outside `scripts/` | 588 | 586 |
| Files with a literal path reference in `.github/workflows/` | 103 | 103 |

The reference counts are sets of existing tracked paths, not counts of calls.
Many matches come from frozen benchmark JSON source-hash maps that list nearly
every script. They are provenance, not evidence that a script is still used.
Excluding data files and the new baseline list, 155 original root files had no
literal path reference outside `scripts/`. Some are still imported, invoked by
basename, used manually, or needed to regenerate committed fixtures. Active CI
and some reproduction instructions use exact paths, so a mass move would break
them.

## Files without a literal path reference outside `scripts/`

Before this cleanup there were nine:

- `freeze_trust_holdouts.py` selects outcome-blind inputs for the committed
  `validation/trust_holdouts/v1` packet. It now lives beside that packet at
  `validation/trust_holdouts/freeze.py`; its output path is unchanged.
- `rdkit_standardization_classification.py` is imported by
  `rdkit_cleanup_parity.py` and `rdkit_cleanup_source_gate.py`. It is retained.
- Seven files in `scripts/tests/` are test modules. Five were already invoked
  by module name in CI. This cleanup adds the other two to the Test job, after
  converting `test_changed_public_crates.py` to `unittest`. A path-only scan is
  not grounds for deleting tests.

## Removed diagnostics

- `analyze_logp_mismatches.py` relied on a local `~/Downloads/SMILES.csv` and
  wrote an untracked TSV under `scripts/`. It had no CI, documentation, or open
  issue reference and no committed result to regenerate.
- `aromaticity_mechanism_probe.py` likewise used a local CSV and printed an
  exploratory report. It had no CI, documentation, or open issue reference and
  no committed result. Later aromaticity gates and records remain in place.
- `check_generated_v3000_roundtrip.mjs` duplicated the SGROUP-preservation and
  SGROUP-before-COLLECTION checks already covered by the WASM unit test and the
  generated browser-artifact reproduction in
  `benchmarks/2026-09-10-generated-wasm-v3000-gate-v1.0.11.md`. It had no
  active CI or reproduction-command reference.

Frozen benchmark JSON source-hash maps are unchanged; Git history preserves
these scripts at the original revisions. Other candidates still need a
purpose-by-purpose review before deletion. The direct-child list is frozen in
`tools/script_root_baseline.list`; `tools/check_script_layout.py` rejects new
unlisted direct children and unknown subdirectories. New scripts go into a
purpose-specific subdirectory. Moving a retained script requires updating its
CI and documentation references, imports, and repository-root resolution.
