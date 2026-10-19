# Maintenance script inventory (2026-10-20)

This is a path and reference audit of `origin/main` at `e842ae41`, followed by
the first layout cleanup. It is not a claim that every historical script is
still executed. `scripts/README.md` describes the retention policy and main
entry points.

| Measure | Before | After this cleanup |
|---|---:|---:|
| Tracked files under `scripts/` | 597 | 596 |
| Direct children of `scripts/` | 516 | 515 |
| Files under `scripts/tests/` | 81 | 81 |
| Files with a literal `scripts/<path>` reference outside `scripts/` | 588 | 588 |
| Files with a literal path reference in `.github/workflows/` | 103 | 103 |

The reference counts are sets of existing tracked paths, not counts of calls.
A literal reference may be historical; absence of one may mean a Python import,
module-name invocation, test discovery, or manual operation. CI and benchmark
references use exact paths, so a mass move would break current checks and
published reproduction commands.

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

No other script was deleted in this pass. The direct-child list is frozen in
`tools/script_root_baseline.list`; `tools/check_script_layout.py` rejects new
unlisted direct children and unknown subdirectories. New scripts go into a
purpose-specific subdirectory. Moving a retained script requires updating its
CI and documentation references, imports, and repository-root resolution.
