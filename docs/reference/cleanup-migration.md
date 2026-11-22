# Upgrade code that uses Cleanup

Rust `chematic_chem::cleanup(&mol)` retains its v1.2.2 return type
`Molecule` and best-effort behavior. Use `cleanup_checked(&mol)` to detect
ambiguous metal transformations. Python and WASM use the checked path.

## Handle refusals

| Binding | Checked call | On refusal |
|---|---|---|
| Rust | `chematic_chem::cleanup_checked(&mol)` | `Err(CleanupError)`; inspect `error.code()` |
| Python | `mol.cleanup()` | `chematic.ChematicInputError`; inspect `exc.code` |
| WASM / Node | `cleanup_mol(mol)` | Thrown string with the code before the colon |

In Rust, propagate errors with `let cleaned = chematic_chem::cleanup_checked(&mol)?;`
from a function with a compatible error type. For batch work, keep the input
and its error together so refusals can be reviewed or counted.

| Code | Meaning |
|---|---|
| `cleanup_metal_hydrogen` | Metal disconnection leaves a hydrogen change the algorithm cannot account for |
| `cleanup_unstable_metal` | Cleanup would change the detached metal's charge on a second pass |
| `cleanup_metal_stereo` | A metal-bound tetrahedral center declares five ligands |

Cleanup can also change accepted metal-containing outputs. Check formula,
formal charge, stereochemistry, and repeat Cleanup on representative inputs
before regenerating stored results. See the
[released and candidate row audits](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-21-cleanup-v122-regression-audit.md).
