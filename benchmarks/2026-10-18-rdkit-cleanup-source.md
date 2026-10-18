# RDKit MolStandardize Cleanup comparison

- Date: 2026-10-18
- Comparator: RDKit 2026.09.1
- CheMatic: source-built 1.1.7 candidate on macOS arm64, CPython 3.13.

The local extension was linked with rust-lld and the macOS 26 SDK; this is
source-candidate evidence, not a published-wheel measurement.

The earlier Cleanup comparison called CheMatic's `Mol.standardize()`. That
profile also neutralizes charges and selects a canonical tautomer. RDKit's
[`Cleanup`](https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MolStandardize.html)
uses hydrogen removal, sanitization, metal disconnection, normalization,
reionization, and stereochemistry assignment. CheMatic now has a separate
`Mol.cleanup()` path for the implemented subset; `standardize()` keeps its
existing behavior.

## Frozen sample

The script selects 256 evenly spaced entries from the 4,072-row unusual-SMILES
corpus and keeps the 254 rows whose parsed RDKit-compatible inputs are
identical. All counts below use that same filter. The comparator version,
corpus hash, source hashes, native-extension hash, excluded rows, and every
remaining mismatch are recorded in the [raw result](2026-10-18-rdkit-cleanup-source.json).

| Comparison | Exact outputs |
|---|---:|
| Previous `standardize()` vs RDKit Cleanup | 185/254 (72.8%) |
| `normalize_groups()` vs RDKit Cleanup | 234/254 (92.1%) |
| New `cleanup()` vs RDKit Cleanup | 240/254 (94.5%) |
| `normalize_groups()` vs RDKit Normalize | 254/254 (100%) |

The new Cleanup path removes removable explicit hydrogens, disconnects the
supported metal-ligand bonds while balancing formal charges, and normalizes
functional groups. It retains aromatic metal-carbon bonds and ordinary Sn-C
bonds. A Li-bound carbonyl donor is treated as dative, avoiding an invalid
oxygen valence after disconnection. If a metal-disconnected graph cannot be
written as RDKit-compatible SMILES, Cleanup retains its original metal bonds
and normalizes the other groups.

Fourteen rows remain different: eleven metal-containing structures and three
reionization cases. The metal group includes non-tetrahedral stereo retained by
RDKit on isolated metal atoms and more complex coordination/charge-placement
cases. RDKit's acid/base-priority reionization is not implemented in this path.
This result is a bounded exact-output comparison, not evidence of complete
RDKit Cleanup parity or replacement accuracy on arbitrary molecules.

As a broader check on the same source corpus, all 4,072 entries were also
evaluated. Ten had different parsed input semantics; RDKit Cleanup rejected
three more. On the remaining 4,059 common-input rows, the new path agreed on
3,826/4,059 (94.3%), versus 3,019/4,059 for the old `standardize()`
comparison. The [full-corpus result](2026-10-18-rdkit-cleanup-full-corpus-source.json)
records 233 residual rows and seven CheMatic Cleanup output errors. This is a
larger check within the same generated corpus, not an independent holdout.

## Reproduce

Build the current CheMatic Python extension, install RDKit 2026.09.1, then run:

```bash
python scripts/rdkit_cleanup_parity.py \
  --output benchmarks/2026-10-18-rdkit-cleanup-source.json
python scripts/rdkit_cleanup_parity.py --sample-size 4072 \
  --output benchmarks/2026-10-18-rdkit-cleanup-full-corpus-source.json
```
