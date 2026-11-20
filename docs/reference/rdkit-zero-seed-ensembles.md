# Zero-seed conformer ensembles

RDKit [#9738](https://github.com/rdkit/rdkit/issues/9738) reports that
`EmbedMultipleConfs(randomSeed=0)` repeats the same conformer. Chematic's Rust
`rdkit_embed_multiple_confs` and Python `Mol.rdkit_embed_multiple_confs` fix
this by using attempt seeds `0, 2, 3, ...`. Seed one is skipped because zero
and one alias in the underlying random generator. The first result preserves
the existing single-conformer seed-zero behavior. Positive base seeds retain
the existing `base * (index + 1)` schedule.

```python
import chematic

mol = chematic.from_smiles("CCCCCCO").add_hydrogens()
conformers = mol.rdkit_embed_multiple_confs(5, random_seed=0)
assert conformers == mol.rdkit_embed_multiple_confs(5, random_seed=0)
assert conformers[0] == mol.rdkit_embed(random_seed=0)
```

Candidate ordering remains independent of native worker count. Corrected
gradient and MMFF ensemble APIs share this schedule; the quality profile also
advances across additional seed-zero batches. Impossible signed seed ranges
are rejected before allocating the candidate ensemble. Failed embeddings are
still omitted. Different seeds do not guarantee different conformational
minima, and no RMS pruning is added. Seed-zero ensembles intentionally differ
from RDKit 2026.09.1 after the first candidate. WASM currently exposes single
RDKit-style embedding, which is unchanged.

The fixed diagnostic in `validation/rdkit-9738-2026-10-24/` records four flexible
SMILES, seeds 0/1/42 and five unoptimized conformers each. With seed zero, the
prior local Chematic runtime and RDKit 2026.09.1 show no change in heavy-atom
internal distances across candidates; the patched runtime varies by more than
0.1 angstrom for all four inputs. All 12 runs preserve repetition and prefix
stability. All 12 first conformers and all eight positive-seed ensembles match
the prior runtime exactly on the measured host. This small regression packet
does not establish ensemble coverage, energy quality or a published release.
