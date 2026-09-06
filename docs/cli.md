# Command-line interface

`chematic-cli` publishes the `chematic` binary for small, composable molecule
conversion workflows. It reads stdin and writes stdout by default, so it can be
used in shell pipelines without a temporary file.

```sh
printf 'CCO\n' | chematic convert \
  --input-format smiles --output-format mol2

chematic convert --input-format mol2 --output-format smiles \
  --input ligand.mol2 --output ligand.smi
```

`convert` enforces a 64 MiB input limit by default. Use
`--max-input-bytes N` to select a stricter or larger explicit bound; inputs
that exceed it are rejected before parsing.

```sh
chematic descriptors 'CC(=O)Oc1ccccc1C(=O)O'
```

The initial CLI contract covers topology-bearing SMILES, MOL/SDF (V2000), MOL
V3000, MOL2, CML, ChemicalJSON, MolJSON, and CDXML. Format aliases such as
`smi`, `.mol2`, and `sdf` are accepted. Coordinates and format-specific
metadata are intentionally outside this first command; use the format-specific
Rust, Python, or WASM APIs when they must be preserved.

Errors are written to stderr and return exit status 2. Unsupported formats and
parse failures are explicit errors; the command never silently forwards the
input unchanged.

`descriptors` emits compact JSON containing canonical SMILES, formula, atom
counts, molecular weight, exact mass, Crippen logP, TPSA, HBD, HBA, and
rotatable-bond count.

`fingerprint <SMILES>` emits JSON with the selected algorithm (`ecfp4`,
`ecfp6`, or `maccs`), width, set-bit indices, popcount, and canonical SMILES.
Set-bit indices
are used instead of an opaque binary string so the output is easy to inspect
and remains language-neutral.

`similarity <SMILES_A> <SMILES_B>` emits the selected algorithm and Tanimoto
similarity. It accepts the same `ecfp4`, `ecfp6`, and `maccs` algorithms as
`fingerprint`.

`substructure <SMILES> <SMARTS>` emits canonical SMILES, the query, match
count, and query-order atom-index mappings as JSON. Invalid SMARTS is reported
as a command error.

`standardize <SMILES>` applies the default standardization pipeline and emits
the canonical input and output SMILES plus an audit report. The report includes
the overall status, whether the structure changed, each pipeline stage's
enabled/changed flags and before/after snapshots, and machine-readable
warnings. This makes normalization decisions inspectable in shell pipelines.

```sh
chematic standardize 'C[NH3+]'
```

`reaction <REACTION_SMILES>` normalizes reaction components and emits the
canonical reaction SMILES, component counts, and a mapped reaction-center
summary. Changed atoms and broken/formed bonds are reported by atom index;
without atom maps, the center is empty rather than inferred.

```sh
chematic reaction '[CH3:1][OH:2]>>[CH3:1][O-:2]'
```

`report <SMILES>` emits the complete single-molecule analysis record used by
the workflow API: canonical identity, formula, Murcko scaffold, descriptor
summary, drug-likeness filters and alert names, functional groups, and named
groups. Invalid input is rejected with a non-zero exit status.

```sh
chematic report 'CC(=O)Oc1ccccc1C(=O)O'
```

`parse <SMILES>` is the lightweight validation path. It emits the original
input, canonical SMILES, formula, atom/bond counts, and formal charge without
running the full descriptor or alert workflow.

```sh
chematic parse 'C[NH3+]'
```

`reaction-match <REACTION_SMILES> <REACTION_SMARTS>` checks whether all
reactant and product patterns match. It emits the normalized reaction, the
original query, and a boolean `matched` field. Invalid reaction input or query
syntax is reported as a command error.

```sh
chematic reaction-match 'CCO>>CCO' '[#6]>>[#6]'
```

`reaction-balance <REACTION_SMILES>` checks element and implicit-hydrogen
counts on the reactant and product sides. It emits both count maps and human-
readable differences; agents are excluded from the balance calculation.

```sh
chematic reaction-balance 'CO.CO>>COC.O'
```

`reaction-fingerprint <REACTION_SMILES>` emits the 2048-bit reaction
fingerprint as set-bit indices. The default `--mode xor` emphasizes structural
transformation differences; `--mode or` reports the composition union. The
output also includes reactant/product popcounts and the normalized reaction.

```sh
chematic reaction-fingerprint 'CCO>>CC=O'
```

`reaction-similarity <REACTION_A> <REACTION_B>` compares the default XOR
reaction ECFP4 fingerprints and emits a Tanimoto similarity together with both
normalized reactions. Invalid reaction input is rejected.

```sh
chematic reaction-similarity 'CCO>>CC=O' 'CCO>>CC=O'
```

`batch-report` reads one SMILES per line from stdin (or `--input FILE`) and
returns a JSON manifest with one record per non-empty, non-comment line.
Successful records contain the complete report; invalid records retain their
input index and an error, so one bad molecule does not discard the batch.

```sh
printf 'CCO\nC1CC\nCCN\n' | chematic batch-report
```

`batch-descriptors` is the lighter batch path: it returns only the compact
descriptor record for each valid SMILES, with per-record errors and aggregate
valid/error counts. It accepts the same stdin or `--input FILE` line contract
as `batch-report`.

```sh
printf 'CCO\nC1CC\nCCN\n' | chematic batch-descriptors
```

`batch-fingerprints` is the corresponding batch path for `ecfp4`, `ecfp6`, or
`maccs`. It emits the selected algorithm, set-bit fingerprint records, and
the same aggregate/per-record error manifest as `batch-descriptors`.

```sh
printf 'CCO\nC1CC\nCCN\n' | chematic batch-fingerprints --algorithm ecfp4
```

`batch-standardize` applies the default auditable standardization pipeline to
each line and returns canonical input/output SMILES, stage reports, warnings,
and per-record errors in input order. It accepts the same stdin or `--input
FILE` contract as the other batch commands.

```sh
printf 'C[NH3+]\nC1CC\nCCO\n' | chematic batch-standardize
```

`batch-similarity` reads one tab-separated pair per line in the form
`SMILES_A<TAB>SMILES_B`, using the selected ecfp4/ecfp6/maccs algorithm. It
returns per-record similarity JSON, input order, and retained errors for bad
pairs or molecules.

```sh
printf 'CCO\tCCO\nCCO\tCCN\n' | chematic batch-similarity --algorithm ecfp4
```

`batch-substructure` reads one `SMILES<TAB>SMARTS` pair per line and returns
match counts and query-order atom mappings in an input-order manifest. Invalid
SMARTS or molecules remain as per-record errors.

```sh
printf 'c1ccccc1\tc\nCCO\t[\n' | chematic batch-substructure
```

`batch-reactions` processes one reaction SMILES per line and returns normalized
reaction records with component counts and reaction-center summaries. Invalid
reaction lines remain as per-record errors in the manifest.

Reaction parsing also enforces bounded input and component sizes through the
public `ReactionParseLimits` contract. The default limits are 16 MiB per
reaction, 10,000 components per side, 100,000 atoms per component, and 200,000
bonds per component; callers can use `parse_reaction_with_limits` for a
stricter policy.

```sh
printf 'CCO>>CCO\nCCO\n' | chematic batch-reactions
```

All batch commands enforce the same resource bounds: `--max-input-bytes`
(default 64 MiB), `--max-records` (default 100,000), and
`--max-line-bytes` (default 1 MiB). An exceeded bound is a command error before
the batch is processed; malformed chemistry within the bounds remains a
per-record error. For example:

```sh
chematic batch-descriptors --max-records 10000 --max-line-bytes 65536 \
  --input molecules.smi
```

Every batch command also includes a versioned partial-result envelope while
retaining its existing top-level operation fields. The envelope contains
`schema_version`, `operation`, `status`, `record_count`, and the effective
`limits`. Records remain in input order, and malformed records remain inline
with their index and error. `status: "complete"` means processing completed
within the declared bounds; cancellation and backpressure are separate future
streaming gates and are not implied by this envelope.
