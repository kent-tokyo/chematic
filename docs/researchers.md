# chematic for researchers

chematic is a pure-Rust cheminformatics runtime with Python, Rust, JavaScript/
WebAssembly, and MCP entry points. The useful research claim is operational:
the same core implementation can be used in a local analysis, a Python
pipeline, and a browser or agent-facing tool while preserving explicit limits
and typed failures.

## Start with one reproducible result

```bash
python -m venv .venv
. .venv/bin/activate
pip install chematic==1.0.12
```

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
print({"mw": mol.mw, "logp": mol.logp, "tpsa": mol.tpsa})
```

Pin the package version and record the input corpus, configuration, and
platform. Do not treat a descriptor agreement or benchmark result as a
universal guarantee: the scope is always the named metric, corpus, runtime,
and failure policy.

## Choose the right evidence

| Research question | Start here | Boundary |
| --- | --- | --- |
| Can I parse and screen a batch safely? | [error and resource limits](error-and-limits.md), [format matrix](format-capabilities.md) | Limits and rejected records are part of the result |
| Can I compare descriptors or fingerprints? | [validation](validation.md), [benchmark index](benchmark.md) | Compare only the pinned metric and corpus |
| Can I reproduce a release? | [release gate](v1.0-local-release-gate.md), [CHANGELOG](changelog.md) | Local evidence is not hosted or independent-review evidence |
| Can I replace an RDKit call? | [RDKit migration guide](rdkit-migration.md) | chematic is not a full RDKit clone; unsupported options fail explicitly |

3D generation, pKa/ADMET screening, IUPAC naming, rich reaction semantics,
and some canonicalization cases remain bounded or experimental. Read the
compatibility and validation documents before using those paths in a paper.

## Citation

The repository contains a machine-readable [`CITATION.cff`](https://github.com/kent-tokyo/chematic/blob/main/CITATION.cff).
Cite the exact release used by the study, including the version and repository
URL. Once the GitHub repository is connected to Zenodo, use the version DOI for
the archived release and keep the Git commit or tag in the methods section.

## How to report a useful issue

Include the smallest reproducible input, chematic version, binding, feature
flags, expected behavior, observed output, and whether the case is a correctness
failure, an unsupported input, or a measurement gap. Remove proprietary or
unlicensed structures before submitting a fixture.
