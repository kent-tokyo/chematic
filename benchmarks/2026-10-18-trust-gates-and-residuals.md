# Trust gates and residual classification — 2026-10-18

This record closes the bounded investigation requested after v1.1.7. It does
not claim universal RDKit compatibility or toolkit-wide speed superiority.

## Identity and protocol

- source commit: `f5d6cd6fae6a369585a72af0c2659c44b3e18a35`
- comparator: RDKit `2026.09.1`
- stress corpus: 4,072 rows, SHA-256
  `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`
- native extension SHA-256:
  `100e6805b16a89b899775587de273d103c3e8d14e8693cbb064d687d0ef02cf4`
- exact row outputs and errors are retained in
  [`validation/results/rdkit-cleanup-current-full-2026-10-18.json`](../validation/results/rdkit-cleanup-current-full-2026-10-18.json)
  and
  [`validation/results/rdkit-molblock-readback-current-2026-10-18.json`](../validation/results/rdkit-molblock-readback-current-2026-10-18.json).

The Python extension was built from the named source on macOS arm64. The local
linker's generated Mach-O string-table offset required alignment and ad-hoc
re-signing before loading. The extension's code and data were otherwise
unchanged. This is source evidence, not a claim about a registry wheel.

## Normalize, Reionize and Cleanup

The historical 256-row selection now has 256 common inputs (two formerly
excluded inputs are accepted). Normalize is 256/256 exact, Cleanup is 242/256,
and Reionize is 205/256. The 14 Cleanup and 51 Reionize residuals are all
classified in the machine record; none is silently dropped.

On all 4,072 stress rows, RDKit Cleanup succeeds on 4,069 rows:

| Operation | Exact | Residual boundary |
|---|---:|---|
| Normalize | 4,055/4,069 | 7 aromatic N-oxide kekulization, 4 functional-group normalization, 3 graph-assignment gaps |
| Cleanup | 3,836/4,069 | 204 not implemented, 26 declared metal-stereo policy differences, 3 sulfur-oxide representation-only differences |
| Reionize | 3,150/4,069 | 883 intentional independent-ionization policy differences; 36 missing hydrogen, resonance or charge assignments |

Cleanup's 204 unimplemented rows comprise 129 metal-disconnection, 33
acid/base transfer, 26 metal-charge, 7 aromatic N-oxide kekulization, 4
functional-group normalization, 3 explicit-H normalization and 2 metal-bond
order rows. This classification is a support contract, not an assertion that
every RDKit output is the only chemically valid answer.

## CTfile readback

A fresh current-source run replaces the earlier 4,000 exact / 57 different /
15 failed snapshot:

| Class | Rows |
|---|---:|
| exact canonical-isomeric readback | 4,045 |
| stereo-only difference | 10 |
| chematic block not readable by RDKit | 1 |
| RDKit reference write/read failure | 16 |

The 16 oracle failures are not charged to chematic. The one chematic readback
failure and ten stereo-only rows remain open chemistry work; the earlier broad
residual count is no longer current-source evidence.

## Outcome-blind holdouts

Candidate `f5d6cd6f` was committed before the source universe was downloaded.
The NCI `NCISMA99` archive and PubChem CID 4,097–20,000 Isomeric SMILES were
then acquired, hashed, and used to select three files by SHA-256 rank without
evaluating any candidate or comparator output:

- MCS comparison: 64 rows;
- SDF attributes: 64 rows;
- stereo/3D: 96 rows.

The NCI source has 249,081 rows (decompressed SHA-256
`031736b438a9596d89cf11802fa97934c3878ebbf021924c417a57e02c7930cd`).
The PubChem bundle has 15,904 records (SHA-256
`5e416ad053da70488f87bdd7422118963e9b92f52342ee45ad55e7fb8ab175ab`).
[`validation/trust_holdouts/v1/manifest.json`](../validation/trust_holdouts/v1/manifest.json)
records candidate/source acquisition times, source and request hashes, and
every selected-file hash. CI verifies the post-candidate ordering, row counts,
hashes and absence of expected outputs. These rows are frozen evaluation data,
have not yet been scored, and must not be used for tuning.

## Performance and compatibility gates

Repeated QED now uses the molecule's derived scalar cache. Graph edits clear
the value and clones start uncached; the descriptor-cache contract tests cover
both rules. First-use chemistry and returned bits are unchanged.

The macOS arm64 source-candidate run used 20 rotating fresh-process blocks on
1,000 fixed molecules. Candidate output matched both public v1.1.7 and RDKit
on every timed operation. Raw samples, confidence intervals, artifact hashes
and process RSS are in
[`2026-10-18-trust-performance-macos-arm64.json`](2026-10-18-trust-performance-macos-arm64.json).

| Lane | Candidate / public v1.1.7 | Candidate / RDKit | Candidate / COSMolKit |
|---|---:|---:|---:|
| QED first-use | 0.990x | 3.254x | 2.762x |
| QED repeated-use | 1,589.759x | 6,552.005x | 5,540.527x |
| QED parse-inclusive | 0.990x | 3.523x | 3.003x |
| chiral Morgan first-use | 0.994x | 1.368x | 0.887x |
| chiral Morgan repeated-use | 1.048x | 5.360x | 3.331x |
| Labute ASA first-use | 0.947x | 0.690x | 0.336x |
| Labute ASA repeated-use | 0.972x | 2.592x | 1.373x |

These are median paired speed ratios. The record's paired confidence rule
confirms the repeated-QED gain and the listed comparator wins, but not a
first-use QED improvement over v1.1.7. First-use chiral Morgan still trails
COSMolKit, and first-use Labute ASA trails both comparators. Peak candidate
RSS remained about 64 MiB versus about 112 MiB for each comparator process.

The hosted Linux x86-64 workflow uses the same 20-block protocol and uploads a
separate raw file. No two-host claim is made until that workflow completes.

Compatibility CI now:

1. runs `cargo-semver-checks` against v1.1.7 for every changed publishable
   crate;
2. rejects removed Python declarations, newly required arguments and narrowed
   dependency contracts;
3. rejects removed generated WASM TypeScript declarations;
4. keeps the outcome-blind holdout manifest immutable.

## Reproduction

```bash
python3 scripts/check_trust_holdouts.py
python3 scripts/check_binding_api_compat.py \
  --baseline v1.1.7 --wasm-dts demo/pkg/chematic_wasm.d.ts
python3 scripts/rdkit_cleanup_parity.py \
  --sample-size 10000 \
  --output validation/results/rdkit-cleanup-current-full-2026-10-18.json
python3 scripts/rdkit_molblock_readback_parity.py \
  --output validation/results/rdkit-molblock-readback-current-2026-10-18.json
```

Rebuilding the frozen selection additionally requires the exact external NCI
source and a PubChem bundle produced by
`scripts/acquire_pubchem_trust_source.py`; the manifest pins both inputs.

The RDKit commands require the pinned comparator environment described by the
JSON records. Performance reproduction uses
`scripts/bench_python_speed_frontier.py`; its raw 20-block files, not this
summary, are the timing evidence.
