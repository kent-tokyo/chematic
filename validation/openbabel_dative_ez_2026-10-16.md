# V3000 dative carbonyl E/Z boundary, 2026-10-16

## Defect and correction

RDKit 2025.09.3 wrote a 2D V3000 MOL for stress corpus line 3945. Its
`C(=[O][Na])N` group has a C=O double bond and an O-to-Na dative bond (V3000
bond type 9). The pre-fix CheMatic 1.1.5 CLI emitted `C(=O/[Na])/N` in that
group. RDKit read those slashes as specified alkene geometry, so the output
did not match the source identity. Open Babel 3.2.1 preserved it.

The five-atom form `CC(N)=[O]->[Na]` reproduces the defect. The 2D E/Z
perceptor counted the dative bond as a sigma substituent on carbonyl oxygen,
then marked both sides of C=O with directional bonds. The correction excludes
dative bonds from E/Z substituent and carrier selection. It retains the dative
bond itself. The source-built candidate now writes `C(=O->[Na])(C)N` for the
minimal MOL; RDKit reads this as the same canonical isomeric SMILES as the
source. The original line 3945 also matches after conversion.

## Paired result

The comparison script gives CheMatic and Open Babel identical RDKit-written
2D MOL inputs, after RDKit has reopened each file and confirmed its identity.
RDKit then reads each output SMILES. Source failures and format substitutions
are excluded before scoring. An ordinary explicit-H representation difference
is accepted only for molecules with no specified stereo, as described in
`openbabel_stereo_2d_2026-10-16.md`.

| Cohort | Format | Comparable | CheMatic before | CheMatic after | Open Babel |
| --- | --- | ---: | ---: | ---: | ---: |
| Last 76 of 1,216 unique stereo-syntax stress molecules | V2000 | 71 | 71 | 71 | 61 |
| Last 76 of 1,216 unique stereo-syntax stress molecules | V3000 | 75 | 74 | 75 | 59 |

The sole changed row is source line 3945. The other 74 comparable V3000 rows
and all 71 comparable V2000 rows retain their result. There are five V2000
and one V3000 source exclusions in this cohort.

The source-built candidate was then rerun on all 1,216 selected unique
stereo-syntax molecules, not just the final cohort. The complete per-row
report is `results/openbabel-2d-stereo-dative-stress-2026-10-16.json.gz`.

| Full corpus | Comparable | CheMatic after | Open Babel | Source exclusions |
| --- | ---: | ---: | ---: | ---: |
| 2D MOL V2000 → SMILES | 1,163 | 1,163 | 979 | 53 |
| 2D MOL V3000 → SMILES | 1,198 | 1,198 | 941 | 18 |

Neither format had a comparable case where both converters missed after the
fix. These are semantic identity counts under the stated oracle, not speed
measurements.

## Provenance and checks

- Corpus: `validation/cosmolkit_comparison/corpora/stress-v1.smi`, SHA-256
  `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`.
- Original line 3945 V3000 MOL SHA-256:
  `d3183976cec72f290133b913792a8421184ab005cbd7077c537036396409b734`.
- Comparison script SHA-256:
  `2a2bc972eaefca0f7fd72723f18f3bc0694d45da2fa69d850e2306edb01339cc`.
- Complete compressed report SHA-256:
  `b1bafa465e13e732135e2a96abe18c1981bf6554b3cab4c881e5886b822da1d3`.
- Source-built CheMatic 1.1.5 candidate CLI SHA-256:
  `70b5f087201ddb2d6155819d70f373795ef890410ba5e28aeef2379c44e183d9`.
  It was built from `57d318c0` plus this change.
- Open Babel 3.2.1 binary SHA-256:
  `733b2c925d60b17a5c9d8ce9ec841e918cc20c7634301e23f9008a187d1d9bc5`.
- `cargo test --offline -p chematic-perception --lib`: 272 passed, 1 ignored.
- `cargo test --offline -p chematic-mol --lib`: 782 passed.
- `cargo fmt --all -- --check`: passed.

To repeat the full corpus comparison from the repository root with the pinned
tool versions installed:

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline \
  -p chematic-cli --bin chematic -j 1
python3 scripts/compare_openbabel_stereo_formats.py \
  --cli target/debug/chematic \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --stress-first-unique 1216 \
  --output /private/tmp/chematic-openbabel-stress-all-fixed.json
```

This verifies the declared MOL-to-SMILES workflow on one host with RDKit as
the comparison oracle. It does not establish experimental stereochemical
truth, published-package behavior, or toolkit-wide superiority over RDKit,
Open Babel, or COSMolKit.
