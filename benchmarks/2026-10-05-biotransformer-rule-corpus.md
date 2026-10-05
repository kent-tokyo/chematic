# BioTransformer rule corpus and v1.0.35 published reruns

> Superseded for the source lane by
> [the follow-up record](2026-10-05-biotransformer-corpus-followups.md)
> (1,051/1,051 implicit-H and 6,601/6,652 explicit-H pairs, no refusals).

#734 asked for BioTransformer's metabolism rules (the reporter used about
1,200) to be applied by RDKit and chematic to the same molecules. The public
BioTransformer repository keeps its rule tables as JSON, so no jar is needed.
One Linux x86-64 host, RDKit 2026.03.6, CPython 3.11. Source lane: code
revision `e2b188a8` (the evidence commit that follows changes documents and
result files only).

## Inputs

| | |
|---|---|
| Rules | `database_metabolicReactions.json`, `database_ENVMICRO_metabolicReactions.json`, `database_standardizationReactions.json` from [Wishartlab-openscience/Biotransformer](https://github.com/Wishartlab-openscience/Biotransformer) `7149f7ec`; 983 SMIRKS, 974 run (3 RDKit cannot parse, 6 have two reactant templates) |
| Reactants | 400 drug-like molecules (`bt_reactants_400.smi`, sha256 `ac68df7b…`), each as parsed and with explicit hydrogens (`AddHs`) |
| Comparison | product sets after RDKit's own `SanitizeMol` (an RDKit product that fails it is not a product) against chematic `run_smirks_checked(..., rdkit_compat=True)` |

The rule tables are LGPL-3.0 with a non-commercial clause (ENVMICRO: CC
BY-NC-SA), so they are not vendored; the summary pins them by SHA-256.
[`scripts/biotransformer_rule_corpus.py`](../scripts/biotransformer_rule_corpus.py)
runs a local copy. Both engines read each reactant as the same canonical
SMILES (RDKit keeps the H count of an input bracket atom such as `[C]` fixed,
which its canonical SMILES drops; 26 input rows have one).

## Results

Of the 400 × 974 pairs per mode, RDKit returns a sanitizable product for
**1,051** (implicit H) and **6,652** (explicit H).

| Pairs where RDKit has a product | v1.0.35 | Source |
|---|---:|---:|
| Implicit H: same product sets | 775 | **1,045** |
| Implicit H: chematic refuses | 241 | 2 |
| Implicit H: different sets | 35 | 4 |
| Explicit H: same product sets | 3,807 | **6,453** |
| Explicit H: chematic refuses | 2,478 | 95 |
| Explicit H: different sets or RDKit truncated | 367 | 104 |

Rules refused on every reactant: v1.0.35 241, source 5. Where RDKit returns
raw products that all fail its own sanitize and chematic returns products:
27 → 0 implicit, 211 → 193 explicit. Run time for the corpus: v1.0.35
2,915 s, source 618 s (v1.0.35 enumerates up to 256 template variants before
refusing).
[Summary](../validation/results/biotransformer-rule-corpus-2026-10-05.json),
[source rows](../validation/results/biotransformer-rule-corpus-rows-2026-10-05.jsonl)
(rule IDs, reactants and products; no SMIRKS).

### What changed in the source

| v1.0.35 behaviour | Rules | Source |
|---|---:|---|
| "atomic-number SMARTS expansion exceeds limit" | 136 | the RDKit profile reads `[#6:n]` as RDKit does (reactant SMARTS, product aliphatic); native falls back to that reading above 256 variants |
| product bond expressions (`!@-`, `@-`, `-;!@`, `=,:`) refused | about 80 | first order primitive; ring-led and new `~` bonds are zero-order |
| `/`, `\` in SMARTS-only reactant templates refused | | match a single or aromatic bond, as in RDKit |
| product OR atoms (`[F,Cl,Br,I;-:7]`, `[N;X3H1+0,X4H2+:1]`) refused | | keep the matched element, apply the first H count/charge |
| products of templates writing Kekulé bonds into an aromatic ring dropped | | judged after RDKit's sanitize steps |
| product `[#6:1]` over reactant `[c:1]`, and product `[#6;a:1]`, aromatic | | aliphatic in the RDKit profile, as RDKit builds them (BTMR0826 ring cleavage matches; BTMR0842 gives RDKit's empty set) |
| product kekulization failed on a thiophene S-oxide `[s+]` and on a ring-fusion N beside `c(=O)` | | RDKit's candidate rule for neutral N/P and `[o+]`/`[s+]` (also native) |
| E/Z stored on an aromatic ring bond (`C/N=c1\c…`) dropped by every rule, even `[Cl:1]>>[Cl:1]` | | carried into products (both profiles) |
| explicit-H ring dearomatized by the template left radicals | | ring carbons refilled as RDKit does (BTMR0655–0657, 454 rows) |

### Remaining differences (source)

| Mode | Class | Rows | Rules |
|---|---|---:|---|
| implicit | unmapped product atom naming no element (RDKit builds a dummy `*`) | 1 | BTMR0596; with BTMR0507/0561/0562/0847 the 5 rules refused everywhere |
| implicit | ambiguous stereo bond order (typed unsupported) | 1 | BTMR0822 |
| implicit | E/Z of a double bond the template creates or rebuilds | 4 | BTMR0709, 0710, 0895 |
| both | `R1,R2` on a quinuclidine N: RDKit's symmetrized ring count gives no match | 2 | BTMR1099 |
| explicit | explicit-H reactant over the atom limit (`reactant_too_large`) | 92 | 58 rules |
| explicit | E/Z only: 18 rows where RDKit drops an E/Z whose reference atom is an explicit H; 36 where chematic's product carries a conjugated double bond's direction on an H atom, which RDKit's `RemoveHs` (and chematic's `remove_hydrogens`) then drop (two macrolides, BTMR0657) | 54 | 31 rules |
| explicit | RDKit keeps H atoms where chematic gives the implicit-H product: `[SH2]` | 15 | BTMR0215, two thiourea standardizations |
| explicit | the same, giving RDKit radicals (`[CH]`) on carried atoms | 20 | BTMR0094, 0666, 0699, 0701, 0868, 1070, 1071, 1134 |
| explicit | RDKit raw products all fail its sanitize; chematic gives the implicit-H product | 193 | 36 rules |
| explicit | RDKit stopped at 1,000 raw products (H permutations); its set is a subset | 10 | BTMR1063 |
| explicit | other product-set differences | 6 | BTMR0072, 0747, 0774, 1031 |
| explicit | product valence refusal | 2 | BTMR0254, 0747 |

The explicit-H classes follow the implicit-H equivalence policy of #734
item 5 ([apply model](../docs/smirks-apply-model.md), step 3), plus the atom
limit. No implicit-H pair has chematic products where RDKit's raw products
fail its sanitize.

## Published v1.0.35 reruns

| Artifact | Gate | Result |
|---|---|---|
| PyPI Linux wheel, CPython 3.9 (sha256 `79eca3a9…`) | 83 reaction fixtures, checked provenance | 80 graph/origin/map, 3 jointly invalid ([json](../validation/results/v1.0.35-published-pypi-linux-cp39-reaction-checked-provenance-83.json)) |
| npm `chematic@1.0.35` (sha512 `X1l4L+…`) | same | 80 + 3 ([json](../validation/results/v1.0.35-published-npm-reaction-checked-provenance-83.json)) |
| crates.io 1.0.35 (`tools/published_rust_checked_gate`) | same | 80 + 3 ([json](../validation/results/v1.0.35-published-crate-reaction-checked-provenance-83.json)) |
| PyPI sdist (sha256 `190b89eb…`), built for CPython 3.11 | xsmarts-autoconf lane (`6ce94c9`) | 75/83, the same 8 expected differences, 0 unexpected ([json](../validation/results/v1.0.35-published-sdist-xsmarts-autoconf-lane.json)) |

The source after this work gives the same on both gates (80 + 3; 75/83 with
0 unexpected). macOS and Windows wheels were not run on this host.
