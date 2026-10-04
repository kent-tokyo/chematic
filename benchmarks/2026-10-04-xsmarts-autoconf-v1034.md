# xsmarts-autoconf on published v1.0.34 (#734, #754)

The contributor of [#734](https://github.com/kent-tokyo/chematic/issues/734)
published [xsmarts-autoconf](https://github.com/swamidasslab/xsmarts-autoconf)
in [#754](https://github.com/kent-tokyo/chematic/issues/754). It probes 83
SMARTS/SMIRKS behaviour flags with small executable cases and runs a
differential sweep. This record pins its result for the **published** PyPI
`chematic==1.0.34` wheel against RDKit 2026.03.6.

## Setup

| Item | Value |
|---|---|
| xsmarts-autoconf | commit `6ce94c945652fa4fc47719387b41be6bf0419402`, rules digest `b2367621cbe9c1d1` |
| chematic | PyPI `chematic==1.0.34` (published wheel, Linux x86-64, CPython 3.11) |
| RDKit | 2026.03.6 (Python wheel) |
| Commands | `xsmarts-autoconf update chematic`; `sweep rdkit chematic` (8 seeds); `tweak rdkit chematic -n 3000 --seed 1` |

The adapter uses the native APIs: `smarts_find` and `run_smirks`.

## Flags

65 of 83 flags give RDKit's value (v1.0.30, from the tool's own config:
42 of 83). The 23 flags aligned since v1.0.30 are the #734 work: SMARTS
reactant templates, element changes, sanitize policy, H counts after charge
changes, explicit-H reactants, inherited charge, mapped `[nH]`, recursive
SMARTS, `z`/`Z`/ranges/dative bonds, carried stereo and order-independent
site choice.

The 18 remaining flags, classified:

| Flag | RDKit | v1.0.34 | Class | Queue |
|---|---|---|---|---|
| `match.unspecified_isotope` | `[12C]` matches only labelled C | isotope ignored | wrong-confident | 1 |
| `match.aromatic_valence` | benzene c is `v4` | aromatic bond counts 1 | wrong-confident | 1 |
| `match.directional_bond` | `C/C` matches single/aromatic | matches nothing | wrong-confident | 1 |
| `match.bond_stereo` | `/` `\` not enforced | enforced, result depends on input spelling | wrong-confident | 1 |
| `smirks.reactant_aromaticity` | Kekulé benzene perceived aromatic | taken as written | wrong-confident | 1 |
| `syntax.insaturation_i` | parse error | parses, matches nothing | wrong-confident | 1 |
| `syntax.hybridization` | `^0`–`^5` | also accepts `^6`, `^7` | lenient parse | 1 |
| `match.ring_count_basis` | symmetrized rings | SSSR | documented profile | 4 |
| `smirks.valence_invalid_raw` | returns unsanitized product | drops at apply | policy (sanitized views agree) | 3 |
| `smirks.hcount_query_explicit_h` | `C[CH6]` | `CC` | policy (RDKit raw is invalid) | 3 |
| `smirks.product_h_pin` | `C[CH5]` with explicit H | `CC` | policy (RDKit raw is invalid) | 3 |
| `smirks.product_h_atom_added` | H atom dropped, `CC[O-]` | no product | policy | 3 |
| `smirks.star_on_explicit_h` | also edits H atoms | heavy atoms only | policy | 3 |
| `smirks.undefined_product_bond` | first alternative of `=,:` | typed error | typed refusal | 3 |
| `smirks.undefined_order_both_sides` | first alternative | typed error | typed refusal | 3 |
| `products.grouped_product_component` | `(A.B)` is one object | parse error | typed refusal | 3 |
| `smirks.fragmented_reactant` | one template with `.` matches a pair | typed error | typed refusal | 3 |
| `smirks.colon_product_bond` | `:` sets aromatic flags | aromatic bond, atoms aliphatic | policy | 3 |

## Sweep and grammar tweaks

The differential sweep (8 seeds) stops at five findings not explained by a
known flag:

| Case | RDKit | v1.0.34 | Class |
|---|---|---|---|
| `[C:1][C:2]>>[C:1].[C:2]` on explicit-H spiro `C1CCC2(C1)CCCCC2` | refilled H | radical `[C]` on the former spiro atom | wrong-confident |
| `[c:1]>>[c+:1]` on benzene; `[c;H1:1]>>[c+:1]` on naphthalene | product kept by sanitize | no product | missing product |
| `[C;$(CO):1]>>[C:1]C(=O)C` on `C[C@H](N)O` | stereo kept | stereo dropped (native) | policy (RDKit profile keeps) |
| `[H:1]>>[H+:1]` on `[H+]` | `[H+]` | no product | edge |

The 3,000-example grammar-tweak run differs in 16 kinds; besides the known
flags it shows `@SP1` constraining a match RDKit ignores, `k0` (RDKit matches
every atom), and the quadruple bond `$` being rejected.

These are measurements of one published artifact on one host; they do not
claim general SMARTS or SMIRKS parity.
