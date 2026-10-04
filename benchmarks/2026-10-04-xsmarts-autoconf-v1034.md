# xsmarts-autoconf on v1.0.34 and the #754 source fixes (#734, #754)

The contributor of [#734](https://github.com/kent-tokyo/chematic/issues/734)
published [xsmarts-autoconf](https://github.com/swamidasslab/xsmarts-autoconf)
in [#754](https://github.com/kent-tokyo/chematic/issues/754). It probes 83
SMARTS/SMIRKS behaviour flags with small executable cases and runs a
differential sweep. This record pins its result for the v1.0.34 release and
for the source branch that follows, against RDKit 2026.03.6, together with
the 83-row reaction gate on the published v1.0.34 binaries.

## Setup

| Item | Value |
|---|---|
| xsmarts-autoconf | commit `6ce94c945652fa4fc47719387b41be6bf0419402`, rules digest `b2367621cbe9c1d1` |
| v1.0.34 column | PyPI sdist `chematic-1.0.34.tar.gz` built locally (maturin 1.15.0, CPython 3.11, Linux x86-64). PyPI publishes no CPython 3.10+ Linux x86-64 wheel for 1.0.34, and xsmarts-autoconf needs Python ≥ 3.10, so this column is the published source, not a published binary |
| Source column | branch `fix/734-754-dialect-parity`, release-profile wheel (CPython 3.11) |
| RDKit | 2026.03.6 (Python wheel) |
| Commands | `xsmarts-autoconf update chematic`; `sweep rdkit chematic --seeds 8`; `tweak rdkit chematic -n 3000 --seed 1` |

The adapter uses the native APIs: `smarts_find` and `run_smirks`.
Flags are compared by value with the tool's own RDKit 2026.03.6 config.

## Flags

| | v1.0.30 (tool config) | v1.0.34 | source |
|---|---|---|---|
| Flags giving RDKit's value | 42 / 83 | 65 / 83 | **72 / 83** |

The 23 flags aligned from v1.0.30 to v1.0.34 are the #734 work: SMARTS
reactant templates, element changes, sanitize policy, H counts after charge
changes, explicit-H reactants, inherited charge, mapped `[nH]`, recursive
SMARTS, `z`/`Z`/ranges/dative bonds, carried stereo and order-independent
site choice.

The 18 flags v1.0.34 differs on, classified, and the source result:

| Flag | RDKit | v1.0.34 | Class | Queue | Source |
|---|---|---|---|---|---|
| `match.unspecified_isotope` | `[12C]` matches only labelled C | isotope ignored | wrong-confident | 1 | aligned |
| `match.aromatic_valence` | benzene c is `v4` | aromatic bond counts 1 | wrong-confident | 1 | aligned |
| `match.directional_bond` | `C/C` matches single/aromatic | matches nothing | wrong-confident | 1 | aligned |
| `match.bond_stereo` | `/` `\` not enforced | enforced, result depends on input spelling | wrong-confident | 1 | aligned |
| `smirks.reactant_aromaticity` | Kekulé benzene perceived aromatic | taken as written | wrong-confident | 1 | aligned |
| `syntax.insaturation_i` | parse error | parses, matches nothing | wrong-confident | 1 | aligned |
| `syntax.hybridization` | `^0`–`^5` | also accepts `^6`, `^7` | lenient parse | 1 | aligned |
| `match.ring_count_basis` | symmetrized rings | SSSR | documented profile | 4 | differs |
| `smirks.valence_invalid_raw` | returns unsanitized product | drops at apply | policy (sanitized views agree) | 3 | differs |
| `smirks.hcount_query_explicit_h` | `C[CH6]` | `CC` | policy (RDKit raw is invalid) | 3 | differs |
| `smirks.product_h_pin` | `C[CH5]` with explicit H | `CC` | policy (RDKit raw is invalid) | 3 | differs |
| `smirks.product_h_atom_added` | H atom dropped, `CC[O-]` | no product | policy | 3 | differs |
| `smirks.star_on_explicit_h` | also edits H atoms | heavy atoms only | policy | 3 | differs |
| `smirks.undefined_product_bond` | first alternative of `=,:` | typed error | typed refusal | 3 | differs |
| `smirks.undefined_order_both_sides` | first alternative | typed error | typed refusal | 3 | differs |
| `products.grouped_product_component` | `(A.B)` is one object | parse error | typed refusal | 3 | differs |
| `smirks.fragmented_reactant` | one template with `.` matches a pair | typed error | typed refusal | 3 | differs |
| `smirks.colon_product_bond` | `:` sets aromatic flags | aromatic bond, atoms aliphatic | policy | 3 | differs |

No wrong-confident flag remains in the source; the 11 that differ are typed
refusals or documented policies, listed with their values in
[`validation/xsmarts_autoconf_expected.json`](../validation/xsmarts_autoconf_expected.json).
[`scripts/xsmarts_autoconf_lane.py`](../scripts/xsmarts_autoconf_lane.py)
reruns the pinned tool in CI (Linux source wheel) and fails on any other
difference, or when an expected one disappears.

## Sweep

The differential sweep (8 seeds) stops at the findings it cannot explain by
a known flag.

| Case | RDKit | v1.0.34 | Source |
|---|---|---|---|
| `[C:1][C:2]>>[C:1].[C:2]` on explicit-H spiro `C1CCC2(C1)CCCCC2` | refilled H | radical `[C]` on the former spiro atom | fixed |
| `[c:1]>>[c+:1]` on benzene; `[c;H1:1]>>[c+:1]` on naphthalene | product kept | no product | fixed (products) |
| `[c:1]>>[c-:1]` on benzene | product kept | no product | fixed (found after the row above) |
| `[C:1]=[O:2]>>[C:1]-[O:2]` on explicit-H acetic acid | `CC(O)O` | radicals `[C]`, `[O]` | fixed (found after the spiro row) |
| `[O:1]>>[O+:1]` on explicit-H dimethyl ether | `C[OH+]C` | `C[O+]C` | fixed (found after the row above) |
| `[C;$(CO):1]>>[C:1]C(=O)C` on `C[C@H](N)O` | stereo kept | stereo dropped (native) | policy: the RDKit profile (`run_smirks_checked(..., rdkit_compat=True)`) keeps it |
| `[H:1]>>[H+:1]` on `[H+]` | `[H+]` | no product | open (edge) |
| `[n:1]>>[n+2:1]` on pyrrole | `c1cc[nH+2]c1` | no product | open (edge) |

The naphthalene cation is non-aromatic after RDKit's sanitize, so its Kekulé
spelling depends on atom order; the source gives the same set of cation
positions, but RDKit may write a different Kekulé form for the same position.

## Grammar tweaks

The 3,000-example tweak run (seed 1) differs in 16 kinds for v1.0.34 and 10
for the source. Besides queue 3 flags, the source still differs on: `k0`
(RDKit matches every atom), the quadruple bond `$` (rejected), `^3` on a
carboxylate O⁻ (`[C,^3]-[C;R0]` on `CC(=O)[O-]`: RDKit 1 match, chematic 2),
`[c;!h1:1]>>[c:1]C` on explicit-H toluene (RDKit also returns the
over-valent ipso product) and `[!c:1]>>[c:1]~C` on pyrrole (RDKit
`C~c1cccc1`, chematic no product).

## Published v1.0.34 reaction gate (83 rows)

The checked RDKit-compatible profile of the published binaries was run on
the 83 pinned reaction rows. The published Linux x86-64 wheel is CPython 3.9
only and RDKit 2026.03.6 ships no 3.9 wheel, so
[`reaction_python_checked_candidates.py`](../scripts/reaction_python_checked_candidates.py)
ran the wheel in CPython 3.9 and
`reaction_python_checked_provenance_gate.py --candidates` classified the rows
against RDKit in CPython 3.11; the npm package ran through
[`run_published_npm_checked_candidates.mjs`](../scripts/run_published_npm_checked_candidates.mjs)
(Node 22).

| Artifact | Digest | Outcome |
|---|---|---|
| PyPI `chematic-1.0.34-cp39-cp39-manylinux_2_17_x86_64.manylinux2014_x86_64.whl` | sha256 `b4e3da6664f0d43c55606258a8beb97b0183ef12e9ee6fc2d6edf067f358f34b` (matches PyPI) | 80 graph/origin/map match, 3 joint invalid input |
| npm `@kent-tokyo/chematic@1.0.34` tarball | sha512 `E2/7tny7jbLriHP3fzD3AykeAbbXTKMYGJBBWLAEFiH/Q1R+WmFFJC21w2lQ32A6IXY5ador0DbHVfBNCQSrPg==` (matches registry integrity) | 80 graph/origin/map match, 3 joint invalid input |
| Source branch (CPython 3.11 release wheel) | — | 80 graph/origin/map match, 3 joint invalid input |

These are measurements of the named artifacts on one Linux x86-64 host; they
do not claim general SMARTS or SMIRKS parity, and the macOS and Windows
wheels were not run.
