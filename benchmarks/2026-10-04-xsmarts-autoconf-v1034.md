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
| Flags giving RDKit's value | 42 / 83 | 65 / 83 | **75 / 83** |

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
| `smirks.undefined_product_bond` | first alternative of `=,:` | typed error | typed refusal | 3 | aligned |
| `smirks.undefined_order_both_sides` | first alternative | typed error | typed refusal | 3 | aligned |
| `products.grouped_product_component` | `(A.B)` is one object | parse error | typed refusal | 3 | aligned |
| `match.ring_count_basis` | symmetrized rings | SSSR | documented profile | 4 | differs |
| `smirks.valence_invalid_raw` | returns unsanitized product | drops at apply | policy: sanitized views agree | 3 | differs |
| `smirks.product_h_atom_added` | H atom dropped, `CC[O-]` | no product | policy: sanitized views agree | 3 | differs |
| `smirks.colon_product_bond` | `:` sets aromatic flags | aromatic bond, atoms aliphatic | policy: sanitized views agree | 3 | differs (atoms now aromatic, product dropped) |
| `smirks.hcount_query_explicit_h` | `C[CH6]` | `CC` | policy: implicit-H equivalence | 3 | differs |
| `smirks.product_h_pin` | `C[CH5]` with explicit H | `CC` | policy: implicit-H equivalence | 3 | differs |
| `smirks.star_on_explicit_h` | also edits H atoms | heavy atoms only | policy: implicit-H equivalence | 3 | differs |
| `smirks.fragmented_reactant` | one template with `.` matches a pair | typed error | adapter difference | 3 | differs |

No wrong-confident flag remains in the source, and each of the 8 that differ
has a decided outcome, listed with its values and a note in
[`validation/xsmarts_autoconf_expected.json`](../validation/xsmarts_autoconf_expected.json):

- **Sanitized views agree.** RDKit returns a raw product that its own
  `SanitizeMol` rejects (`C=OC`, `CC[OH-]`, aromatic atoms outside a ring);
  chematic returns only products that pass that check, so after sanitizing
  both give nothing.
- **Implicit-H equivalence.** An explicit-H reactant gives the products of
  its implicit-H form, as #734 item 5 asked; every RDKit product in these
  cases (`C[CH6]`, `C[CH5]`, `C[H]Cl`) fails its own `SanitizeMol`.
- **Adapter difference.** The tool's RDKit adapter passes the molecule once
  per reactant template; `run_smirks` given the same two reactants returns
  RDKit's `COCO`.
- **Documented profile.** Native ring counts stay SSSR.
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
| `[H:1]>>[H+:1]` on `[H+]` | `[H+]` | product removed by `remove_hydrogens` | fixed: isolated H atoms are kept, as by RDKit's `RemoveHs` |
| `[n:1]>>[n+2:1]` on pyrrole | `c1cc[nH+2]c1` | no product | decided: RDKit's own SMILES parser rejects this spelling; typed refusal in the checked API |
| `[C:1]>>[N:1]` on `[13CH3]C` | `CN`, `[13CH3]N` | `[13CH3]N` only | fixed: an element change re-derives H |
| `[c:1]>>[c:1]` / `[n:1]>>[n:1]` on explicit-H pyrrole | `c1cc[nH]c1` | `c1ccnc1` / no product | fixed: `remove_hydrogens` keeps the N's H; aromatic atoms keep their H count |
| `[N+:1]>>[N+:1]~C` on `C[NH3+]`, `[N:1]>>[N+:1]~C` on `CN(C)C` | `C[NH3+]~C`, `C[NH+](C)(C)~C` | `C[NH2+]~C`, `C[N+](C)(C)~C` | policy: RDKit's raw products are five-valent N and fail its sanitize |

The naphthalene cation is non-aromatic after RDKit's sanitize, so its Kekulé
spelling depends on atom order; the source gives the same set of cation
positions, but RDKit may write a different Kekulé form for the same position.

The final source sweep (8 seeds) stops at 6 findings: the two `~` rows and
the two native-stereo rows above, the `[nH+2]` row and the naphthalene
cation spelling.

## Grammar tweaks

The 3,000-example tweak run (seed 1) differs in 16 kinds for v1.0.34 and 4
for the final source: one queue 3 policy case, `k0` (RDKit matches every
atom), `[D4:1][H]>>[O:1]O` on explicit-H decalin (RDKit returns over-valent
products), and `[N+,^3:1]>>[N:1]` on explicit-H methylammonium (RDKit's
raw products keep the H atoms).

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
| crates.io `chematic` 1.0.34 ([`tools/published_rust_checked_gate`](../tools/published_rust_checked_gate)) | registry checksum `0ca1c7b6797ea106bdc44574f350e9c6561a38b75517085b6ecefaa2e07ce974` (lockfile) | 80 graph/origin/map match, 3 joint invalid input |
| Source branch (CPython 3.11 release wheel) | — | 80 graph/origin/map match, 3 joint invalid input |

Row-level reports:
[PyPI](../validation/results/v1.0.34-published-pypi-linux-cp39-reaction-checked-provenance-83.json),
[npm](../validation/results/v1.0.34-published-npm-reaction-checked-provenance-83.json),
[crate](../validation/results/v1.0.34-published-crate-reaction-checked-provenance-83.json).
The three published artifacts return the same rows apart from one invalid
reactant's error text.

## Published v1.0.34 npm formula and E/Z JSON

[`check_published_npm_formula_ez_json.mjs`](../scripts/check_published_npm_formula_ez_json.mjs)
reruns the two v1.0.31 WASM fixes on the npm 1.0.34 tarball. `formula()` on
the first 5,000 rows of `scripts/chembl_accuracy_corpus_4999.smi` (sha256
`1c47371d…`) equals published Python `Mol.formula` (PyPI Linux wheel,
CPython 3.9, via [`emit_python_formulas.py`](../scripts/emit_python_formulas.py))
on all 5,000 rows; published v1.0.30 npm spelled 1,575 differently. Five E/Z
cases through `run_reactants` and `enumerate_library_2way` all return
parseable JSON whose products keep their `/`/`\` markers, and RDKit reads
each product as the expected E/Z isomer.
[Report](../validation/results/v1.0.34-published-npm-formula-ez-json.json).

## Stereo in `random_smiles` and `write` (found while checking reaction products)

On the 1,687 rows of the exposed 10k corpus with `@`, `/` or `\`, three
`random_smiles` seeds each (5,061 strings) were parsed by RDKit 2026.03.6 and
compared with the input's canonical SMILES: 1,999 differed with v1.0.34
source, none with the fix. Writing each parsed random SMILES (and its
`add_hydrogens`/`remove_hydrogens` round trip) again with `write` and
parsing it back gave a different canonical SMILES for 3,965 of 10,122
molecules before the `write` fix and none after.

These are measurements of the named artifacts on one Linux x86-64 host; they
do not claim general SMARTS or SMIRKS parity, and the macOS and Windows
wheels were not run.
