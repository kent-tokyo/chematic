# How chematic applies a SMIRKS

This page answers the question asked in
[#734](https://github.com/kent-tokyo/chematic/issues/734): what chematic does
when it applies a reaction template, and where that differs from RDKit's
`RunReactants` on purpose. It describes the source after v1.0.35.

There are two entry points:

| Entry point | Profile |
|---|---|
| `run_smirks`, Rust `run_reactants` | native |
| `run_smirks_checked(..., rdkit_compat=True)`, Rust `run_reactants_traced_rdkit_2026_03_6` | RDKit 2026.03.6 profile: results RDKit would give, or a typed refusal |

Both share the steps below; the RDKit profile differs where noted.

## 1. Reading the template

**Reactant side** is SMARTS. A component that is also valid SMILES and carries
no stereo is matched with its SMARTS reading (`[*:1]` is any atom, the bond
between two aromatic atoms is single or aromatic); one with `@`/`@@` or
`/`/`\` keeps its SMILES reading in the native profile, which checks that
stereo on the matched atoms. In the RDKit profile tetrahedral and
double-bond stereo in a reactant template are not matched (RDKit does not
match them either); `/` and `\` match a single or aromatic bond, as in plain
SMARTS. The native profile refuses `/`/`\` in a component that is not SMILES,
since it cannot check them there.

Reactants are matched on RDKit's aromaticity view: a Kekulé benzene matches
`c`.

**Product side** is a specification, read as RDKit reads a product SMARTS:

| Written | Applied |
|---|---|
| element, isotope, charge, `H<n>`, `@`/`@@` | as written |
| `X`, `D`, `R`, `r`, `x`, `v`, `h`, `^`, `$()` | nothing (query-only) |
| a list of elements (`[c,n:1]`, `[F,Cl,Br,I;-:7]`) on a mapped atom | the matched atom's element; a charge, H count or isotope spelled beside it still applies |
| an H count or charge inside a list (`[N;X3H1+0,X4H2+:1]`) | the first one written (RDKit applies it with a warning) |
| a list of elements on an unmapped atom (`[#6,#7]`) | typed refusal (RDKit builds a dummy `*` atom) |
| a bond list (`=,:`, `-;!@`) | the first order primitive: `=,:` is double, `!=` is double (RDKit skips the `!`) |
| a bond starting with a ring primitive (`!@-`, `@=`, `@`) | a zero-order bond (RDKit's `UNSPECIFIED`): it joins the atoms but adds nothing to valence |
| `~` | between two atoms already bonded in the reactant, that bond; otherwise a zero-order bond |
| `(A.B)` | one product object |

`[#6:1]`-style atoms: the native profile enumerates aliphatic and aromatic
spellings and pairs each mapped product atom with its reactant atom's
spelling (#679), also when the reactant spells `[c:1]`; it reads a product
`[#6;a:1]` as aromatic. The RDKit profile, and the native profile when the
combinations exceed 256 (rules spelling dozens of such atoms), read the
reactant atoms as SMARTS and the product atoms as aliphatic, as RDKit does
(in the RDKit profile `a`/`A` beside no element symbol is query-only);
step 4 then perceives aromaticity again.

## 2. Building a product

For each match:

1. Each mapped atom is copied from the reactant, then the template's element
   (when it differs from the reactant template's), charge, isotope, H count
   and aromaticity are applied. An unspecified charge keeps the reactant's;
   an unspecified H count keeps the reactant's when the atom's template
   degree is unchanged and is re-derived from valence otherwise.
2. Product-only atoms are created; unmapped reactant-template atoms are
   deleted.
3. Bonds between template atoms come from the product template; a bond
   between two mapped atoms that the reactant template does not spell is
   kept (RDKit #1387). A template single bond over a reactant `/`/`\` bond
   keeps the direction, so the E/Z of an adjacent double bond survives; so
   does a direction read on an aromatic ring bond (`C/N=c1\c(O)c(O)c1`).
4. Atoms outside the template are carried with their bonds (native
   `run_reactants_strict` does not carry them).

## 3. Explicit hydrogens

A reactant with explicit H atoms (`add_hydrogens`) gives the products of its
implicit-H form: the H atoms on mapped atoms are folded into H counts and
re-derived after the edit. RDKit keeps them as atoms, so some RDKit products
are over-valent and fail its own sanitize (a CH₃ that gains an O while its
three mapped H atoms stay) where chematic gives the sanitizable product. An
H atom that carried a double bond's only `/`/`\` passes the direction to
another substituent of that end. Templates that spell `[H]` atoms match the
explicit H atoms, as in RDKit.

## 4. Sanitizing

RDKit returns raw products and `SanitizeMol` rejects the invalid ones;
chematic returns only products that pass that check. When a template has
written Kekulé bonds into an aromatic ring, the product is first treated as
RDKit's sanitize treats it: the bonds still aromatic are kekulized, aromatic
flags are cleared and aromaticity is perceived again. Then each atom's
valence is checked against RDKit's allowed valences, aromatic atoms outside
rings are rejected, and RDKit's clean-up spellings (`[N+][O-]` for a neutral
five-valent N=O) are accepted.

Kekulization follows RDKit's candidate rule (an atom takes a double bond
when its bond orders, aromatic as 1, and H count come to one less than its
default valence) for atoms that already have a non-aromatic double bond, for
neutral N and P, and for `[o+]`/`[s+]`: a thiophene S-oxide `[s+]` and a
ring-fusion N beside a `c(=O)` take no double bond. Other atoms use
chematic's ordinary rule.

## 5. Stereo

Carried tetrahedral centres keep their configuration. In the native profile
a centre written `@`/`@@` in both templates follows the template; the RDKit
profile copies RDKit's raw tag handling (inversion flags, neighbour
permutation). A centre that gains a neighbour in place of its implicit H
keeps its geometry in the native profile; RDKit's profile gives RDKit's
opposite arrangement.

## Known differences from RDKit

These are policies, not bugs, and each is measured in the BioTransformer
rule corpus record
([benchmarks/2026-10-05-biotransformer-rule-corpus.md](../benchmarks/2026-10-05-biotransformer-rule-corpus.md)):

- explicit-H reactants give the implicit-H products (step 3): RDKit's
  products that keep H atoms over-valent (`[SH2]`), as radicals or failing
  its own sanitize differ;
- products RDKit's sanitize rejects are not returned (step 4);
- unmapped product atoms naming no element are refused, not dummy atoms;
- `R<n>` in a reactant template counts SSSR rings; RDKit counts its
  symmetrized ring set (a quinuclidine N is `R3` there);
- where RDKit drops the E/Z of a template-matched double bond whose stereo
  reference atoms are explicit H atoms, chematic keeps it;
- a template that matches `[H]` atoms leaves them as atoms; when such an H
  carries the direction of a conjugated double bond, removing the H atoms
  afterwards (RDKit `RemoveHs`, chematic `remove_hydrogens`) can drop that
  E/Z.
