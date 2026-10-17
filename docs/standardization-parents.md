# Parent transformations

Parent methods reduce selected kinds of molecular variability for grouping or
deduplication. They already run in chematic's Rust core. Their Python signatures
and generated reference are now published together, addressing the discovery
problem described in [RDKit #9651](https://github.com/rdkit/rdkit/issues/9651).

The names below help locate related RDKit convenience functions. They do not
imply identical chemistry, cleanup order, defaults or output strings.

| Related RDKit name | chematic method | Native scope / Python return |
| --- | --- | --- |
| FragmentParent | [`fragment_parent`](api/chematic.md#chematic.Mol.fragment_parent) | Select one classified/ranked fragment → new `Mol` |
| ChargeParent | [`charge_parent`](api/chematic.md#chematic.Mol.charge_parent) | Fragment selection, then local charge neutralization → new `Mol` |
| IsotopeParent | [`isotope_parent`](api/chematic.md#chematic.Mol.isotope_parent) | Remove isotope labels across all fragments → new `Mol` |
| StereoParent | [`stereo_parent`](api/chematic.md#chematic.Mol.stereo_parent) | Remove stored stereo annotations across all fragments → new `Mol` |
| TautomerParent | [`tautomer_parent`](api/chematic.md#chematic.Mol.tautomer_parent) | Native bounded tautomer selection → `(Mol, str)` |
| SuperParent | [`super_parent`](api/chematic.md#chematic.Mol.super_parent) | Compose the five steps below → `(Mol, str)` |
| Inspect intermediate stages | [`super_parent_report`](api/chematic.md#chematic.Mol.super_parent_report) | Final SMILES/status and five named intermediate SMILES → `SuperParentReport` |

## Choose the axis you intend to remove

`fragment_parent` uses the default native classification of water, ions and
structural salts. It chooses among retained candidates by heavy-atom count,
then carbon preference and canonical-SMILES tie breaking. When every fragment
is classified for removal, it returns a size-ranked fallback. The Python method
returns only the selected molecule; fragment selection reasons/warnings are
available through the Rust transformation record, not this Python return.

`charge_parent` first selects a fragment, then applies the native local rules
for carbon-bound O- and protonated N+/O+ motifs. A quaternary ammonium ion may
remain charged. This is not a pH model or a promise of zero net charge. Isotope
and stereo parent steps do not run here.

`isotope_parent` and `stereo_parent` operate on the entire input, retaining all
fragments. Removing an isotope label does not delete an explicit graph atom:
deuterium becomes ordinary explicit H. Isotope removal does not run the stereo
parent step. Stereo removal clears atom chirality, directional single bonds and
stored stereo metadata while retaining isotope labels and bond orders.

Each returned `Mol` is a new molecule. The input graph and Python properties
are unchanged; the returned molecule's Python property store starts empty.

```python
import chematic

salt = chematic.from_smiles("CC(=O)[O-].[Na+]")
salt.SetProp("source", "example")
fragment = salt.fragment_parent()
assert "." not in fragment.smiles
assert "[O-]" in fragment.smiles
assert "[O-]" not in salt.charge_parent().smiles
assert fragment.GetPropNames() == []
assert salt.GetProp("source") == "example"

labeled = chematic.from_smiles("[13CH3]O.[Na+]").isotope_parent()
assert "[13" not in labeled.smiles and "." in labeled.smiles
unstereo = chematic.from_smiles("N[C@@H](C)C(=O)O.[Na+]").stereo_parent()
assert "@" not in unstereo.smiles and "." in unstereo.smiles
```

## Compose and inspect the stages

The native composition order is **fragment → charge → isotope → stereo →
tautomer**. Every stage consumes the previous stage's output. This order is
specific to the Parent workflow and differs from `standardize`'s pipeline.

```python
import chematic

mol = chematic.from_smiles("[NH3+][C@@H]([2H])C(=O)[O-].Cl")
original = mol.smiles
report: chematic.SuperParentReport = mol.super_parent_report()
assert report["status"] == "Completed"
assert [stage["name"] for stage in report["stages"]] == [
    "fragment", "charge", "isotope", "stereo", "tautomer",
]
assert report["stages"][0]["smiles"] == "[2H][C@@H](C(=O)[O-])[NH3+]"
assert report["stages"][1]["smiles"] == "[2H][C@@H](C(O)=O)N"
assert report["stages"][2]["smiles"] == "[H][C@@H](C(O)=O)N"
assert report["stages"][3]["smiles"] == "[H]C(C(O)=O)N"
assert report["smiles"] == report["stages"][-1]["smiles"]
assert mol.smiles == original
```

[`ParentStage`](api/chematic.md#chematic.ParentStage) is a `TypedDict` containing
`name` (one of the five literal stage names) and `smiles` (`str`).
[`SuperParentReport`](api/chematic.md#chematic.SuperParentReport) contains
`smiles: str`, `status: str`, and `stages: list[ParentStage]`. These describe the
existing ordinary mutable dictionaries; no new report computation is added.
Intermediate SMILES are not the full Rust audit: the report does not expose
applied-rule records or fragment-selection reasons.

## Handle bounded results

`tautomer_parent`, `super_parent` and `super_parent_report` accept the same
optional native limits:

| Argument | Default | Meaning |
| --- | --- | --- |
| `max_transforms` | `16` | Nonnegative native transformation budget |
| `max_tautomers` | `32` | Native candidate-comparison budget; zero skips the direct-aromatic-shift comparison, rather than all tautomer work |
| `timeout_ms` | `None` | Optional nonnegative wall-clock budget in milliseconds, for the tautomer stage only |

`Completed` means the native bounded algorithm completed. Check the returned
status before treating the result as an accepted grouping key. A limit or
timeout can return a partial result. Status strings observed by these APIs
include `Completed`, `MaxTransformsReached`, `MaxTautomersReached`, `TimedOut`,
and `InvalidInput(EmptyMolecule)`. The type remains `str`: the native enum can
be extended and reason payloads are represented in the string.

```python
mol = chematic.from_smiles("OC=C")
limited, status = mol.tautomer_parent(max_transforms=0)
assert status == "MaxTransformsReached"
assert limited.smiles == mol.smiles

timed, status = mol.super_parent(timeout_ms=0)
assert status == "TimedOut"

empty = chematic.RWMol().GetMol()
_, status = empty.super_parent()
assert status == "InvalidInput(EmptyMolecule)"
```

Unsigned native budgets reject negative/out-of-range integers with
`OverflowError`, and non-integer values with `TypeError`. A wall-clock timeout
can depend on machine load; it is not a deterministic identity guarantee.
The native Parent methods expose no cancellation callback. Their algorithms
also differ from `rdkit_canonical_tautomer` / `rdkit_tautomers`, the separate
RDKit-profile APIs. Completion is not proof that every chemical tautomer or
RDKit standardization option is covered.

The generated [API reference](api/chematic.md) includes all seven signatures,
defaults, arguments and returns. A source gate checks native/stub coverage; a
build gate checks the actual generated HTML anchors and this guide. Both run
in CI, and the HTML gate also runs before the existing Pages deployment.
