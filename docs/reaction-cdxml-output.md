# Checked single-step reaction CDXML output

[Open Babel #2703](https://github.com/openbabel/openbabel/issues/2703) reports
reaction SMILES losing reactant/product roles and the arrow during CDXML export.
Open Babel 3.2.1 reproduces that loss on three original examples on this host:
one merged fragment, no reaction step and no arrow. chematic's existing writer
handles individual molecules; this addition provides an explicit reaction writer.

## Use

```python
import chematic

# Coordinates follow reactants (CCO, O), agent (Cl), then product (CC=O).
coords = [
    [[0, 0], [30, 15], [60, 0]],
    [[90, 0]],
    [[150, -45]],
    [[240, 0], [270, 15], [300, 0]],
]
text = chematic.write_reaction_cdxml_checked(
    "CCO.O>Cl>CC=O", coords, [[100, 0], [200, 0]]
)
```

Coordinates use the CDXML Y-down convention; arrow endpoints are **tail, head**.
Each component needs exactly one pair per atom, in native parsed atom order.
There is no automatic layout or position-based role inference. Component order
and duplicates are preserved, including duplicate molecules on the same side.
The examples test data transport, not chemical reaction feasibility or balance.

Rust exports `chematic_mol::write_reaction_cdxml_checked(&Reaction,
&[Vec<[f64; 2]>], [[f64; 2]; 2]) -> Result<String, ReactionCdxmlError>`.
Python parses reaction SMILES with bounded native parsing and consumes coordinate
iterables incrementally. It accepts nested numeric sequences and lazy outer/row
iterables, consuming at most one excess item to detect a shape mismatch.
Original iterator exceptions propagate; arbitrary caller callbacks have no timeout.

## Supported profile

One page contains distinct fragments and globally unique object IDs, one arrow,
and one scheme/step. The step records reactants, products, objects above the arrow
(agents), and the arrow's ID. Atom element, isotope, formal charge, bracket H
count, indexed connectivity and accepted bond orders are written explicitly.
The role attributes follow the published
[CDXML reaction-step specification](https://iupac.github.io/IUPAC-FAIRSpec/cdx_sdk/ReactionStep.htm).

The initial graph profile is non-aromatic, unmapped, without query/pseudoatoms,
R groups, stereo or enhanced stereo groups. Accepted bonds are single, triple
and **terminal double** bonds (at least one degree-one endpoint). Internal double
bonds are refused to avoid deriving unspecified E/Z stereo from supplied layout.
Bracket H counts must equal native normal-valence inference. That model infers
implicit H only for organic-subset atoms, so bracket H on other elements (for
example `[SiH4]`) is outside this profile. Isolated or bridging
H atoms are refused; these checks prevent known reader changes such as `[O]`
becoming water or `[H]` acquiring another H. Unsupported graphs return errors
instead of silently dropping annotations. This is not a general valence sanitizer.

| Bound | Maximum |
| --- | ---: |
| Total components | 32 |
| Atoms per component | 256 |
| Total atoms | 1,024 |
| Total bonds | 4,096 |
| Coordinate/arrow component magnitude | 1,000,000 |
| Python reaction SMILES bytes | 16,384 |
| Python parsed bonds per component | 512 |

At least one reactant and product are required. Coordinate pairs must be finite;
arrow endpoints must differ. The native output buffer is reserved fallibly after
validation. Invalid shape/profile/limits raise Python `ValueError`, incompatible
types raise `TypeError`, and native output reservation failure raises `MemoryError`.

This writer takes a graph-based `Reaction`, not a `ReactionDocument`; source
metadata, conditions, labels and atom-map correspondences are outside its input
contract. Multi-step writing, automatic depiction, CDX output and WASM bindings
are separate work. Existing molecular readers still project molecular fragments
and do not provide a new reaction reader. ChemDraw application acceptance has
not been tested.

## Local evidence

`validation/openbabel_issues/reaction-cdxml/probe.py` supplies eight original
symbolic recipes, all role-order reversals and four coordinate transforms:
256 documents. Independent XML checks verify roles, globally unique IDs, local
bond endpoints, exact coordinates and arrow endpoints. Indigo 1.46.0 reads all
three roles and their molecules correctly in all 256 cases. RDKit 2026.09.1
reads the reactants/products correctly in all 256 cases, but omits the agent role
in the 160 cases containing agents on this host. That reader difference remains
in the evidence; full RDKit reaction round-trip compatibility is not claimed.

The dated record at `validation/results/reaction-cdxml-2026-10-17/` retains
compressed documents/outputs, Open Babel diagnostics, versions, binary/source
hashes and test logs. These are local candidate checks, not a published package.
