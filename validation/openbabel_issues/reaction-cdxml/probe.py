#!/usr/bin/env python3
"""Original single-step role/layout corpus with independent external readers.

PYTHONPATH=<candidate>:<rdkit> python3 probe.py <output-directory>
Requires Indigo, RDKit, and local Open Babel CLI. This does not test chemical
reaction feasibility, automatic layout or ChemDraw application acceptance.
"""
import gzip
import hashlib
import itertools
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import chematic
from indigo import Indigo
from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions

RECIPES = [
    (["CCO", "O"], ["Cl"], ["CC=O"]),
    (["CCO", "O"], [], ["CC=O", "O"]),
    (["C#N"], ["O", "Cl"], ["C(=O)N"]),
    (["[13CH3]O", "[Na+]"], ["[Cl-]"], ["[13CH3]O"]),
    (["O", "O"], [], ["O"]),
    (["[2H]O"], [], ["[2H]O"]),
    (["CC(=O)O"], ["N"], ["CC(=O)N", "O"]),
    (["N#N", "C"], ["[Na+]", "[Cl-]"], ["CO"]),
]


def layout(groups, mode):
    coordinates = []
    for role, molecules in enumerate(groups):
        for i, smiles in enumerate(molecules):
            mol = Chem.MolFromSmiles(smiles)
            positions = []
            for j in range(mol.GetNumAtoms()):
                x, y = role * 350 + i * 110 + j * 30, (j % 2) * 15 - (
                    60 if role == 1 else 0
                )
                if mode == 1:
                    x, y = x + 17, y - 33
                elif mode == 2:
                    x, y = 100 - y, 200 + x
                elif mode == 3:
                    x, y = -x * 0.5 - 70, y * 0.5 + 13
                positions.append([x, y])
            coordinates.append(positions)
    arrow = [[300, 0], [600, 0]]
    if mode == 1:
        arrow = [[x + 17, y - 33] for x, y in arrow]
    elif mode == 2:
        arrow = [[100 - y, 200 + x] for x, y in arrow]
    elif mode == 3:
        arrow = [[-x * 0.5 - 70, y * 0.5 + 13] for x, y in arrow]
    return coordinates, arrow


def main():
    output = Path(sys.argv[1])
    output.mkdir(parents=True, exist_ok=True)
    session = Indigo()
    records, failures, agent_diagnostics = [], [], []
    for recipe, reverse, mode in itertools.product(
        range(len(RECIPES)), itertools.product([False, True], repeat=3), range(4)
    ):
        groups = [
            list(reversed(g)) if rev else list(g)
            for g, rev in zip(RECIPES[recipe], reverse)
        ]
        smiles = ">".join(".".join(g) for g in groups)
        coordinates, arrow = layout(groups, mode)
        text = chematic.write_reaction_cdxml_checked(smiles, coordinates, arrow)
        root = ET.fromstring(text)
        fragments, step = root.findall(".//fragment"), root.find(".//step")
        object_ids = [node.attrib["id"] for node in root.iter() if "id" in node.attrib]
        assert len(object_ids) == len(set(object_ids))
        assert len(fragments) == sum(map(len, groups))
        start = 0
        for group, key in zip(
            groups,
            [
                "ReactionStepReactants",
                "ReactionStepObjectsAboveArrow",
                "ReactionStepProducts",
            ],
        ):
            assert step.attrib[key].split() == [
                fragment.attrib["id"]
                for fragment in fragments[start : start + len(group)]
            ]
            start += len(group)
        for fragment, points in zip(fragments, coordinates):
            assert [
                list(map(float, n.attrib["p"].split())) for n in fragment.findall("n")
            ] == points
            nodes = {n.attrib["id"] for n in fragment.findall("n")}
            assert all(
                b.attrib["B"] in nodes and b.attrib["E"] in nodes
                for b in fragment.findall("b")
            )
        arrow_node = root.find(".//arrow")
        assert step.attrib["ReactionStepArrows"] == arrow_node.attrib["id"]
        assert list(map(float, arrow_node.attrib["Tail3D"].split()))[:2] == arrow[0]
        assert list(map(float, arrow_node.attrib["Head3D"].split()))[:2] == arrow[1]
        reaction = session.loadReaction(text)
        indigo_roles = [
            [m.canonicalSmiles() for m in molecules]
            for molecules in [
                reaction.iterateReactants(),
                reaction.iterateCatalysts(),
                reaction.iterateProducts(),
            ]
        ]
        indigo_reference = [
            [session.loadMolecule(s).canonicalSmiles() for s in g] for g in groups
        ]
        if indigo_roles != indigo_reference:
            failures.append(
                dict(
                    case=len(records),
                    reader="indigo",
                    expected=indigo_reference,
                    actual=indigo_roles,
                )
            )
        read = rdChemReactions.ReactionsFromCDXMLBlock(text, True, False)
        assert len(read) == 1
        rdkit_roles = [
            [Chem.MolToSmiles(m) for m in g]
            for g in [
                read[0].GetReactants(),
                read[0].GetAgents(),
                read[0].GetProducts(),
            ]
        ]
        rdkit_reference = [
            [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in g] for g in groups
        ]
        if [rdkit_roles[0], rdkit_roles[2]] != [rdkit_reference[0], rdkit_reference[2]]:
            failures.append(
                dict(
                    case=len(records),
                    reader="rdkit_reactants_products",
                    expected=rdkit_reference,
                    actual=rdkit_roles,
                )
            )
        if rdkit_roles[1] != rdkit_reference[1]:
            agent_diagnostics.append(
                dict(
                    case=len(records),
                    expected=rdkit_reference[1],
                    actual=rdkit_roles[1],
                )
            )
        records.append(
            dict(
                recipe=recipe,
                reversed=list(reverse),
                layout=mode,
                reaction=smiles,
                coordinates=coordinates,
                arrow=arrow,
                cdxml=text,
                indigo_roles=indigo_roles,
                indigo_reference=indigo_reference,
                rdkit_roles=rdkit_roles,
                rdkit_reference=rdkit_reference,
            )
        )
    atom_states = []
    for smiles in ["[O]", "[N]", "[CH3]", "[H]"]:
        legacy = chematic.from_smiles(smiles).to_cdxml([[0, 0]])
        diagnostic = dict(
            input=smiles,
            indigo_reference=session.loadMolecule(smiles).canonicalSmiles(),
            indigo_legacy_read=session.loadMolecule(legacy).canonicalSmiles(),
            legacy_cdxml=legacy,
        )
        try:
            chematic.write_reaction_cdxml_checked(
                smiles + ">>O", [[[0, 0]], [[100, 0]]], [[40, 0], [80, 0]]
            )
            failures.append(
                dict(
                    reader="checked_atom_state",
                    input=smiles,
                    reason="expected profile refusal",
                )
            )
        except ValueError as error:
            diagnostic["checked_refusal"] = str(error)
        atom_states.append(diagnostic)
    baseline = []
    for i in [0, 1, 6]:
        groups = RECIPES[i]
        smiles = ">".join(".".join(g) for g in groups)
        result = subprocess.run(
            ["/opt/homebrew/bin/obabel", "-:" + smiles, "-ocdxml", "--gen2D"],
            text=True,
            capture_output=True,
            timeout=20,
        )
        (output / f"openbabel-{i}.cdxml").write_text(result.stdout)
        (output / f"openbabel-{i}.stderr").write_text(result.stderr)
        root = ET.fromstring(result.stdout)
        baseline.append(
            dict(
                reaction=smiles,
                exit_code=result.returncode,
                fragments=len(root.findall(".//fragment")),
                arrows=len(root.findall(".//arrow")),
                schemes=len(root.findall(".//scheme")),
                steps=len(root.findall(".//step")),
            )
        )
    summary = dict(
        schema_version=1,
        scope="bounded supplied single-step roles and layout; no reaction-feasibility or full ChemDraw claim",
        datasets=len(records),
        xml_role_layout_checks=len(records),
        indigo_cross_reads=len(records),
        rdkit_reactant_product_cross_reads=len(records),
        failures=failures,
        failure_count=len(failures),
        rdkit_agent_role_differences=agent_diagnostics,
        rdkit_agent_role_difference_count=len(agent_diagnostics),
        rdkit_version=rdBase.rdkitVersion,
        indigo_version=session.version(),
        openbabel_version=subprocess.check_output(
            ["/opt/homebrew/bin/obabel", "-V"], text=True
        ).strip(),
        openbabel_diagnostics=baseline,
        legacy_atom_state_diagnostics=atom_states,
        native_binary_sha256=hashlib.sha256(
            Path(chematic.chematic.__file__).read_bytes()
        ).hexdigest(),
        python_version=sys.version,
    )
    with gzip.open(output / "corpus-results.json.gz", "wt") as stream:
        json.dump(records, stream, separators=(",", ":"))
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k != "rdkit_agent_role_differences"},
            indent=2,
        )
    )
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
