"""RDKit oracle for explicitly supplied, bounded concrete-core matches."""
import json
from pathlib import Path

import chematic as ch
from rdkit import Chem, rdBase

SOURCES = [
    "CC", "CCC", "CCO", "CCOC", "CC(C)O", "CC(C)(C)C", "CCCCCC",
    "CC(=O)O", "CC(=O)NCC", "CCOC(=O)C", "CC[NH3+]", "CC(=O)[O-]",
    "[13CH3:7]CO", "[2H]OC", "[H]OC", "CC#N", "CCC=C", "c1ccccc1CO",
    "c1ccccc1CC(=O)O", "Cc1ncccc1", "CC1CCCCC1", "C1CC1CCO", "CCO.CN",
    "NCCSCCP", "ClCCBr", "CCS(=O)(=O)C", "CCOP(=O)(O)O",
    "CC[Si](C)(C)O", "COc1ccc(C(=O)NCC)cc1", "CCN(CC)CC",
    "CC(=O)OC1=CC=CC=C1C(=O)O",
    "CC(=O)Oc1ccccc1C(=O)O",
]
CORES = ["C", "CC", "CO", "O", "N", "CN", "C(=O)O", "c1ccccc1",
         "C1CC1", "C#N", "[13CH3]", "[NH3+]", "[2H]", "C.C"]


def atom_key(atom):
    return (atom.GetAtomicNum(), atom.GetIsotope(), atom.GetFormalCharge(),
            atom.GetAtomMapNum(), atom.GetTotalNumHs())


def eligible(source, core, matching):
    # Independent scope filtering; native refusals are not silently skipped.
    for i, target in enumerate(matching):
        a, b = source.GetAtomWithIdx(target), core.GetAtomWithIdx(i)
        if (a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge(), a.GetIsAromatic()) != (
            b.GetAtomicNum(), b.GetIsotope(), b.GetFormalCharge(), b.GetIsAromatic()
        ):
            return False
    selected = set(matching)
    return all(not bond.IsInRing() and bond.GetBondType() == Chem.BondType.SINGLE
               for bond in source.GetBonds()
               if (bond.GetBeginAtomIdx() in selected) != (bond.GetEndAtomIdx() in selected))


def main():
    params = Chem.SmilesParserParams()
    params.removeHs = False
    rows, refusals, empty, indexed_cases = [], [], 0, 0
    for text in SOURCES:
        source, original = ch.from_smiles(text), Chem.MolFromSmiles(text, params)
        assert original is not None
        native_edges = json.loads(source.to_moljson())["bonds"]
        native_atoms = json.loads(source.to_moljson())["atoms"]
        for core_text in CORES:
            core, rd_core = ch.from_smiles(core_text), Chem.MolFromSmiles(core_text, params)
            core_atoms = json.loads(core.to_moljson())["atoms"]
            for matching in original.GetSubstructMatches(rd_core, uniquify=False, maxMatches=256):
                if not eligible(original, rd_core, matching):
                    continue
                fields = ("element", "isotope", "charge", "aromatic")
                if any(tuple(native_atoms[target][key] for key in fields) !=
                       tuple(core_atoms[i][key] for key in fields) for i, target in enumerate(matching)):
                    # Native strict representation differs from sanitized RDKit
                    # perception for Kekule/aromatic inputs. Verify refusal;
                    # retain this boundary explicitly rather than hiding it.
                    try:
                        ch.replace_core_checked(source, core, matching)
                    except ValueError:
                        refusals.append(dict(source=text, core=core_text, matching=matching,
                                             reason="native representation mismatch"))
                        continue
                    raise AssertionError("native representation mismatch was accepted")
                try:
                    result = ch.replace_core_checked(source, core, matching)
                except ValueError as error:
                    raise AssertionError((text, core_text, matching)) from error
                expected = Chem.ReplaceCore(original, rd_core, matching, labelByIndex=True)
                assert expected is not None
                Chem.SanitizeMol(expected)
                parsed = Chem.MolFromSmiles(result.molecule.smiles, params)
                assert parsed is not None
                assert Chem.MolToSmiles(parsed) == Chem.MolToSmiles(expected), (text, core_text, matching)
                retained = [i for i in range(source.atom_count) if i not in matching]
                wanted_map = tuple(retained.index(i) if i in retained else None for i in range(source.atom_count))
                assert result.original_to_output == wanted_map
                reference_attachments = sorted((a.GetIsotope(), retained[a.GetNeighbors()[0].GetIdx()])
                                               for a in expected.GetAtoms() if a.GetAtomicNum() == 0)
                assert sorted((a.core_atom, a.source_anchor) for a in result.attachments) == reference_attachments
                for attachment in result.attachments:
                    assert matching[attachment.core_atom] == attachment.source_core_atom
                    assert wanted_map[attachment.source_anchor] == attachment.anchor
                    bond = native_edges[attachment.source_bond]
                    ends = {int(bond["source_id"][1:]) - 1, int(bond["target_id"][1:]) - 1}
                    assert ends == {attachment.source_core_atom, attachment.source_anchor}
                if result.molecule.atom_count:
                    block = result.molecule.to_mol_v3000([[0., 0.]] * result.molecule.atom_count)
                    indexed = Chem.MolFromMolBlock(block, removeHs=False)
                    assert indexed is not None
                    for old, new in enumerate(wanted_map):
                        if new is not None:
                            assert atom_key(indexed.GetAtomWithIdx(new)) == atom_key(original.GetAtomWithIdx(old))
                    for a in result.attachments:
                        dummy = indexed.GetAtomWithIdx(a.dummy)
                        assert atom_key(dummy) == (0, a.core_atom, 0, 0, 0)
                        assert [n.GetIdx() for n in dummy.GetNeighbors()] == [a.anchor]
                    indexed_cases += 1
                else:
                    assert expected.GetNumAtoms() == 0
                    empty += 1
                rows.append(dict(source=text, core=core_text, matching=matching,
                                 canonical=Chem.MolToSmiles(expected),
                                 original_to_output=wanted_map, attachments=reference_attachments))
    directory = Path(__file__).parent
    report = dict(rdkit=rdBase.rdkitVersion, source_molecules=len(SOURCES), core_graphs=len(CORES),
                  cases=len(rows), passed=len(rows), empty_outputs=empty,
                  indexed_v3000_cases=indexed_cases,
                  representation_refusals=len(refusals), refusal_details=refusals, failures=0)
    (directory / "oracle.json").write_text(json.dumps(report, indent=2) + "\n")
    (directory / "oracle.jsonl").write_text("".join(json.dumps(r, separators=(",", ":")) + "\n" for r in rows))
    print(json.dumps(report))


if __name__ == "__main__":
    main()
