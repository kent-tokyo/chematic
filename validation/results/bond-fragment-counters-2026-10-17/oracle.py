"""Run against the candidate native package and RDKit 2026.09.1.

MolJSON is used only to map concrete source bond indices; fragment graph
identity is independently read from SMILES and V3000 (never MolJSON/V2000).
"""
import itertools
import json
from pathlib import Path

import chematic as ch
from rdkit import Chem, rdBase

SAMPLES = [
    "CC", "CCO", "CCOC", "CC(C)O", "CC(C)(C)C", "CCCCCC", "CC(=O)O",
    "CC(=O)NCC", "CCOC(=O)C", "CC[NH3+]", "CC(=O)[O-]", "[13CH3:7]CO",
    "[2H]OC", "[H]OC", "CC#N", "CCC=C", "c1ccccc1CO", "c1ccccc1CC(=O)O",
    "Cc1ncccc1", "CC1CCCCC1", "C1CC1CCO", "CCO.CN", "NCCSCCP",
    "ClCCBr", "CCS(=O)(=O)C", "CCOP(=O)(O)O", "CC[Si](C)(C)O",
    "COc1ccc(C(=O)NCC)cc1", "CCN(CC)CC", "CC(=O)OC1=CC=CC=C1C(=O)O",
]
params = Chem.SmilesParserParams()
params.removeHs = False


def atom_key(atom):
    return (atom.GetAtomicNum(), atom.GetIsotope(), atom.GetFormalCharge(),
            atom.GetAtomMapNum(), atom.GetTotalNumHs())


def bond_keys(mol):
    return sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                   max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                   b.GetBondTypeAsDouble()) for b in mol.GetBonds())


def main():
    rows = []
    for smiles in SAMPLES:
        source = ch.from_smiles(smiles)
        reference = Chem.MolFromSmiles(smiles, params)
        assert reference is not None, smiles
        graph = json.loads(source.to_moljson())
        edges = [(int(b["source_id"][1:]) - 1, int(b["target_id"][1:]) - 1)
                 for b in graph["bonds"]]
        eligible = [i for i, (a, b) in enumerate(edges)
                    if reference.GetBondBetweenAtoms(a, b).GetBondType()
                    == Chem.BondType.SINGLE
                    and not reference.GetBondBetweenAtoms(a, b).IsInRing()]
        selections = [(), *[(i,) for i in eligible]]
        selections += list(itertools.combinations(eligible, 2))
        if len(eligible) > 2:
            selections += [tuple(eligible), tuple(reversed(eligible))]
        for cuts in dict.fromkeys(selections):
            for mode in ("indices", "custom", "zero-high"):
                labels = None if mode == "indices" else [
                    (7 + i, 23 + i) if mode == "custom" else (0, 65535)
                    for i in range(len(cuts))]
                rd_indices, rd_labels = [], []
                for j, cut in enumerate(cuts):
                    a, b = edges[cut]
                    bond = reference.GetBondBetweenAtoms(a, b)
                    pair = labels[j] if labels is not None else (a, b)
                    rd_indices.append(bond.GetIdx())
                    rd_labels.append(pair if bond.GetBeginAtomIdx() == a else pair[::-1])
                expected_counts = [99] * source.atom_count
                if cuts:
                    expected = Chem.FragmentOnBonds(
                        reference, rd_indices, dummyLabels=rd_labels,
                        cutsPerAtom=expected_counts)
                else:
                    # RDKit rejects an empty cut list; chematic defines a clone.
                    expected = Chem.Mol(reference)
                    expected_counts[:] = [0] * source.atom_count
                counts = [-99] * source.atom_count
                result = ch.fragment_on_bonds_checked(
                    source, cuts, dummy_labels=labels, cuts_per_atom=counts)
                assert counts == expected_counts == list(result.cuts_per_atom)
                canonical = Chem.MolToSmiles(expected)
                parsed = Chem.MolFromSmiles(result.molecule.smiles, params)
                assert parsed is not None
                assert Chem.MolToSmiles(parsed) == canonical, (smiles, cuts, mode)
                # V3000 preserves atom order, so compare every indexed atom and
                # edge too. Query metadata of R dummy records is not compared.
                block = result.molecule.to_mol_v3000(
                    [[0.0, 0.0]] * result.molecule.atom_count)
                indexed = Chem.MolFromMolBlock(block, removeHs=False)
                assert indexed is not None, (smiles, cuts, mode)
                assert [atom_key(a) for a in indexed.GetAtoms()] == [
                    atom_key(a) for a in expected.GetAtoms()], (smiles, cuts, mode)
                assert bond_keys(indexed) == bond_keys(expected), (smiles, cuts, mode)
                for j, attachment in enumerate(result.attachments):
                    cut = cuts[j // 2]
                    a, b = edges[cut]
                    pair = labels[j // 2] if labels is not None else (a, b)
                    assert attachment.source_bond == cut
                    assert attachment.anchor == (b if j % 2 == 0 else a)
                    assert attachment.dummy == source.atom_count + j
                    assert attachment.isotope == pair[j % 2]
                    dummy = indexed.GetAtomWithIdx(attachment.dummy)
                    assert dummy.GetAtomicNum() == 0
                    assert dummy.GetIsotope() == attachment.isotope
                    assert [a.GetIdx() for a in dummy.GetNeighbors()] == [attachment.anchor]
                rows.append(dict(smiles=smiles, cuts=cuts, labels=mode,
                                 counts=counts, canonical=canonical))
    output = dict(rdkit=rdBase.rdkitVersion, molecules=len(SAMPLES), cases=len(rows),
                  passed=len(rows), failures=0, rows=rows)
    Path(__file__).with_name("oracle.json").write_text(json.dumps(output, indent=2) + "\n")
    print(f"RDKit {rdBase.rdkitVersion}: {len(rows)} cases / {len(SAMPLES)} molecules passed")


if __name__ == "__main__":
    main()
