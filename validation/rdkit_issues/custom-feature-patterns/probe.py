"""RDKit #9601/#9615 wrapper reproduction and direct SMARTS mask oracle."""
import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator as FG


def requested_masks(mol, patterns):
    masks = [0] * mol.GetNumAtoms()
    for index, text in enumerate(patterns):
        query = Chem.MolFromSmarts(text)
        for match in mol.GetSubstructMatches(query):
            for atom in match:
                masks[atom] |= 1 << index
    return masks


def wrapper_probe():
    mol = Chem.MolFromSmiles("CCO")
    cases = []
    for patterns in [["O"], ["C"], [], ["O"] * 32, ["O"] * 33]:
        generator = FG.GetMorganFeatureAtomInvGen(patterns=[Chem.MolFromSmarts(p) for p in patterns])
        output = FG.AdditionalOutput()
        output.AllocateAtomToBits()
        FG.GetMorganGenerator(radius=0, atomInvariantsGenerator=generator).GetSparseCountFingerprint(mol, additionalOutput=output)
        cases.append({"patterns": patterns, "actual_atom_masks": [list(v) for v in output.GetAtomToBits()],
                      "requested_atom_masks": requested_masks(mol, patterns) if len(patterns) <= 32 else None})
    return cases


def compare():
    import chematic
    smiles = ["CCO", "NCC(=O)O", "COC", "C=O", "O=C=O", "CCC", "c1ccccc1",
              "C1=CC=CC=C1", "n1ccccc1", "N1=CC=CC=C1", "[13CH3]CO", "Oc1ccccc1",
              "C1CCCCC1", "NCCN", "CC(=O)N"]
    definitions = ["[#8]", "[#6]", "C=O", "[a]", "[N,O]", "[13C]", "[$([N;H2])]", "[R]"]
    sets = [definitions[:n] for n in range(9)] + [["O"] * 32]
    cases = []
    for smi in smiles:
        native = chematic.from_smiles(smi)
        oracle = Chem.MolFromSmiles(smi)
        for patterns in sets:
            expected = requested_masks(oracle, patterns)
            actual = chematic.FeaturePatternSet(iter(patterns)).invariants(native)
            assert actual == expected, (smi, patterns, actual, expected)
            cases.append({"smiles": smi, "patterns": patterns, "masks": actual})
    module = Path(__import__("chematic.chematic", fromlist=["__file__"]).__file__)
    return {"passed": len(cases), "cases": cases,
            "module_sha256": hashlib.sha256(module.read_bytes()).hexdigest()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare", action="store_true")
    args = parser.parse_args()
    result = {"rdkit_version": rdBase.rdkitVersion, "wrapper_cases": wrapper_probe()}
    if args.compare:
        result["direct_smarts_mask_comparison"] = compare()
    print(json.dumps(result, indent=2))
