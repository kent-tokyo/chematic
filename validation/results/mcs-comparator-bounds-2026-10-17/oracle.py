"""Compare MCS atom/bond counts, not SMARTS or returned Mol identity."""
import itertools
import json
import random
from pathlib import Path

import chematic as ch
from rdkit import Chem, rdBase
from rdkit.Chem import rdFMCS

SAMPLES = [
    "C", "N", "O", "CC", "NN", "CO", "CN", "NO", "CCC", "NNN", "OOO",
    "COC", "CCN", "NCCN", "CCCO", "CCCC", "NCCO", "COCO", "CC(C)C",
    "NC(N)N", "CC(C)O", "CCCCCC", "NNNNNN", "CCNCCO", "CCC=C", "CC=NN",
    "CC#N", "CCOC(=O)C", "CC(=O)NCC", "[H][H]", "[2H]O[2H]", "[H]OC",
    "CC[NH3+]", "[13CH3]CO", "ClCCBr",
]
ATOM_MODES = [
    (ch.McsAtomCompare.ELEMENTS, rdFMCS.AtomCompare.CompareElements),
    (ch.McsAtomCompare.ANY, rdFMCS.AtomCompare.CompareAny),
    (ch.McsAtomCompare.ANY_HEAVY_ATOM, rdFMCS.AtomCompare.CompareAnyHeavyAtom),
]
BOND_MODES = [
    (ch.McsBondCompare.ORDER_OR_AROMATIC, rdFMCS.BondCompare.CompareOrderExact),
    (ch.McsBondCompare.ANY, rdFMCS.BondCompare.CompareAny),
]


def main():
    pairs = list(itertools.combinations(SAMPLES, 2))
    random.Random(9601).shuffle(pairs)
    groups = pairs[:160] + [("CCC", "NNN"), ("CCCO", "NCCN"), ("COC", "NNN")]
    groups += [("CCCC", "NNN", "OOO"), ("CCCO", "NCCN", "COCO")]
    parser = Chem.SmilesParserParams()
    parser.removeHs = False
    rows = []
    for group in groups:
        for inputs in [group, group[::-1]]:
            mols = [ch.from_smiles(s) for s in inputs]
            reference = [Chem.MolFromSmiles(s, parser) for s in inputs]
            assert all(m is not None for m in reference)
            for atom, rd_atom in ATOM_MODES:
                for bond, rd_bond in BOND_MODES:
                    actual, timed_out = ch.find_mcs_checked(
                        mols, atom_compare=atom, bond_compare=bond
                    )
                    expected = rdFMCS.FindMCS(
                        reference, atomCompare=rd_atom, bondCompare=rd_bond, timeout=5
                    )
                    counts = [actual.atom_count, len(json.loads(actual.to_moljson())["bonds"])] if actual else [0, 0]
                    wanted = [expected.numAtoms, expected.numBonds]
                    profile_wanted = wanted
                    if atom is ch.McsAtomCompare.ANY_HEAVY_ATOM:
                        # chematic's existing mode excludes H altogether; RDKit
                        # CompareAnyHeavyAtom permits H-to-H. Check the native
                        # definition using the independent heavy induced graph,
                        # while retaining raw RDKit differences in every row.
                        heavy = []
                        for mol in reference:
                            edited = Chem.RWMol(mol)
                            for index in reversed(range(mol.GetNumAtoms())):
                                if mol.GetAtomWithIdx(index).GetAtomicNum() == 1:
                                    edited.RemoveAtom(index)
                            heavy.append(edited.GetMol())
                        if any(m.GetNumAtoms() == 0 for m in heavy):
                            profile_wanted = [0, 0]
                        else:
                            profile = rdFMCS.FindMCS(
                                heavy, atomCompare=rd_atom, bondCompare=rd_bond, timeout=5
                            )
                            assert not profile.canceled
                            profile_wanted = [profile.numAtoms, profile.numBonds]
                    row = dict(inputs=inputs, atom_mode=atom.value, bond_mode=bond.value,
                               actual=counts, expected=wanted,
                               native_profile_expected=profile_wanted,
                               native_timed_out=timed_out, rdkit_canceled=expected.canceled)
                    row["passed"] = counts == wanted and not timed_out and not expected.canceled
                    row["native_profile_passed"] = counts == profile_wanted and not timed_out and not expected.canceled
                    rows.append(row)
    failures = [row for row in rows if not row["passed"]]
    report = dict(rdkit=rdBase.rdkitVersion, samples=len(SAMPLES), groups=len(groups),
                  comparisons=len(rows), passed=len(rows) - len(failures),
                  failures=failures,
                  native_profile_passed=sum(r["native_profile_passed"] for r in rows))
    Path(__file__).with_name("oracle.json").write_text(json.dumps(report, indent=2) + "\n")
    Path(__file__).with_name("oracle.jsonl").write_text(
        "".join(json.dumps(row, separators=(",", ":")) + "\n" for row in rows)
    )
    print(f"RDKit {rdBase.rdkitVersion}: {report['passed']}/{len(rows)} MCS count comparisons")
    for failure in failures[:10]:
        print(failure)
    print(f"Native mode definitions (heavy-only induced reference): {report['native_profile_passed']}/{len(rows)}")
    if report["native_profile_passed"] != len(rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
