#!/usr/bin/env python3
"""IO-3 dedicated kekulize/stereo option check, independent of the main
206-fixture oracle gate (`mrv_io_parity.py`). Verifies:

- kekulize=True/False: the written bond-order token shape matches the
  option (order="A" absent when True, present when False), and RDKit reads
  both variants back to the SAME canonical structure (kekulize is a
  representation choice, not a structural change).
- include_stereo=False: given a real RDKit-generated MRV fixture that
  already encodes stereo via a native wedge/dash bond, turning the option
  off correctly DROPS the stereo assignment when RDKit re-reads chematic's
  output (documented, expected loss -- not a bug). `include_stereo=True`'s
  round trip is already covered by `mrv_io_parity.py`'s phase2 check on
  the main fixture pool.

  This only holds for TETRAHEDRAL centers (wedge/dash on a single bond,
  gated by the option). E/Z double-bond stereo is encoded geometrically
  via the 2D coordinates of the four substituent atoms, which are always
  written regardless of `include_stereo` -- RDKit (or chematic's own
  `assign_ez_from_2d`) perceives E/Z directly from that layout, so
  `include_stereo=False` has NO effect on E/Z bonds by design, not a bug.
  E/Z cases are checked only for connectivity + stereo-presence-on-write,
  not for a "without_stereo drops it" expectation.

Usage:
    python scripts/mrv_kekulize_stereo_check.py --dump <mrv_kekulize_stereo.json>
    python scripts/mrv_kekulize_stereo_check.py --self-test
"""

from __future__ import annotations

import argparse
import json
import sys

try:
    from rdkit import Chem
except ImportError:
    print("rdkit is required", file=sys.stderr)
    raise


def mol_from_mrv(block):
    try:
        mol = Chem.MolFromMrvBlock(block)
    except Exception:
        return None
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    return mol


def has_defined_stereo(mol):
    if mol is None:
        return False
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    ez = [b for b in mol.GetBonds() if b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)]
    return bool(centers) or bool(ez)


def check_kekulize(case):
    mol_true = mol_from_mrv(case["kekulize_true_mrv"])
    mol_false = mol_from_mrv(case["kekulize_false_mrv"])
    ok = (
        case["kekulize_true_has_aromatic_token"] is False
        and case["kekulize_false_has_aromatic_token"] is True
        and mol_true is not None
        and mol_false is not None
        and Chem.MolToSmiles(mol_true) == Chem.MolToSmiles(mol_false)
    )
    return {
        "id": case["id"],
        "kind": "kekulize",
        "ok": ok,
        "canonical_true": Chem.MolToSmiles(mol_true) if mol_true else None,
        "canonical_false": Chem.MolToSmiles(mol_false) if mol_false else None,
    }


def is_ez_case(mol):
    """True if the molecule's only defined stereo is double-bond E/Z (no
    tetrahedral centers) -- `include_stereo` has no effect on this kind,
    since E/Z is geometry-derived from always-written 2D coordinates."""
    if mol is None:
        return False
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    return not centers


def check_stereo(case):
    mol_original = mol_from_mrv(case["original_mrv"])
    mol_with = mol_from_mrv(case["with_stereo_mrv"])
    mol_without = mol_from_mrv(case["without_stereo_mrv"])

    original_has_stereo = has_defined_stereo(mol_original)
    with_has_stereo = has_defined_stereo(mol_with)
    without_has_stereo = has_defined_stereo(mol_without)

    # connectivity (ignoring stereo) must match in all three
    def flat_smiles(m):
        if m is None:
            return None
        m2 = Chem.Mol(m)
        Chem.RemoveStereochemistry(m2)
        return Chem.MolToSmiles(m2)

    connectivity_matches = flat_smiles(mol_original) == flat_smiles(mol_with) == flat_smiles(mol_without)
    ez_only = is_ez_case(mol_original)

    if ez_only:
        # include_stereo has no effect on geometry-derived E/Z -- only
        # check that stereo survives on write and connectivity is intact.
        ok = original_has_stereo and with_has_stereo and connectivity_matches
    else:
        ok = original_has_stereo and with_has_stereo and not without_has_stereo and connectivity_matches

    return {
        "id": case["id"],
        "kind": "stereo",
        "ez_only": ez_only,
        "ok": ok,
        "original_has_stereo": original_has_stereo,
        "with_stereo_has_stereo": with_has_stereo,
        "without_stereo_has_stereo": without_has_stereo,
        "connectivity_matches": connectivity_matches,
    }


def run(cases):
    results = []
    for case in cases:
        if case["kind"] == "kekulize":
            results.append(check_kekulize(case))
        else:
            results.append(check_stereo(case))
    all_ok = all(r["ok"] for r in results)
    return results, all_ok


def run_self_test():
    checks = []

    def mrv_for(smi, wedge=False):
        mol = Chem.MolFromSmiles(smi)
        return Chem.MolToMrvBlock(mol)

    # kekulize: both variants read back to the same structure -> ok
    benzene_true = mrv_for("c1ccccc1")  # kekulize=True default in RDKit's own writer
    case = {
        "id": "self_test_benzene",
        "kind": "kekulize",
        "kekulize_true_has_aromatic_token": False,
        "kekulize_false_has_aromatic_token": True,
        "kekulize_true_mrv": benzene_true,
        "kekulize_false_mrv": benzene_true.replace('order="1"', 'order="A"').replace('order="2"', 'order="A"'),
    }
    results, ok = run([case])
    checks.append(("kekulize_matching_structures_pass", ok))

    case_bad = dict(case)
    case_bad["kekulize_true_has_aromatic_token"] = True  # wrong -- should fail
    results, ok = run([case_bad])
    checks.append(("kekulize_wrong_token_shape_fails", ok is False))

    # stereo: a real chiral fixture
    chiral_mol = Chem.MolFromSmiles("N[C@@H](C)C(=O)O")
    chiral_mrv = Chem.MolToMrvBlock(chiral_mol)
    flat_mol = Chem.MolFromSmiles("NC(C)C(=O)O")
    flat_mrv = Chem.MolToMrvBlock(flat_mol)
    stereo_case_ok = {
        "id": "self_test_stereo_ok",
        "kind": "stereo",
        "original_mrv": chiral_mrv,
        "with_stereo_mrv": chiral_mrv,
        "without_stereo_mrv": flat_mrv,
    }
    results, ok = run([stereo_case_ok])
    checks.append(("stereo_correct_drop_passes", ok))

    stereo_case_bad = {
        "id": "self_test_stereo_bad",
        "kind": "stereo",
        "original_mrv": chiral_mrv,
        "with_stereo_mrv": flat_mrv,  # stereo missing even though include_stereo=True -- should fail
        "without_stereo_mrv": flat_mrv,
    }
    results, ok = run([stereo_case_bad])
    checks.append(("stereo_missing_when_expected_fails", ok is False))

    ok_all = True
    for name, passed in checks:
        status = "OK" if passed else "FAIL"
        print(f"  self-test {name}: {status}")
        ok_all = ok_all and passed
    return ok_all


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dump")
    p.add_argument("--summary-out", default=None)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()

    if args.self_test:
        sys.exit(0 if run_self_test() else 1)

    if not args.dump:
        p.error("--dump is required unless --self-test")

    with open(args.dump) as f:
        cases = json.load(f)

    results, all_ok = run(cases)
    print(json.dumps(results, indent=2))
    if args.summary_out:
        with open(args.summary_out, "w") as f:
            json.dump({"all_ok": all_ok, "results": results}, f, indent=2)
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
