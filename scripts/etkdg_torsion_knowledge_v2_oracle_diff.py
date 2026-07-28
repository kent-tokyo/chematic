#!/usr/bin/env python3
"""RDKit oracle differential for the torsion-knowledge v2 layer (3D
Breakthrough Program, Wave 2, Agent E, spec section 12).

Run against an ISOLATED venv with RDKit installed -- never the repo's shared
`.venv` (this program's standing rule; see CLAUDE.md-adjacent process notes).
Example setup:

    python3 -m venv /tmp/oracle_venv
    /tmp/oracle_venv/bin/pip install rdkit
    /tmp/oracle_venv/bin/python3 scripts/etkdg_torsion_knowledge_v2_oracle_diff.py

What this script compares, and what it honestly cannot:

1. RING CLASSIFICATION (comparable): reads
   `validation/etkdg_torsion_knowledge_v2_chematic_side.json` (written by
   `cargo run --release -p chematic-3d --example
   torsion_knowledge_v2_gap_check`), independently computes RDKit's own
   per-atom element sequence and per-bond SSSR ring-size list for the SAME
   fixture SMILES, and reports mismatches. Atom-index correspondence between
   chematic-smiles's parser and `Chem.MolFromSmiles` is VERIFIED here (the
   element-sequence check), not assumed.

2. TORSION DISTRIBUTION / MINIMA (comparable, but only empirically): for a
   couple of representative bonds, generates an RDKit ETKDG conformer
   ensemble with useExpTorsionAnglePrefs=True (and, where relevant,
   useSmallRingTorsions/useMacrocycleTorsions=True) and reports the resulting
   dihedral-angle histogram, compared qualitatively against this PR's own
   predicted Fourier-term minima for the same bond.

3. MATCHED RULE FAMILY / 1-4 PAIR SELECTION (NOT achievable via RDKit's
   public Python API in any version checked): RDKit's internal
   ExperimentalTorsionAngle / BoundsMatrixBuilder machinery that performs
   this matching has no public accessor exposing which SMARTS rule matched
   which bond, or which 1-4 pairs it adjusted. This is a genuine, disclosed
   limitation of what spec section 12 asks for -- not silently narrowed.
   Reading RDKit's C++/.in source (already done, see
   validation/manifests/etkdg_torsion_knowledge_sources.json) is a
   translation-provenance activity, not a live differential, and does not
   substitute for this missing piece.
"""

import json
import statistics
import sys
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom
except ImportError:
    print(
        "RDKit not importable. Run this in an ISOLATED venv with `pip install rdkit`, "
        "never the repo's shared .venv.",
        file=sys.stderr,
    )
    sys.exit(1)

REPO_ROOT = Path(__file__).resolve().parent.parent
CHEMATIC_SIDE_PATH = REPO_ROOT / "validation" / "etkdg_torsion_knowledge_v2_chematic_side.json"
RESULT_PATH = REPO_ROOT / "validation" / "etkdg_torsion_knowledge_v2_rdkit_oracle_diff.json"


def rdkit_ring_sizes_per_bond(mol):
    ri = mol.GetRingInfo()
    out = {}
    for b in mol.GetBonds():
        idx = b.GetIdx()
        if ri.NumBondRings(idx) > 0:
            sizes = sorted(ri.BondRingSizes(idx))
        else:
            sizes = []
        out[(b.GetBeginAtomIdx(), b.GetEndAtomIdx())] = sizes
    return out


def compare_ring_classification(fixtures):
    print("--- 1. Ring classification (RDKit SSSR vs. chematic RingMembershipIndex) ---")
    results = []
    total_bonds = 0
    matching_bonds = 0
    atom_order_mismatches = []
    for fx in fixtures:
        name, smiles = fx["name"], fx["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  {name}: RDKit failed to parse {smiles!r}, skipped")
            continue
        rdkit_atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        if rdkit_atoms != fx["atoms"]:
            atom_order_mismatches.append(name)
            print(f"  {name}: ATOM-ORDER MISMATCH chematic={fx['atoms']} rdkit={rdkit_atoms}")
            continue  # cannot compare per-bond without verified atom correspondence

        rdkit_bonds = rdkit_ring_sizes_per_bond(mol)
        bond_diffs = []
        for b in fx["bonds"]:
            key = (b["a"], b["b"])
            rev_key = (b["b"], b["a"])
            chematic_sizes = sorted(b["ring_sizes"])
            rdkit_sizes = rdkit_bonds.get(key, rdkit_bonds.get(rev_key))
            if rdkit_sizes is None:
                bond_diffs.append((key, chematic_sizes, "NO_SUCH_BOND_IN_RDKIT"))
                continue
            total_bonds += 1
            if chematic_sizes == rdkit_sizes:
                matching_bonds += 1
            else:
                bond_diffs.append((key, chematic_sizes, rdkit_sizes))
        status = "MATCH" if not bond_diffs else f"{len(bond_diffs)} bond(s) differ"
        print(f"  {name}: {status}")
        for key, c_sizes, r_sizes in bond_diffs:
            print(f"    bond {key}: chematic={c_sizes} rdkit={r_sizes}")
        results.append(
            {
                "name": name,
                "smiles": smiles,
                "atom_order_verified": True,
                "bond_diffs": [
                    {"bond": list(k), "chematic_ring_sizes": c, "rdkit_ring_sizes": r}
                    for k, c, r in bond_diffs
                ],
            }
        )
    print(
        f"  TOTAL: {matching_bonds}/{total_bonds} ring bonds agree on ring-size set "
        f"across {len(fixtures) - len(atom_order_mismatches)} fixtures with verified atom order "
        f"({len(atom_order_mismatches)} fixture(s) had unverifiable atom order, excluded from the bond count)"
    )
    return {
        "total_bonds": total_bonds,
        "matching_bonds": matching_bonds,
        "atom_order_mismatches": atom_order_mismatches,
        "per_fixture": results,
    }


def torsion_distribution_for_bond(smiles, atom_indices, n_confs=50, **etkdg_flags):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = rdDistGeom.ETKDGv3()
    for k, v in etkdg_flags.items():
        setattr(params, k, v)
    params.randomSeed = 42
    cids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    angles = []
    for cid in cids:
        conf = mol.GetConformer(cid)
        angle = AllChem.GetDihedralDeg(conf, *atom_indices)
        angles.append(angle)
    return angles


def compare_torsion_distributions():
    print("\n--- 2. Torsion distribution (RDKit ETKDG conformer ensemble, empirical) ---")
    cases = [
        {
            "name": "n_methylacetamide_amide_bond",
            "smiles": "CC(=O)NC",
            # atoms after AddHs keep heavy-atom indices 0..4 first (RDKit
            # convention: AddHs appends H atoms after existing heavy atoms,
            # verified by inspection): C0-C1(=O2)-N3-C4. Central bond is C1-N3.
            "atoms": (0, 1, 3, 4),
            "chematic_prediction": "two-well (0 deg / 180 deg): rule "
            "standard:secondary_amide, term (n=1, s=-1, V=100.0) -> minimum at phi=180 deg",
            "flags": {"useExpTorsionAnglePrefs": True},
        },
        {
            "name": "butane_cc_bond",
            "smiles": "CCCC",
            "atoms": (0, 1, 2, 3),
            "chematic_prediction": "no standard-tier rule matches a plain "
            "alkane C-C-C-C bond in this PR's curated subset (reported as "
            "unmatched, not a gauche/anti claim) -- this case is included "
            "to show what RDKit's OWN unprefixed distribution looks like "
            "for the same bond, as a sanity check, not a rule comparison.",
            "flags": {"useExpTorsionAnglePrefs": True},
        },
    ]
    results = []
    for case in cases:
        angles = torsion_distribution_for_bond(case["smiles"], case["atoms"], **case["flags"])
        angles_mod = [a % 360 for a in angles]
        print(f"  {case['name']} ({case['smiles']}):")
        print(f"    chematic prediction: {case['chematic_prediction']}")
        print(
            f"    RDKit empirical distribution (n={len(angles)}): "
            f"min={min(angles):.1f} max={max(angles):.1f} "
            f"mean_abs={statistics.mean(abs(a) for a in angles):.1f}"
        )
        # Coarse bucket into 6 60-degree bins to show multi-modality
        buckets = [0] * 6
        for a in angles_mod:
            buckets[int(a // 60) % 6] += 1
        print(f"    60-deg-bucket histogram [0,60,120,180,240,300): {buckets}")
        results.append(
            {
                "name": case["name"],
                "smiles": case["smiles"],
                "chematic_prediction": case["chematic_prediction"],
                "n_conformers": len(angles),
                "angles_deg": angles,
                "bucket_histogram_60deg": buckets,
            }
        )
    return results


def main():
    if not CHEMATIC_SIDE_PATH.exists():
        print(
            f"Missing {CHEMATIC_SIDE_PATH}. Run "
            "`cargo run --release -p chematic-3d --example torsion_knowledge_v2_gap_check` "
            "from the repo root first.",
            file=sys.stderr,
        )
        sys.exit(1)

    with open(CHEMATIC_SIDE_PATH) as f:
        chematic_side = json.load(f)

    ring_result = compare_ring_classification(chematic_side["molecules"])
    torsion_result = compare_torsion_distributions()

    print(
        "\n--- 3. Known, disclosed limitation (not measured here) ---\n"
        "  Matched-rule-family and macrocycle-1-4-pair-selection differentials against\n"
        "  RDKit's own internal choices are NOT achievable via RDKit's public Python API\n"
        "  in any version checked (rdkit==2026.03.4 in this isolated venv) -- the\n"
        "  ExperimentalTorsionAngle/BoundsMatrixBuilder C++ machinery that performs this\n"
        "  matching has no public accessor. This PR's rule-to-fixture matching is instead\n"
        "  verified via chematic's own SMARTS-parse + unit/integration tests, and its\n"
        "  translation provenance via direct source citation (source manifest), not a\n"
        "  live differential."
    )

    out = {
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "ring_classification": ring_result,
        "torsion_distribution": torsion_result,
        "known_unmeasured": (
            "matched-rule-family and macrocycle-1-4-pair-selection vs. RDKit's own "
            "internal choices: no public API access found"
        ),
    }
    RESULT_PATH.write_text(json.dumps(out, indent=2))
    print(f"\nWrote full results to {RESULT_PATH}")


if __name__ == "__main__":
    main()
