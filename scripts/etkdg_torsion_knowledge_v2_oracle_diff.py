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

3. MATCHED RULE FAMILY / CENTRAL BOND SELECTION (comparable -- corrected
   after independent review found the public accessor an earlier draft of
   this script claimed didn't exist): `rdkit.Chem.rdDistGeom.
   GetExperimentalTorsions(mol, useExpTorsionAnglePrefs, useSmallRingTorsions,
   useMacrocycleTorsions, useBasicKnowledge, ETversion)` returns, per bond,
   the matched SMARTS, coefficients, and the exact atom quadruple RDKit
   bound -- reads
   `validation/etkdg_torsion_knowledge_v2_chematic_torsions.json` (written by
   the same gap-check example run) and compares: which bonds each engine
   assigns ANY torsion to (central-bond-selection agreement), and, for
   shared bonds, whether the SAME atom quadruple was chosen (a strong,
   exact, non-fuzzy proxy for rule-family agreement -- two engines landing
   on the identical 4 atoms essentially always means they identified the
   same real substructure, without needing a subjective SMARTS-text-to-
   SMARTS-text semantic classifier).

4. MACROCYCLE 1-4 PAIR SELECTION (comparable): `rdDistGeom.
   GetMoleculeBoundsMatrix(mol, set15bounds, scaleVDW, doTriangleSmoothing,
   useMacrocycle14config)` takes the macrocycle-1-4 flag directly --
   diffing the bounds matrix with/without it (same molecule, everything else
   held fixed) directly answers "which pairs did RDKit's own macrocycle-1-4
   logic touch, and by how much", compared against this PR's own
   `macrocycle_14_bound_adjustments`' proposed pairs for the same molecule.
"""

import json
import statistics
import sys
from pathlib import Path

try:
    import numpy as np
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
CHEMATIC_TORSIONS_PATH = (
    REPO_ROOT / "validation" / "etkdg_torsion_knowledge_v2_chematic_torsions.json"
)
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


def compare_rule_family_and_central_bond(fixtures):
    """Compares which bonds get ANY torsion assigned, and (for shared
    bonds) whether the exact atom quadruple matches, between chematic's
    `full_config()` (use_exp_torsions + use_small_ring_torsions +
    use_macrocycle_torsions, matching `useBasicKnowledge=True` implicitly
    via chematic's tier 5 being gated on the same flag as tier 4 -- see
    matcher.rs) and RDKit's `GetExperimentalTorsions` with the same 4 flags.
    """
    print(
        "\n--- 3. Rule-family / central-bond-selection "
        "(RDKit GetExperimentalTorsions vs. chematic) ---"
    )
    total_chematic_bonds = 0
    total_rdkit_bonds = 0
    shared_bonds = 0
    quadruple_matches = 0
    quadruple_mismatches = []
    atom_order_mismatches = []
    per_fixture = []

    for fx in fixtures:
        name, smiles = fx["name"], fx["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        rdkit_atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        if rdkit_atoms != fx["atoms"]:
            atom_order_mismatches.append(name)
            continue  # cannot compare per-bond without verified atom correspondence

        rdkit_result = rdDistGeom.GetExperimentalTorsions(mol, True, True, True, True, 2, False)
        rdkit_bonds = {}
        for item in rdkit_result:
            b = mol.GetBondWithIdx(item["bondIndex"])
            key = tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
            rdkit_bonds[key] = (item["smarts"], tuple(item["atomIndices"]))

        chematic_bonds = {}
        for tb in fx["torsion_bonds"]:
            key = tuple(sorted((tb["a"], tb["b"])))
            chematic_bonds[key] = (tb["rule_id"], tuple(tb["atoms"]))

        total_chematic_bonds += len(chematic_bonds)
        total_rdkit_bonds += len(rdkit_bonds)
        shared = set(chematic_bonds) & set(rdkit_bonds)
        shared_bonds += len(shared)

        fixture_mismatches = []
        for key in shared:
            c_rule, c_atoms = chematic_bonds[key]
            r_smarts, r_atoms = rdkit_bonds[key]
            # A quadruple match counts whether it's read forwards (A,B,C,D)
            # or backwards (D,C,B,A) -- both describe the identical
            # dihedral (see matcher.rs's own doc note on this symmetry).
            if c_atoms == r_atoms or c_atoms == tuple(reversed(r_atoms)):
                quadruple_matches += 1
            else:
                fixture_mismatches.append(
                    {
                        "bond": list(key),
                        "chematic_rule": c_rule,
                        "chematic_atoms": list(c_atoms),
                        "rdkit_smarts": r_smarts,
                        "rdkit_atoms": list(r_atoms),
                    }
                )
        if fixture_mismatches:
            quadruple_mismatches.extend(
                {"name": name, **m} for m in fixture_mismatches
            )
        per_fixture.append(
            {
                "name": name,
                "chematic_bond_count": len(chematic_bonds),
                "rdkit_bond_count": len(rdkit_bonds),
                "shared_bond_count": len(shared),
            }
        )

    print(
        f"  chematic assigns a torsion to {total_chematic_bonds} bonds, "
        f"RDKit to {total_rdkit_bonds} bonds, across {len(fixtures) - len(atom_order_mismatches)} "
        f"fixtures with verified atom order ({len(atom_order_mismatches)} unverifiable, excluded)"
    )
    print(f"  {shared_bonds} bonds get a torsion from BOTH engines")
    if shared_bonds:
        print(
            f"  of those, {quadruple_matches}/{shared_bonds} "
            f"({100 * quadruple_matches / shared_bonds:.1f}%) chose the IDENTICAL atom quadruple"
        )
    for m in quadruple_mismatches:
        print(
            f"    MISMATCH {m['name']}: bond {m['bond']} chematic={m['chematic_rule']}"
            f"{m['chematic_atoms']} vs rdkit={m['rdkit_smarts']}{m['rdkit_atoms']}"
        )

    return {
        "total_chematic_bonds": total_chematic_bonds,
        "total_rdkit_bonds": total_rdkit_bonds,
        "shared_bonds": shared_bonds,
        "quadruple_matches": quadruple_matches,
        "quadruple_mismatches": quadruple_mismatches,
        "atom_order_mismatches": atom_order_mismatches,
        "per_fixture": per_fixture,
    }


def compare_macrocycle_14_pairs(fixtures):
    """Diffs RDKit's own bounds matrix with/without `useMacrocycle14config`
    (same molecule, everything else held fixed) to find which atom pairs
    RDKit's own macrocycle-1-4 logic actually touches, compared against this
    PR's `macrocycle_14_bound_adjustments` proposed pairs for the same
    fixtures.
    """
    print("\n--- 4. Macrocycle 1-4 pair selection (RDKit GetMoleculeBoundsMatrix diff) ---")
    results = []
    for fx in fixtures:
        chematic_pairs = fx.get("macrocycle_14_pairs", [])
        if not chematic_pairs:
            continue
        name, smiles = fx["name"], fx["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        n = mol.GetNumAtoms()
        bm_off = np.array(rdDistGeom.GetMoleculeBoundsMatrix(mol, True, False, True, False))
        bm_on = np.array(rdDistGeom.GetMoleculeBoundsMatrix(mol, True, False, True, True))
        rdkit_changed_pairs = set()
        for i in range(n):
            for j in range(i + 1, n):
                lo_off, hi_off = bm_off[j, i], bm_off[i, j]
                lo_on, hi_on = bm_on[j, i], bm_on[i, j]
                if abs(lo_off - lo_on) > 1e-3 or abs(hi_off - hi_on) > 1e-3:
                    rdkit_changed_pairs.add((i, j))

        chematic_pair_set = {tuple(sorted((p["a"], p["b"]))) for p in chematic_pairs}
        overlap = chematic_pair_set & rdkit_changed_pairs
        print(
            f"  {name}: chematic proposes {len(chematic_pair_set)} pairs, "
            f"RDKit's useMacrocycle14config flag changes {len(rdkit_changed_pairs)} pairs, "
            f"{len(overlap)} in common"
        )
        results.append(
            {
                "name": name,
                "chematic_pairs": sorted(list(p) for p in chematic_pair_set),
                "rdkit_changed_pairs": sorted(list(p) for p in rdkit_changed_pairs),
                "overlap": sorted(list(p) for p in overlap),
            }
        )
    return results


def main():
    if not CHEMATIC_SIDE_PATH.exists() or not CHEMATIC_TORSIONS_PATH.exists():
        print(
            f"Missing {CHEMATIC_SIDE_PATH} and/or {CHEMATIC_TORSIONS_PATH}. Run "
            "`cargo run --release -p chematic-3d --example torsion_knowledge_v2_gap_check` "
            "from the repo root first.",
            file=sys.stderr,
        )
        sys.exit(1)

    with open(CHEMATIC_SIDE_PATH) as f:
        chematic_side = json.load(f)
    with open(CHEMATIC_TORSIONS_PATH) as f:
        chematic_torsions = json.load(f)

    ring_result = compare_ring_classification(chematic_side["molecules"])
    torsion_result = compare_torsion_distributions()
    family_result = compare_rule_family_and_central_bond(chematic_torsions["molecules"])
    pair_result = compare_macrocycle_14_pairs(chematic_torsions["molecules"])

    out = {
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "ring_classification": ring_result,
        "torsion_distribution": torsion_result,
        "rule_family_and_central_bond": family_result,
        "macrocycle_14_pairs": pair_result,
    }
    RESULT_PATH.write_text(json.dumps(out, indent=2))
    print(f"\nWrote full results to {RESULT_PATH}")


if __name__ == "__main__":
    main()
