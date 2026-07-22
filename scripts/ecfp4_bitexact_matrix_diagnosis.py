#!/usr/bin/env python3
"""diag/ecfp4-bitexact-api: full parameter-matrix bit-exactness diagnosis for
`rdkit_morgan_ecfp4_experimental` (crates/chematic-fp/src/rdkit_morgan_ecfp4.rs)
against a live RDKit oracle. Diagnosis only -- no production code touched.

Companion Rust dump: `crates/chematic-fp/examples/rdkit_ecfp4_bitexact_matrix_dump.rs`
Fixture corpus (single source of truth for both sides): `scripts/ecfp4_bitexact_matrix_fixtures.csv`
RFC write-up: `docs/ecfp4_bitexact_api_rfc.md`

The production API exposes exactly ONE point in the full RDKit Morgan
parameter space: radius=2, fpSize=2048, includeRedundantEnvironments=false,
useChirality=false, useBondTypes=true, default (GetConnectivityInvariants)
atom invariant. Every matrix cell below is classified into one of:

  verified_bit_exact                 -- driven end-to-end against the oracle
                                         at the production API's one real
                                         config, matches.
  verified_reachable_via_internal_math -- the port's underlying hash math
                                         (expand_one_pass/rdkit_morgan_raw_trace,
                                         diagnostics feature) generalizes to
                                         this cell (e.g. radius != 2) and
                                         matches the oracle, even though no
                                         production entry point exposes it.
  verified_via_postfold_of_public_data -- derived by folding the production
                                         API's already-public `sparse_counts`
                                         field in Python (no source touched)
                                         and matching the oracle's real fold.
  architecturally_unimplemented      -- no code path in
                                         rdkit_morgan_hash.rs/rdkit_morgan_ecfp4.rs
                                         can express this option at all (confirmed
                                         by reading connectivity_invariant/
                                         checked_bond_invariant); the oracle
                                         comparison below only demonstrates that
                                         real RDKit *does* distinguish the option
                                         (so the gap is consequential, not just
                                         theoretical).
  mismatch_aromaticity_kekulize_hardfail -- traced to the diag/aromaticity-
                                         rdkit-parity PR's already-identified
                                         kekulize() hard-fail classes.
  mismatch_unclassified               -- anything else. Must be empty for
                                         this diagnosis to be considered done.

Usage:
    .venv/bin/python scripts/ecfp4_bitexact_matrix_diagnosis.py \\
        --fixtures scripts/ecfp4_bitexact_matrix_fixtures.csv \\
        --chematic-dump validation/results/ecfp4_bitexact_matrix_dump.jsonl \\
        --summary-out validation/results/ecfp4_bitexact_matrix_summary.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

NBITS_SWEEP = [128, 256, 512, 1024, 2048]
RADIUS_SWEEP = [0, 1, 2, 3]
PRODUCTION_RADIUS = 2
PRODUCTION_NBITS = 2048

_PARSE_PARAMS = Chem.SmilesParserParams()
_PARSE_PARAMS.removeHs = False  # keep isotope-labeled/explicit atoms as real graph atoms


def rd_parse(smi):
    return Chem.MolFromSmiles(smi, _PARSE_PARAMS)


def make_gen(radius, fp_size=2048, include_chirality=False, use_bond_types=True,
             atom_invariants_generator=None, include_redundant=False):
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=radius,
        fpSize=fp_size,
        includeChirality=include_chirality,
        useBondTypes=use_bond_types,
        atomInvariantsGenerator=atom_invariants_generator,
        includeRedundantEnvironments=include_redundant,
    )


def rd_default_pairs(mol, max_radius):
    """(atom, radius) -> raw_id, RDKit's real includeRedundantEnvironments=False
    lifecycle, unfolded (sparse) -- ground truth for the radius axis."""
    gen = make_gen(radius=max_radius, include_redundant=False)
    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateBitInfoMap()
    gen.GetSparseFingerprint(mol, additionalOutput=ao)
    pairs = {}
    for raw_id, envs in ao.GetBitInfoMap().items():
        rid = raw_id & 0xFFFFFFFF
        for atom_idx, radius in envs:
            pairs[(atom_idx, radius)] = rid
    return pairs


def rd_sparse_counts(mol, radius):
    gen = make_gen(radius=radius, include_redundant=False)
    fp = gen.GetSparseCountFingerprint(mol)
    return {rid & 0xFFFFFFFF: c for rid, c in fp.GetNonzeroElements().items()}


def rd_folded_on_bits(mol, radius, fp_size):
    gen = make_gen(radius=radius, fp_size=fp_size, include_redundant=False)
    return set(gen.GetFingerprint(mol).GetOnBits())


def rd_folded_counts(mol, radius, fp_size):
    gen = make_gen(radius=radius, fp_size=fp_size, include_redundant=False)
    fp = gen.GetCountFingerprint(mol)
    return dict(fp.GetNonzeroElements())


def load_fixtures(path):
    fixtures = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("id|tags|smiles"):
                continue
            fid, tags, smi = line.split("|", 2)
            fixtures.append({"id": fid, "tags": tags.split(","), "smiles": smi})
    return fixtures


def load_chematic_dump(path):
    rows = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            rows[r["fixture_id"]] = r
    return rows


def classify_fixture_tags(tags):
    buckets = set()
    for t in tags:
        if t.startswith("charged_kekulize_fail"):
            buckets.add("charged_kekulize_fail")
        elif t.startswith("charged_ok"):
            buckets.add("charged_ok")
        elif t == "isotope":
            buckets.add("isotope")
        elif t == "disconnected":
            buckets.add("disconnected")
        elif t.startswith("aromatic_form") or t.startswith("kekule_form"):
            buckets.add("aromatic_vs_kekule")
        elif t.startswith("stereo_tetrahedral_pair"):
            buckets.add("stereo_tetrahedral")
        elif t.startswith("stereo_ez_pair"):
            buckets.add("stereo_ez")
        elif t == "baseline":
            buckets.add("baseline")
    return buckets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixtures", required=True)
    ap.add_argument("--chematic-dump", required=True)
    ap.add_argument("--summary-out", required=True)
    args = ap.parse_args()

    fixtures = load_fixtures(args.fixtures)
    chem_rows = load_chematic_dump(args.chematic_dump)

    if set(f["id"] for f in fixtures) != set(chem_rows.keys()):
        print("PIPELINE ERROR: fixture id set mismatch between CSV and chematic dump",
              file=sys.stderr)
        sys.exit(1)

    per_fixture = []
    molecule_shape_bucket = defaultdict(lambda: {"total": 0, "verified_bit_exact": 0,
                                                   "mismatch_aromaticity_kekulize_hardfail": 0,
                                                   "mismatch_unclassified": 0})
    radius_bucket = {r: {"total": 0, "match": 0, "n_a_no_aromaticity": 0} for r in RADIUS_SWEEP}
    nbits_binary_bucket = {n: {"total": 0, "match": 0} for n in NBITS_SWEEP}
    nbits_count_bucket = {n: {"total": 0, "match": 0} for n in NBITS_SWEEP}

    unclassified = []

    for fx in fixtures:
        fid, smi, tags = fx["id"], fx["smiles"], fx["tags"]
        chem = chem_rows[fid]
        rd_mol = rd_parse(smi)
        shape_buckets = classify_fixture_tags(tags)

        row = {"fixture_id": fid, "smiles": smi, "tags": tags,
               "chematic_status": chem["status"]}

        if rd_mol is None:
            row["oracle_parse_failed"] = True
            unclassified.append(row)
            per_fixture.append(row)
            continue

        # ---- one-cell production config: r=2, 2048 bits, binary ----
        is_kekulize_fail_fixture = "aromaticity_rfc_kekulize_hardfail" in tags
        if chem["status"] == "success":
            rd_pairs_r2 = rd_default_pairs(rd_mol, PRODUCTION_RADIUS)
            chem_pairs_r2 = {(a, r): rid for a, r, rid in chem["radius_default_pairs"]
                             if r <= PRODUCTION_RADIUS}
            rd_counts_r2 = rd_sparse_counts(rd_mol, PRODUCTION_RADIUS)
            chem_counts_r2 = {rid: c for rid, c in chem["sparse_counts"]}
            rd_folded_r2_2048 = rd_folded_on_bits(rd_mol, PRODUCTION_RADIUS, PRODUCTION_NBITS)
            chem_folded = set(chem["folded_on_bits_2048"])

            one_cell_match = (
                chem_pairs_r2 == rd_pairs_r2
                and chem_counts_r2 == rd_counts_r2
                and chem_folded == rd_folded_r2_2048
            )
            row["one_cell_config_match"] = one_cell_match
            if one_cell_match:
                cell_status = "verified_bit_exact"
            else:
                cell_status = "mismatch_unclassified"
                unclassified.append(row)
        elif chem["status"] == "error_kekulization_failed" and is_kekulize_fail_fixture:
            cell_status = "mismatch_aromaticity_kekulize_hardfail"
            row["one_cell_config_match"] = None
        elif chem["status"].startswith("error_"):
            cell_status = "mismatch_unclassified"
            row["one_cell_config_match"] = None
            unclassified.append(row)
        else:
            cell_status = "mismatch_unclassified"
            unclassified.append(row)

        row["one_cell_cell_status"] = cell_status
        for b in shape_buckets:
            molecule_shape_bucket[b]["total"] += 1
            if cell_status in molecule_shape_bucket[b]:
                molecule_shape_bucket[b][cell_status] += 1

        # ---- radius axis (only meaningful where aromaticity succeeded) ----
        if chem["radius_default_pairs"]:
            for r in RADIUS_SWEEP:
                rd_pairs_r = rd_default_pairs(rd_mol, r)
                rd_pairs_r = {k: v for k, v in rd_pairs_r.items() if k[1] == r}
                chem_pairs_r = {(a, rad): rid for a, rad, rid in chem["radius_default_pairs"]
                                if rad == r}
                radius_bucket[r]["total"] += 1
                if chem_pairs_r == rd_pairs_r:
                    radius_bucket[r]["match"] += 1
                else:
                    # Pre-existing latent bug, fixed here: (atom_idx, radius) tuple
                    # dict keys aren't JSON-serializable. Never exercised before K1
                    # since no fixture reached this branch with a mismatch. Stringify
                    # for the diagnostic sample only -- the match/no-match decision
                    # above (the actual comparator) is untouched.
                    unclassified.append({**row, "axis": "radius", "radius": r,
                                          "chem": {str(k): v for k, v in chem_pairs_r.items()},
                                          "rdkit": {str(k): v for k, v in rd_pairs_r.items()}})
        else:
            for r in RADIUS_SWEEP:
                radius_bucket[r]["n_a_no_aromaticity"] += 1

        # ---- nBits / count-vs-binary axes: fold production's public sparse_counts ----
        if chem["status"] == "success":
            chem_counts_r2 = {rid: c for rid, c in chem["sparse_counts"]}
            for n in NBITS_SWEEP:
                manual_bits = {rid % n for rid in chem_counts_r2}
                rd_bits = rd_folded_on_bits(rd_mol, PRODUCTION_RADIUS, n)
                nbits_binary_bucket[n]["total"] += 1
                if manual_bits == rd_bits:
                    nbits_binary_bucket[n]["match"] += 1
                else:
                    unclassified.append({**row, "axis": "nbits_binary", "nbits": n})

                manual_counts = defaultdict(int)
                for rid, c in chem_counts_r2.items():
                    manual_counts[rid % n] += c
                rd_counts = rd_folded_counts(rd_mol, PRODUCTION_RADIUS, n)
                nbits_count_bucket[n]["total"] += 1
                if dict(manual_counts) == rd_counts:
                    nbits_count_bucket[n]["match"] += 1
                else:
                    unclassified.append({**row, "axis": "nbits_count", "nbits": n})

        per_fixture.append(row)

    # ---- architecturally-unimplemented capability checks (oracle-only, real RDKit) ----
    capability_checks = {}

    alanine_l, alanine_d = rd_parse("C[C@H](N)C(=O)O"), rd_parse("C[C@@H](N)C(=O)O")
    gen_no_chir = make_gen(radius=2, include_chirality=False)
    gen_chir = make_gen(radius=2, include_chirality=True)
    capability_checks["use_chirality"] = {
        "status": "architecturally_unimplemented",
        "reason": ("rdkit_morgan_hash.rs's connectivity_invariant/checked_bond_invariant have "
                   "no chirality byte or branch at all -- useChirality=true cannot be requested."),
        "rdkit_default_useChirality_false_L_eq_D": (
            gen_no_chir.GetFingerprint(alanine_l).ToBitString()
            == gen_no_chir.GetFingerprint(alanine_d).ToBitString()
        ),
        "rdkit_useChirality_true_L_eq_D": (
            gen_chir.GetFingerprint(alanine_l).ToBitString()
            == gen_chir.GetFingerprint(alanine_d).ToBitString()
        ),
    }

    pyridine = rd_parse("c1ccncc1")
    gen_bt_true = make_gen(radius=2, use_bond_types=True)
    gen_bt_false = make_gen(radius=2, use_bond_types=False)
    capability_checks["use_bond_types_false"] = {
        "status": "architecturally_unimplemented",
        "reason": ("checked_bond_invariant only implements the useBondTypes=true branch "
                   "(module docs: 'the chirality-perturbed branch is out of scope' -- likewise "
                   "no useBondTypes=false branch exists)."),
        "rdkit_true_vs_false_differ_same_molecule": (
            gen_bt_true.GetFingerprint(pyridine).ToBitString()
            != gen_bt_false.GetFingerprint(pyridine).ToBitString()
        ),
    }

    phenol = rd_parse("c1ccccc1O")
    gen_default_inv = make_gen(radius=2)
    gen_feature_inv = make_gen(radius=2,
                                atom_invariants_generator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen())
    capability_checks["alternative_atom_invariants_fcfp_style"] = {
        "status": "architecturally_unimplemented",
        "reason": ("connectivity_invariant hardcodes RDKit's default GetConnectivityInvariants "
                   "component set; no FCFP-style feature-invariant path exists in "
                   "rdkit_morgan_hash.rs."),
        "rdkit_default_vs_feature_invariant_differ": (
            gen_default_inv.GetFingerprint(phenol).ToBitString()
            != gen_feature_inv.GetFingerprint(phenol).ToBitString()
        ),
    }

    summary = {
        "schema_version": "1",
        "production_api": "chematic_fp::rdkit_morgan_ecfp4_experimental",
        "production_fixed_config": {
            "radius": PRODUCTION_RADIUS, "nbits": PRODUCTION_NBITS,
            "include_redundant_environments": False, "use_chirality": False,
            "use_bond_types": True, "atom_invariant": "default (GetConnectivityInvariants)",
        },
        "fixture_count": len(fixtures),
        "molecule_shape_matrix": dict(molecule_shape_bucket),
        "radius_axis_matrix": radius_bucket,
        "nbits_binary_axis_matrix": nbits_binary_bucket,
        "nbits_count_axis_matrix": nbits_count_bucket,
        "capability_checks_architecturally_unimplemented": capability_checks,
        "unclassified_count": len(unclassified),
        "unclassified_sample": unclassified[:10],
    }

    with open(args.summary_out, "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    # ---- console headline ----
    print("=== ECFP4 bit-exact API matrix diagnosis ===")
    print(f"fixtures: {len(fixtures)}, unclassified mismatches: {len(unclassified)}")
    print("\n-- molecule-shape axis (one-cell production config) --")
    for bucket, counts in sorted(molecule_shape_bucket.items()):
        print(f"  {bucket:28s} total={counts['total']:3d} "
              f"verified={counts['verified_bit_exact']:3d} "
              f"kekulize_hardfail={counts['mismatch_aromaticity_kekulize_hardfail']:3d} "
              f"unclassified={counts['mismatch_unclassified']:3d}")
    print("\n-- radius axis (internal hash math, default lifecycle) --")
    for r, c in radius_bucket.items():
        print(f"  radius={r}  total={c['total']:3d} match={c['match']:3d} "
              f"n/a(no_aromaticity)={c['n_a_no_aromaticity']:3d}")
    print("\n-- nBits axis, binary (post-fold of public sparse_counts) --")
    for n, c in nbits_binary_bucket.items():
        print(f"  nbits={n:5d}  total={c['total']:3d} match={c['match']:3d}")
    print("\n-- nBits axis, count (post-fold of public sparse_counts) --")
    for n, c in nbits_count_bucket.items():
        print(f"  nbits={n:5d}  total={c['total']:3d} match={c['match']:3d}")
    print("\n-- architecturally unimplemented capabilities (oracle-only demonstration) --")
    for name, c in capability_checks.items():
        print(f"  {name}: {json.dumps(c)}")

    if unclassified:
        print(f"\nWARNING: {len(unclassified)} unclassified mismatches -- see "
              f"{args.summary_out}'s unclassified_sample", file=sys.stderr)
        sys.exit(1)

    print(f"\nwrote {args.summary_out}")


if __name__ == "__main__":
    main()
