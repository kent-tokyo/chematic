"""Live-RDKit-oracle validation for `mmff94_stbn_equivalence_diagnostic_227.rs`
(issue #227, diagnostic-only, no chematic-ff production code changed).

Reads the diagnostic example's own JSONL output (one row per StretchBend
`routing_bug_candidate`, `crates/chematic-3d/examples/
mmff94_stbn_equivalence_diagnostic_227.rs`) and, for every row, calls a live
RDKit `MMFFMolProperties.GetMMFFStretchBendParams(mol, i, j, k)` on the SAME
(molecule, atom-triple) -- not a re-derivation from source, a direct call
into the real library (`rdkit==2026.3.3`, matching the pinned commit
`e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f` release tag used throughout
`scripts/mmff94_provenance/`).

Reports two things per candidate:
  1. Does the diagnostic's own derived `stretch_bend_type`/parameter match
     the oracle exactly?
  2. Separately: does chematic's OWN bond-order/aromaticity perception
     (`bond_type_ij`/`bond_type_jk`, as computed by chematic-ff's already
     production `bond_type_for` fed chematic's own parsed bond orders) agree
     with what RDKit's real (post-sanitization) bond typing gives for the
     same two bonds (`GetMMFFBondStretchParams`)? Any stretch-bend-type
     mismatch on a row where bond typing DISAGREES is downstream of an
     upstream aromaticity-perception gap between the two engines, not a bug
     in this diagnostic's RDKit-formula port -- this split is what lets the
     two be told apart instead of conflated into one "mismatch" number.

Usage:
    .venv/bin/python scripts/mmff94_stbn_oracle_validate_227.py \\
        validation/results/mmff94_stbn_equivalence_diagnostic_227.jsonl
"""

import json
import sys

from rdkit import Chem
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties

rows = [json.loads(line) for line in open(sys.argv[1]) if line.strip()]
print(f"total candidate rows: {len(rows)}", file=sys.stderr)

n_index_mismatch = 0
n_props_invalid = 0
n_sbt_match = 0
n_sbt_mismatch = 0
n_outcome_match = 0
n_outcome_mismatch = 0
n_bond_typing_agrees = 0
n_bond_typing_disagrees = 0
n_sbt_match_clean = 0  # bond typing agreed AND sbt matched
n_sbt_mismatch_clean = 0  # bond typing agreed BUT sbt still mismatched (would be a real formula-port bug)
mismatch_examples = []
clean_mismatch_examples = []

mol_cache = {}

for row in rows:
    smiles = row["smiles"]
    a_idx, b_idx, c_idx = row["atoms"]
    za, zb, zc = row["atomic_numbers"]

    if smiles not in mol_cache:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            mol_cache[smiles] = None
        else:
            try:
                props = MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
            except Exception:
                props = None
            mol_cache[smiles] = (mol, props)
    cached = mol_cache[smiles]
    if cached is None or cached[1] is None:
        n_props_invalid += 1
        continue
    mol, props = cached

    # Sanity: atomic numbers at the same index must match (rules out any
    # atom-ordering mismatch between chematic_smiles::parse and RDKit's
    # MolFromSmiles before trusting the comparison at all).
    try:
        rd_z = [
            mol.GetAtomWithIdx(a_idx).GetAtomicNum(),
            mol.GetAtomWithIdx(b_idx).GetAtomicNum(),
            mol.GetAtomWithIdx(c_idx).GetAtomicNum(),
        ]
    except RuntimeError:
        n_index_mismatch += 1
        continue
    if rd_z != [za, zb, zc]:
        n_index_mismatch += 1
        continue

    oracle = props.GetMMFFStretchBendParams(mol, a_idx, b_idx, c_idx)

    oracle_bt_ij = props.GetMMFFBondStretchParams(mol, a_idx, b_idx)
    oracle_bt_jk = props.GetMMFFBondStretchParams(mol, b_idx, c_idx)
    bond_typing_agrees = (
        oracle_bt_ij is not None
        and oracle_bt_jk is not None
        and oracle_bt_ij[0] == row["bond_type_ij"]
        and oracle_bt_jk[0] == row["bond_type_jk"]
    )
    if bond_typing_agrees:
        n_bond_typing_agrees += 1
    else:
        n_bond_typing_disagrees += 1

    predicted_outcome = row["selected_parameter_kind"]
    predicted_sbt = row["rdkit_classification"]

    if oracle is None:
        oracle_sbt = None
    else:
        oracle_sbt, kba_ijk, kba_kji = oracle

    if oracle_sbt is not None:
        if oracle_sbt == predicted_sbt:
            n_sbt_match += 1
            if bond_typing_agrees:
                n_sbt_match_clean += 1
        else:
            n_sbt_mismatch += 1
            example = {
                "kind": "sbt",
                "molecule_id": row["molecule_id"],
                "atoms": row["atoms"],
                "predicted_sbt": predicted_sbt,
                "oracle_sbt": oracle_sbt,
                "predicted_outcome": predicted_outcome,
                "oracle": oracle,
                "bond_typing_agrees": bond_typing_agrees,
                "oracle_bt": [
                    oracle_bt_ij[0] if oracle_bt_ij else None,
                    oracle_bt_jk[0] if oracle_bt_jk else None,
                ],
                "chematic_bt": [row["bond_type_ij"], row["bond_type_jk"]],
            }
            if len(mismatch_examples) < 15:
                mismatch_examples.append(example)
            if bond_typing_agrees:
                n_sbt_mismatch_clean += 1
                if len(clean_mismatch_examples) < 15:
                    clean_mismatch_examples.append(example)

    # outcome comparison: "found_but_zero_dropped"/"unresolved" both predict
    # oracle is None; "exact"/"dfsb_fallback" both predict oracle is found
    # (RDKit's Python binding doesn't distinguish exact-vs-Dfsb in its
    # return value, so "found vs none" + sbt/value is the finest check
    # available through this API).
    predicted_none = predicted_outcome in ("found_but_zero_dropped", "unresolved")
    outcome_ok = (oracle is None) == predicted_none
    if outcome_ok:
        n_outcome_match += 1
    else:
        n_outcome_mismatch += 1
        if len(mismatch_examples) < 15:
            mismatch_examples.append(
                {
                    "kind": "outcome",
                    "molecule_id": row["molecule_id"],
                    "atoms": row["atoms"],
                    "predicted_outcome": predicted_outcome,
                    "oracle": oracle,
                }
            )

    if oracle is not None and predicted_outcome in ("exact", "dfsb_fallback"):
        pred_val = row["selected_parameter_value"]
        if pred_val is not None:
            ok = abs(pred_val[0] - kba_ijk) < 1e-6 and abs(pred_val[1] - kba_kji) < 1e-6
            if not ok and len(mismatch_examples) < 15:
                mismatch_examples.append(
                    {
                        "kind": "value",
                        "molecule_id": row["molecule_id"],
                        "atoms": row["atoms"],
                        "predicted_value": pred_val,
                        "oracle_value": [kba_ijk, kba_kji],
                    }
                )

print(f"n_index_mismatch (atom reindexing between parsers, excluded) = {n_index_mismatch}")
print(f"n_props_invalid (RDKit couldn't build MMFF props at all, excluded) = {n_props_invalid}")
n_compared = len(rows) - n_index_mismatch - n_props_invalid
print(f"n_compared = {n_compared}")
print(f"n_sbt_match = {n_sbt_match}, n_sbt_mismatch = {n_sbt_mismatch}")
print(f"n_outcome_match = {n_outcome_match}, n_outcome_mismatch = {n_outcome_mismatch}")
print(
    "n_bond_typing_agrees (chematic bond_type_ij/jk == oracle's real bond type, "
    f"i.e. upstream aromaticity/bond-order perception matches) = {n_bond_typing_agrees}"
)
print(f"n_bond_typing_disagrees = {n_bond_typing_disagrees}")
print(f"n_sbt_match_clean (bond typing agreed AND sbt matched) = {n_sbt_match_clean}")
print(
    "n_sbt_mismatch_clean (bond typing agreed but sbt STILL mismatched -- "
    f"would be a real formula-port bug, not an upstream confound) = {n_sbt_mismatch_clean}"
)
print()
print("clean mismatch examples (bond typing agreed, sbt still wrong; up to 15):")
for ex in clean_mismatch_examples:
    print(" ", ex)
print()
print("mismatch examples (up to 15):")
for ex in mismatch_examples:
    print(" ", ex)
