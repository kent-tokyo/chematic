"""Live-RDKit-oracle validation for `mmff94_torsion_equivalence_diagnostic_227.rs`
(issue #227, diagnostic-only, no chematic-ff production code changed).

Reads the diagnostic example's own JSONL output (one row per Torsion
`routing_bug_candidate`, `crates/chematic-3d/examples/
mmff94_torsion_equivalence_diagnostic_227.rs`) and, for every row, calls a
live RDKit `MMFFMolProperties.GetMMFFTorsionParams(mol, i, j, k, l)` on the
SAME (molecule, atom-quadruple) -- not a re-derivation from source, a direct
call into the real library (`rdkit==2026.3.3`, matching the pinned commit
`e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f` release tag used throughout
`scripts/mmff94_provenance/`).

`GetMMFFTorsionParams` returns `(torsionType, V1, V2, V3)` or `None`. Two
API quirks, both handled explicitly below (verified from RDKit's own C++
source, `AtomTyper.cpp:3650`, not guessed):
  - If the FINAL resolved table torsion type is 0 (the type-agnostic generic
    code), RDKit's own `torsionType = (mmffTorPair.first ? mmffTorPair.first
    : torTypePair.first)` reports the *original pre-fallback* classification
    instead of 0 (C++ `0 ? ... : ...` truthiness quirk) -- so a reported
    `torsionType` of, e.g., 1 does not by itself prove the match was NOT a
    type-0 generic fallback. This script therefore treats `(V1,V2,V3)` as
    the primary ground truth and `torsionType` as a secondary, looser check.
  - A table hit whose V1=V2=V3=0.0 is treated by RDKit as "no term"
    (`res = !(isDoubleZero(V1)&&isDoubleZero(V2)&&isDoubleZero(V3))`) and
    also returns `None`, indistinguishable at the API level from "no row
    found at all, empirical rule also gave zero".

Writes the final, oracle-enriched JSONL (adds `oracle_torsion_type`,
`oracle_value`, `oracle_kind`, and the fully-resolved
`used_exact`/`used_equivalence`/`used_empirical`/`used_unresolved` booleans
requested by the task spec) to `<input>_oracle_enriched.jsonl`, and prints
the headline resolution counts to stdout.

Usage:
    .venv/bin/python scripts/mmff94_torsion_oracle_validate_227.py \\
        validation/results/mmff94_torsion_equivalence_diagnostic_227.jsonl
"""

import json
import sys

from rdkit import Chem
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties

in_path = sys.argv[1]
out_path = in_path.replace(".jsonl", "_oracle_enriched.jsonl")
rows = [json.loads(line) for line in open(in_path) if line.strip()]
print(f"total candidate rows: {len(rows)}", file=sys.stderr)

n_index_mismatch = 0
n_props_invalid = 0
n_value_match = 0
n_value_mismatch = 0
n_oracle_none = 0
n_oracle_present = 0

# Headline resolution buckets (final answer to "how many of the 1,107
# resolve correctly with a proper eqLevel ladder vs RDKit's real behavior").
oracle_kind_counts = {}
mismatch_examples = []

mol_cache = {}
enriched_rows = []

for row in rows:
    smiles = row["smiles"]
    i_idx, j_idx, k_idx, l_idx = row["atoms"]
    zi, zj, zk, zl = row["atomic_numbers"]

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
            mol.GetAtomWithIdx(i_idx).GetAtomicNum(),
            mol.GetAtomWithIdx(j_idx).GetAtomicNum(),
            mol.GetAtomWithIdx(k_idx).GetAtomicNum(),
            mol.GetAtomWithIdx(l_idx).GetAtomicNum(),
        ]
    except RuntimeError:
        n_index_mismatch += 1
        continue
    if rd_z != [zi, zj, zk, zl]:
        n_index_mismatch += 1
        continue

    oracle = props.GetMMFFTorsionParams(mol, i_idx, j_idx, k_idx, l_idx)

    predicted_kind = row["selected_parameter_kind"]
    predicted_value = row["selected_parameter_value"]

    if oracle is None:
        n_oracle_none += 1
        oracle_torsion_type = None
        oracle_value = None
    else:
        n_oracle_present += 1
        oracle_torsion_type, ov1, ov2, ov3 = oracle
        oracle_value = [ov1, ov2, ov3]

    # Determine the FINAL oracle-informed kind for this row, per the
    # methodology documented in the Rust file's own doc comment:
    #   - predicted_kind != "table_unresolved": this file's own table+ladder
    #     port claims a specific row. Confirmed if oracle_value matches.
    #   - predicted_kind == "table_unresolved": this file's port found
    #     nothing in the table at all. If the oracle STILL returns a nonzero
    #     value, that is direct evidence RDKit's empirical rule produced it
    #     ("empirical_rule"). If the oracle also returns None, genuinely
    #     "unresolved" even under RDKit's complete algorithm.
    is_found_but_zero_dropped = False
    if predicted_kind == "table_unresolved":
        if oracle is not None:
            oracle_kind = "empirical_rule"
        else:
            oracle_kind = "unresolved"
        value_ok = None  # not applicable -- no self-predicted value to check
    elif (
        oracle is None
        and predicted_value is not None
        and all(abs(v) < 1e-9 for v in predicted_value)
    ):
        # RDKit's real getMMFFTorsionParams (AtomTyper.cpp:3627-3668) does
        # NOT skip an all-zero table row and keep searching the ladder for a
        # nonzero one -- it accepts the first row the ladder finds (zero or
        # not) and only THEN applies one final isDoubleZero(V1)&&(V2)&&(V3)
        # gate over whatever was found, dropping to `None` if it's all zero
        # (bypassing the empirical rule entirely, since that branch only
        # triggers when the ladder found NO row at all, not when it found a
        # zero one). Our own ladder found the identical zero row (same
        # mechanism, verified: this is what "our predicted_value is exactly
        # [0,0,0]" means), so this is RDKit's real, documented behavior, not
        # a port bug -- same "found_but_zero_dropped" pattern the sibling
        # stretch-bend diagnostic (PR #273) reported. Net physical effect
        # (no torsion energy contribution) matches "unresolved" for the
        # required schema fields; `is_found_but_zero_dropped` records the
        # distinction separately.
        oracle_kind = "unresolved"
        value_ok = True
        n_value_match += 1
        is_found_but_zero_dropped = True
    else:
        is_found_but_zero_dropped = False
        oracle_kind = predicted_kind
        if oracle_value is not None and predicted_value is not None:
            value_ok = (
                abs(predicted_value[0] - oracle_value[0]) < 1e-6
                and abs(predicted_value[1] - oracle_value[1]) < 1e-6
                and abs(predicted_value[2] - oracle_value[2]) < 1e-6
            )
        else:
            # Our port claims a table hit but the oracle returned None (or
            # vice versa) -- a real port bug, not just a numeric rounding
            # difference.
            value_ok = False
        if value_ok:
            n_value_match += 1
        else:
            n_value_mismatch += 1
            if len(mismatch_examples) < 20:
                mismatch_examples.append(
                    {
                        "molecule_id": row["molecule_id"],
                        "atoms": row["atoms"],
                        "predicted_kind": predicted_kind,
                        "predicted_value": predicted_value,
                        "oracle_torsion_type": oracle_torsion_type,
                        "oracle_value": oracle_value,
                    }
                )

    oracle_kind_counts[oracle_kind] = oracle_kind_counts.get(oracle_kind, 0) + 1

    used_exact = oracle_kind == "exact"
    used_equivalence = oracle_kind.startswith("equivalence_level_")
    used_empirical = oracle_kind == "empirical_rule"
    used_unresolved = oracle_kind == "unresolved"

    enriched = dict(row)
    enriched["oracle_torsion_type"] = oracle_torsion_type
    enriched["oracle_value"] = oracle_value
    enriched["oracle_kind"] = oracle_kind
    enriched["oracle_value_matches_prediction"] = value_ok
    enriched["is_found_but_zero_dropped"] = is_found_but_zero_dropped
    enriched["used_exact"] = used_exact
    enriched["used_equivalence"] = used_equivalence
    enriched["used_empirical"] = used_empirical
    enriched["used_unresolved"] = used_unresolved
    enriched_rows.append(enriched)

with open(out_path, "w") as f:
    for r in enriched_rows:
        f.write(json.dumps(r) + "\n")

n_compared = len(rows) - n_index_mismatch - n_props_invalid
print(f"n_index_mismatch (atom reindexing between parsers, excluded) = {n_index_mismatch}")
print(f"n_props_invalid (RDKit couldn't build MMFF props at all, excluded) = {n_props_invalid}")
print(f"n_compared = {n_compared}")
print(f"n_oracle_none (oracle found no torsion term at all) = {n_oracle_none}")
print(f"n_oracle_present (oracle found a real, nonzero torsion term) = {n_oracle_present}")
print()
print("=== HEADLINE: final oracle-informed selected_parameter_kind breakdown ===")
for k, v in sorted(oracle_kind_counts.items(), key=lambda kv: -kv[1]):
    print(f"  {k} = {v}  ({100 * v / n_compared:.1f}%)")
resolved = sum(v for k, v in oracle_kind_counts.items() if k != "unresolved")
print(f"resolved by SOME mechanism (exact/equivalence/empirical) = {resolved}/{n_compared} ({100 * resolved / n_compared:.1f}%)")
ladder_resolved = sum(
    v for k, v in oracle_kind_counts.items() if k == "exact" or k.startswith("equivalence_level_")
)
print(f"resolved specifically by the table+eqLevel-ladder (exact/equivalence_level_N, i.e. what a ladder fix alone would close) = {ladder_resolved}/{n_compared} ({100 * ladder_resolved / n_compared:.1f}%)")
print()
print(f"n_value_match (of non-table_unresolved rows, our self-predicted value matches oracle exactly, incl. found_but_zero_dropped) = {n_value_match}")
print(f"n_value_mismatch (real, unexplained port discrepancies) = {n_value_mismatch}")
n_found_but_zero_dropped = sum(1 for r in enriched_rows if r["is_found_but_zero_dropped"])
print(f"n_found_but_zero_dropped (our ladder found the identical row RDKit's real algorithm finds, but it's an explicit all-zero row RDKit's own isDoubleZero gate drops to None -- expected RDKit behavior, not a port bug, folded into oracle_kind='unresolved' above) = {n_found_but_zero_dropped}")
print()
print("value mismatch examples (up to 20):")
for ex in mismatch_examples:
    print(" ", ex)
print()
print(f"enriched JSONL written to {out_path}")
