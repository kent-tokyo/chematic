#!/usr/bin/env python3
"""Aromaticity-A1-0 step 5: join the Rust trace against the RDKit oracle,
and check that the false_positive/false_negative/negative_control corpus
buckets are polarized the way their labels claim.

Reads validation/results/aromaticity_a1_0_trace.jsonl (from
`cargo run -p chematic-perception --release --example aromaticity_a1_0_report`),
adds a `rdkit_atom_aromatic` column per (smiles, atom_idx), writes
validation/results/aromaticity_a1_0_diagnosis.jsonl, and prints a report.

Diagnostic only -- reads two already-frozen inputs, writes one joined
output; does not change any production code or behavior.

Run:
    cargo run -p chematic-perception --release --example aromaticity_a1_0_report \
        -- validation/aromaticity_a1_0_corpus.jsonl \
        > validation/results/aromaticity_a1_0_trace.jsonl
    python3 scripts/aromaticity_a1_0_diagnosis.py
"""
import json
import sys
from pathlib import Path

try:
    from rdkit import Chem as RDChem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    print("RDKit not installed -- this script requires RDKit. Install with:")
    print("  pip install rdkit")
    sys.exit(1)

REPO = Path(__file__).resolve().parent.parent
TRACE_IN = REPO / "validation" / "results" / "aromaticity_a1_0_trace.jsonl"
OUT = REPO / "validation" / "results" / "aromaticity_a1_0_diagnosis.jsonl"


def rdkit_aromatic_atoms(smiles):
    mol = RDChem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {atom.GetIdx(): atom.GetIsAromatic() for atom in mol.GetAtoms()}


def main():
    rows = [json.loads(line) for line in open(TRACE_IN) if line.strip()]
    print(f"loaded {len(rows)} trace rows from {TRACE_IN}")

    rdkit_cache = {}
    unparseable = set()
    for row in rows:
        smi = row["smiles"]
        if smi not in rdkit_cache and smi not in unparseable:
            result = rdkit_aromatic_atoms(smi)
            if result is None:
                unparseable.add(smi)
                print(f"WARNING: RDKit could not parse {smi!r}")
            else:
                rdkit_cache[smi] = result

    joined = []
    for row in rows:
        smi = row["smiles"]
        if smi in unparseable:
            row["rdkit_atom_aromatic"] = None
        else:
            row["rdkit_atom_aromatic"] = rdkit_cache[smi].get(row["atom_idx"])
        joined.append(row)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        for row in joined:
            f.write(json.dumps(row) + "\n")
    print(f"wrote {len(joined)} joined rows to {OUT}")

    # ---- overall agreement ----
    valid = [r for r in joined if r["rdkit_atom_aromatic"] is not None]
    agree = sum(1 for r in valid if r["current_engine_atom_aromatic"] == r["rdkit_atom_aromatic"])
    print(f"\ncurrent_engine vs rdkit per-atom-row agreement: {agree}/{len(valid)} ({100*agree/len(valid):.2f}%)")

    # ---- Aromaticity-A1-1a: exhaustive oracle vs rdkit, and vs current_engine ----
    # The oracle is a discovery tool (candidates built from the SAME per-atom
    # rules that are wrong for the false-positive family), NOT a corrected
    # engine -- do not read "oracle agrees with RDKit more often" as "the
    # oracle is more correct" in general; it's only informative case-by-case.
    oracle_agree = sum(1 for r in valid if r["oracle_atom_aromatic"] == r["rdkit_atom_aromatic"])
    print(f"oracle vs rdkit per-atom-row agreement:        {oracle_agree}/{len(valid)} ({100*oracle_agree/len(valid):.2f}%)")
    oracle_vs_engine = sum(1 for r in valid if r["oracle_atom_aromatic"] == r["current_engine_atom_aromatic"])
    print(f"oracle vs current_engine per-atom-row agreement: {oracle_vs_engine}/{len(valid)} ({100*oracle_vs_engine/len(valid):.2f}%)")

    # ---- polarization check, per molecule (case_id), not per row ----
    by_case = {}
    for r in joined:
        by_case.setdefault(r["case_id"], {"bucket": r["bucket"], "smiles": r["smiles"], "rows": []})
        by_case[r["case_id"]]["rows"].append(r)

    print("\n--- polarization check ---")
    fp_ok, fp_total = 0, 0
    fn_ok, fn_total = 0, 0
    nc_ok, nc_total = 0, 0
    problems = []

    for case_id, case in by_case.items():
        bucket = case["bucket"]
        rows_ = [r for r in case["rows"] if r["rdkit_atom_aromatic"] is not None]
        if not rows_:
            continue
        over = any(r["current_engine_atom_aromatic"] and not r["rdkit_atom_aromatic"] for r in rows_)
        under = any(not r["current_engine_atom_aromatic"] and r["rdkit_atom_aromatic"] for r in rows_)
        all_agree = all(r["current_engine_atom_aromatic"] == r["rdkit_atom_aromatic"] for r in rows_)

        if bucket == "false_positive":
            fp_total += 1
            if over:
                fp_ok += 1
            else:
                problems.append(f"false_positive {case_id} ({case['smiles']}): expected an over-aromatized atom, found none")
        elif bucket == "false_negative":
            fn_total += 1
            if under:
                fn_ok += 1
            else:
                problems.append(f"false_negative {case_id} ({case['smiles']}): expected an under-aromatized atom, found none")
        elif bucket == "negative_control":
            nc_total += 1
            if all_agree:
                nc_ok += 1
            else:
                mismatched = [r["atom_idx"] for r in rows_ if r["current_engine_atom_aromatic"] != r["rdkit_atom_aromatic"]]
                problems.append(f"negative_control {case_id} ({case['smiles']}): expected full agreement, mismatched atoms {mismatched}")

    print(f"false_positive bucket polarized correctly:   {fp_ok}/{fp_total}")
    print(f"false_negative bucket polarized correctly:   {fn_ok}/{fn_total}")
    print(f"negative_control bucket fully agrees:         {nc_ok}/{nc_total}")

    if problems:
        print(f"\n{len(problems)} POLARIZATION PROBLEMS:")
        for p in problems:
            print(f"  - {p}")
    else:
        print("\nAll buckets polarized as labeled. 0 problems.")

    # ---- A1-1a spot-check: oracle vs RDKit on the newly-added confirmatory cases ----
    spot = {
        "O=c1cccccc1": "tropone",
        "O=c1cccc[nH]1": "2-pyridone",
        "O=c1ccocc1": "4-pyranone",
        "c1ccn2ccccc12": "indolizine",
        "c1ccc2cc3ccccc3cc2c1": "anthracene",
        "c1cnc2[nH]cnc2n1": "purine",
        "C1=CC2=CC=CC=CC2=C1": "azulene",
    }
    print("\n--- A1-1a oracle vs RDKit, confirmatory cases ---")
    for smi, name in spot.items():
        rows_ = [r for r in by_case.values() if r["smiles"] == smi]
        if not rows_:
            continue
        rows_ = [r for r in rows_[0]["rows"] if r["rdkit_atom_aromatic"] is not None]
        oracle_match = all(r["oracle_atom_aromatic"] == r["rdkit_atom_aromatic"] for r in rows_)
        engine_match = all(r["current_engine_atom_aromatic"] == r["rdkit_atom_aromatic"] for r in rows_)
        print(f"  {name}: oracle matches RDKit={oracle_match}, current_engine matches RDKit={engine_match}")

    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
