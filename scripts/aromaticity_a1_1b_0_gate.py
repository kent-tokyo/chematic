#!/usr/bin/env python3
"""Aromaticity-A1-1b-0 gate: join rdkit_parity_aromaticity's output against
real RDKit and measure the user's exact gate criteria.

Reads validation/results/aromaticity_a1_1b_0_trace.jsonl (from
`cargo run -p chematic-perception --release --example aromaticity_a1_1b_0_report`),
joins against real RDKit per (smiles, atom_idx), and reports:
  - 33 false-positive corpus molecules: all fixed?
  - 5 false-negative corpus molecules: all fixed?
  - 17 negative-control molecules: all maintained?
  - RDKit atom-flag agreement: 100%?
  - unexplained differences: 0?

Diagnostic only -- reads frozen inputs, writes one joined output file; does
not change any production code or behavior.

Run:
    cargo run -p chematic-perception --release --example aromaticity_a1_1b_0_report \
        -- validation/aromaticity_a1_0_corpus.jsonl \
        > validation/results/aromaticity_a1_1b_0_trace.jsonl
    python3 scripts/aromaticity_a1_1b_0_gate.py
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
TRACE_IN = REPO / "validation" / "results" / "aromaticity_a1_1b_0_trace.jsonl"
OUT = REPO / "validation" / "results" / "aromaticity_a1_1b_0_gate.jsonl"


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
        row["rdkit_atom_aromatic"] = None if smi in unparseable else rdkit_cache[smi].get(row["atom_idx"])
        joined.append(row)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        for row in joined:
            f.write(json.dumps(row) + "\n")
    print(f"wrote {len(joined)} joined rows to {OUT}\n")

    valid = [r for r in joined if r["rdkit_atom_aromatic"] is not None]
    agree = sum(1 for r in valid if r["rdkit_parity_atom_aromatic"] == r["rdkit_atom_aromatic"])
    print(f"rdkit_parity_aromaticity vs real RDKit, per-atom-row agreement: {agree}/{len(valid)} ({100*agree/len(valid):.2f}%)")

    by_case = {}
    for r in joined:
        by_case.setdefault(r["case_id"], {"bucket": r["bucket"], "smiles": r["smiles"], "rows": []})
        by_case[r["case_id"]]["rows"].append(r)

    fp_fixed, fp_total = 0, 0
    fn_fixed, fn_total = 0, 0
    nc_ok, nc_total = 0, 0
    problems = []

    for case_id, case in by_case.items():
        rows_ = [r for r in case["rows"] if r["rdkit_atom_aromatic"] is not None]
        if not rows_:
            continue
        all_agree = all(r["rdkit_parity_atom_aromatic"] == r["rdkit_atom_aromatic"] for r in rows_)
        bucket = case["bucket"]
        if bucket == "false_positive":
            fp_total += 1
            if all_agree:
                fp_fixed += 1
            else:
                mismatched = [r["atom_idx"] for r in rows_ if r["rdkit_parity_atom_aromatic"] != r["rdkit_atom_aromatic"]]
                problems.append(f"false_positive {case_id} ({case['smiles']}): still wrong, mismatched atoms {mismatched}")
        elif bucket == "false_negative":
            fn_total += 1
            if all_agree:
                fn_fixed += 1
            else:
                mismatched = [r["atom_idx"] for r in rows_ if r["rdkit_parity_atom_aromatic"] != r["rdkit_atom_aromatic"]]
                problems.append(f"false_negative {case_id} ({case['smiles']}): still wrong, mismatched atoms {mismatched}")
        elif bucket == "negative_control":
            nc_total += 1
            if all_agree:
                nc_ok += 1
            else:
                mismatched = [r["atom_idx"] for r in rows_ if r["rdkit_parity_atom_aromatic"] != r["rdkit_atom_aromatic"]]
                problems.append(f"negative_control {case_id} ({case['smiles']}): REGRESSED, mismatched atoms {mismatched}")

    print(f"\nfalse_positive corpus fixed:  {fp_fixed}/{fp_total}")
    print(f"false_negative corpus fixed: {fn_fixed}/{fn_total}")
    print(f"negative_control maintained: {nc_ok}/{nc_total}")

    gate_pass = (
        fp_fixed == fp_total
        and fn_fixed == fn_total
        and nc_ok == nc_total
        and agree == len(valid)
    )
    print(f"\n=== GATE: {'PASS' if gate_pass else 'NOT MET'} ===")
    print(f"RDKit atom-flag agreement: {100*agree/len(valid):.2f}% (target 100.00%)")
    print(f"Unexplained differences: {len(problems)} (target 0)")

    if problems:
        print(f"\n{len(problems)} case(s) not yet matching RDKit:")
        for p in problems:
            print(f"  - {p}")

    return 0 if gate_pass else 1


if __name__ == "__main__":
    sys.exit(main())
