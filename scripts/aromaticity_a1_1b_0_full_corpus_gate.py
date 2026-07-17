#!/usr/bin/env python3
"""Aromaticity-A1-1b-0 full-corpus gate: SET-level (not count-level) join of
`rdkit_parity_aromaticity` against real RDKit across the full 5000-molecule
benchmark corpus.

A prior count-only comparison (aromatic_atom_count == rd_atom_count) is blind
to same-count/different-atoms mismatches. This joins per (smiles, atom_idx)
and per (smiles, bond_atoms) instead, same join key style as the 55-molecule
diagnosis-corpus gate (scripts/aromaticity_a1_1b_0_gate.py).

Run:
    cargo run -p chematic-perception --release --example rdkit_parity_full_corpus \
        -- ~/Downloads/SMILES.csv > /tmp/rdkit_parity_full_corpus_trace.jsonl
    python3 scripts/aromaticity_a1_1b_0_full_corpus_gate.py /tmp/rdkit_parity_full_corpus_trace.jsonl
"""
import json
import sys

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def main():
    trace_path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/rdkit_parity_full_corpus_trace.jsonl"
    rows = [json.loads(line) for line in open(trace_path) if line.strip()]
    print(f"loaded {len(rows)} rows from {trace_path}")

    by_smiles = {}
    for row in rows:
        by_smiles.setdefault(row["smiles"], []).append(row)

    atom_total = atom_agree = 0
    bond_total = bond_agree = 0
    atom_mismatches = []
    bond_mismatches = []
    unparseable = 0

    for smi, smi_rows in by_smiles.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            unparseable += 1
            continue
        rd_atom_arom = {a.GetIdx(): a.GetIsAromatic() for a in mol.GetAtoms()}
        rd_bond_arom = {}
        for b in mol.GetBonds():
            key = tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
            rd_bond_arom[key] = b.GetIsAromatic()

        for row in smi_rows:
            if row["kind"] == "atom":
                atom_total += 1
                rd = rd_atom_arom.get(row["atom_idx"])
                ours = row["rdkit_parity_atom_aromatic"]
                if rd == ours:
                    atom_agree += 1
                else:
                    atom_mismatches.append((smi, row["atom_idx"], ours, rd))
            else:
                bond_total += 1
                key = tuple(sorted(row["bond_atoms"]))
                rd = rd_bond_arom.get(key)
                ours = row["rdkit_parity_bond_aromatic"]
                if rd == ours:
                    bond_agree += 1
                else:
                    bond_mismatches.append((smi, key, ours, rd))

    print(f"\nunparseable-by-RDKit molecules skipped: {unparseable}")
    print(f"atom-level set agreement: {atom_agree}/{atom_total} ({100*atom_agree/atom_total:.4f}%)")
    print(f"bond-level set agreement: {bond_agree}/{bond_total} ({100*bond_agree/bond_total:.4f}%)")

    if atom_mismatches:
        print(f"\n{len(atom_mismatches)} atom mismatch(es) (smiles, atom_idx, ours, rdkit):")
        for m in atom_mismatches[:20]:
            print(f"  {m}")
    if bond_mismatches:
        print(f"\n{len(bond_mismatches)} bond mismatch(es) (smiles, bond_atoms, ours, rdkit):")
        for m in bond_mismatches[:20]:
            print(f"  {m}")

    gate_pass = not atom_mismatches and not bond_mismatches
    print(f"\n=== SET-LEVEL GATE: {'PASS' if gate_pass else 'FAIL'} ===")
    return 0 if gate_pass else 1


if __name__ == "__main__":
    sys.exit(main())
