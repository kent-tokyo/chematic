#!/usr/bin/env python3
"""RDKit oracle side of the platinum coordination-chemistry benchmark.

Reads the same corpus (validation/platinum/pt_corpus.jsonl) chematic's own
`cargo run -p chematic-mol --example platinum_benchmark` reads, and produces
a comparable JSONL of RDKit's own parse/formula/mass/charge/canonicalization/
MOL-V3000-round-trip results -- an independent process, not fed by or
feeding chematic's own run.

Also runs a secondary, RDKit-only check (not part of the main per-row JSONL):
whether RDKit's own extended square-planar stereo descriptors (@SP1/@SP2/
@SP3) actually distinguish cis/trans when explicitly present in the input,
to separate "chematic cannot express this at all" from "no one wrote it into
this corpus's plain SMILES".

Usage: python scripts/platinum_rdkit_oracle.py [validation/results/platinum_baseline_rdkit.jsonl]
"""

import json
import sys

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

CORPUS_PATH = "validation/platinum/pt_corpus.jsonl"


def formula_counts_from_mol(mol):
    counts = {}
    for atom in mol.GetAtoms():
        an = atom.GetAtomicNum()
        counts[an] = counts.get(an, 0) + 1
        h = atom.GetTotalNumHs()
        if h:
            counts[1] = counts.get(1, 0) + h
    return counts


def parse_expected_formula(formula):
    """Very small Hill-notation parser matching chematic_chem::formula::parse_formula's
    shape (element symbol + optional count), just enough for this corpus's formulas."""
    import re

    counts = {}
    for sym, num in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if not sym:
            continue
        n = int(num) if num else 1
        z = Chem.GetPeriodicTable().GetAtomicNumber(sym)
        counts[z] = counts.get(z, 0) + n
    return counts


def net_charge(mol):
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


def pt_coordination_number(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "Pt":
            return atom.GetDegree()
    return None


def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else None
    rows = []
    with open(CORPUS_PATH) as f:
        lines = [line for line in f if line.strip()]

    for line in lines:
        entry = json.loads(line)
        cid = entry["id"]
        smiles = entry["smiles_dative"]
        formula_expected = entry["formula_expected"]
        charge_expected = entry["charge_expected"]
        coordination_expected = entry.get("pt_coordination_number")

        row = {
            "id": cid,
            "smiles_input": smiles,
            "formula_expected": formula_expected,
            "charge_expected": charge_expected,
        }

        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            row["parse_ok"] = False
            rows.append(row)
            continue
        try:
            Chem.SanitizeMol(mol)
            row["sanitize_ok"] = True
        except Exception as e:
            row["sanitize_ok"] = False
            row["sanitize_error"] = str(e)
        row["parse_ok"] = True
        row["atom_count"] = mol.GetNumAtoms()
        row["bond_count"] = mol.GetNumBonds()

        counts = formula_counts_from_mol(mol)
        expected_counts = parse_expected_formula(formula_expected)
        row["formula_counts"] = {str(k): v for k, v in counts.items()}
        row["formula_matches_expected"] = counts == expected_counts
        row["net_charge"] = net_charge(mol)
        row["charge_matches_expected"] = net_charge(mol) == charge_expected
        row["molecular_weight"] = rdMolDescriptors.CalcExactMolWt(mol)  # placeholder, overwritten below
        try:
            from rdkit.Chem import Descriptors

            row["molecular_weight"] = Descriptors.MolWt(mol)
            row["exact_mass"] = Descriptors.ExactMolWt(mol)
        except Exception as e:
            row["mass_error"] = str(e)

        cn = pt_coordination_number(mol)
        row["pt_coordination_number_observed"] = cn
        row["pt_coordination_number_matches_expected"] = (
            True if coordination_expected is None else coordination_expected == cn
        )
        row["connected_components"] = len(Chem.GetMolFrags(mol))

        try:
            row["canonical_smiles"] = Chem.MolToSmiles(mol)
        except Exception as e:
            row["canonical_smiles_error"] = str(e)

        # MOL V3000 round-trip
        try:
            molblock = Chem.MolToMolBlock(mol, forceV3000=True)
            rt_mol = Chem.MolFromMolBlock(molblock, sanitize=False)
            if rt_mol is not None:
                Chem.SanitizeMol(rt_mol, catchErrors=True)
                rt_counts = formula_counts_from_mol(rt_mol)
                row["mol_v3000_roundtrip_ok"] = True
                row["mol_v3000_roundtrip_formula_preserved"] = rt_counts == counts
                row["mol_v3000_roundtrip_charge_preserved"] = net_charge(rt_mol) == net_charge(mol)
                row["mol_v3000_roundtrip_coordination_preserved"] = (
                    pt_coordination_number(rt_mol) == cn
                )
                dative_before = sum(
                    1 for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DATIVE
                )
                dative_after = sum(
                    1 for b in rt_mol.GetBonds() if b.GetBondType() == Chem.BondType.DATIVE
                )
                row["mol_v3000_roundtrip_dative_bonds_before"] = dative_before
                row["mol_v3000_roundtrip_dative_bonds_after"] = dative_after
                row["mol_v3000_roundtrip_dative_preserved"] = dative_before == dative_after
            else:
                row["mol_v3000_roundtrip_ok"] = False
        except Exception as e:
            row["mol_v3000_roundtrip_ok"] = False
            row["mol_v3000_roundtrip_error"] = str(e)

        # Fingerprint smoke check (Morgan/ECFP4-equivalent)
        try:
            from rdkit.Chem import AllChem

            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            row["ecfp4_ok"] = True
            row["ecfp4_deterministic"] = fp1.ToBitString() == fp2.ToBitString()
        except Exception:
            row["ecfp4_ok"] = False

        rows.append(row)

    cisplatin_canon = next((r.get("canonical_smiles") for r in rows if r["id"] == "cisplatin"), None)
    transplatin_canon = next(
        (r.get("canonical_smiles") for r in rows if r["id"] == "transplatin"), None
    )
    killer_benchmark = {
        "cisplatin_canonical": cisplatin_canon,
        "transplatin_canonical": transplatin_canon,
        "cis_trans_distinguished": bool(
            cisplatin_canon and transplatin_canon and cisplatin_canon != transplatin_canon
        ),
    }

    # Secondary check: does RDKit distinguish cis/trans when @SP1/@SP2/@SP3
    # IS present in the input (unlike this corpus's plain SMILES)? Answers
    # "could RDKit do this with better input" vs "chematic cannot do this at
    # all" (chematic rejects @SP1/@SP2/@SP3 syntax outright -- see
    # FEASIBILITY.md).
    sp_variants = {}
    for tag in ("SP1", "SP2", "SP3"):
        smi = f"N[Pt@{tag}](N)(Cl)Cl"
        m = Chem.MolFromSmiles(smi)
        sp_variants[tag] = Chem.MolToSmiles(m) if m else None
    distinct_sp_canonicals = len(set(v for v in sp_variants.values() if v is not None))
    square_planar_stereo_check = {
        "sp_variant_canonical_smiles": sp_variants,
        "distinct_canonical_forms": distinct_sp_canonicals,
        "rdkit_can_distinguish_with_explicit_tags": distinct_sp_canonicals > 1,
    }

    if out_path:
        with open(out_path, "w") as f:
            for row in rows:
                f.write(json.dumps(row) + "\n")
        summary_path = out_path.replace(".jsonl", "_summary.json")
        with open(summary_path, "w") as f:
            json.dump(
                {
                    "killer_benchmark": killer_benchmark,
                    "square_planar_stereo_check": square_planar_stereo_check,
                },
                f,
                indent=2,
            )
        print(f"wrote {len(rows)} rows to {out_path}", file=sys.stderr)
        print(f"summary written to {summary_path}", file=sys.stderr)
    else:
        print(json.dumps({"rows": rows, "killer_benchmark": killer_benchmark}, indent=2))


if __name__ == "__main__":
    main()
