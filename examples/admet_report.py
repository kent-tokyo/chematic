"""
ADMET profile report → CSV
===========================
Reads a SMILES list, computes ADMET predictions for each molecule,
and writes a filtered CSV. Shows the typical lead-optimisation screening
workflow: load → filter → profile → export.

Run:
    python examples/admet_report.py                      # uses built-in demo set
    python examples/admet_report.py smiles.txt           # one SMILES per line
    python examples/admet_report.py smiles.txt out.csv   # custom output path

Dependencies:
    pip install chematic
"""
import csv
import sys
import chematic

DEMO_SMILES = [
    ("aspirin",     "CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen",   "CC(C)Cc1ccc(CC(C)C(=O)O)cc1"),
    ("caffeine",    "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("verapamil",   "COc1ccc(CCN(C)CCCC(C#N)(c2ccc(OC)c(OC)c2)C(C)C)cc1OC"),
    ("erythromycin","CCC1OC(=O)C(CC(CC(C(C(C(OC2CC(CC(O2)C)N(C)C)C)OC3OC(C)CC(C3O)OC)C)O)OC)C(C1O)(C)O"),
]

FIELDNAMES = [
    "name", "smiles", "mw", "logp", "tpsa", "hbd", "hba", "qed",
    "lipinski_passes", "pains_passes", "brenk_passes",
    "bbb_penetrant", "caco2_high", "herg_risk", "cyp3a4_risk",
]


def process(name: str, smiles: str) -> dict | None:
    mol = chematic.from_smiles(smiles)
    if mol is None:
        return None
    admet = mol.admet()
    return {
        "name":            name,
        "smiles":          mol.smiles,
        "mw":              round(mol.mw, 2),
        "logp":            round(mol.logp, 2),
        "tpsa":            round(mol.tpsa, 1),
        "hbd":             mol.hbd,
        "hba":             mol.hba,
        "qed":             round(mol.qed, 3),
        "lipinski_passes": mol.lipinski_passes,
        "pains_passes":    mol.pains_passes,
        "brenk_passes":    mol.brenk_passes,
        "bbb_penetrant":   admet.get("bbb_penetrant"),
        "caco2_high":      admet.get("caco2_high"),
        "herg_risk":       admet.get("herg_risk"),
        "cyp3a4_risk":     admet.get("cyp3a4_risk"),
    }


def main() -> None:
    smiles_arg = sys.argv[1] if len(sys.argv) > 1 else None
    out_path   = sys.argv[2] if len(sys.argv) > 2 else "admet_report.csv"

    if smiles_arg:
        with open(smiles_arg) as f:
            pairs = [(f"mol_{i}", line.strip()) for i, line in enumerate(f, 1) if line.strip()]
    else:
        pairs = DEMO_SMILES

    rows = [r for name, smi in pairs if (r := process(name, smi)) is not None]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} molecules → {out_path}")
    for r in rows:
        flag = "" if r["lipinski_passes"] and r["pains_passes"] else "  ⚠"
        print(f"  {r['name']:15s}  MW={r['mw']:6.1f}  LogP={r['logp']:5.2f}  QED={r['qed']:.3f}  BBB={r['bbb_penetrant']}{flag}")


if __name__ == "__main__":
    main()
