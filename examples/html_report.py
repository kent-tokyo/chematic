"""
HTML compound report
====================
Generates a self-contained HTML file with a grid of molecule cards.
Cards display 2D structure, MW/LogP/TPSA/QED, and Lipinski/PAINS/Brenk flags.
Sorted by QED descending (most drug-like first).

Run:
    python examples/html_report.py                      # demo set → report.html
    python examples/html_report.py smiles.txt           # one SMILES per line
    python examples/html_report.py smiles.txt out.html  # custom output path

Dependencies:
    pip install chematic
"""
import sys
import chematic

DEMO = [
    ("aspirin",     "CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen",   "CC(C)Cc1ccc(CC(C)C(=O)O)cc1"),
    ("caffeine",    "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("naproxen",    "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
    ("celecoxib",   "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1"),
]


def main() -> None:
    out_path = sys.argv[2] if len(sys.argv) > 2 else "report.html"

    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            lines = [l.strip() for l in f if l.strip()]
        pairs = [(f"mol_{i}", smi) for i, smi in enumerate(lines, 1)]
    else:
        pairs = DEMO

    mols = []
    names = []
    for name, smi in pairs:
        try:
            mol = chematic.from_smiles(smi)
        except ValueError:
            continue
        if mol is not None:
            mols.append(mol)
            names.append(name)

    chematic.report(
        mols,
        names=names,
        title=f"chematic report — {len(mols)} molecules",
        output=out_path,
    )
    print(f"Wrote {len(mols)} molecules → {out_path}")
    print(f"Open in browser: open {out_path}")


if __name__ == "__main__":
    main()
