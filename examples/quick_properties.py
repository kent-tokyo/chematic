"""
Quick property summary for one or more SMILES strings
=====================================================
Shows mol.describe() — the shortest path from SMILES to a human-readable
property report.

Run:
    python examples/quick_properties.py
    python examples/quick_properties.py "CC(=O)Oc1ccccc1C(=O)O"
    python examples/quick_properties.py "CC(=O)Oc1ccccc1C(=O)O" "CC(C)Cc1ccc(CC(C)C(=O)O)cc1"

Dependencies:
    pip install chematic
"""
import sys
import chematic

DEFAULT_SMILES = [
    ("aspirin",    "CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen",  "CC(C)Cc1ccc(CC(C)C(=O)O)cc1"),
    ("caffeine",   "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("paracetamol","CC(=O)Nc1ccc(O)cc1"),
]


def report(name: str, smiles: str) -> None:
    mol = chematic.from_smiles(smiles)
    if mol is None:
        print(f"[{name}] invalid SMILES: {smiles}\n")
        return
    print(f"=== {name} ===")
    print(f"SMILES: {mol.smiles}")
    print(mol.describe())
    print()


def main() -> None:
    if len(sys.argv) > 1:
        for i, smi in enumerate(sys.argv[1:], 1):
            report(f"compound_{i}", smi)
    else:
        for name, smi in DEFAULT_SMILES:
            report(name, smi)


if __name__ == "__main__":
    main()
