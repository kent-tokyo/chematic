"""
chematic — built-in ML solubility model demo
=============================================
Shows how to use `mol.ml_solubility` — the ECFP4-based MLP predictor
embedded in the chematic binary (WASM-compatible, zero C++ dependencies).

In placeholder mode (default), ml_solubility falls back to the Delaney
ESOL linear regression.  To activate the full neural network:

  1. Download training data:
       wget -O data/delaney.csv "https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/delaney-processed.csv"

  2. Train and emit Rust weights:
       python scripts/train_solubility_mlp.py --csv data/delaney.csv

  3. Paste the printed constants into crates/chematic-chem/src/mlp.rs
     and set MLP_SOLUBILITY_TRAINED = true, then rebuild:
       maturin develop --release -m crates/chematic-py/Cargo.toml

Dependencies:
    pip install chematic

Run:
    python examples/ml_builtin_model.py
"""

import chematic

MOLECULES = [
    ("water",        "O"),
    ("ethanol",      "CCO"),
    ("benzene",      "c1ccccc1"),
    ("naphthalene",  "c1ccc2ccccc2c1"),
    ("caffeine",     "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
    ("ibuprofen",    "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"),
    ("aspirin",      "CC(=O)Oc1ccccc1C(=O)O"),
    ("glucose",      "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
    ("cholesterol",  "CC(C)CCCC(C)[C@@H]1CC[C@@H]2[C@H]1CC=C1CC(O)CC[C@]21C"),
    ("hexadecane",   "CCCCCCCCCCCCCCCC"),
]

print("chematic built-in ML solubility predictor")
print(f"  Mode: {'trained MLP' if chematic.from_smiles('C').ml_solubility != chematic.from_smiles('C').esol else 'placeholder (ESOL fallback)'}")
print()
print(f"{'Compound':<15}  {'SMILES':<45}  {'logS (ML)':>10}  {'logS (ESOL)':>12}")
print("-" * 90)

for name, smi in MOLECULES:
    mol = chematic.from_smiles(smi)
    ml   = mol.ml_solubility   # ECFP4 MLP (or ESOL fallback)
    esol = mol.esol            # Delaney ESOL linear regression
    diff = ml - esol
    diff_str = f"({diff:+.2f})" if abs(diff) > 0.01 else ""
    print(f"{name:<15}  {smi:<45}  {ml:>10.3f}  {esol:>12.3f}  {diff_str}")

print()
print("Interpretation: logS -2 = moderately soluble, -5 = poorly, -7 = very poorly")
print()

# ---------------------------------------------------------------------------
# Comparison with sklearn model (external, using chematic ECFP4 as features)
# ---------------------------------------------------------------------------
try:
    import numpy as np
    from sklearn.linear_model import Ridge

    print("Training a quick sklearn ridge regression on the same molecules ...")
    smiles  = [smi for _, smi in MOLECULES]
    targets = [chematic.from_smiles(smi).esol for smi in smiles]  # use ESOL as pseudo-labels

    fps = chematic.bulk.ecfp4(smiles)   # (N, 2048) uint8, sklearn-ready
    clf = Ridge(alpha=1.0).fit(fps, targets)
    preds = clf.predict(fps)

    print(f"{'Compound':<15}  {'sklearn ridge':>14}  {'chematic MLP':>13}")
    print("-" * 46)
    for (name, _), sk, ch in zip(MOLECULES, preds, [chematic.from_smiles(smi).ml_solubility for _, smi in MOLECULES]):
        print(f"{name:<15}  {sk:>14.3f}  {ch:>13.3f}")
except ImportError:
    print("(install scikit-learn to see the sklearn comparison)")
