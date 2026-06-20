#!/usr/bin/env python3
"""
Train an ECFP4-based MLP for aqueous solubility prediction and emit
Rust constant arrays ready to paste into mlp.rs.

Usage
-----
1. Download the Delaney ESOL dataset (CSV with 'smiles' and 'measured log(solubility:mol/L)'):
       wget -O data/delaney.csv "https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/delaney-processed.csv"

2. (Optional) Use the larger AqSolDB instead:
       # 9982 molecules, column: 'Solubility'
       # https://github.com/PatWalters/solubility/blob/master/aqsoldb.csv

3. Run:
       python scripts/train_solubility_mlp.py --csv data/delaney.csv
       # or with AqSolDB:
       python scripts/train_solubility_mlp.py --csv data/aqsoldb.csv --col Solubility --smiles SMILES

4. Copy the printed Rust constants into:
       crates/chematic-chem/src/mlp.rs

5. Set:  pub const MLP_SOLUBILITY_TRAINED: bool = true;

Dependencies
------------
    pip install chematic scikit-learn pandas numpy
"""

import argparse
import sys
import numpy as np
import pandas as pd

try:
    import chematic
except ImportError:
    sys.exit("Install chematic first:  pip install chematic  (or: maturin develop)")

try:
    from sklearn.neural_network import MLPRegressor
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import cross_val_score
    from sklearn.metrics import mean_squared_error, r2_score
except ImportError:
    sys.exit("Install scikit-learn:  pip install scikit-learn")


def main():
    parser = argparse.ArgumentParser(description="Train ECFP4-MLP solubility model")
    parser.add_argument("--csv",    default="data/delaney.csv",
                        help="Path to CSV file with SMILES and logS columns")
    parser.add_argument("--smiles", default="smiles",
                        help="Column name for SMILES (default: 'smiles')")
    parser.add_argument("--col",    default="measured log(solubility:mol/L)",
                        help="Column name for logS values")
    parser.add_argument("--hidden", default=64, type=int,
                        help="Hidden-layer size (default: 64)")
    parser.add_argument("--seed",   default=42, type=int)
    args = parser.parse_args()

    # ── 1. Load data ──────────────────────────────────────────────────────
    print(f"Loading {args.csv} ...", flush=True)
    df = pd.read_csv(args.csv)
    smiles_col = args.smiles
    label_col  = args.col

    # Validate columns
    missing = [c for c in [smiles_col, label_col] if c not in df.columns]
    if missing:
        print(f"ERROR: columns not found: {missing}")
        print(f"Available: {list(df.columns)}")
        sys.exit(1)

    df = df[[smiles_col, label_col]].dropna()
    print(f"  {len(df)} molecules after dropna")

    # ── 2. ECFP4 fingerprints ─────────────────────────────────────────────
    print("Computing ECFP4 fingerprints ...", flush=True)
    smiles_list = df[smiles_col].tolist()
    labels      = df[label_col].to_numpy(dtype=np.float32)

    # Filter invalid SMILES
    valid_idx   = []
    valid_smiles = []
    for i, smi in enumerate(smiles_list):
        try:
            chematic.from_smiles(smi)
            valid_idx.append(i)
            valid_smiles.append(smi)
        except Exception:
            pass

    if len(valid_idx) < len(smiles_list):
        n_invalid = len(smiles_list) - len(valid_idx)
        print(f"  Skipped {n_invalid} invalid SMILES")
    labels = labels[valid_idx]

    fps = chematic.bulk.ecfp4(valid_smiles)   # (N, 2048) uint8
    X   = fps.astype(np.float32)
    y   = labels

    print(f"  X.shape={X.shape}  y.shape={y.shape}")

    # ── 3. Train MLP ──────────────────────────────────────────────────────
    print(f"Training MLP (hidden={args.hidden}) ...", flush=True)
    mlp = MLPRegressor(
        hidden_layer_sizes=(args.hidden,),
        activation="relu",
        max_iter=1000,
        random_state=args.seed,
        early_stopping=True,
        validation_fraction=0.1,
    )
    mlp.fit(X, y)

    # ── 4. Evaluate ───────────────────────────────────────────────────────
    y_pred = mlp.predict(X)
    rmse   = mean_squared_error(y, y_pred) ** 0.5
    r2     = r2_score(y, y_pred)
    print(f"Train RMSE={rmse:.3f}  R²={r2:.3f}")

    cv_r2 = cross_val_score(mlp, X, y, cv=5, scoring="r2")
    print(f"5-fold CV R²: {cv_r2.mean():.3f} ± {cv_r2.std():.3f}")

    # ── 5. Emit Rust constants ────────────────────────────────────────────
    W0, B0 = mlp.coefs_[0], mlp.intercepts_[0]   # (2048 → hidden)
    W1, B1 = mlp.coefs_[1], mlp.intercepts_[1]   # (hidden → 1)

    assert W0.shape == (2048, args.hidden), f"unexpected shape {W0.shape}"
    assert W1.shape == (args.hidden, 1),    f"unexpected shape {W1.shape}"

    # sklearn stores W as (in, out); Rust expects row-major (out, in)
    W0_rust = W0.T.flatten().astype(np.float32)   # (hidden, 2048)
    B0_rust = B0.flatten().astype(np.float32)
    W1_rust = W1.T.flatten().astype(np.float32)   # (1, hidden)
    B1_rust = B1.flatten().astype(np.float32)

    def fmt_array(name: str, arr: np.ndarray) -> str:
        vals = ", ".join(f"{v:.6e}_f32" for v in arr)
        return f"const {name}: &[f32] = &[{vals}];"

    print()
    print("=" * 72)
    print("// Copy the constants below into crates/chematic-chem/src/mlp.rs")
    print(f"// and set: pub const MLP_SOLUBILITY_TRAINED: bool = true;")
    print(f"// Architecture: ECFP4(2048) → Dense({args.hidden}, ReLU) → Dense(1, linear)")
    print(f"// Train RMSE={rmse:.3f}  5-fold CV R²={cv_r2.mean():.3f}")
    print(f"const HIDDEN1: usize = {args.hidden};")
    print(fmt_array("W0", W0_rust))
    print(fmt_array("B0", B0_rust))
    print(fmt_array("W1", W1_rust))
    print(fmt_array("B1", B1_rust))
    print("=" * 72)


if __name__ == "__main__":
    main()
