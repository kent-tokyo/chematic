"""
Molecular descriptors → pandas → ML-ready feature matrix
=========================================================
Shows how bulk.descriptors() produces a pandas DataFrame in one line,
ready for exploratory analysis, Lipinski filtering, and sklearn pipelines.

Dependencies:
    pip install chematic pandas numpy scikit-learn

Run:
    python examples/descriptors_pandas.py
"""

import pandas as pd
import numpy as np
import chematic

# ---------------------------------------------------------------------------
# 1. Molecules of interest (drug-like diversity set)
# ---------------------------------------------------------------------------
SMILES = [
    "c1ccccc1",                              # benzene
    "CC(=O)Oc1ccccc1C(=O)O",               # aspirin
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",          # caffeine
    "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",     # ibuprofen
    "OCC(O)C(O)C(O)C(O)CO",                # sorbitol
    "c1ccc2ccccc2c1",                       # naphthalene
    "c1ccncc1",                             # pyridine
    "CC(=O)Nc1ccc(O)cc1",                  # paracetamol
    "NCC(=O)O",                             # glycine
    "CN(C)C(=N)NC(=N)N",                   # metformin
]
NAMES = [
    "benzene", "aspirin", "caffeine", "ibuprofen", "sorbitol",
    "naphthalene", "pyridine", "paracetamol", "glycine", "metformin",
]

# ---------------------------------------------------------------------------
# 2. Compute 55+ descriptors in one call → list[dict] → DataFrame
#    bulk.descriptors() runs in parallel (Rayon), returns one dict per molecule.
# ---------------------------------------------------------------------------
descs = chematic.bulk.descriptors(SMILES)   # list[dict], ≈1 ms for 10 mols
df = pd.DataFrame(descs)
df.insert(0, "name", NAMES)
df.insert(1, "smiles", SMILES)

# ---------------------------------------------------------------------------
# 3. Key physicochemical properties
# ---------------------------------------------------------------------------
COLS = ["name", "mw", "logp", "tpsa", "hbd", "hba", "rotatable_bonds",
        "fsp3", "aromatic_ring_count", "qed"]
print("── Physicochemical Properties ──────────────────────────────────────")
print(df[COLS].to_string(index=False))

# ---------------------------------------------------------------------------
# 4. Lipinski Ro5 filter
# ---------------------------------------------------------------------------
ro5 = df[
    (df["mw"]  <= 500) &
    (df["logp"] <= 5)  &
    (df["hbd"]  <= 5)  &
    (df["hba"]  <= 10)
]
print(f"\nLipinski Ro5 pass : {len(ro5)}/{len(df)}")
print("  Failing:", [n for n in NAMES if n not in ro5["name"].values])

# ---------------------------------------------------------------------------
# 5. Drug-likeness ranking by QED
# ---------------------------------------------------------------------------
print("\n── QED Ranking (higher = more drug-like) ───────────────────────────")
print(df[["name", "qed", "mw", "logp", "tpsa"]].sort_values("qed", ascending=False).to_string(index=False))

# ---------------------------------------------------------------------------
# 6. Feature matrix for sklearn
#    Select numeric columns, drop non-numeric, fill NaN → ready for ML
# ---------------------------------------------------------------------------
feature_cols = [c for c in df.columns if c not in ("name", "smiles")]
X = df[feature_cols].select_dtypes(include="number").fillna(0).values
print(f"\nFeature matrix for sklearn : shape={X.shape}")

# Example: PCA to visualise chemical space
try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2).fit_transform(X_scaled)
    print("PCA projections (PC1, PC2):")
    for name, (pc1, pc2) in zip(NAMES, pca):
        print(f"  {name:12s}  {pc1:+.2f}  {pc2:+.2f}")
except ImportError:
    print("(install scikit-learn to run the PCA section)")

# ---------------------------------------------------------------------------
# 7. Fingerprint diversity — pairwise Tanimoto with tanimoto_matrix
# ---------------------------------------------------------------------------
fps = [chematic.from_smiles(s).ecfp4() for s in SMILES]
sim_matrix = chematic.tanimoto_matrix(fps, fps)   # (N, N) list[list[float]]
sim_arr = np.array(sim_matrix)
print(f"\nTanimoto matrix    : shape={sim_arr.shape}, "
      f"mean_offdiag={sim_arr[~np.eye(len(SMILES), dtype=bool)].mean():.3f}")
