"""
QSAR with chematic + scikit-learn
==================================
Demonstrates that chematic's bulk.ecfp4() integrates directly with sklearn
— no conversion, no wrapper, no glue code required.

Dependencies:
    pip install chematic scikit-learn numpy

Run:
    python examples/qsar_sklearn.py
"""

import numpy as np
import chematic
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score

# ---------------------------------------------------------------------------
# 1. Dataset
#    Two classes: BBB-penetrant (1) vs non-penetrant (0) toy example.
#    In a real project, load from a CSV/SDF file.
# ---------------------------------------------------------------------------
SMILES = [
    # BBB penetrant (logP high, TPSA low, MW small)
    "c1ccccc1",                             # benzene
    "c1ccc2ccccc2c1",                       # naphthalene
    "c1ccncc1",                             # pyridine
    "Cc1ccccc1",                            # toluene
    "c1ccc(cc1)C",                          # toluene (canonical)
    "CCc1ccccc1",                           # ethylbenzene
    "CN1CCCCC1",                            # N-methylpiperidine
    "C1CCNCC1",                             # piperidine
    # Non-penetrant (polar, high TPSA, or large MW)
    "OCC(O)C(O)C(O)C(O)CO",               # sorbitol
    "CC(=O)Oc1ccccc1C(=O)O",              # aspirin
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",         # caffeine
    "OC(=O)c1ccccc1O",                     # salicylic acid
    "NC(=O)c1ccccc1",                      # benzamide
    "NCC(=O)O",                            # glycine
    "CC(N)C(=O)O",                         # alanine
    "OC(=O)CC(=O)O",                       # malonic acid
]
LABELS = [1, 1, 1, 1, 1, 1, 1, 1,
          0, 0, 0, 0, 0, 0, 0, 0]

# ---------------------------------------------------------------------------
# 2. Fingerprint generation — bulk.ecfp4() returns (N, 2048) uint8 ndarray
#    No conversion needed: sklearn accepts uint8 directly.
# ---------------------------------------------------------------------------
fps = chematic.bulk.ecfp4(SMILES)
print(f"Fingerprint matrix : shape={fps.shape}, dtype={fps.dtype}")
# → shape=(16, 2048), dtype=uint8

# ---------------------------------------------------------------------------
# 3. RandomForest QSAR model with cross-validation
# ---------------------------------------------------------------------------
clf = RandomForestClassifier(n_estimators=200, random_state=42)
cv = StratifiedKFold(n_splits=4, shuffle=True, random_state=0)
scores = cross_val_score(clf, fps, LABELS, cv=cv, scoring="roc_auc")
print(f"CV ROC-AUC         : {scores.mean():.3f} ± {scores.std():.3f}")

# ---------------------------------------------------------------------------
# 4. Train on full set, predict a new compound
# ---------------------------------------------------------------------------
clf.fit(fps, LABELS)

query_smiles = [
    "c1ccc(cc1)Cl",   # chlorobenzene  (expect BBB=1)
    "OC(=O)c1ccc(N)cc1",  # 4-aminobenzoic acid (expect BBB=0)
]
q_fps = chematic.bulk.ecfp4(query_smiles)
proba = clf.predict_proba(q_fps)[:, 1]
for smi, p in zip(query_smiles, proba):
    print(f"  {smi:40s}  BBB prob = {p:.2f}")

# ---------------------------------------------------------------------------
# 5. MACCS keys — (N, 166) uint8 — also plugs straight into sklearn
# ---------------------------------------------------------------------------
maccs = chematic.bulk.maccs(SMILES)
clf_maccs = RandomForestClassifier(n_estimators=100, random_state=0)
scores_maccs = cross_val_score(clf_maccs, maccs, LABELS, cv=cv, scoring="roc_auc")
print(f"MACCS CV ROC-AUC   : {scores_maccs.mean():.3f} ± {scores_maccs.std():.3f}")

# ---------------------------------------------------------------------------
# 6. Tanimoto similarity screen (no sklearn needed)
# ---------------------------------------------------------------------------
query_mol = chematic.from_smiles("c1cccnc1")  # pyridine
ranked = sorted(
    [(smi, chematic.tanimoto(query_mol.ecfp4(), chematic.from_smiles(smi).ecfp4()))
     for smi in SMILES],
    key=lambda x: -x[1],
)
print("\nTanimoto nearest neighbours (query = pyridine):")
for smi, sim in ranked[:5]:
    print(f"  {smi:45s}  sim={sim:.3f}")
