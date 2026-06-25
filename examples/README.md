# chematic — Python Examples

Runnable Python scripts. Each script works standalone with `pip install chematic`.

## Setup

```bash
pip install chematic                          # required for all scripts
pip install chematic scikit-learn pandas numpy  # for ML examples
```

## Scripts

### `quick_properties.py` — SMILES → human-readable property report

The simplest way to see what chematic gives you. Calls `mol.describe()` on one or more SMILES strings.

```
python examples/quick_properties.py
python examples/quick_properties.py "CC(=O)Oc1ccccc1C(=O)O"
```

Output:
```
=== aspirin ===
SMILES: OC(=O)c1ccccc1OC(=O)C
Molecular weight 180.2 Da, formula C9H8O4.
LogP 1.31 (mildly lipophilic), TPSA 63.6 Å².
HBD 1, HBA 3, 3 rotatable bond(s), 1 aromatic ring(s).
Drug-likeness: no Lipinski rule-of-5 violations. likely orally bioavailable (passes Veber criteria).
QED 0.56 (0 = non-drug-like, 1 = ideal).
Structural alerts: Brenk alert.
```

### `admet_report.py` — compound library → ADMET CSV

Reads a SMILES list (built-in demo set or a text file), computes ADMET properties, and writes `admet_report.csv`.

```
python examples/admet_report.py
python examples/admet_report.py smiles.txt
python examples/admet_report.py smiles.txt out.csv
```

Output:
```
Wrote 6 molecules → admet_report.csv
  aspirin          MW= 180.2  LogP= 1.31  QED=0.563  BBB=True
  ibuprofen        MW= 220.3  LogP= 3.15  QED=0.813  BBB=True
  caffeine         MW= 194.2  LogP=-1.03  QED=0.538  BBB=False
  paracetamol      MW= 151.2  LogP= 1.35  QED=0.595  BBB=True
  verapamil        MW= 454.6  LogP= 5.09  QED=0.420  BBB=True  ⚠
  erythromycin     MW= 663.9  LogP= 1.60  QED=0.187  BBB=False  ⚠
```

### `html_report.py` — self-contained HTML compound report

Generates a browser-ready HTML file with a grid of molecule cards (2D structure + descriptors + drug-likeness flags). Sorted by QED. No external CSS or JavaScript.

```
python examples/html_report.py                  # demo set → report.html
python examples/html_report.py smiles.txt        # one SMILES per line
python examples/html_report.py smiles.txt out.html
```

### `ai_agent_demo.py` — LLM / AI agent context builder

Demonstrates `mol.describe()` and `mol.diff()` as building blocks for LLM prompts. Prints a ready-to-paste comparison prompt for Claude / GPT-4.

```
python examples/ai_agent_demo.py
python examples/ai_agent_demo.py "CC(=O)Oc1ccccc1C(=O)O" "CC(C)Cc1ccc(CC(C)C(=O)O)cc1"
```

### `qsar_sklearn.py` — QSAR with scikit-learn

Builds a binary classifier (BBB penetration toy dataset) using ECFP4 fingerprints and RandomForest.
Demonstrates that `bulk.ecfp4()` → `(N, 2048) uint8 ndarray` plugs directly into sklearn with zero conversion.

```
python examples/qsar_sklearn.py
```

Key patterns:

```python
import chematic
from sklearn.ensemble import RandomForestClassifier

fps = chematic.bulk.ecfp4(smiles_list)   # (N, 2048) uint8 ndarray — sklearn-ready
clf = RandomForestClassifier().fit(fps, labels)

maccs = chematic.bulk.maccs(smiles_list)  # (N, 166) uint8 ndarray
```

### `descriptors_pandas.py` — Descriptors → pandas → ML feature matrix

Computes 55+ molecular descriptors in one call, converts to a pandas DataFrame,
applies Lipinski filtering, ranks by QED, and feeds descriptors into sklearn PCA.

```
python examples/descriptors_pandas.py
```

Key patterns:

```python
import pandas as pd
import chematic

descs = chematic.bulk.descriptors(smiles_list)  # list[dict]
df = pd.DataFrame(descs)                        # one-liner → DataFrame

# Lipinski filter
ro5 = df[(df["mw"] <= 500) & (df["logp"] <= 5) & (df["hbd"] <= 5) & (df["hba"] <= 10)]

# Pairwise Tanimoto matrix
fps = [chematic.from_smiles(s).ecfp4() for s in smiles_list]
sim = chematic.tanimoto_matrix(fps, fps)  # (N, N) list[list[float]]
```

### `aizynthfinder_integration.py` — AiZynthFinder + chematic retrosynthesis pipeline

End-to-end retrosynthetic planning tutorial showing how to use chematic as a backend
for AiZynthFinder.  Runs fully without AiZynthFinder installed (chematic-only mode
with mock routes); activates the ML multi-step sections when AiZynthFinder is present.

```
python examples/aizynthfinder_integration.py
python examples/aizynthfinder_integration.py --smiles "O=C(O)c1ccc(N)cc1"
python examples/aizynthfinder_integration.py --config aizynthfinder_data/config.yml
```

Sections:

1. **Target preparation** — parse, validate, drug-likeness profile, ADMET, SA score
2. **BRICS one-step retrosynthesis** — instant rule-based disconnection via `mol.brics_fragments()`
3. **AiZynthFinder multi-step** — ML-based route search (real or mock)
4. **Building block scoring** — SA score, Lipinski, Tanimoto to known BB library
5. **Route ranking** — composite feasibility score combining chematic + AiZynthFinder metrics

Key patterns:

```python
import chematic

# 1. Prepare target
mol = chematic.from_smiles("O=C(Nc1ccc(S(N)(=O)=O)cc1)c1ccc(N)cc1")
d = mol.descriptors()
print(f"SA score: {d['sa_score']:.2f}  QED: {d['qed']:.3f}")

# 2. BRICS one-step retrosynthesis
fragments = mol.brics_fragments()        # list[Mol] with dummy attachment points
for frag in fragments:
    fd = frag.descriptors()
    print(f"  SA={fd['sa_score']:.2f}  MW={fd['mw']:.1f}  {frag.smiles}")

# 3. Score building blocks from AiZynthFinder routes
bb = chematic.from_smiles("c1ccc(N)cc1")
sim = chematic.tanimoto(mol.ecfp4(), bb.ecfp4())
print(f"Tanimoto to aniline: {sim:.3f}")

# 4. ADMET filter on proposed building blocks
admet = bb.admet()
print(f"BBB: {admet['bbb_penetrant']}  hERG risk: {admet['herg_risk']:.2f}")
```

---

## What chematic returns (no conversion needed)

| API | Return type | sklearn-compatible? |
|-----|-------------|---------------------|
| `bulk.ecfp4(smiles)` | `(N, 2048)` `uint8` ndarray | Yes |
| `bulk.maccs(smiles)` | `(N, 166)` `uint8` ndarray | Yes |
| `bulk.descriptors(smiles)` | `list[dict]` → `pd.DataFrame` | Yes (after DataFrame) |
| `mol.ecfp4()` | `bytes` (256 bytes = 2048 bits) | Via `np.frombuffer` |
| `tanimoto_matrix(fps_a, fps_b)` | `list[list[float]]` | Via `np.array` |
