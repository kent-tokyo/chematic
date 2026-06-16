# chematic cookbook — よくある20タスク

各タスクはコピペで動くコードです。出力例は aspirin (`CC(=O)Oc1ccccc1C(=O)O`) で示しています。

```python
import chematic
```

---

## 1. SMILES から基本物性を取得する

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin

print(mol.mw)              # 180.16
print(mol.logp)            # 1.31
print(mol.tpsa)            # 63.6
print(mol.hbd, mol.hba)   # 1 4
print(mol.qed)             # 0.55
print(mol.formula)         # C9H8O4
```

## 2. ドラッグライクネスフィルターを確認する

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

print(mol.lipinski_passes)  # True
print(mol.veber_passes)     # True
print(mol.ghose_passes)     # True
print(mol.pains_passes)     # True
print(mol.brenk_passes)     # True
```

## 3. pKa を予測する

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

pka = mol.pka()
print(pka["most_acidic"])   # 3.49  (カルボン酸)
print(pka["most_basic"])    # None  (塩基性サイトなし)
```

## 4. ADMET プロファイルを取得する

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

profile = mol.admet()
# {
#   "bbb": False,
#   "bbb_score": -1.3,
#   "caco2": -4.9,         # log Papp (Caco-2)
#   "herg_risk": 0.12,     # 0–1 リスクスコア
#   "cyp3a4_risk": 0.21,
# }
```

## 5. フィンガープリントで類似度を計算する

```python
aspirin  = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")

sim = chematic.tanimoto(aspirin.ecfp4(), ibuprofen.ecfp4())
print(f"Tanimoto (ECFP4): {sim:.3f}")  # 0.148
```

## 6. SDF ファイルを読んで DataFrame を作る

```python
import pandas as pd

records = list(chematic.iter_sdf("library.sdf"))
df = pd.DataFrame({
    "name":  [r.name for r in records],
    "smiles": [r.mol.smiles for r in records],
    "mw":    [r.mol.mw for r in records],
    "logp":  [r.mol.logp for r in records],
    "qed":   [r.mol.qed for r in records],
})
print(df.head())
```

## 7. 大量の SMILES を並列で記述子計算する

```python
import pandas as pd

smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "c1cccnc1", "CCCCCCCC"]

# CPU コア数に応じて自動並列化
df = pd.DataFrame(chematic.bulk.descriptors(smiles_list))
print(df[["mw", "logp", "tpsa", "qed"]].round(2))
```

## 8. ECFP4 行列を numpy で取得する（ML 用途）

```python
import numpy as np

smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]

# shape: (N, 2048), dtype: uint8
X = chematic.bulk.ecfp4(smiles_list)

# scikit-learn で直接使える
from sklearn.ensemble import RandomForestClassifier
# clf.fit(X, y)
```

## 9. Tanimoto 類似度行列を計算する

```python
smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "c1cccnc1"]

# shape: (N, N), dtype: float32
matrix = chematic.bulk.tanimoto(smiles_list, smiles_list)
print(matrix.round(2))
```

## 10. 仮想スクリーニング（類似度サーチ）

```python
library_smiles = [...]  # 数万件でも OK

# ターゲット分子に似た上位候補を探す
query = "CC(=O)Oc1ccccc1C(=O)O"
scores = chematic.bulk.tanimoto_search(query, library_smiles)  # (N,) float32

# 上位10件
top_indices = scores.argsort()[::-1][:10]
for i in top_indices:
    print(f"{library_smiles[i]}: {scores[i]:.3f}")
```

## 11. LSH インデックスで高速近傍探索

```python
# 大規模ライブラリ（数十万件）向け
idx = chematic.SimilarityIndex.from_smiles(library_smiles)

hits = idx.search("CC(=O)Oc1ccccc1C(=O)O", threshold=0.5, k=20)
for mol_idx, score in hits:
    print(f"{library_smiles[mol_idx]}: {score:.3f}")
```

## 12. SMARTS でサブ構造検索

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

# マッチするかどうか
if chematic.smarts_match("[CX3](=O)[OX2H1]", mol):
    print("カルボン酸あり")

# マッチした原子インデックスを取得
matches = chematic.smarts_find("[CX3](=O)[OX2H1]", mol)
print(matches)  # [[7, 8, 9], ...]
```

## 13. 分子を標準化する

```python
# 塩・溶媒の除去、電荷の中和、互変異性体の正規化
mol = chematic.from_smiles("[Na+].[O-]c1ccccc1")
clean = mol.standardize()
print(clean.smiles)  # Oc1ccccc1 （フェノール）

# 個別操作
mol.largest_fragment()   # 最大フラグメントのみ
mol.neutralize()         # 電荷を中和
mol.remove_isotopes()    # 同位体ラベル除去
mol.remove_stereo()      # 立体化学を削除
```

## 14. Murcko スカフォールドを抽出する

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

scaffold = mol.scaffold()
print(scaffold.smiles)         # c1ccc(CC(=O)O)cc1 (Murcko)

generic = mol.generic_scaffold()
print(generic.smiles)          # C1ccc(CC(C)C)cc1 (ジェネリック)
```

## 15. BRICS フラグメント分解

```python
mol = chematic.from_smiles("CC(=O)Nc1ccc(O)cc1")  # paracetamol

frags = mol.brics_fragments()
for f in frags:
    print(f.smiles)
# [NH2+]=[C@@H](C)C  など
```

## 16. タウトマーを列挙・正規化する

```python
mol = chematic.from_smiles("OC1=CC=CC=N1")  # 2-pyridinol

canonical = mol.canonical_tautomer()
print(canonical.smiles)       # O=C1CC=CC=N1

all_tautomers = mol.enumerate_tautomers()
print(len(all_tautomers))     # 複数のタウトマー
```

## 17. 最大共通部分構造 (MCS) を求める

```python
mol1 = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")   # aspirin
mol2 = chematic.from_smiles("CC(=O)Nc1ccc(O)cc1")       # paracetamol

mcs = chematic.find_mcs([mol1, mol2])
if mcs:
    print(mcs.smiles)   # 共通部分構造
```

## 18. SMIRKS で反応を適用する

```python
# フェノールを脱プロトン化
phenol = chematic.from_smiles("c1ccccc1O")
products = chematic.run_smirks("[OH:1]>>[O-:1]", [phenol])

for product_set in products:
    for p in product_set:
        print(p.smiles)  # [O-]c1ccccc1
```

## 19. Jupyter で SVG 描画する

```python
from IPython.display import SVG, display

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

# 単分子
display(SVG(mol.svg()))

# 原子ハイライト（SMARTS マッチ部位）
matches = chematic.smarts_find("[CX3](=O)[OX2H1]", mol)
atoms = [i for match in matches for i in match]
display(SVG(mol.svg_highlighted(atoms, color="#FF6B6B")))

# グリッド表示
mols = [chematic.from_smiles(s) for s in ["CCO", "c1ccccc1", "CC(=O)O", "CCCC"]]
display(SVG(chematic.depict_grid(mols, cols=2)))
```

## 20. 全記述子 + フィルター結果を一括出力する

```python
import pandas as pd

smiles_list = [
    "CC(=O)Oc1ccccc1C(=O)O",   # aspirin
    "CC(C)Cc1ccc(CC(C)C(=O)O)cc1",  # ibuprofen
    "c1ccc2ccccc2c1",            # naphthalene
    "CCCCCCCCCCCCCCCCCC(=O)O",   # stearic acid
]

df = pd.DataFrame(chematic.bulk.descriptors(smiles_list))
# lipinski_passes, pains_passes, brenk_passes なども含む
print(df[["mw", "logp", "tpsa", "qed", "lipinski_passes", "pains_passes"]].round(2))
```

---

## 補足：InChI / InChIKey

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

print(mol.inchi)     # InChI=1S/C9H8O4/...
print(mol.inchikey)  # BSYNRYMUTXBXSQ-UHFFFAOYSA-N

# InChI → Mol
mol2 = chematic.from_inchi(mol.inchi)
assert mol2.smiles == mol.smiles
```
