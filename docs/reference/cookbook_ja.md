# chematic cookbook

Pythonでよく使う処理を、短い実行例にまとめています。`pip install chematic`でインストールしてから始めます。

```python
import chematic

aspirin = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
```

## 分子を調べる

```python
print(aspirin.smiles)
print(aspirin.formula, aspirin.mw)
print(aspirin.logp, aspirin.tpsa)
print(aspirin.hbd, aspirin.hba, aspirin.qed)

print(aspirin.lipinski_passes)
print(aspirin.veber_passes)
print(aspirin.pains_passes)
print(aspirin.brenk_passes)
```

`pka()`と`admet()`はルールベースのスクリーニング機能です。実験値や臨床結果を予測するものではありません。

```python
print(aspirin.pka())
print(aspirin.admet())
print(aspirin.describe())
```

## 記述子とフィンガープリントを一括計算する

```python
import pandas as pd

smiles = ["CCO", "c1ccccc1", "CC(=O)O", "c1cccnc1"]
df = pd.DataFrame(chematic.bulk.descriptors(smiles))
print(df[["mw", "logp", "tpsa", "qed"]])

# shape (N, 2048)、dtype uint8
X = chematic.bulk.ecfp4(smiles)

# RDKitとのビット一致が必要な場合は、互換APIを明示して使う
morgan = aspirin.rdkit_ecfp4_numpy()
avalon = aspirin.rdkit_avalon_fp_numpy()
```

## 類似度検索を行う

```python
ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")
score = chematic.tanimoto(aspirin.ecfp4(), ibuprofen.ecfp4())

matrix = chematic.bulk.tanimoto(smiles, smiles)
scores = chematic.bulk.tanimoto_search(aspirin.smiles, smiles)

index = chematic.SimilarityIndex.from_smiles(smiles)
hits = index.search(aspirin.smiles, threshold=0.3, k=10)
```

## SMARTS検索と構造操作

```python
pattern = "[CX3](=O)[OX2H1]"
print(chematic.smarts_match(pattern, aspirin))
print(chematic.smarts_find(pattern, aspirin))

salt = chematic.from_smiles("[Na+].[O-]c1ccccc1")
print(salt.standardize().smiles)
print(aspirin.scaffold().smiles)
print(aspirin.generic_scaffold().smiles)

tautomers = chematic.from_smiles("OC1=CC=CC=N1").enumerate_tautomers()
canonical = chematic.from_smiles("OC1=CC=CC=N1").canonical_tautomer()

for fragment in chematic.from_smiles("CC(=O)Nc1ccc(O)cc1").brics_fragments():
    print(fragment.smiles)
```

`standardize()`は、文書化されたchematicの処理順を使います。RDKitとの置換互換性が必要な場合は、すべて同じと仮定せず、[互換性の範囲](compatibility-scope.md)を確認してください。

## MCS・R-group・反応

```python
paracetamol = chematic.from_smiles("CC(=O)Nc1ccc(O)cc1")
mcs = chematic.find_mcs([aspirin, paracetamol])
print(mcs.smiles if mcs else None)

mols = [
    chematic.from_smiles("Cc1ccc(CC)cc1"),
    chematic.from_smiles("CCc1ccc(N)cc1"),
]
columns = chematic.rgroup_decompose_columns(
    "c1cc([*:1])ccc1[*:2]", mols
)
print(columns["R1"], columns["R2"])

phenol = chematic.from_smiles("c1ccccc1O")
products = chematic.run_smirks("[OH:1]>>[O-:1]", [phenol])
print([[mol.smiles for mol in product_set] for product_set in products])
```

R-group分解は、保証範囲を限定しています。曖昧なラベルや未対応の指定を推測で補わず、明示的なエラーとして返します。

## SDFを逐次処理し、レポートを作る

```python
records = chematic.iter_sdf("library.sdf")
rows = [
    {"name": r.name, "smiles": r.mol.smiles, "mw": r.mol.mw}
    for r in records
]
df = pd.DataFrame(rows)

report = chematic.report(
    [aspirin, ibuprofen],
    names=["aspirin", "ibuprofen"],
    title="Compound review",
)
report.save("report.html")
```

`iter_sdf()`は、ファイル全体を一度にメモリへ読み込みません。情報損失の報告が必要な場合は、checked APIまたはsource-preserving APIを使います。詳しくは[エラーと制限](error-and-limits.md)を参照してください。

## Jupyterで描画する

```python
from IPython.display import SVG, display

display(SVG(aspirin.svg()))

matches = chematic.smarts_find("[CX3](=O)[OX2H1]", aspirin)
atoms = [atom for match in matches for atom in match]
display(SVG(aspirin.svg_highlighted(atoms, color="#FF6B6B")))

display(SVG(chematic.depict_grid([aspirin, ibuprofen], cols=2)))
```

## InChIを変換する

```python
print(aspirin.inchi)
print(aspirin.inchikey)

roundtrip = chematic.from_inchi(aspirin.inchi)
assert roundtrip.smiles == aspirin.smiles
```

API全体は[Python API](language-bindings.md)、RDKitからの移行は[RDKit対応表](rdkit_cheatsheet.md)を参照してください。
