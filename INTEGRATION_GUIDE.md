# chematic v0.1.32 統合ガイド

## 概要

このガイドでは、chematic v0.1.32 の新機能（Section 4 + Step 2&3）をアプリケーションに統合するための手順を説明します。

---

## Section 4: WASM & API 改善

### 1. fastrand js feature による WASM RNG シード修正

**問題**: v0.1.30 では MD シミュレーションの初期速度が WASM で毎回同じシード（0x4d595df4d0f33173）で初期化されていました。

**修正内容**: `crates/chematic-3d/Cargo.toml` に target 条件付き依存を追加し、WASM ターゲット向けに `fastrand` の `js` feature を有効化。

**効果**:
- WASM MD シミュレーション が暗号学的ランダム性を使用
- 物理的に意味のある非決定的軌跡を生成
- ネイティブビルドは影響なし

---

### 2. V3000 座標復元機能

**新関数**: `parse_mol_v3000_with_coords()`

```rust
use chematic_mol::parse_mol_v3000_with_coords;

let v3000_block = "..."; // V3000 MOL string
let (mol, metadata, coords) = parse_mol_v3000_with_coords(v3000_block)?;

// coords: Vec<(f64, f64)> — atom[i] の 2D 座標（x, y）
println!("Atom 0: ({:.2}, {:.2})", coords[0].0, coords[0].1);
```

**メリット**:
- V2000 の `parse_mol_with_coords()` と同じ API
- 2D 座標のラウンドトリップ保持
- SDF レンダリング・レイアウト作成に活用可能

---

### 3. Y座標系の仕様明確化

**重要**: 座標系がファイル形式・処理エンジンで異なります。

| 機能 | 座標系 | 説明 |
|------|--------|------|
| `compute_layout()` | SVG Y-down | ブラウザ・SVG 標準（Y は下方向） |
| `parse_cml()` | 化学 Y-up | CML 標準（Y は上方向） |
| `parse_cdxml()` | ChemDraw Y-down | ChemDraw 標準（Y は下方向） |

**変換ルール**:
```rust
// CML → SVG レンダリング
let (mol, cml_coords) = parse_cml(cml_str)?;
let svg_coords: Vec<(f64, f64)> = cml_coords.iter()
    .map(|(x, y)| (x, -y))  // Y を反転
    .collect();
```

---

### 4. エラー型の Display + Error trait

**対応済み（13 型）**:
- SmartsError, ValenceError, StereoError (Display + Error)
- CmlError, CdxmlError, Mol2Error, RxnParseError
- MolError, IupacError, ConformerError, RxnError, TransformError

```rust
use std::error::Error;

fn process() -> Result<(), Box<dyn Error>> {
    let mol = parse_smiles("invalid")?;
    Ok(())
}
```

---

## Step 2: 3D 制約充足 (Distance Geometry Constraints)

### API

**`build_constraints(mol) -> ConstraintSet`**
```rust
let constraints = build_constraints(&mol);
```

**`satisfy_constraints(coords, mol, constraints, max_iterations) -> Coords3D`**
```rust
let projected = satisfy_constraints(&coords, &mol, &constraints, 20);
```

**`generate_and_minimize_constrained(mol) -> Coords3D`**
```rust
let optimized = generate_and_minimize_constrained(&mol);
```

### パフォーマンス

| 分子 | 原子数 | 実行時間 |
|------|--------|--------|
| Benzene | 6 | ~150 µs |
| Naphthalene | 10 | ~400 µs |
| Caffeine | 14 | ~700 µs |

---

## Step 3: 芳香族性モデル厳密化

### API

**`ring_classifications(mol) -> Vec<RingAromaticity>`**
```rust
let model = assign_aromaticity(&mol);
let classifications = model.ring_classifications(&mol);
```

**`antiaromatic_rings(mol) -> Vec<Vec<AtomIdx>>`**
```rust
let antiaromatic = model.antiaromatic_rings(&mol);
```

**`has_antiaromaticity(mol) -> bool`**
```rust
if model.has_antiaromaticity(&mol) {
    println!("Molecule contains antiaromatic rings");
}
```

---

## WASM 統合例

```javascript
import init, {
    molecule_report_json,
    generate_3d_optimized_pdb
} from '@kent-tokyo/chematic';

await init();

const report = JSON.parse(
    molecule_report_json("CC(=O)Oc1ccccc1C(=O)O")
);
console.log(report.descriptors.logp);

const pdb = generate_3d_optimized_pdb("c1ccccc1");
console.log(pdb);
```

---

## 参考資料

- **GitHub**: https://github.com/kent-tokyo/chematic
- **npm**: https://www.npmjs.com/package/chematic-wasm@0.1.32
- **Docs**: https://docs.rs/chematic/
- **Demo**: https://kent-tokyo.github.io/chematic/

