# chematic

[English](README.md) | [日本語](README_ja.md)

纯 Rust 实现的化学信息学库，目标是与 RDKit 功能对等，零 C/C++ FFI。

---

## 设计目标

**纯 Rust，零 C/C++ FFI**
不依赖 rdkit-sys，不使用 openbabel 绑定。所有算法均以安全 Rust 实现。

**兼容 WASM，体积轻量**
核心 crate 无需修改即可编译至 `wasm32-unknown-unknown`。二进制体积仅数百 KB，而 C++ FFI 封装通常达数十 MB。

**化学领域专用算法**
不封装通用图形库，而是直接实现化学专用算法：Kekulization（Kekulé化）、Hückel 芳香性、CIP 立体化学、SSSR 环感知。

**可重现性与确定性**
指纹采用固定不变量排序的 FNV-1a 哈希。相同的 SMILES 输入始终产生相同的位串。无随机数，无平台相关行为。

---

## 当前状态

所有阶段已完成。544 个测试，全部通过。

| Crate                 | 说明                                                                               | 测试数 |
|-----------------------|------------------------------------------------------------------------------------|--------|
| `chematic-core`       | Atom、Bond、Molecule、Element、Kekulization（无依赖）                              | 30     |
| `chematic-smiles`     | OpenSMILES 解析器、写入器、规范 SMILES                                             | 52     |
| `chematic-perception` | SSSR（Balducci-Pearlman）、Hückel 芳香性                                           | 14     |
| `chematic-mol`        | MOL/SDF V2000+V3000 解析器与写入器                                                 | 37     |
| `chematic-depict`     | 2D SVG 绘制，CPK 配色，原子/键高亮                                                 | 15     |
| `chematic-chem`       | 分子描述符、BRICS 碎片化、QED、标准化、Murcko 骨架、CIP                            | 216    |
| `chematic-fp`         | ECFP4/6、MACCS 166位、拓扑路径、AtomPair、Torsion FP、Tanimoto/Dice               | 44     |
| `chematic-smarts`     | SMARTS 解析器（递归、价键、杂化），VF2 子图同构，MCS                               | 76     |
| `chematic-3d`         | 3D 坐标生成，PDB/XYZ 文件格式                                                      | 25     |
| `chematic-rxn`        | 反应 SMILES 解析器与写入器                                                         | 15     |
| `chematic-wasm`       | WebAssembly 绑定 — npm：`@kent-tokyo/chematic`                                     | 18     |
| `chematic`            | 带功能标志的伞形 crate（含所有子 crate）                                           | 1      |

```
cargo test --workspace   # 544 个测试，全部通过
```

---

## 快速开始

### 使用伞形 crate

```toml
# Cargo.toml
[dependencies]
chematic = { git = "https://github.com/kent-tokyo/chematic", features = ["smiles", "fp"] }
```

```rust
use chematic::smiles::{parse, canonical_smiles};
use chematic::fp::ecfp4;
```

### 使用单独 crate

```toml
# Cargo.toml
[dependencies]
chematic-smiles     = { git = "https://github.com/kent-tokyo/chematic" }
chematic-perception = { git = "https://github.com/kent-tokyo/chematic" }
chematic-fp         = { git = "https://github.com/kent-tokyo/chematic" }
```

```rust
use chematic_smiles::{parse, canonical_smiles};
use chematic_perception::{find_sssr, assign_aromaticity};
use chematic_fp::{ecfp4, tanimoto_ecfp4};

fn main() {
    let benzene = parse("c1ccccc1").unwrap();
    let toluene = parse("Cc1ccccc1").unwrap();

    // 环感知与芳香性
    let rings = find_sssr(&benzene);
    println!("rings: {}", rings.ring_count()); // 1

    // 指纹相似度
    let sim = tanimoto_ecfp4(&benzene, &toluene);
    println!("Tanimoto(benzene, toluene): {sim:.3}"); // ~0.5

    // 规范 SMILES
    println!("{}", canonical_smiles(&benzene)); // c1ccccc1
}
```

---

## SMARTS 子结构搜索

```rust
use chematic_smiles::parse;
use chematic_smarts::{parse_smarts, find_matches};

let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap(); // 阿司匹林
let query = parse_smarts("[$(C(=O)O)]").unwrap();   // 羧基 / 酯基 C
let matches = find_matches(&query, &mol);
println!("C(=O)O groups: {}", matches.len()); // 2
```

---

## 分子描述符

```rust
use chematic_smiles::parse;
use chematic_chem::{molecular_weight, tpsa, logp_crippen, fsp3, qed, lipinski_passes};

let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
println!("MW:       {:.2}", molecular_weight(&aspirin)); // ~180.16
println!("TPSA:     {:.2}", tpsa(&aspirin));             // ~63.6
println!("LogP:     {:.2}", logp_crippen(&aspirin));     // ~1.2
println!("Fsp3:     {:.3}", fsp3(&aspirin));             // ~0.111
println!("QED:      {:.3}", qed(&aspirin));              // 类药性评分
println!("Lipinski: {}", lipinski_passes(&aspirin));     // true
```

---

## BRICS 碎片化

```rust
use chematic_smiles::parse;
use chematic_chem::brics_fragments;

let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
let frags = brics_fragments(&aspirin);
println!("fragments: {}", frags.len()); // ≥ 2
```

---

## 分子指纹

```rust
use chematic_smiles::parse;
use chematic_fp::{ecfp4, atom_pair_fp, torsion_fp};

let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
let caffeine = parse("Cn1cnc2c1c(=O)n(c(=O)n2C)C").unwrap();

let sim_ecfp4    = ecfp4(&aspirin).tanimoto(&ecfp4(&caffeine));
let sim_atompair = atom_pair_fp(&aspirin).tanimoto(&atom_pair_fp(&caffeine));
let sim_torsion  = torsion_fp(&aspirin).tanimoto(&torsion_fp(&caffeine));
```

---

## 2D 绘制

```rust
use chematic_smiles::parse;
use chematic_depict::depict_svg;

let caffeine = parse("Cn1cnc2c1c(=O)n(c(=O)n2C)C").unwrap();
let svg = depict_svg(&caffeine);
std::fs::write("caffeine.svg", svg).unwrap();
```

### 高亮绘制

```rust
use std::collections::HashSet;
use chematic_smiles::parse;
use chematic_depict::depict_svg_highlighted;

let mol = parse("c1ccncc1").unwrap(); // 吡啶
let n_idx = mol.atoms().find(|(_, a)| a.element.atomic_number() == 7)
               .map(|(i, _)| i).unwrap();
let svg = depict_svg_highlighted(&mol, &HashSet::from([n_idx]), &HashSet::new());
```

---

## JavaScript / TypeScript（WebAssembly）

```sh
npm install @kent-tokyo/chematic
```

```js
import init, { parse_smiles, tanimoto_ecfp4, tanimoto_atom_pair, brics_fragment_count } from '@kent-tokyo/chematic';

await init();

const mol = parse_smiles('CC(=O)Oc1ccccc1C(=O)O'); // 阿司匹林
console.log(mol.molecular_weight()); // ~180.16
console.log(mol.logp_crippen());     // ~1.2
console.log(mol.qed());              // 类药性 [0,1]
console.log(mol.fsp3());             // sp3 碳比例
console.log(brics_fragment_count(mol)); // BRICS 碎片数

const caffeine = parse_smiles('Cn1cnc2c1c(=O)n(c(=O)n2C)C');
console.log(tanimoto_ecfp4(mol, caffeine));    // ECFP4 相似度
console.log(tanimoto_atom_pair(mol, caffeine)); // AtomPair 相似度
```

---

## 与其他化学信息学库的比较

| 功能                             | chematic                | RDKit (rdkit-sys)  | OpenBabel FFI  | chemcore / purr   |
|----------------------------------|-------------------------|--------------------|----------------|-------------------|
| 语言                             | 纯 Rust                 | Rust + C++ FFI     | Rust + C++ FFI | 纯 Rust           |
| WASM 目标                        | 支持                    | 不支持             | 不支持         | 部分支持          |
| 二进制体积（核心）               | ~500 KB                 | ~50 MB             | ~20 MB         | ~200 KB           |
| OpenSMILES 解析器                | 完整                    | 完整               | 完整           | 部分              |
| SMILES 写入 / 规范化             | 支持                    | 支持               | 支持           | 不支持            |
| Kekulization                     | 支持                    | 支持               | 支持           | 不支持            |
| 芳香性感知                       | 支持（Hückel）          | 支持               | 支持           | 部分支持          |
| 环感知（SSSR）                   | 支持                    | 支持               | 支持           | 不支持            |
| SDF/MOL V2000+V3000              | 支持                    | 支持               | 支持           | 不支持            |
| 2D 绘制（SVG，CPK 配色）         | 支持                    | 支持               | 支持           | 不支持            |
| ECFP 指纹                        | 支持（ECFP4/6）         | 支持               | 支持           | 不支持            |
| AtomPair / Torsion 指纹          | 支持                    | 支持               | 支持           | 不支持            |
| MACCS 指纹                       | 支持（166位）           | 支持               | 支持           | 不支持            |
| SMARTS / 子结构搜索              | 支持（VF2 + 递归）      | 支持               | 支持           | 不支持            |
| 分子描述符                       | 支持（MW/LogP/TPSA/Fsp3/QED/…）| 支持      | 支持           | 不支持            |
| BRICS 碎片化                     | 支持                    | 支持               | 不支持         | 不支持            |
| 3D 坐标生成                      | 支持（规则驱动）        | 支持（ETKDG）      | 支持           | 不支持            |
| PDB/XYZ 文件格式                 | 支持                    | 支持               | 支持           | 不支持            |
| CIP 立体化学（R/S、E/Z）         | 支持                    | 支持               | 支持           | 不支持            |
| 力场能量最小化                   | 支持（规则驱动）        | 支持（UFF/MMFF）   | 支持           | 不支持            |
| 反应 SMILES/SMIRKS               | 支持                    | 支持               | 支持           | 不支持            |
| Unsafe Rust                      | 无                      | 大量使用           | 大量使用       | 无                |
| 维护状态（2026）                 | 活跃                    | 活跃               | 最低限度       | 已归档            |

注：
- 二进制体积为估算值，实际取决于启用的功能。
- chemcore 和 purr 已归档；chematic 在功能范围上超越了它们。

---

## 路线图

### 第一阶段 — 基础（已完成）
核心类型、OpenSMILES 解析/写入、Kekulization、规范 SMILES。

### 第二阶段 — 分子感知（已完成）
SSSR、Hückel 芳香性、SDF/MOL V2000+V3000、2D SVG 绘制。

### 第三阶段 — 化学智能（已完成）
描述符（MW、LogP、TPSA、Fsp3、Lipinski）、QED、BRICS 碎片化、
ECFP4/6 指纹、SMARTS+VF2（递归 SMARTS、价键、杂化），
分子标准化、Murcko 骨架、CIP R/S 和 E/Z。

### 第四阶段 — 相似性与搜索（已完成）
MACCS 166位键、拓扑路径 FP、AtomPair FP、Topological Torsion FP、
MCS、互变异构体规范化。

### 第五阶段 — 3D 化学（已完成）
基于规则的 3D 坐标生成、PDB/XYZ 格式、类 UFF 能量最小化。

### 第六阶段 — RDKit 对等（已完成）
反应 SMILES/SMIRKS ✓、带功能标志的伞形 crate ✓、
WASM npm 包 `@kent-tokyo/chematic` ✓、CPK 配色 + 高亮绘制 ✓、
ChEMBL 37 全集验证（2,897,819 个分子，100.000%）✓。

---

## 仓库结构

```
chematic/
├── Cargo.toml               工作区根目录
├── CHANGELOG.md             版本历史
├── crates/
│   ├── chematic-core/       Atom, Bond, Molecule, Element, kekulization
│   ├── chematic-smiles/     OpenSMILES 解析器、写入器、规范 SMILES
│   ├── chematic-perception/ SSSR 环感知、Hückel 芳香性
│   ├── chematic-mol/        MOL/SDF V2000+V3000 解析器与写入器
│   ├── chematic-depict/     2D SVG 绘制引擎（CPK 配色，高亮）
│   ├── chematic-chem/       描述符、BRICS、QED、标准化、骨架
│   ├── chematic-fp/         ECFP4/6、MACCS、路径、AtomPair、Torsion FP
│   ├── chematic-smarts/     SMARTS 解析器 + VF2 子图同构，MCS
│   ├── chematic-3d/         3D 坐标生成，PDB/XYZ 格式
│   ├── chematic-rxn/        反应 SMILES 解析器与写入器
│   └── chematic/            带功能标志的伞形 crate
└── tasks/
    ├── todo.md              详细路线图清单（日语）
    └── lessons.md           开发经验总结
```

---

## 开发命令

```bash
cargo build --workspace      # 构建所有 crate
cargo test --workspace       # 运行所有测试（544 个）
cargo check --workspace      # 仅类型检查，不构建
cargo clippy --workspace     # 代码检查
```

---

## 许可证

可选择 Apache License 2.0 或 MIT License 中的任意一种。
