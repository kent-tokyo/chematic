# chematic

[English](README.md) | [日本語](README_ja.md)

纯 Rust 实现的化学信息学库，目标是与 RDKit 功能对等，**零 C/C++ FFI**。

> **为什么零 C/C++ 如此重要？**
> RDKit.js、Indigo WASM 和 OpenBabel 均使用 Emscripten 编译 C++ 代码，
> 这意味着 **30〜50 MB 的 WASM 包**、复杂的构建工具链和平台相关的构建错误。
> chematic 只需 `wasm-pack build` 即可生成 **〜550 KB 的 WASM 包**，
> 整个依赖树中没有任何 `-sys` crate、`cc` 构建依赖或 `build.rs` C 编译。

---

## 在线演示

**[https://kent-tokyo.github.io/chematic/](https://kent-tokyo.github.io/chematic/)** — 可在浏览器中通过 WebAssembly 运行的交互式演示：描述符计算、类药性规则、相似度比较、3D 查看器、反应方案、SAR 分析。

---

## 设计目标

**纯 Rust，零 C/C++ FFI — 已验证**
不使用 `rdkit-sys`、`openbabel-sys`、`cc` 构建依赖或 `bindgen`。从 SSSR 环感知到 ECFP 指纹再到力场最小化，所有算法均以 100% 安全 Rust 实现。已验证整个依赖树无 FFI。

**WASM 兼容，体积轻量**
所有 crate 无需修改即可编译至 `wasm32-unknown-unknown`。npm 包 `@kent-tokyo/chematic` **〜550 KB**，远小于 C++ FFI 替代方案的 30〜50 MB。无需 cmake、emcc 或 Emscripten 工具链。

**100+ WebAssembly API**
WASM 层提供 100 余个函数，涵盖描述符、指纹、骨架分析、立体异构体枚举、3D 几何、多样性选择、MMP 分析、R 基团分解及分子编辑。带完整 TypeScript 类型定义。

**化学领域专用算法**
不封装通用图形库，而是直接实现化学专用算法：Kekulization、Hückel 芳香性、CIP 立体化学、SSSR 环感知、Gasteiger 电荷、MaxMin/Butina 多样性筛选。

**可重现性与确定性**
指纹采用固定不变量排序的 FNV-1a 哈希。相同的 SMILES 输入始终产生相同的位串。无随机数，无平台相关行为。

---

## 当前状态

所有阶段已完成 + 第 4 部分（WASM、API 改进）+ Sprint v0.1.33（CXSMILES/CXSMARTS + 审计）+ Sprint v0.1.34（InChI 环闭合 + 立体层）+ Sprint v0.1.35（wasmBridge 支持）+ Sprint v0.1.36（Issue #1 审计：BUG-2/3/4 修复）+ Sprint v0.1.37（mol_transforms API + 随机 SMILES）+ **Sprint v0.1.69–v0.1.74（RDKit 功能缺口分析：6 项功能实现）+ v0.1.88–v0.1.89（功能缺口分析 89% 完成：A1–A6、B1–B2 实现）+ v0.1.91–v0.1.94（功能缺口分析 100% 完成：A1–A5、B3 实现）+ v0.1.95–v0.1.99（真实算法 + MMFF94 BCI + MinHash LSH + WASM API 扩展 + LogP 烯酮）。** 1,681 个库测试，全部通过。零 C/C++ 依赖。

最新版本：**v0.1.99**（2026-06-13）— LogP 烯酮系 vinyl C 区分

| Crate                 | 说明                                                                                                   | 测试数 |
|-----------------------|--------------------------------------------------------------------------------------------------------|--------|
| `chematic-core`       | Atom、Bond、Molecule、Element、Kekulization（无依赖）；可变 API、`fragments`、`validate_valence`、`formula_with_isotopes`；`StereoGroup`/`StereoGroupKind` | 48     |
| `chematic-smiles`     | OpenSMILES 解析器、写入器、规范 SMILES、**CXSMILES 元数据支持**                                      | 57     |
| `chematic-perception` | SSSR、Hückel 芳香性 + 反芳香性（4n+2 规则）、`apply_aromaticity`/`aromatize`/`kekulize_inplace`、`assign_stereo_from_2d`、`assign_ez_from_2d`、`cip_ez_descriptor` | 34     |
| `chematic-mol`        | MOL/SDF V2000+V3000（读写含 2D 坐标）、CML（读写）、CDXML（读）；`SdfRecord`（含坐标+属性）、MDL RXN V2000 读写；V3000 立体基团 COLLECTION 读写 | 63     |
| `chematic-depict`     | 2D SVG 绘制（CPK 配色、高亮、网格）、`detect_crossings`/`render_svg_with_metadata`、反应 SVG；Y 坐标系文档已更新 | 43     |
| `chematic-chem`       | 70+ 描述符、互变异构体评分、骨架网络、BRICS、QED、标准化、分子哈希、立体化学、`parse_condensed`、CIP、IFG、Gasteiger、VSA（EState+Labute）、`isotope_distribution`、`num_amide_bonds`、`num_ester_bonds` | 375    |
| `chematic-fp`         | ECFP2/4/6、FCFP4/6、MACCS 166位、TopoPF、AtomPair、Torsion FP — bitvec + Tanimoto/Dice               | 50     |
| `chematic-smarts`     | SMARTS、VF2、MCS（支持 `match_chiral_tag` 手性匹配、原子/键比较模式）                                 | 87     |
| `chematic-3d`         | 3D 坐标生成、力场最小化、形状描述符、ConformerEnsemble（含 RMSD 剪枝）、PDB/XYZ 格式                 | 147    |
| `chematic-rxn`        | 反应 SMILES/SMIRKS、`find_reaction_center` — `run_reactants`（含生成物价键验证）                      | 30     |
| `chematic-inchi`      | InChI/InChIKey 生成；formula/connectivity/hydrogen/stereo/charge/isotope 层；环闭合支持               | 28     |
| `chematic-wasm`       | **110+ WASM 导出** — npm：`@kent-tokyo/chematic` v0.1.95；InChI API + 立体反演                         | 175    |
| `chematic-iupac`      | 本地 IUPAC 命名（纯 Rust·离线）— 烷烃/环烷烃、烯烃/炔烃、醇、胺、卤代烃、醛、酮、羧酸、酯、酰胺、苯、芳香杂环 | 14     |
| `chematic`            | 带功能标志的伞形 crate（含所有子 crate，含 `iupac`、`inchi`）                                         | 1      |

```
cargo test --workspace --lib   # 1,681 个库测试，全部通过
```

---

## 快速开始

### 使用伞形 crate

```toml
# Cargo.toml
[dependencies]
chematic = { version = "0.1.94", features = ["smiles", "fp", "chem", "mol", "depict"] }
```

### 使用单独 crate

```toml
# Cargo.toml
[dependencies]
chematic-smiles     = "0.1.94"
chematic-perception = "0.1.94"
chematic-fp         = "0.1.94"
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

## CXSMILES 保留元数据

```rust
use chematic_smiles::parse_cxsmiles;

let cx = parse_cxsmiles("CCO |$ethanol$,atomProp:1.role.acceptor,^2:0|").unwrap();
// cx.atom_labels: ["ethanol"]
// cx.atom_props: [(atom: 1, key: "role", value: "acceptor")]
// cx.atom_radicals: [None, 2, None]

// 保持 CX 信息进行规范化
let canonical = chematic_smiles::write_cxsmiles(&cx);
println!("{}", canonical); // CCO |$ethanol$,atomProp:1.role.acceptor,^2:0|
```

---

## 标准化管道与审计报告

```rust
use chematic_chem::{StandardizationPipeline, StandardizeOptions};

let opts = StandardizeOptions {
    largest_fragment_only: true,
    neutralize_charges: true,
    remove_explicit_h: false,
    canonical_tautomer: false,
};

let pipeline = StandardizationPipeline::new(opts);
let (standardized, report) = pipeline.run(&mol);

// 检查整体状态（Unchanged / Modified / CompletedWithWarnings）
println!("Status: {:?}", report.status);

// 跟踪每一步的变化
for step in &report.steps {
    if step.changed {
        println!("  {}: {} atoms → {} atoms",
            step.step.as_str(),
            step.before.atoms,
            step.after.atoms
        );
    }
}

// 检查警告（金属键、价电子错误等）
for warning in &report.warnings {
    println!("⚠️  {}: {}", warning.code, warning.message);
}
```

---

## JavaScript / TypeScript（WebAssembly）

> **〜550 KB，零 C/C++ 依赖。** 支持浏览器和 Node.js。
> 相比之下，RDKit.js 通过 Emscripten 构建约 30 MB。

```sh
npm install @kent-tokyo/chematic
```

```js
import init, {
  parse_smiles, canonical_tautomer, murcko_scaffold,
  tanimoto_ecfp4, tanimoto_ecfp6, tanimoto_maccs,
  brics_fragments_json, mcs_smiles_json,
  get_descriptors_json, enumerate_stereo_isomers_json,
  mmp_pairs_json, rgroup_decompose_json,
  mol_with_atom_added, mol_with_atom_charge, mol_with_atom_element,
  depict_data_json, cpk_color,
} from '@kent-tokyo/chematic';

await init();

// ── 解析与描述符 ─────────────────────────────────────────────
const mol = parse_smiles('CC(=O)Oc1ccccc1C(=O)O'); // 阿司匹林
console.log(mol.molecular_weight()); // ~180.16
console.log(mol.qed());              // 类药性 [0,1]

// 一次性获取所有描述符（JSON 对象）
const desc = JSON.parse(get_descriptors_json(mol));
console.log(desc.mw, desc.tpsa, desc.logP);

// ── 立体异构体枚举 ───────────────────────────────────────────
const isomers = JSON.parse(enumerate_stereo_isomers_json(parse_smiles('C(F)(Cl)Br')));

// ── SAR 分析 ─────────────────────────────────────────────────
const smiles_json = '["CCc1ccccc1","CCCc1ccccc1","CCCCc1ccccc1"]';
const pairs = JSON.parse(mmp_pairs_json(smiles_json));
const rgroups = JSON.parse(rgroup_decompose_json(smiles_json, 'c1ccc(*)cc1'));

// ── 分子编辑 API ─────────────────────────────────────────────
const mol2 = mol_with_atom_added(mol, 'N');
const mol3 = mol_with_atom_charge(mol, 0, 1);    // 将原子 0 的电荷设为 +1
const mol4 = mol_with_atom_element(mol, 0, 'O'); // 将原子 0 的元素改为 O
```

---

## 与其他化学信息学库的比较

| 功能                                      | **chematic**             | RDKit.js (WASM)   | OCL.js | Indigo WASM |
|-------------------------------------------|--------------------------|-------------------|--------|-------------|
| **C/C++ 依赖**                            | **零 — 纯 Rust**         | C++（Emscripten）| △      | C++（Emscripten）|
| **WASM 二进制体积**                       | **〜550 KB**             | 〜30 MB           | 〜5 MB | 〜10 MB     |
| 描述符丰富度                              | **◎ 40+**                | ○ 〜30            | △      | △           |
| 指纹种类与设置自由度                      | **◎ 7 种 bitvec + 相似度**| ◎                | ○      | △           |
| 立体化学（CIP + 枚举）                    | **◎**                    | ○                 | ○      | △           |
| 3D + 构象管理                             | **◎**                    | ○                 | △      | △           |
| 多样性筛选（MaxMin/Butina）               | **◎**                    | ○                 | ✗      | ✗           |
| MMP 分析                                  | ✓                        | ✓                 | ✗      | ✗           |
| R 基团分解                                | ✓                        | ✓                 | ✗      | ✗           |
| 分子编辑 API                              | **◎ with_atom_* 系列**   | ○                 | ○      | ○           |
| CML 读写                                  | ✓                        | ✓                 | ✓      | ✓           |
| CDXML 读取（多分子片段 + 立体化学）       | ✓                        | ✓                 | ✓      | ✓           |
| InChI / InChIKey                          | ✗（依赖 C 库）           | ✓                 | ✓      | ✓           |
| Unsafe Rust                               | **无**                   | —                 | —      | —           |

---

## 路线图

### 第一阶段〜第六阶段（已完成）
基础、分子感知、化学智能、相似性搜索、3D 化学、生态系统。

### 第七阶段（已完成）
扩展描述符、多样性、SA 评分、EState、IFG、Gasteiger、VSA。

### 第八阶段（v0.1.20〜v0.1.22，已完成）
100+ WASM 导出、CML/CDXML、Mutable Molecule API、DepictData、MMP、R 基团、ConformerEnsemble、SDF/V3000 写入、MCS 环感知约束。

### 第十五阶段（v0.1.29〜32，已完成）
可变 `Molecule`（`add/remove_atom/bond`、`fragments`、`is_connected`），
`assign_stereo_from_2d`（楔形键→R/S），`aromatize`/`kekulize_inplace`，
`depict_reaction_svg`，`SdfRecord`（含坐标+属性），MDL RXN V2000 读写，
`expand_abbreviation`（30 种缩写），`formula_with_isotopes`。

### 第十四阶段（v0.1.28，已完成）
`xlogp3()`（Cheng 2007 原子类型），`chematic-iupac`（纯 Rust 离线 IUPAC 命名），
`BricsConfig { min_fragment_size }`，`MatchConfig { max_matches }`，
`McsConfig { atom_compare, bond_compare }` 支持异环 scaffold hopping。

### 第十六阶段（v0.1.27，已完成）
`assign_ez_from_2d(mol, coords)` / `cip_ez_descriptor(mol, bond_idx, coords)` — 通过 2D 坐标叉积计算 E/Z 双键立体化学（无需楔形键，采用 1-sphere CIP 优先级）。
`StereoGroup` / `StereoGroupKind`（Absolute / Or / And）加入 `chematic-core`；`Molecule` 新增 `stereo_groups` 字段；
V3000 MOL 解析器与写入器支持 `BEGIN COLLECTION / MDLV30/STEABS / MDLV30/STEOR<n> / MDLV30/STEAND<n>`。
`isotope_distribution(mol, resolution) -> Vec<(f64, f64)>` — 基于卷积的同位素包络计算（支持显式同位素标记优先，H/C/N/O/S/Cl/Br 等 14+ 元素）。

### 第十三阶段（v0.1.27，已完成）
`MolMetadata::default().with_name("阿司匹林").with_comment("...")` — MOL/SDF 元数据 fluent builder。

### 第十二阶段（v0.1.26，已完成）
`atom_color_rgb(atomic_number: u8) -> [u8; 3]` — 无需解析十六进制字符串，直接获取 CPK 颜色 RGB 字节三元组。

### 当前 main（v0.1.94 之后）
`chematic-iupac`：支持带位置编号的酮、羧酸、酯、第一/第二酰胺、苯及常见芳香杂环。CI 中 `cargo clippy --workspace -- -D warnings` 已通过当前 stable Clippy。

### 第十一阶段（v0.1.94，已完成）
`suggest_bond_direction(mol, atom, layout) -> f64`（弧度）：基于 sp2/sp3 角偏移 + 最大最小角分离选择的化学自然新键方向建议。导出 `BOND_LEN` 常量。

### 第十阶段（v0.1.24，已完成）
`validate_valence(mol) -> Vec<ValenceError>` 公开 API（chematic-core + chematic-perception 重导出），`run_reactants` 自动过滤含过价原子的生成物集合。

### 第九阶段（v0.1.23，已完成）
`Element::vdw_radius()` / `covalent_radius()`（Bondi/Alvarez 表，118 个元素），
`Molecule::implicit_hydrogen_count()` / `total_formula()`（含隐式 H 的 Hill 式），
`apply_aromaticity()`（Kekulize 分子 → 芳香性标志应用），
`with_atom_aromatic()` / `with_bond_order()` 不可变更新 API，
`minimize_uff()` 别名（UFF 力场最小化可发现性提升）。

---

## 仓库结构

```
chematic/
├── Cargo.toml               工作区根目录
├── CHANGELOG.md             版本历史
├── crates/
│   ├── chematic-core/       Atom, Bond, Molecule, Element, Kekulization
│   ├── chematic-smiles/     OpenSMILES 解析器、写入器、规范 SMILES
│   ├── chematic-perception/ SSSR 环感知、Hückel 芳香性
│   ├── chematic-mol/        MOL/SDF V2000+V3000、CML、CDXML
│   ├── chematic-depict/     2D SVG 绘制引擎（CPK 配色，DepictData）
│   ├── chematic-chem/       描述符、BRICS、QED、MMP、标准化、CIP
│   ├── chematic-fp/         ECFP4/6、FCFP4/6、MACCS、AtomPair、Torsion FP
│   ├── chematic-smarts/     SMARTS 解析器 + VF2 子图同构，MCS
│   ├── chematic-3d/         3D 坐标生成、ConformerEnsemble、PDB/XYZ 格式
│   ├── chematic-rxn/        反应 SMILES/SMIRKS
│   └── chematic/            带功能标志的伞形 crate
└── tasks/
    ├── todo.md              详细路线图清单（日语）
    └── lessons.md           开发经验总结
```

---

## 许可证

可选择 Apache License 2.0 或 MIT License 中的任意一种。
