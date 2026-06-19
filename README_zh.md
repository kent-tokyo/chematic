# chematic

[English](README.md) | [日本語](README_ja.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/chematic.svg)](https://crates.io/crates/chematic)
[![PyPI](https://img.shields.io/pypi/v/chematic.svg)](https://pypi.org/project/chematic/)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic.svg)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![ライセンス](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](LICENSE-MIT)
[![Docs](https://img.shields.io/badge/docs-site-blue)](https://kent-tokyo.github.io/chematic/getting_started/installation/)
[![演示](https://img.shields.io/badge/demo-live-brightgreen)](https://kent-tokyo.github.io/chematic/playground/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/kent-tokyo/chematic/blob/main/notebooks/quickstart.ipynb)

纯 Rust 实现的化学信息学库，目标是与 RDKit 功能对等，**默认零 C/C++ FFI**。

> **为什么零 C/C++ 如此重要？**
> RDKit.js、Indigo WASM 和 OpenBabel 均使用 Emscripten 编译 C++ 代码，
> 这意味着 **30〜50 MB 的 WASM 包**、复杂的构建工具链和平台相关的构建错误。
> chematic 只需 `wasm-pack build` 即可生成 **〜550 KB 的 WASM 包**，
> 默认依赖树中没有任何 `-sys` crate、`cc` 构建依赖或 `build.rs` C 编译。
> *（例外：`native-inchi` feature 是唯一的可选项，需要 C 编译器，WASM 构建不受影响。）*

---

## 在线演示

**[https://kent-tokyo.github.io/chematic/playground/](https://kent-tokyo.github.io/chematic/playground/)** — 可在浏览器中通过 WebAssembly 运行的交互式演示：描述符计算、类药性规则、相似度比较、3D 查看器、反应方案、SAR 分析。

---

## 设计目标

**纯 Rust，零 C/C++ FFI — 默认构建已验证**
不使用 `rdkit-sys`、`openbabel-sys`、`bindgen`。从 SSSR 环感知到 ECFP 指纹再到力场最小化，所有算法均以 100% 安全 Rust 实现。已验证整个默认依赖树无 FFI。

> **可选例外**：`chematic-inchi` 的 `native-inchi` feature 可链接 vendored IUPAC InChI C 库 (v1.07.5)，生成与 IUPAC 参考实现逐位一致的标准 InChI/InChIKey。需要 C 编译器，完全 opt-in，默认构建仍保持零 FFI。

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

所有阶段已完成 + **v0.3.x 系列（超越所有主要竞争库）**：MCP 服务器（AI 代理集成）、pKa 预测（15 条 SMARTS 规则）、ADMET 概况（BBB/Caco-2/hERG/CYP3A4）、IUPAC 25+ 类、WASM pKa/ADMET 绑定、criterion 性能基准测试。**1,991 个测试，全部通过。零 C/C++ 依赖（默认构建）。**

最新版本：**v0.3.2**（2026-06-15）— v0.3.0: MCP+pKa+ADMET | v0.3.1: WASM 绑定 | v0.3.2: criterion 基准

| Crate                 | 说明                                                                                                   | 测试数 |
|-----------------------|--------------------------------------------------------------------------------------------------------|--------|
| `chematic-core`       | Atom、Bond、Molecule、Element、Kekulization（无依赖）；可变 API、`fragments`、`validate_valence`、`formula_with_isotopes`；`StereoGroup`/`StereoGroupKind` | 48     |
| `chematic-smiles`     | OpenSMILES 解析器、写入器、规范 SMILES、**CXSMILES 元数据支持**                                      | 57     |
| `chematic-perception` | SSSR、Hückel 芳香性 + 反芳香性（4n+2 规则）、`apply_aromaticity`/`aromatize`/`kekulize_inplace`、`assign_stereo_from_2d`、`assign_ez_from_2d`、`cip_ez_descriptor` | 34     |
| `chematic-mol`        | MOL/SDF V2000+V3000（读写含 2D 坐标）、CML（读写）、CDXML（读）；`SdfRecord`（含坐标+属性）、MDL RXN V2000 读写；V3000 立体基团 COLLECTION 读写 | 63     |
| `chematic-depict`     | 2D SVG 绘制（CPK 配色、高亮、网格）、`detect_crossings`/`render_svg_with_metadata`、反应 SVG；Y 坐标系文档已更新 | 43     |
| `chematic-chem`       | 70+ 描述符、互变异构体、骨架、BRICS、QED、标准化；**pKa 预测**（15 条 SMARTS 规则）；**ADMET 概况**（BBB/Caco-2/hERG/CYP3A4）；**HBA 与 RDKit 一致率 99.98%**（5,000 分子基准） | 496    |
| `chematic-fp`         | ECFP2/4/6、FCFP4/6、MACCS、TopoPF、AtomPair、Torsion、Layered、Pattern、Pharmacophore、Reaction、**MAP4** — Tanimoto/Dice | 55     |
| `chematic-ff`         | **MMFF94 全 7 能量项**（Halgren 1996）：OOP（117 条）+ STRE-BEN（282 条）；L-BFGS；DREIDING | 98     |
| `chematic-smarts`     | SMARTS、VF2、MCS；**SmartsCache**（LRU 5–20×）；**named_pattern()** 库（20 种模式） | 87     |
| `chematic-3d`         | 3D 坐标生成、ETKDG KB（20+ 模式）、力场最小化、形状描述符、ConformerEnsemble、PDB/XYZ | 147    |
| `chematic-rxn`        | 反应 SMILES/SMIRKS、`find_reaction_center`、`run_reactants` | 30     |
| `chematic-inchi`      | InChI/InChIKey：纯 Rust 近似（WASM 兼容）**+ `native-inchi` feature 提供 IUPAC 标准**（vendored C 库 1.07.5，逐位一致）；**parse_inchi** 读取 | 28 (+14*)   |
| `chematic-wasm`       | **130+ WASM 导出** — npm：`@kent-tokyo/chematic` v0.3.2；**pKa/ADMET/BBB/Caco-2/hERG/CYP3A4** WASM API | 209    |
| `chematic-iupac`      | 本地 IUPAC 命名（纯 Rust·离线）— **25+ 化合物类**：烷烃、环烷烃、醇、胺、卤代烃、酮、酸、酯、酰胺、**哌啶、吗啉、哌嗪、萘、硫醚** | 45     |
| `chematic-mcp`        | **MCP（模型上下文协议）服务器** — AI 代理集成；**15 个工具**：parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, name_to_smiles | 28     |
| `chematic`            | 带功能标志的伞形 crate                                                                                   | 1      |

```
cargo test --workspace --lib --quiet                                               # 1,991 个库测试，全部通过
cargo test -p chematic-inchi --features native-inchi --test standard_inchi         # +14 IUPAC 标准 InChI 集成测试
```

---

## 快速开始

### 使用伞形 crate

```toml
# Cargo.toml
[dependencies]
chematic = { version = "0.2.11", features = ["smiles", "fp", "chem", "mol", "depict"] }
```

### 使用单独 crate

```toml
# Cargo.toml
[dependencies]
chematic-smiles     = "0.2.11"
chematic-perception = "0.2.11"
chematic-fp         = "0.2.11"
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
| **C/C++ 依赖**                            | **零（默认）**†          | C++（Emscripten）| △      | C++（Emscripten）|
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
| InChI / InChIKey                          | ✓ `native-inchi` feature（IUPAC C 库 1.07.5）| ✓       | ✓      | ✓           |
| Unsafe Rust                               | **无**                   | —                 | —      | —           |

† 仅限默认构建。`native-inchi` feature 是可选例外，需要 C 编译器（vendored IUPAC InChI C 库 1.07.5）。其余所有 crate 保持零 FFI。

---

## 近期更新

**v0.4.5**（2026-06-19）：Kekulization blossom 算法、E/Z 立体化学、6 个新 MCP 工具、BOILED-Egg
- **Kekulization 4-pass + Edmonds blossom**：5,000 分子语料库中**仅 2 个**失败（硼芳香环、纯 H₂）。
- **E/Z 立体化学**：SMILES 解析器精确读写 E/Z 双键。
- **MCP 新增 6 个工具**：`pains_check`, `brenk_check`, `sa_score`, `admet_profile`, `boiled_egg`, `lipinski_check`，共计 15 个工具。
- **BOILED-Egg**：在 LogP vs TPSA 空间可视化 BBB 渗透性 / GI 吸收的过滤器实现。

**v0.3.2–v0.3.0**：criterion 基准测试、WASM pKa/ADMET 绑定、MCP 服务器 + pKa + ADMET

**v0.2.x**：MMFF94 全 7 项、MAP4 指纹、SMARTS 缓存

**v0.1.x**：核心基础 — SSSR、Kekulization、CIP、3D 几何、WASM API

---

## 已知限制

### Kekulization（5,000 分子中**仅 2 个**失败 — 基本解决）

`chematic-core` 的 Kekulé 赋值采用 4-pass 策略：

- **Pass 1/2**：BFS 增广路径（升序 / 降序）。
- **Pass 3**：桥头 N 排除 — 位于环连接处的 N 原子（芳香度 ≥ 3）提供孤对电子而非占据双键，剩余 C 原子在二部子图上匹配。修复 indolizine 类系统（语料库约 109 例）。
- **Pass 4**：Edmonds' blossom 算法（O(n²m)）— 处理含奇数环的非二部 C 芳香子图（如 corannulene C₂₀H₁₀）。修复剩余的复杂多环系统。

在 5,000 分子语料库（issue #11）中，经上述修复后 Kekulization 仍失败的**仅 2 个**：

| 类别 | 数量 | 示例 |
|---|---|---|
| 硼芳香环 | 1 | `b1ccccn1` |
| 纯 H₂（无重原子） | 1 | `[H][H]` |

**影响**：明确返回 `KekuleError`，不产生无声的错误输出。

### 芳香性模型（Hückel vs RDKit）

chematic 使用 **Hückel 4n+2 规则独立应用于每个 SSSR 环**，而 RDKit 使用更复杂的稠合环电子离域模型。差异在 N-杂环（吡啶酮、喹啉酮、吲哚嗪）中最为明显。

**5,000 分子语料库（issue #12）当前状态：**

| 特征 | issue #12 关闭时 | 当前 | 状态 |
|---|---|---|---|
| `[nH]` SMARTS 匹配 | 67% | **100% recall / 99.8% precision** | 已解决 |
| HBA 计数 | 87.7% | **99.98%**（4,999 / 5,000） | 已解决 |
| 芳香环计数 | 92.6% | **95.6%**（4,778 / 5,000） | 已改善 |

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
```

---

## 许可证

可选择 Apache License 2.0 或 MIT License 中的任意一种。
