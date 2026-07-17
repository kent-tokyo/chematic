# chematic

[English](README.md) | [日本語](README_ja.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![docs.rs](https://docs.rs/chematic/badge.svg)](https://docs.rs/chematic)

![Pure Rust](https://img.shields.io/badge/Pure%20Rust-zero%20C%2B%2B-orange?logo=rust)
![WASM](https://img.shields.io/badge/WASM-719%20KB-blueviolet?logo=webassembly)
![MCP](https://img.shields.io/badge/MCP-agent%20ready-purple)
[![许可证](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)
[![演示](https://img.shields.io/badge/demo-live-brightgreen)](https://kent-tokyo.github.io/chematic/playground/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/kent-tokyo/chematic/blob/main/notebooks/quickstart.ipynb)

面向 Python、Rust 和浏览器的化学信息学库。

**默认快速，设计安全的化学信息学库。**  
纯 Rust · 零 C/C++ · Python · WebAssembly · [在线演示](https://kent-tokyo.github.io/chematic/playground/)

| | chematic | RDKit (Python) | RDKit.js (WASM) |
|---|---|---|---|
| **快速上手** | `pip install chematic` | 需要 conda / cmake | 无 Python 绑定 |
| **浏览器包体积** | **719 KB** | 不支持 | ~30 MB（大约 42 倍） |
| **批量指纹速度** | **~78 µs/mol**（快 2–3 倍） | ~160–235 µs/mol | — |
| **内存安全性** | 编译器保证（Rust） | C++ | C++ |
| **源码构建** | 仅需 `cargo build` | cmake + clang + Boost | Emscripten SDK |

所有数据均可复现 — 参阅[基准测试详情](https://kent-tokyo.github.io/chematic/benchmark/)。  
WASM 包体积对比：chematic **719 KB** · RDKit.js ~30 MB · Indigo WASM ~40 MB

**功能成熟度一览：**

| 功能 | 状态 |
|---|---|
| SMILES / SMARTS / 指纹 / 描述符 | 稳定 |
| 3D 构象生成（DG + MMFF94） | 实验性 |
| pKa / ADMET | 基于规则的筛选（不适用于临床） |
| IUPAC 命名生成 | 部分实现（25+ 类别） |
| 纯 Rust InChI | 近似值（精确值需启用 `native-inchi` feature） |

---

## 何时使用 chematic

**适合使用 chematic 的场景：**

- 需要在浏览器中运行化学计算（WASM，719 KB，无需服务器）
- 需要纯 Rust 技术栈，不依赖 C++ 工具链
- 部署到 `pip install rdkit` 不可行的环境（Cloudflare Workers、Lambda、嵌入式设备）
- 构建 AI 代理并需要原生 MCP 工具集成
- 需要批量高吞吐量处理分子（ECFP4：比 RDKit 快 2–3 倍，Rayon 并行）
- 希望 `pip install chematic` 在任何环境都能直接使用，无需编译器

**适合使用 RDKit 的场景：**

- 需要最大的生态兼容性和 20 年以上的生产验证
- 需要带 ML 辅助扭转角修正的出版级 3D 结构（RDKit 的 ETKDGv3）
- 需要在不启用 `native-inchi` feature 的情况下获得逐位精确的标准 InChI
- 依赖基于 RDKit Python API 编写的社区插件

---

## 快速开始

### 安装

```bash
pip install chematic  # 无需 C/C++ 编译器
```

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # 阿司匹林

# 在 Jupyter 中直接写 mol 即可自动渲染 2D 结构
mol

# 访问 190+ 描述符值（属性形式）
print(mol.mw, mol.logp, mol.tpsa)           # 180.16  1.31  63.6
print(mol.lipinski_passes, mol.pains_passes) # True   True

# 子结构搜索
mol.has_substructure("[OH]")   # True
mol.find_matches("[CX3](=O)O") # → [[1, 2, 3], [7, 8, 9]]

# 自然语言摘要（适用于 LLM / MCP 代理）
print(mol.describe())
# → "Molecular weight 180.2 Da, formula C9H8O4. LogP 1.31 (mildly lipophilic)..."

# 两分子结构差异
ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")
d = mol.diff(ibuprofen)  # {"summary": "+C7, -O2. ΔLogP +2.75 ...", "delta_mw": 66.1, ...}

# 批量处理（并行，支持 numpy）
fps = chematic.bulk.ecfp4(["CCO", "c1ccccc1"])  # (2, 2048) uint8

# 一行生成 DataFrame
df = chematic.descriptors_df(["CCO", "c1ccccc1", "CC(=O)O"])
df[["mw", "logp", "tpsa", "qed"]]
```

Rust 和 JavaScript/TypeScript 的详细示例请参阅[文档](https://kent-tokyo.github.io/chematic/)。

---

## AI / LLM 开发者

chematic 是首个内置 **MCP（模型上下文协议）服务器**的化学信息学库，可直接集成 AI 代理。

```json
// Claude Desktop (~/.config/claude/claude_desktop_config.json)
{
  "mcpServers": {
    "chematic": { "command": "chematic-mcp" }
  }
}
```

15 个化学工具可由任意 MCP 兼容代理调用：

| 工具 | 功能 |
|---|---|
| `name_to_smiles` | 通过 PubChem 将化合物名（"阿司匹林"…）解析为 SMILES |
| `calc_properties` | MW、LogP、TPSA、HBA/HBD、QED、SA Score、pKa、ADMET |
| `smarts_match` | 子结构搜索 |
| `pains_check` / `brenk_check` | 筛查检测干扰或活性基团 |
| `generate_3d` | 3D 坐标生成（ETKDG + MMFF94） |
| `find_mcs` | 最大公共子结构 |
| 其余 9 个 | `ecfp4`、`tanimoto`、`canonical_smiles`、`admet_profile`、`boiled_egg`、`sa_score`、`lipinski_check`… |

---

## 为什么选择纯 Rust？

### 快速

Rust 的零成本抽象和所有权模型从源头消除开销。
chematic 的 ECFP4 指纹批处理在多样化分子语料库上达到 **~78 µs/mol** — 在同等硬件上比
RDKit Python API 快 2–3×（通过 Rayon 在所有 CPU 核心上并行）。无 GIL，无解释器开销，无 `_sys` crate 中隐藏的 FFI 调用开销。

### 安全

整个默认依赖树在 15,000 余行 Rust 代码中仅含 **~6 个 `unsafe` 块**。
无 C++ 堆损坏，无因畸形 SMILES 输入导致的段错误，无 `-sys` crate 引起的平台相关构建失败。
编译器在每个调用点强制执行内存安全。

> `native-inchi` feature 是唯一的可选例外 — 它 vendors IUPAC InChI C 库 (v1.07.5)
> 以获得逐位精确的标准 InChI。其余所有 crate 保持零 FFI。

### 随处可用

纯 Rust 无需 Emscripten、`cmake`、`clang` 即可原生编译至 `wasm32-unknown-unknown`。
npm 包 `@kent-tokyo/chematic` 为 **719 KB gzip** — 比 RDKit.js 小约 42 倍。
一套代码库在 Linux、macOS、Windows 及任意浏览器中运行。

---

## 与其他化学信息学库的比较

| 功能                                         | **chematic**                                 | RDKit (rdkit-sys)  | OpenBabel FFI | RDKit.js (WASM)   |
|----------------------------------------------|----------------------------------------------|--------------------|---------------|-------------------|
| **C/C++ 依赖**                               | **零（默认）**†                              | 大量 C++           | 大量 C++      | C++（Emscripten） |
| **WASM 二进制体积**                          | **〜550 KB**                                 | N/A（不支持 WASM） | N/A           | 〜30 MB           |
| **构建要求**                                 | 仅需 `cargo build`                           | cmake + clang      | cmake + clang | Emscripten SDK    |
| **Python 绑定**                              | **有** (`pip install chematic`, PyO3)        | 有（rdkit-sys）    | 有            | 无                |
| Unsafe Rust                                  | **无**                                       | 大量               | 大量          | N/A               |
| Kekulization                                 | **4-pass（含 Edmonds' blossom）**            | 有                 | 有            | 有                |
| SDF/MOL V2000+V3000                          | 有                                           | 有                 | 有            | 有                |
| Tripos MOL2 格式                             | **有**（读写 + Python）                      | 有                 | 有            | 无                |
| 分子描述符                                   | **190+ 描述符值**（71 个函数；MQN×42、BCUT2D、autocorr2d 返回多值数组） | 〜30               | 〜20          | 〜30              |
| **MAP4 指纹**                                | **有**（Minervini 2020）                     | 无（外部包）       | 无            | 无                |
| MMFF94 全 7 能量项                           | **有**                                       | 有                 | 有            | 无                |
| 3D 坐标生成                                  | 有（DG + MMFF94/DREIDING + L-BFGS）          | 有（ETKDG）        | 有            | 有                |
| 多样性筛选（MaxMin/Butina）                  | **有**                                       | 有                 | 无            | 无                |
| InChI / InChIKey                             | **有** — 纯 Rust（默认）+ **IUPAC 标准**（`native-inchi`）| 需 C 库 | 需 C 库 | 需 C 库 |
| **pKa 预测**                                 | **有（15 条 SMARTS 规则）**                  | 无                 | 无            | 无                |
| **ADMET 简况 + BOILED-Egg**                  | **有**                                       | 部分               | 无            | 部分              |
| **MCP 服务器（AI Agent API）**               | **有 — 20 个工具（含 Name→SMILES）**        | 无                 | 无            | 无                |
| IUPAC 命名生成                               | **有（25+ 化合物类）**                       | 无                 | 无            | 部分              |
| 维护状态（2026）                             | 活跃                                         | 活跃               | 最小维护      | 活跃              |

† 仅限默认构建。`native-inchi` feature 需要 C 编译器，为可选例外。其余所有 crate 保持零 FFI。

---

## JavaScript / TypeScript（WebAssembly）

**719 KB gzip — 比 RDKit.js 小约 42 倍。** 无需 Emscripten 或 cmake，可直接在浏览器和 Node.js 中使用。

```sh
npm install @kent-tokyo/chematic
```

```js
import init, { parse_smiles, get_descriptors_json, tanimoto_ecfp4,
               generate_3d_minimized_pdb, enumerate_stereo_isomers_json,
               maxmin_picks_ecfp4_json } from '@kent-tokyo/chematic';

await init();

const mol = parse_smiles('CC(=O)Oc1ccccc1C(=O)O'); // 阿司匹林
console.log(mol.molecular_weight(), mol.qed(), mol.lipinski_passes());

const desc = JSON.parse(get_descriptors_json(mol));  // 一次获取所有描述符
const caffeine = parse_smiles('Cn1cnc2c1c(=O)n(c(=O)n2C)C');
console.log(tanimoto_ecfp4(mol, caffeine));           // ECFP4 相似度: 0.26

const pdb = generate_3d_minimized_pdb(mol);           // 3D 坐标生成
const isomers = JSON.parse(enumerate_stereo_isomers_json(parse_smiles('C(F)(Cl)Br')));
const picks = JSON.parse(maxmin_picks_ecfp4_json('["CC","c1ccccc1","CCO","CCCC"]', 2));
```

130+ 个导出函数涵盖描述符、指纹、3D 几何、反应、多样性筛选和 SDF 处理。
完整 API 请参阅 [WASM API 参考文档](https://kent-tokyo.github.io/chematic/)。
---

## Crate 列表

| Crate                 | 说明                                                                                                   | 测试数 |
|-----------------------|--------------------------------------------------------------------------------------------------------|--------|
| `chematic-core`       | Atom、Bond、Molecule、Element、Kekulization（无依赖）；可变 API、`fragments`、`validate_valence`、`formula_with_isotopes`；`StereoGroup`/`StereoGroupKind` | 71     |
| `chematic-smiles`     | OpenSMILES 解析器、写入器、规范 SMILES、**CXSMILES 元数据支持**                                      | 109     |
| `chematic-perception` | SSSR、Hückel 芳香性 + 反芳香性（4n+2 规则）、`apply_aromaticity`/`aromatize`/`kekulize_inplace`、`assign_stereo_from_2d`、`assign_ez_from_2d`、`cip_ez_descriptor` | 101     |
| `chematic-mol`        | MOL/SDF V2000+V3000（读写含 2D 坐标）、CML（读写）、CDXML（读）；`SdfRecord`（含坐标+属性）、MDL RXN V2000 读写；V3000 立体基团 COLLECTION 读写 | 130     |
| `chematic-depict`     | 2D SVG 绘制（CPK 配色、高亮、网格）、`detect_crossings`/`render_svg_with_metadata`、反应 SVG；Y 坐标系文档已更新 | 64     |
| `chematic-chem`       | 190+ 描述符值（71 个函数）、互变异构体、骨架、BRICS、QED、标准化；**pKa 预测**（15 条 SMARTS 规则）；**ADMET 概况**（BBB/Caco-2/hERG/CYP3A4）；**HBA 与 RDKit 一致率 99.98%**（4,999 分子 ChEMBL 基准）；**TPSA ±0.1 Å² 98.1% / LogP ±0.01 96.5% / HBD 100%** | 662    |
| `chematic-fp`         | ECFP2/4/6、FCFP4/6、MACCS、TopoPF、AtomPair、Torsion、Layered、Pattern、Pharmacophore、Reaction、**MAP4** — Tanimoto/Dice | 185     |
| `chematic-ff`         | **MMFF94 全 7 能量项**（Halgren 1996）：OOP（117 条）+ STRE-BEN（282 条）；L-BFGS；DREIDING | 98     |
| `chematic-smarts`     | SMARTS、VF2、MCS；**SmartsCache**（LRU 5–20×）；**named_pattern()** 库（20 种模式）；**SMARTS 原子映射 `:N`**（`[O;D1;H0:3]` — 作为元数据存储，不用于匹配） | 142    |
| `chematic-3d`         | 3D 坐标生成、ETKDG KB（40 种模式，自适应噪声）、力场最小化、形状描述符、ConformerEnsemble、PDB/XYZ | 265    |
| `chematic-rxn`        | 反应 SMILES/SMIRKS、`run_reactants`/`run_reactants_strict`；**`retro_disconnect()`** — 60 个 retro-SMIRKS 模板（AmideBond/Ester/Ether/CNBond/CCBond/CSBond）+ SA 分数排序 | 137     |
| `chematic-inchi`      | InChI/InChIKey：纯 Rust 近似（WASM 兼容）**+ `native-inchi` feature 提供 IUPAC 标准**（vendored C 库 1.07.5，逐位一致）；**parse_inchi** 读取 | 96 (+16*)   |
| `chematic-wasm`       | **130+ WASM 导出** — npm：`@kent-tokyo/chematic` v0.4.30（719 KB gzip）；**pKa/ADMET/BBB/Caco-2/hERG/CYP3A4** WASM API | 211    |
| `chematic-iupac`      | 本地 IUPAC 命名（纯 Rust·离线）— **25+ 化合物类**：烷烃、环烷烃、醇、胺、卤代烃、酮、酸、酯、酰胺、**哌啶、吗啉、哌嗪、萘、硫醚** | 47     |
| `chematic-mcp`        | **MCP（模型上下文协议）服务器** — AI 代理集成；**20 个工具**：parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, name_to_smiles, retrosynthesis, smiles_to_moljson, moljson_to_smiles, representation_router, molecule_context_pack | 31     |
| `chematic`            | 带功能标志的伞形 crate                                                                                   | 1      |

```
cargo test --workspace --lib --quiet                                               # 2,366 个库测试，全部通过
cargo test -p chematic-inchi --features native-inchi --test standard_inchi         # +16 IUPAC 标准 InChI 集成测试
```

---

## 近期开发（v0.4.x）

**v0.4.19**（2026-06-23）：**PDF/EPS 输出、ChemicalJSON、新描述符、WASM −38.5%**
- `chematic-depict`：`depict_pdf()` / `depict_eps()` — PDF 和 EPS 输出；纯 Rust，无需外部工具
- `chematic-mol`：**ChemicalJSON** — `parse_cjson()` / `write_cjson()` 支持 Avogadro2 / MolSSI 互操作
- `chematic-chem`：4 个新描述符 — `schultz_mti()`, `gutman_mti()`, `vabc()`（Bondi vdW 体积）, `gravitational_index()`
- `chematic-3d`：**Spectrophores** 3D 指纹（药效团壳编码）
- `chematic-py`：`mol.to_pdf()`, `mol.to_eps()`, `mol.to_cjson()`, `from_cjson()`；`bulk.substructure_match(smarts, mols)` 并行 VF2；`estate_all()` / `ring_bundle` 加入 bulk
- **WASM 包：819 → 504 KB gzip（−38.5%）** — `tiny_skia` 改为可选、内联 SHA-256、`opt-level="z" lto=true codegen-units=1`

**v0.4.18**（2026-06-23）：**Python API 扩展 + 基准文档**
- `chematic-py`：**Jupyter 自动显示** — 在单元格中写 `mol` 即可渲染 2D 结构（`_repr_svg_()` 钩子）；`mol.has_substructure(smarts)`, `mol.find_matches(smarts)`；`from_smiles_list()`, `descriptors_df()`
- `chematic-chem`：`chi_all()` — 单次遍历计算全部 10 个 Hall-Kier 连接性指数；`cns_mpo_from_parts()`；`pains_passes_and_matches()` / `brenk_passes_and_matches()` — 单次扫描同时返回标志和名称
- 文档：新增基准页面（ECFP4 比 RDKit 快 5–14×，4,999 分子 ChEMBL 语料库描述符 100% 准确）

**v0.4.16–v0.4.17**（2026-06-22–23）：**SSSR 共享性能冲刺**
- `chematic-smarts`：`find_matches_with_rings()` — 批量模式下共享一次预计算的 `RingSet`
- `chematic-chem`：Crippen 117 SSSR → 1 次/调用；PAINS ~480 → 1；QED 113 → 1；pKa 42 → 1；新增 `logp_and_mr()`, `logd_from_logp()`, `pka_both()`
- `chematic-fp`：MHFP 增量 BFS — 每分子 3N → N 次 BFS（radius=2 时）

**v0.4.15**（2026-06-21）：**TPSA 校准 + 反应 E/Z 立体**
- `chematic-chem`：TPSA ±0.1 Å² 校准冲刺 — **HBA 100%、HBD 100%、芳香环计数 100%**（4,999 分子 ChEMBL 语料库）；TPSA 86.7% → 93.3%（4,999 分子），175 分子药物样集 100%
- `chematic-rxn`：`run_reactants` 新增 E/Z 几何过滤 — 通过 `smirks_ez_stereo_ok()` / `ez_stereo_outward()` 进行 SMIRKS `/`/`\` 几何匹配

**v0.4.14**（2026-06-21）：**拓扑描述符 + 立体正确性**
- `chematic-chem`：8 个拓扑描述符 — `petitjean_index()`, `graph_eccentricities()`, `graph_diameter()`, `graph_radius()`, `eccentric_connectivity_index()`, `hosoya_index()`, `moran_autocorr()`, `geary_autocorr()`
- `chematic-3d`：GETAWAY HATS-matrix（19 维）；`whim_getaway_combined()` 扩展至 29 维
- `chematic-smiles`：累积双键立体化学 `C=C=C` `@`/`@@` — 往返稳定
- `chematic-smarts`：`[kN]` 环大小原语；当查询原子数超过目标时 VF2 提前退出
- `chematic-rxn`：奇偶校验感知的 SMIRKS 手性匹配；产物括号原子清理

**v0.4.13**（2026-06-21）：**模板逆合成 + 描述符修复**
- `chematic-rxn`：`retro_disconnect()` — 60 个 retro-SMIRKS 模板（AmideBond / Ester / Ether / CNBond / CCBond / CSBond）附 SA 分数排序；Python `mol.retro_disconnect(reaction_class=...)`
- `chematic-3d`：ETKDG 扭转知识库 28 → 40 种模式；自适应噪声
- `chematic-chem`：`hbd_count()` 新增 S-H（硫醇）；TPSA 硝基-N / 芳香氧桥 / Kekulé-N 修复

**v0.4.9–v0.4.12**（2026-06-19–21）：**AutoDock、UFF、SMARTS 原子映射、环感知**
- `chematic-mol`：AutoDock PDBQT 读写；`write_sdf_with_charges`（部分电荷）
- `chematic-ff`：金属/有机金属 UFF 力场（Zn、Fe、Cu…）
- `chematic-smarts`：SMARTS 原子映射 `:N`（`[O;D1;H0:3]` 格式，作为元数据存储）
- `chematic-perception`：稠合多环芳香环计数的迭代 `augmented_ring_set`（修复 bench5k 222 个失败案例）
- MCP：第 15 个工具 `name_to_smiles`（PubChem REST 代理）

**v0.4.5–v0.4.7**（2026-06-19）：**Kekulization Blossom + BOILED-Egg + InChI E/Z**
- Edmonds' blossom 算法（128 → 2 失败）；InChI `/b` E/Z 层；BOILED-Egg + Python/WASM 绑定

**v0.4.0–v0.4.4**（2026-06-17–18）：**PyO3 Python 绑定 + native-inchi**
- `chematic-py`：PyO3/maturin 绑定 — `from_smiles()`, `Mol.aromatic_ring_count`, `Mol.descriptors()`
- `native-inchi` feature：IUPAC 标准 InChI（vendored C 库 v1.07.5）
- HBA 重写：与 RDKit 一致率 99.98%（4,999 分子 ChEMBL 基准）

---

## 已知限制

- **Kekulization**：4,999 分子中仅 2 个失败 — 硼芳香环（`b1ccccn1`）和 `[H][H]`。明确返回 `KekuleError`，不产生无声错误输出。
- **芳香性模型**：Hückel 4n+2 规则独立应用于每个 SSSR 环（RDKit 使用稠合环电子离域模型）。N-杂环中存在差异。4,999 分子 ChEMBL 语料库当前状态：HBA/HBD/芳香环计数 **100%**，TPSA **98.1%**（±0.1 Å²）。

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

---

如果 chematic 对您有帮助，欢迎给项目一个 [GitHub star](https://github.com/kent-tokyo/chematic)，这将帮助更多人发现它。
