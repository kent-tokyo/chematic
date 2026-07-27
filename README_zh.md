# chematic

[English](README.md) | [日本語](README_ja.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![docs.rs](https://docs.rs/chematic/badge.svg)](https://docs.rs/chematic)
[![许可证](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)

[Open in Colab](https://colab.research.google.com/github/kent-tokyo/chematic/blob/main/notebooks/quickstart.ipynb)

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

20 个化学工具可由任意 MCP 兼容代理调用（完整列表见 [`chematic-mcp` README](crates/chematic-mcp/README.md)）：

| 工具 | 功能 |
|---|---|
| `name_to_smiles` | 通过 PubChem 将化合物名（"阿司匹林"…）解析为 SMILES（唯一进行网络通信的工具） |
| `calc_properties` | MW、exact mass、Crippen LogP、TPSA、HBD、HBA、rotatable bonds、QED |
| `smarts_match` | 子结构搜索 |
| `pains_check` / `brenk_check` | 筛查检测干扰或活性基团 |
| `generate_3d` | 基于规则的坐标放置 + DREIDING 力场最小化生成 3D 坐标 |
| `find_mcs` | 最大公共子结构 |
| 其余 13 个 | `ecfp4`、`tanimoto`、`canonical_smiles`、`admet_profile`、`boiled_egg`、`sa_score`、`lipinski_check`、`retrosynthesis`、`smiles_to_moljson`、`moljson_to_smiles`、`representation_router`、`molecule_context_pack`、`parse_smiles` |

**Transport**：仅 stdio（通过标准输入输出的 JSON-RPC 2.0）。以本地进程方式运行，目前没有已托管的 Remote MCP 端点、身份验证或公开服务 SLA；remote 化的重构正在考虑中，尚未实现。

**Protocol**：在同一个 stdio 连接上同时支持旧版（`2024-11-05` 风格的 `initialize` 握手）和 MCP `2026-07-28` 无状态方言（`server/discover`、逐请求 `_meta`、可缓存的 `tools/list`、`structuredContent`）。详见 [`chematic-mcp` README](crates/chematic-mcp/README.md#protocol-eras) 与 [`docs/mcp/2026-07-28-implementation-rfc.md`](docs/mcp/2026-07-28-implementation-rfc.md)。Remote HTTP、OAuth、Tasks 扩展与 MCP Apps 仍不支持。

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
| **MCP 服务器（AI Agent API）**               | **有 — 20 个工具（含 Name→SMILES，仅 stdio）** | 无                 | 无            | 无                |
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
| `chematic-mol`        | MOL/SDF V2000+V3000（读写含 2D 坐标）、CML（读写）、CDXML（读）；`SdfRecord`（含坐标+属性）、MDL RXN V2000 读写；V3000 立体基团 COLLECTION 读写；**读取时自动识别 2D 楔形/虚线四面体 parity 与 E/Z 双键方向**（`read_mol_with_diagnostics`/`read_mol_v3000_with_diagnostics`，类型化 opt-in 诊断） | 130+     |
| `chematic-depict`     | 2D SVG 绘制（CPK 配色、高亮、网格）、`detect_crossings`/`render_svg_with_metadata`、反应 SVG；Y 坐标系文档已更新 | 64     |
| `chematic-chem`       | 190+ 描述符值（71 个函数）、互变异构体、骨架、BRICS、QED、标准化；**pKa 预测**（15 条 SMARTS 规则）；**ADMET 概况**（BBB/Caco-2/hERG/CYP3A4）；**HBA 与 RDKit 一致率 99.98%**（4,999 分子 ChEMBL 基准）；**TPSA ±0.1 Å² 98.1% / LogP ±0.01 96.5% / HBD 100%** | 662    |
| `chematic-fp`         | ECFP2/4/6、FCFP4/6、MACCS、TopoPF、AtomPair、Torsion、Layered、Pattern、Pharmacophore、Reaction、**MAP4** — Tanimoto/Dice | 185     |
| `chematic-ff`         | **MMFF94 全 7 能量项**（Halgren 1996）：OOP（117 条）+ STRE-BEN（282 条）；L-BFGS；DREIDING | 98     |
| `chematic-smarts`     | SMARTS、VF2、MCS；**SmartsCache**（LRU 5–20×）；**named_pattern()** 库（20 种模式）；**SMARTS 原子映射 `:N`**（`[O;D1;H0:3]` — 作为元数据存储，不用于匹配） | 142    |
| `chematic-3d`         | 3D 坐标生成、ETKDG KB（40 种模式，自适应噪声）、力场最小化、形状描述符、ConformerEnsemble、PDB/XYZ | 265    |
| `chematic-rxn`        | 反应 SMILES/SMIRKS、`run_reactants`/`run_reactants_strict`；**`retro_disconnect()`** — 60 个 retro-SMIRKS 模板（AmideBond/Ester/Ether/CNBond/CCBond/CSBond）+ SA 分数排序 | 137     |
| `chematic-inchi`      | InChI/InChIKey：纯 Rust 近似（WASM 兼容）**+ `native-inchi` feature 提供 IUPAC 标准**（vendored C 库 1.07.5，逐位一致）；**parse_inchi** 读取；**带验证的 canonical SMILES 去重**（`dedup` 模块，遇到 legacy CIP 无法解析的指定四面体立体中心时安全失败） | 96 (+16*)   |
| `chematic-cip`        | opt-in 高精度 CIP 引擎（`assign_cip_accurate_experimental`，层次化 digraph，Rules 1a/1b/2/4b/5，RDKit 兼容 MANCUDE 分数原子序数）— 默认的 `assign_cip()`/`CipMode::LegacyFast` 未变更 | —    |
| `chematic-wasm`       | **130+ WASM 导出** — npm：`@kent-tokyo/chematic`（已发布 `0.7.0`，与 crates.io/PyPI 同步）；**pKa/ADMET/BBB/Caco-2/hERG/CYP3A4** WASM API | 211    |
| `chematic-iupac`      | 本地 IUPAC 命名（纯 Rust·离线）— **25+ 化合物类**：烷烃、环烷烃、醇、胺、卤代烃、酮、酸、酯、酰胺、**哌啶、吗啉、哌嗪、萘、硫醚** | 47     |
| `chematic-mcp`        | **MCP（模型上下文协议）服务器** — AI 代理集成；**20 个工具**：parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, name_to_smiles, retrosynthesis, smiles_to_moljson, moljson_to_smiles, representation_router, molecule_context_pack | 31     |
| `chematic`            | 带功能标志的伞形 crate                                                                                   | 1      |

```
cargo test --workspace --lib --quiet                                               # 2,746 个库测试，全部通过（截至 2026-07-26）
cargo test -p chematic-inchi --features native-inchi --test standard_inchi         # +16 IUPAC 标准 InChI 集成测试
```

---

## 近期开发

**v0.7.0**（2026-07-26）：**MOL/SDF 自动识别 2D 楔形/虚线 + E/Z 立体化学、带验证的 canonical SMILES 去重、CIP Rule-5 磷修复、native InChI 显式氢/同位素修复**
- `chematic-mol`/`chematic-perception`：MOL V2000/V3000/SDF 读取器现在读取时自动识别四面体楔形/虚线 parity（PR #154）与 E/Z 双键方向（PR #162），与 CIP 无关。类型化 opt-in 诊断（`StereoDiagnostic`/`EzDirectionDiagnostic`）对畸形/歧义输入从不猜测。大规模语料验证（4,999 分子，对比 RDKit 2026.03.3）：E/Z — 622 个 RDKit 可解析双键，语义反转 0 例，误报 0 例。构建过程中还修复了 V2000 MDL 代码 4 缺陷、两处 V3000 `CFG` 缺陷及一处 V2000 写入器缺陷。PR #162 自身的语料核对发现了**尚未修复的新缺口**：部分分子的 `canonical_smiles()` 会丢失已经在 `write()` 中正确编码的 E/Z 标记（因 aromaticity 无关的 carrier 分组导致）— 尚未建立 issue
- `chematic-inchi`：新增带验证的 canonical SMILES 去重（`dedup` 模块）— 快速 canonical SMILES 候选分桶，并与 native InChI 验证后的同一性进行核对。当指定立体中心的 legacy CIP 排序无法解析时安全失败（`VerificationUnavailable`），而非冒错误合并的风险（修复了在 5,000 分子语料验证中发现的真实 false `VerifiedDuplicate`）。后续追踪为 [#161](https://github.com/kent-tokyo/chematic/issues/161)（accurate CIP preflight 有望恢复大部分保守案例）
- `chematic-cip`：CIP Rule-5 拟不对称 r/s 磷修复
- `chematic-inchi`：native InChI（`native-inchi` feature）显式氢/同位素转换修复
- 详见 `CHANGELOG.md` 的 `[Unreleased]` 部分

**v0.6.0**（2026-07-25）：**RDKit bit-exact ECFP4 跨语言稳定 API、canonical SMILES 的 E/Z 标记一致性、opt-in 芳香性标志权威降级**
- `chematic-fp`/Python/WASM：将 RDKit bit-exact 的 Morgan/ECFP4 路径提升为文档化的跨语言 opt-in API（Python `Mol.rdkit_ecfp4()`，WASM `rdkit_ecfp4_bitvec()`），并推广为独立验证过的 `(radius, fpSize)` 矩阵（4×5=20 个组合，均为不可构造非法值的封闭枚举）。Rust/Python/WASM 三端均实际构建并运行验证 bit-exact 一致，而非仅"设计上应当一致"。`ecfp4()` 的默认行为不变
- `chematic-smiles`：修复了 canonical 输出中 E/Z 方向标记依输入原子顺序落在不同取代基键上的问题——标记放置的置换不变性从 **93.0% 提升到 98.1%**（282 个已知发散分子中 264 个现已收敛；剩余 18 例为避免破坏共享候选键的边界情况而故意保留未解决，追踪于 [#149](https://github.com/kent-tokyo/chematic/issues/149)）。同时排查并**排除**了此前怀疑的桥环双环 canonicalization 缺陷——经独立 RDKit InChI 验证，涉及的两个 SMILES 实际上是不同的分子，并非 chematic 的缺陷
- `chematic-perception`：新增 opt-in 的 `apply_aromaticity_authoritative_experimental`——使芳香性标志的升级/降级双向完全服从 Hückel 模型的实际计算结果（默认的 `apply_aromaticity`/`apply_aromaticity_ex` 保持不变，已验证 byte-identical）。在此过程中修复了影响并环二氮杂环（quinazoline/quinoxaline/purine 型环）的环并键误分类问题，附带解决了 32 个已有的假阳性回归 pin
- 完整细节和已知限制见 `CHANGELOG.md`

**v0.5.0**（2026-07-23）：**CIP 无关的 2D 楔形键 local parity、charge 感知 kekulization（此前 6 类分子解析失败）、PAINS/Brenk 的 budget 化非静默匹配**
- `chematic-perception`：新增 `local_parity_from_wedges`/`apply_local_parity_from_wedges`——完全不依赖 CIP 排序，直接从楔形/虚线键与 2D 坐标计算 `Atom.chirality`/`stereo_neighbor_order`，因此 CIP 排序打平也不会抹去已知的 local parity；符号约定针对 RDKit 原始 chiral tag 实测确定，而非类推得出。尚未接入任何 reader 的默认解析路径
- `chematic-core`：`kekulize()` 的原子匹配规则此前对电荷不敏感且缺失 Tellurium 处理——tropylium、imidazolium、pyridinium、pyrylium、tellurophene、phosphole 现均可成功 kekulize，且与 RDKit 逐键一致，零回归；`Element::normal_valences()` 新增了经源码验证的 Tellurium 条目，修复了由此导致的 ECFP4 芳香性不一致
- `chematic-perception`：charge 感知的 Hückel π 电子计数——tropylium、imidazolium、pyridinium、pyrylium 现已与 RDKit 的芳香原子/键标志完全一致（tellurophene/phosphole 及更广泛的权威降级修复仍待解决，单独追踪）
- `chematic-smiles`：修复两个 writer 缺陷——bracket 强制原子（isotope/charge/atom-map）此前会静默丢弃隐式氢（`[NH4+]` → `[N+]`）；无相邻双键的独立楔形键此前被错误写成无意义的 SMILES `/`、`\` 标记
- `chematic-smarts`/`chematic-chem`：对称目标上的 PAINS/Brenk 子结构匹配此前可能挂起数分钟；VF2 现引入显式访问 budget，返回真正的三态结果（`Found`/`NotFound`/`BudgetExhausted`），探索耗尽绝不会被静默折叠为假阴性——后续的对称感知候选排序（可让现有 budget 正确解决部分保守标记案例）追踪于 [#139](https://github.com/kent-tokyo/chematic/issues/139)
- 完整细节、语料库层面的前后对比数字及已知限制见 `CHANGELOG.md`

**v0.4.30**（2026-07-17）：**`chematic-cip` opt-in 接入全部接口、SMARTS `[rN]` 修复、新 RDKit-parity 芳香性引擎、5 处立体元数据缺陷修复**
- `chematic-smarts`：修复 `[rN]`（环大小 SMARTS，如 `[r5]`/`[r6]`）被错误地等同于 `[kN]` 的"任意环"语义——RDKit 真正的 `[rN]` 含义是"该原子所在的*最小*环恰好为 N 大小"，是完全不同的谓词（实测确认：在共享 5 元环与 6 元环的并环原子上，RDKit 的 `[k6]` 匹配但 `[r6]` 不匹配）。未改动环模型本身——`[rN]` 现拥有基于 chematic 现有 SSSR 计算的独立 `MinRingSize` 原语。SMARTS 匹配集合与 RDKit 的一致率在 5,000 分子语料库上从 **96.9% 提升到 99.93%**，零回归；详见 `docs/rdkit_compat.md` 的 "SMARTS-R0"/"SMARTS-R1" 条目
- Milestone 5A：从所有公开接口 opt-in 访问精确引擎——`chematic_chem::assign_cip_with_mode(mol, CipMode::Accurate)`（Rust）、`Mol.cip_stereo(mode="accurate")` + `Mol.cip_stereo_unresolved()`（Python）、`cip_assignments_accurate_json`/`cip_unresolved_json`（WASM）。所有默认接口（`assign_cip()`、`cip_stereo()`、无后缀的 WASM 函数）保持不变——纯增量，非默认切换；合并语义（精确引擎的四面体 R/S + legacy 的 E/Z，因精确引擎不计算键立体）及"绝不猜测"契约（打平/超预算情况显式暴露，绝不静默回退）详见 `docs/cip_accurate_rfc.md` 的 Milestone 5A 条目
- Milestone 4 关口达成：全语料库、按表示稳定性分层的 oracle-stable 一致率达 99.64%（原始一致率 99.38%，4160/4186）——最后剩余的 11 行磷残差经查明为 oracle 不稳定（RDKit 自身标签在化学中性的 Kekulé 重新拼写下会变化），并非 chematic 缺陷；15 行 Rule 5 笼状家族仍延后处理，不受此关口影响
- 精确引擎（携带溯源信息、逐层递归的 digraph 比较器——Rule 1a/1b/2，外加芳香环立体中心的 RDKit 兼容 MANCUDE 分数原子序数）已通过上述 opt-in API 提供，但尚未成为 `assign_cip()` 背后的默认实现
- 发现并修复一个真实的约 10-14 倍性能回归（SSSR 被误用于布尔型环键判断，替换为 O(V+E) 桥边 DFS）；CI Criterion 关口的引导脚本修复；一项 Criterion 关口可靠性发现（伪重复采样问题，[#70](https://github.com/kent-tokyo/chematic/issues/70)）——已落地流程级重新设计（独立进程运行观测、两阶段筛选、同一二进制空对照）。后续修复（[#117](https://github.com/kent-tokyo/chematic/pull/117)）解决了通过真实 CI 运行发现的两个更具体的缺陷：Stage 1 的 3 组符号检验在结构上永远无法失败（替换为纯幅度路由筛选，阈值来自对 28 次历史无操作运行的离线评估），Stage 2 的符号检验对较小的构建/代码生成差异同样不敏感（改为符号检验 + 实际效应阈值的组合关口）。该关口仍为非必需项；同一二进制空对照对跨二进制代码生成差异的盲区是已知、已报告但尚未修复的缺口
- Milestone 4A：`CipCode::LowerR`/`LowerS`——Rule 5（伪不对称性），仅限 2 行已验证独立案例；发现一个三臂对称笼状家族（15 行）在此两两配对架构下被证明无法触及，延后为 Milestone 4A-2（需要对称性/自同构检测）
- Milestone 4A-0：从零重新冻结残差为 34 行并对其 100% 机械分类（0 项未解释）——15 行 Rule 5/伪不对称性（4A-2 笼状家族）、8 行 Rule 4 候选（通过结构同一性检查正面确认，非排除法推断）、11 行磷（9 行比较器缺陷导致"错误"+2 行确实打平）
- `chematic-perception`：新增 opt-in 的 `assign_aromaticity_rdkit_parity_experimental`/`apply_aromaticity_rdkit_parity_experimental`——对 RDKit 真实芳香性算法的源码级验证移植，在 4,999/5,000 个可比分子上实现 **100.0000% 原子/键一致率**。未接入默认路径（`RdkitLike`/`Huckel` 保持不变）；默认提升被一个与此引擎无关、已存在的 canonical-SMILES-writer 敏感性问题所阻塞
- 修复了同一"元数据未复制"缺陷的 5 个实例（`MoleculeBuilder` 重建时未调用 `copy_stereo_groups_from`/`copy_stereo_from`/`copy_bond_directions_from`），各自静默丢弃 `stereo_neighbor_order` 或更严重：`apply_kekule`（P0）、`enumerate_stereoisomers`（可能静默翻转新分配立体中心的 CIP 编码）、`transfer_hydrogen_aromatic`/`clone_mol`（后者已删除，改用 `Molecule::clone`）、`transfer_hydrogen`、以及 `invert_stereocenter`（结果发现它对纯 `@`/`@@` SMILES 输入是功能性空操作，属于另一个更严重的缺陷）
- `chematic-smiles`：将芳香键方向暂存逻辑统一整合到 3 条 parser 键创建路径（链边/环闭合/分支连接）共用的辅助函数中——修复了一处 canonical round-trip 表示不稳定问题（4,994/5,000 → 5,000/5,000 稳定）；未修复 `assign_ez` 对该旁路通道本身存在的既有盲区（作为后续任务追踪）
- 实验性 CIP 引擎的全语料库准确率从 96.68% 提升到 99.38%（原始）/ 99.64%（oracle-stable），对比现代 RDKit `rdCIPLabeler`（零回归）——完整里程碑历史见 `docs/cip_accurate_rfc.md`
- 基准测试已刷新（`benchmarks/2026-07-17.md`，Apple M4）：此前的 ECFP4 吞吐量头条数字（3.6 µs/mol，比 RDKit 快 5–14 倍）在干净重新测量下未能复现——本 README 及 `docs/` 中已更新为今日实测数字（多样化语料库上约 78 µs/mol / 2–3 倍）；描述符准确率数字复现良好

**v0.4.29**（2026-07-10）：**Kabsch 旋转缺陷修复 + SDF V3000/CDXML 写入、Avalon 指纹、O3A**
- `chematic-3d`：修复 `align_coords` 的 Kabsch 旋转计算方向错误——对任何非纯平移的对齐都会给出严重偏大的 RMSD（在此修复前已随 v0.4.28 发布到 crates.io/PyPI/npm）；新增用于 O3A 原子对应关系的 `correspondence_search`
- `chematic-mol`：SDF V3000 写入接线；CDXML 写入
- `chematic-fp`：Avalon 指纹

**v0.4.28**（2026-07-09）：**SMARTS 性能优化、注册表重新同步**
- `chematic-smarts`：存在性检查短路——`bulk.substructure_search` 比 RDKit 快 2.2 倍
- v0.4.23–v0.4.27 期间未推送 git tag（crates.io 通过手动 `cargo publish` 保持最新，但 PyPI/npm/GitHub Releases 落后）——本版本重新同步三个注册表

**v0.4.27**（2026-07-04）：**描述符修复、RWMol/FCFP、veridict CI 关口**
- `chematic-chem`：修复 `kappa1-3`、`balaban_j`、`labute_asa`、`bcut2d`、`hall_kier_alpha` 描述符
- `chematic-fp`：`useFeatures=True` 的 FCFP
- `chematic-mol`：RWMol 原地编辑
- CI：基于 veridict 的性能/Criterion/准确率漂移回归关口；修复集成测试 CI 覆盖缺口

**v0.4.26**（2026-06-29）：**反应中的 E/Z 立体转移 + 验证 Sprint 6/7**
- `chematic-rxn`：反应产物现在会在 `run_reactants()` 中保留来自反应物的 `/`/`\` 双键几何信息（此前在转化过程中丢失）
- 验证：canonical SMILES 与 RDKit 的差异化验证（Sprint 6）；SMARTS/芳香性差异化测试及 I/O 兼容性（rdkit_compat Sprint 7）；将剩余的 RDKit canonical 差异根因定位到芳香性 round-trip，而非 Morgan 排序

**v0.4.25**（2026-06-29）：**`chematic.rdkit_compat` 层**
- `chematic-py`：RDKit API 兼容层（Sprint 1–5）——Morgan `bitInfo`、Fingerprint/Mol/Atom/Bond/RingInfo 兼容性、针对 RDKit 的差异化测试；流式 `SDMolSupplier`/`SDWriter`/`Mol.GetProp`
- `chematic-perception`：`AromaticityAlgorithm::RdkitLike`——匹配 RDKit 模型的 Se/Te 硫族芳香性

**v0.4.24**（2026-06-29）：**CIP Rule 5、桥头原子/可旋转键/TPSA/摩尔折射率达到 100%、HDF 指纹**
- `chematic-chem`：CIP Rule 5 立体打破平局（立体中心一致率 99.8% → 99.98%，对比 RDKit）；桥头原子检测 98.5% → 100%；可旋转键 99.1% → 100%；TPSA 100%；摩尔折射率 97.5% → 100%（3 环 XOR 增强）——均在 5,000 分子 ChEMBL 语料库上测得
- `chematic-py`：`bulk.descriptors_array()` 列式 numpy 输出；真正流式 SDF（`SdfFileReader`/`iter_sdf_batched`）；`screen()` 化合物筛选工作流
- LLM/RAG：表示路由器（`to_llm_text`、`best_representation`）、分子上下文包、**超维指纹（HDF）**——无需训练的稠密分子向量

**v0.4.23**（2026-06-26）：**LogP 96.5% → 99.7%**
- `chematic-chem`：修复 `crippen_anchor_sets` 使用 `uniquify: false`，使对称三键（内部炔）能得到两种 VF2 匹配方向，而不是其中一种回退到通用的 `[#6]` 值

**v0.4.22**（2026-06-26）：**CITATION.cff + `chematic.doctor()`**
- `chematic-py`：`doctor()` 自诊断；README 新增按功能可靠性矩阵

**v0.4.21**（2026-06-25）：**面向 LLM/Jupyter 的 HTML/Markdown 报告**
- `chematic-py`：`chematic.report()` 自包含 HTML 化合物网格、`chematic.compare()`、`mol.review()` Markdown 分析
- 文档：`benchmarks/`/`validation/` 可复现准确率历史

**v0.4.20**（2026-06-25）：**ETKDG 扭转知识库 44 → 80 条规则、`mol.describe()`/`diff()`**
- `chematic-3d`：为 6/5 元脂肪环新增椅式/信封式环构象；SMARTS 匹配的扭转规则作为高精度预检层
- `chematic-py`：面向 LLM/MCP agent 的 `mol.describe()`/`mol.diff(other)`；`bulk.generate_3d`/`tanimoto_matrix`/`standardize`

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
