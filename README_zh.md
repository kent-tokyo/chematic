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
默认纯 Rust · 可选原生 InChI C FFI · Python · WebAssembly · [官方网站](https://chematic.io/) · [在线演示](https://kent-tokyo.github.io/chematic/playground/)

### v1.0.8 范围

v1.0.8 保持 v1.0.0 的 bounded 兼容性边界，并加入 typed reaction document、
document-level CDXML 编辑、显式且 bounded 的 Markush/polymer 展开、晶体组成汇总、
更安全的 UFF rescue，以及 canonical/SDF 热路径优化。完整任意结构 CDXML 编辑、
复杂拓扑 expansion、完整 RDKit `RWMol` 与完整 ETKDG/MMFF94 parity 仍不支持。
Spectrophores 在 patent/FTO 状态得到独立确认前已从 Rust/Python API 中移除。
基准数字仅适用于记录的语料、操作、硬件与配置；完整方法见[基准文档](docs/benchmark.md)。
算法与第三方来源边界记录在[实现来源文档](docs/implementation-provenance.md)中。

| | chematic | RDKit (Python) | RDKit.js (WASM) |
|---|---|---|---|
| **快速上手** | `pip install chematic` | `pip install rdkit`（官方预编译 wheel）或 conda | `npm install @rdkit/rdkit`，无 Python 绑定 |
| **浏览器包体积** | **raw 3.58 MB / gzip 1.31 MB** | 不适用（Python/C++ 库） | raw 6.91 MB* |
| **ECFP4批量** | **54.7 µs/mol** | 94.3 µs/mol | — |
| **Canonical SMILES** | **24.95 / 18.27 µs/mol** | 25.58 / 26.82 µs/mol | — |
| **SDF graph read / serialization-only write** | **9.48 / 7.62 µs/mol** | 99.96 / 79.54 µs/mol | — |
| **内存安全性** | 编译器保证（Rust） | C++ | C++ |
| **源码构建** | 仅需 `cargo build` | cmake + clang + Boost | Emscripten SDK |

\* RDKit.js 的 gzip 传输体积未独立测量，此处以 raw 体积做同口径比较。RDKit.js 目前处于
维护者交接阶段（详见其仓库）。

canonical/SDF 行是 2026-09-04 macOS arm64 的中位数，仅适用于所记录的语料和
操作边界。参阅[基准测试详情](https://kent-tokyo.github.io/chematic/benchmark/)。
chematic WASM 包体积于 2026-09-06 使用 `wasm-pack 0.13.1` + `wasm-opt 130`
从 v1.0.8 release candidate 构建并测量：**raw 3.58 MB**（**gzip 1.31 MB**）。固定的历史比较项为
RDKit.js **6.91 MB**
（`@rdkit/rdkit@2025.3.4-1.0.0` 的 `RDKit_minimal.wasm`，经 unpkg.com 确认）· Indigo（Ketcher
构建版）**11.24 MB**（`indigo-ketcher@1.45.1` 的主 `.wasm`，经 jsDelivr 确认）—— 以 raw 对 raw
比较，chematic 目前比 RDKit.js 小约 2.1 倍，比 Indigo 的 Ketcher 构建版小约 3.8 倍。
详见 [v1.0.8 artifact 记录](benchmarks/2026-09-06-wasm-size-v1.0.8.md)。

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

- 需要在浏览器中运行化学计算（WASM，1.29 MB gzip，无需服务器）
- 需要纯 Rust 技术栈，不依赖 C++ 工具链
- 部署到难以安装或不支持 RDKit 的环境（Cloudflare Workers、Lambda、嵌入式设备 ——
  RDKit 本身已提供官方 `pip install rdkit` wheel，但仍需标准 CPython 环境）
- 构建 AI 代理并需要原生 MCP 工具集成
- 需要批量高吞吐量处理分子（ECFP4：比 RDKit 快 2–3 倍，Rayon 并行）
- 希望 `pip install chematic` 在任何环境都能直接使用，无需编译器

**适合使用 RDKit 的场景：**

- 需要最大的生态兼容性和 20 年以上的生产验证
- 需要带 ML 辅助扭转角修正的出版级 3D 结构（RDKit 的 ETKDGv3）
- 需要在不启用 `native-inchi` feature 的情况下获得逐位精确的标准 InChI
- 依赖基于 RDKit Python API 编写的社区插件

---

## 选择你的接口

- [Rust](#快速开始)
- [Python](#快速开始)
- [WebAssembly / Node.js](README.md#javascript--typescript-webassembly)
- [材料与模拟格式](docs/format-capabilities.md) — mmCIF、PQR、QCSchema、ORCA、Gaussian Cube、OpenDX、LAMMPS
- [从 RDKit 迁移](docs/rdkit-migration.md) — 按功能划分的 Supported / Partial / Not-supported 对照表

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

chematic 提供用于本地AI代理集成的 **MCP（模型上下文协议）服务器**。

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

**Protocol**：在同一个 stdio 连接上同时支持旧版（`2024-11-05` 风格的 `initialize` 握手）和 MCP `2026-07-28` 无状态方言（`server/discover`、逐请求 `_meta`、可缓存的 `tools/list`、`structuredContent`）。详见 [`chematic-mcp` README](crates/chematic-mcp/README.md#protocol-eras)。Remote HTTP、OAuth、Tasks 扩展与 MCP Apps 仍不支持。

---

## 为什么选择纯 Rust？

常用化学核心采用safe Rust，公开的非可信输入路径具有有限默认值与typed error。`native-inchi`是使用IUPAC InChI C library的显式opt-in FFI例外；依赖crate自身的unsafe code属于单独边界。详见[SECURITY](SECURITY.md)。

历史v0.18.0 ECFP4批处理中位数，在同一5,000分子语料与Apple M4环境下为54.7 µs/mol，v1.0.8候选WASM artifact为3.58 MB raw / 1.31 MB gzip。两者均为固定日期与条件的测量，不是普遍性能保证。

## 比较原则

| 方面 | chematic | RDKit |
|---|---|---|
| 部署 | pure-Rust core、Python wheel、WASM/Node | 最广泛的功能与生态系统 |
| 浏览器 | native WASM | RDKit.js是独立distribution |
| 3D/MMFF94 | Experimental | 成熟的ETKDG/force-field workflow |
| 兼容性 | 命名subset与typed failure | 参考实现 |

功能差异与不支持范围见[RDKit比较](docs/rdkit-comparison.md)和[兼容性范围](docs/compatibility-scope.md)，速度见[benchmark](docs/benchmark.md)。

---

## JavaScript / TypeScript（WebAssembly）

**1.31 MB gzip — 与 RDKit.js 的 raw WASM 相比约小 2.0 倍。** 无需 Emscripten 或 cmake，可直接在浏览器和 Node.js 中使用。

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

218+ 个导出函数（另加 `MolHandle`/`DepictOptions` 类方法，2026-08-21 测量）涵盖描述符、指纹、3D 几何、反应、多样性筛选和 SDF 处理。
完整 API 请参阅 [WASM API 参考文档](https://kent-tokyo.github.io/chematic/)。
---

## Crate 列表

| 领域 | crate |
|---|---|
| 分子图与结构标识 | `chematic-core`, `chematic-smiles`, `chematic-perception`, `chematic-cip` |
| 查询、描述符与指纹 | `chematic-smarts`, `chematic-chem`, `chematic-fp` |
| 文件与反应 | `chematic-mol`, `chematic-rxn`, `chematic-inchi`, `chematic-iupac` |
| 2D、3D 与材料 | `chematic-depict`, `chematic-3d`, `chematic-ff`, `chematic-crystal`, `chematic-ewald` |
| 用户接口 | `chematic`, `chematic-py`, `chematic-wasm`, `chematic-cli`, `chematic-mcp` |

详细信息请参阅[格式支持](docs/format-capabilities.md)、[语言绑定](docs/language-bindings.md)和各crate的README。

---

## 近期开发

**未发布:** 解决 #210 的legacy UFF stereo rescue残差，并继续优化canonical SMILES与SDF hot path。固定条件的测量记录在[benchmarks](benchmarks/)中。

**v1.0.8（2026-09-06）:** 延续 v1.0.7 的 descriptor provenance 与跨绑定契约，增加 ECFP4/MACCS 形状契约以及 `PeriodicStructure::identity_bytes()` 的确定性 identity 序列化。#149/#337 残差仍固定为 fail-closed 或诊断专用。Spectrophores 在完成独立 patent/FTO 审查前不属于公共 API。

公开发布摘要见[CHANGELOG](CHANGELOG.md)；详细开发记录保存在其中链接的archive。

---

## 已知限制

- `canonical_smiles()` 是一种表示形式，不应直接用作dedup/cache key；请使用可返回`None`的`canonical_smiles_stable_key()`。
- Aromaticity/CIP明确区分default与opt-in model，不宣称全面RDKit parity。
- 3D生成和MMFF94仍为Experimental。
- Python `RWMol`、CDXML编辑与Markush/polymer expansion是有意限制的bounded subset。
- pure-Rust InChI为近似实现；标准IUPAC InChI需启用`native-inchi`。

精确契约请参阅[兼容性范围](docs/compatibility-scope.md)、[验证](docs/validation.md)与[错误和资源限制](docs/error-and-limits.md)。

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
│   ├── chematic-crystal/    周期晶体结构：晶格、PBC、近邻、supercell（不依赖 Molecule）
│   └── chematic/            带功能标志的伞形 crate
```

---

## 许可证

可选择 Apache License 2.0 或 MIT License 中的任意一种。
版权归属：Kentaro Tanabe (kent-tokyo)。再分发时的归属声明请参阅 [`NOTICE`](NOTICE)。

---

如果 chematic 对您有帮助，欢迎给项目一个 [GitHub star](https://github.com/kent-tokyo/chematic)，这将帮助更多人发现它。
