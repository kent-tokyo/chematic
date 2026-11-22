# chematic

[English](README.md) | [日本語](README_ja.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg?branch=main&event=push)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/github/deployments/kent-tokyo/chematic/github-pages?label=docs)](https://kent-tokyo.github.io/chematic/)
[![Coverage](https://github.com/kent-tokyo/chematic/actions/workflows/coverage.yml/badge.svg?branch=main&event=push)](https://github.com/kent-tokyo/chematic/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/kent-tokyo/chematic/branch/main/graph/badge.svg)](https://app.codecov.io/github/kent-tokyo/chematic)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![许可证](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)
[![Ask DeepWiki](docs/assets/ask-deepwiki.svg)](https://deepwiki.com/kent-tokyo/chematic)

chematic 是用 Rust 编写的化学信息学库，可用于分子结构的搜索、分析和绘制。
Python、Rust、Node.js 和浏览器应用共用同一个 Rust 内核。

[官方网站](https://chematic.io/) · [文档](https://kent-tokyo.github.io/chematic/) · [Playground](https://kent-tokyo.github.io/chematic/playground/)

## 安装

```bash
pip install chematic
cargo add chematic --features "smiles,perception,chem,3d,fp"
npm install @kent-tokyo/chematic
```

安装 Python wheel 无需 C/C++ 编译器。在浏览器中，WebAssembly 可直接在用户设备上处理分子数据。

### v1.2.3 范围

当前版本为 **v1.2.3**。参见[更新记录](https://github.com/kent-tokyo/chematic/releases/tag/v1.2.3)和[发布验证](docs/reference/validation.md)。下一版本的变更见[未发布条目](CHANGELOG.md#unreleased)。

## 分子搜索与分析

- 从 SMILES 或 SDF 读取分子，计算分子量、LogP 和拓扑极性表面积（TPSA）。
- 用 SMARTS 搜索指定的子结构，或用 Morgan 指纹和 Tanimoto 系数查找结构相似的分子。
- 求出多个分子的最大公共子结构（MCS）。
- 将分子结构绘制为 SVG，用于网页和报告。

使用已发布的 Python 包对同一组 1,000 个分子进行比较时，chematic 1.2.2 与 RDKit 2026.9.2
对每个分子给出的环数和 Morgan 指纹（半径 2，2,048 位）均一致，所有 TPSA 值的差异也都在指定的容差范围内。[比较方法与结果](benchmarks/2026-11-21-v122-published-comparison.md)

## 使用

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
print(mol.mw, mol.logp, mol.tpsa)
print(mol.has_substructure("[OH]"))
```

```rust
use chematic::smiles::parse;

let mol = parse("c1ccccc1").expect("valid SMILES");
println!("{}", mol.atom_count());
```

```js
import init, { parse_smiles, get_descriptors_json } from "@kent-tokyo/chematic";

await init();
console.log(JSON.parse(get_descriptors_json(parse_smiles("CCO"))));
```

`rdkit_compat` 仅覆盖已文档化的子集；不支持的选项会明确失败而非静默近似。

## 范围

SMILES/SMARTS、描述符、指纹、相似度和子结构搜索、常用 MOL/SDF I/O 以及已文档化的
Rust/Python/Node/WASM 绑定属于稳定范围。3D、pKa/ADMET、IUPAC、Markush/polymer、
复杂 CDXML 编辑和 RDKit 风格 API 仍为实验性或有边界的功能。

`canonical_smiles()` 不是通用的身份键；需要时只应在文档化范围内使用 fail-closed 的
`canonical_smiles_stable_key()`。

## 指南

- [选择指南](docs/reference/choosing-chematic.md)
- [入门与 Cookbook](https://kent-tokyo.github.io/chematic/)
- [兼容性范围](docs/reference/compatibility-scope.md) 与 [RDKit 迁移](docs/reference/rdkit-migration.md)
- [验证](docs/reference/validation.md) 与 [基准方法](docs/reference/benchmark.md)
- [格式和绑定](docs/reference/format-capabilities.md)
- [MCP 服务](crates/chematic-mcp/README.md)

## 开发

```bash
cargo test --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

## 许可证

采用 Apache-2.0 或 MIT。归属信息见 [NOTICE](NOTICE)。
