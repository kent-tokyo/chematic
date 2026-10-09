# chematic

[English](README.md) | [日本語](README_ja.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![Coverage](https://github.com/kent-tokyo/chematic/actions/workflows/coverage.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/kent-tokyo/chematic/branch/main/graph/badge.svg)](https://app.codecov.io/github/kent-tokyo/chematic)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![许可证](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)
[![Ask DeepWiki](docs/assets/ask-deepwiki.svg)](https://deepwiki.com/kent-tokyo/chematic)

面向 Python、Rust 和浏览器的纯 Rust 化学信息学工具包：本地优先、输入受限，并为不支持的
化学表达提供明确错误。

[官方网站](https://chematic.io/) · [文档](https://kent-tokyo.github.io/chematic/) · [Playground](https://kent-tokyo.github.io/chematic/playground/)

## 安装

```bash
pip install chematic
cargo add chematic --features "smiles,perception,chem,3d,fp"
npm install @kent-tokyo/chematic
```

Python wheel 无需 C/C++ 编译器；各绑定共享同一 Rust 内核。

### v1.0.39 范围

v1.0.39 增加了默认仅绑定回环地址的 MCP Streamable HTTP，并改进了 RDKit
兼容配置下的 SMARTS／SMIRKS 解析、元素替换和芳香生成物的高价态处理。同时新增
SMIRKS property-based fuzzing、按平台固定的比较策略和可提交给上游的探针。
在 v1.0.39 发布包重新测量前，v1.0.38 的发布包结果仍作为历史记录。详情见
[验证报告](docs/validation.md)，变更见 [CHANGELOG](CHANGELOG.md)。

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

- [入门与 Cookbook](https://kent-tokyo.github.io/chematic/)
- [兼容性范围](docs/compatibility-scope.md) 与 [RDKit 迁移](docs/rdkit-migration.md)
- [验证](docs/validation.md) 与 [基准方法](docs/benchmark.md)
- [格式和绑定](docs/format-capabilities.md)
- [MCP 服务](crates/chematic-mcp/README.md)

## 开发

```bash
cargo test --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

## 许可证

采用 Apache-2.0 或 MIT。归属信息见 [NOTICE](NOTICE)。
