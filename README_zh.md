# chematic

[English](README.md) | [日本語](README_ja.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![许可证](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)

面向 Python、Rust 和浏览器的化学信息学库。核心采用纯 Rust，提供输入限制、
类型化错误和确定性的批处理 API。

[官方网站](https://chematic.io/) · [文档](https://kent-tokyo.github.io/chematic/) · [在线演示](https://kent-tokyo.github.io/chematic/playground/)

## 安装

```bash
pip install chematic
cargo add chematic --features "smiles,perception,chem,3d,fp"
npm install @kent-tokyo/chematic
```

Python 版本无需 C/C++ 编译器。

## Python

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
print(mol.mw, mol.logp, mol.tpsa)
print(mol.has_substructure("[OH]"))

report = chematic.report([mol], names=["aspirin"])
report.save("report.html")
```

也可以使用有限的 RDKit 兼容接口：

```python
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import Descriptors

mol = Chem.MolFromSmiles("CCO")
print(Descriptors.MolWt(mol))
```

这不是完整的 RDKit 克隆；不支持的选项会明确报错。详见
[RDKit 迁移指南](docs/rdkit-migration.md)。

## Rust

```rust
use chematic::smiles::parse;

let mol = parse("c1ccccc1")?;
println!("{}", mol.atom_count());
# Ok::<(), Box<dyn std::error::Error>>(())
```

SMILES/SMARTS、描述符、指纹、反应、SDF/MOL/CDXML、2D/3D及晶体格式等功能均提供为
独立 crate。参阅[格式支持矩阵](docs/format-capabilities.md)。

## JavaScript / WebAssembly

```js
import init, { parse_smiles, get_descriptors_json } from "@kent-tokyo/chematic";

await init();
const mol = parse_smiles("CCO");
console.log(mol.atom_count());
console.log(JSON.parse(get_descriptors_json(mol)));
```

WASM 版本支持分子解析、描述符、指纹、反应、部分 2D/3D 操作和格式转换。
详见 [WASM README](crates/chematic-wasm/README.md)。

## 稳定性和限制

稳定功能包括 SMILES/SMARTS、描述符、指纹、相似度与子结构搜索、SDF/MOL I/O，
以及 Rust/Python/Node/WASM 绑定。

3D、pKa/ADMET、IUPAC 名称、Markush/polymer 展开、CDXML 编辑和 RDKit 兼容接口属于
实验性功能或 bounded subset。`canonical_smiles()` 不一定适合作为去重或缓存键；
需要时请使用会 fail-closed 的 `canonical_smiles_stable_key()`。

精确契约请参阅[兼容性范围](docs/compatibility-scope.md)、[验证](docs/validation.md)
以及[错误和资源限制](docs/error-and-limits.md)。

## MCP 服务

`chematic-mcp` 为支持 MCP 的代理提供本地 stdio 服务。详见
[chematic-mcp README](crates/chematic-mcp/README.md)。

## 开发

```bash
cargo build --workspace
cargo test --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

基准测试见 [benchmark guide](docs/benchmark.md)，版本历史见 [CHANGELOG](CHANGELOG.md)。

## 许可证

可选择 Apache License 2.0 或 MIT License。归属信息请参阅 [NOTICE](NOTICE)。
