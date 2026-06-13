# 变更日志（中文）

All notable changes to chematic will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

v0.1.8 之前的变更历史，请参考 [CHANGELOG.md](CHANGELOG.md)。

---

## [Unreleased]

---

## [0.1.96] — 2026-06-13

### Improved — MMFF94 BCI 部分电荷

- 将 `mmff94_charges()` 从电负性近似（±0.5e）替换为基于 Bond Charge Increment (BCI) 表的模型（±0.1e）。
- 新模块 `chematic-chem/src/mmff94_bci.rs`：包含 25 个 BCI 条目，涵盖 C=O (+0.47)、C–O (+0.04)、C–N (+0.10)、C–N(酰胺) (+0.31)、N–H (-0.16)、O–H (-0.30)、C–F (+0.22) 等（Halgren 1996）。
- 保证电荷守恒：`sum(q) == sum(形式电荷)`。
- 新增测试 +9（乙醇、丙酮、甲胺、乙酸、铵离子、乙酸根、氯甲烷、咪唑、酰胺 vs 胺 BCI 对比）。
- 测试数量：1,657 → 1,666 (+9)，全部通过。

---

## [0.1.95] — 2026-06-13

### Added — `chematic-iupac` 本地命名范围扩展

- 支持带位置编号的酮：`propan-2-one`、`butan-2-one`、`pentan-3-one`。
- 支持羧酸、酯和第一/第二酰胺：`ethanoic acid`、`methyl ethanoate`、`ethanamide`。
- 支持未取代苯和常见芳香杂环：吡啶、呋喃、噻吩、吡咯、咪唑、嘧啶。
- IUPAC 单元测试从 8 个扩展到 14 个。

### Fixed — CI Clippy 兼容性

- `cargo clippy --workspace -- -D warnings` 已通过当前 stable Clippy。
- 保留 crate root 的 deprecated `total_hcount` 兼容 re-export，同时只在该 re-export 上抑制警告。
- 按新 Clippy lint 更新 ECFP 迭代循环、condensed formula guard 和 DG 坐标生成器文档注释。

### Improved — 真指纹算法 (A4/A5/A6)

- **MHFP 规范化哈希**：将依赖原子索引的字节签名替换为 Morgan 式循环片段哈希。对同一分子的不同 SMILES 顺序输入，均生成相同指纹。新增测试 +3。
- **ERG 药效团节点类型**：新增 `assign_pharmacophore_features()`。正确赋予 DONOR (N-H/O-H)、ACCEPTOR (无 H 的吡啶型 N、O、F)、POSITIVE/NEGATIVE（形式电荷）、HYDROPHOBIC（纯 C/H 基团）。区分吡啶 N（受体）和吡咯 N（供体）。新增测试 +5。
- **Reaction FP**：确认 `use_xor: true`（XOR 差分）为默认值，比较表从 ⚠️ 更新为 ✅。
- 测试数量：1,649 → 1,657 (+8)，全部通过。

---

## [0.1.37] — 2026-06-08

### Added — mol_transforms API + 随机 SMILES 生成

#### `chematic-3d` — 分子几何操作

**NEW**: 结合长/角度/二面角的测量和变换的公开 API:
- `get_bond_length(coords, a, b) -> f64` — 结合长（Ångströms）
- `get_bond_angle(coords, a, center, b) -> f64` — 角度（弧度）
- `get_bond_angle_deg(...)` — 角度（度）
- `get_dihedral(coords, a, b, c, d) -> Option<f64>` — 二面角（弧度）
- `get_dihedral_deg(...)` — 二面角（度）
- `set_dihedral(coords, mol, a, b, c, d, angle_rad) -> Coords3D` — 旋转 D 侧子树
- `compute_centroid(coords) -> [f64; 3]` — 所有原子的重心
- `center_on_origin(coords) -> Coords3D` — 平移到原点
- `transform_conformer(coords, 4x4_matrix) -> Coords3D` — 应用 4×4 齐次变换

**公开内部函数**:
- `dihedral()`、`compute_angle()`、`rotate_around_axis()` （先前为 private）
- 新文件: `crates/chematic-3d/src/mol_transforms.rs`

#### `chematic-smiles` — SMILES 多样性生成

**NEW**: 用于 ML 数据扩增的随机 SMILES 生成:
- `random_smiles(mol, seed) -> String` — 使用 xorshift64 RNG permute 原子顺序
- `random_smiles_vect(mol, count, seed) -> Vec<String>` — 生成 N 个唯一变体
- 算法：Fisher-Yates 洗牌
- 用例：数据扩增（同一分子的多种 SMILES 表示）

#### `chematic-wasm` — mol_transforms + 随机 SMILES 的 WASM 绑定

**NEW**: 4 个新 WASM 导出函数:
- `get_bond_length_json(smiles, a, b) -> f64` — 从 SMILES 获取结合长
- `get_dihedral_json(smiles, a, b, c, d) -> JsValue` — 从 SMILES 获取二面角（度）
- `set_dihedral_json(smiles, a, b, c, d, angle_deg) -> Result<String>` — 返回 PDB 块
- `random_smiles_json(smiles, count, seed) -> Result<String>` — SMILES 的 JSON 数组

### 测试覆盖

- **chematic-3d**: +5 mol_transforms 测试（结合长、角度、重心、矩阵变换）
- **chematic-smiles**: +6 random_smiles 测试（确定性、唯一性、往返）
- **总计**: 1,120 → 1,151 (+31 新测试)
- 所有测试通过 ✅

---

## [0.1.36] — 2026-06-08

### Fixed — Issue #1 审计：拓扑正确但化学无意义的结果

#### `chematic-smarts` — VF2 和 MCS 正确性修复

**BUG-2: SMARTS `[h]` 原语（隐式氢计数）**
- **修复**: 解析器现在能够正确区分 `[H]`/`[H2]`（总氢数）和 `[h]`/`[h2]`（仅隐式氢）
- 向 query.rs 添加 `AtomPrimitive::ImplicitHCount(u8)` 变量
- 修复 parser.rs，在元素回退之前处理小写 `h`
- 向 `eval_atom_primitive()` 添加使用 `implicit_hcount()` 函数的匹配分支
- **影响**: 防止芳香族 H 被错误地匹配到显式原子的无声不匹配

**BUG-3: MCS `maximize_bonds` 平局破解**
- **修复**: 修改 mcs.rs 中的 `grow()` 函数，当原子数相等时使用结合数作为平局破解条件
- 添加条件：`|| (maximize_bonds && mapping.size == best.size && mapping.bond_count > best.bond_count)`
- 默认 `maximize_bonds=true` 以匹配 RDKit 行为
- **影响**: 当存在多个同等大小的匹配时，MCS 现在返回一致的结果

**BUG-4: SMARTS `/\` 几何立体键（E/Z）**
- **修复**: 向 `BondPrimitive` enum 添加 `Up` 和 `Down` 变量以支持几何立体化学
- 修改 parser.rs：`is_bond_token()` 现在识别 `/` 和 `\` 字符
- 修改 `consume_bond_prim()` 以将 `/` 解析为 `Up`，`\` 解析为 `Down`
- 向 `eval_bond_primitive()` 添加 `BondPrimitive::Up` 和 `Down` 匹配分支
- **影响**: 像 `/C=C\` 这样的 SMARTS 查询现在可以正确匹配 E/Z 构型的双键

### 测试覆盖

- **chematic-smarts**: 124 个测试全部通过（包括 BUG-2/3/4 验证）
- **工作区**: 1,120+ 个测试全部通过
- clippy 无警告

### 说明

- **Issue #1 模式**: 审计发现 3 个错误，其中算法拓扑正确但返回化学无意义的结果
  - 根本原因：RDKit 有约束选项，但 chematic 中没有实现
  - 示例：VF2 `uniquify`（v0.1.33 修复）、MCS 环感知（v0.1.22 修复）
  - 本 Sprint 添加了 3 个缺失的正确性约束

---

## [0.1.32] — 2026-06-07

### Added — 3D 几何、坐标处理、WASM 稳定性增强

#### `chematic-3d` — 距离几何约束满足

- **NEW**: `build_constraints()` & `satisfy_constraints()` 用于迭代约束投影
  - 强制理想键长（±0.05 Å 容差）和原子价角（±5° 容差）
  - 典型分子在 5-10 次迭代内收敛
  - 每次迭代 O(n²)；适用于小到中等大小分子（< 1000 原子）
- **NEW**: `generate_and_minimize_constrained()` 高级 API：DG → 约束 → DREIDING
  - 改进受应力/问题分子的几何质量
- **NEW**: BondConstraint & AngleConstraint 结构体及强制方法

#### `chematic-mol` — V3000 坐标恢复

- **NEW**: `parse_mol_v3000_with_coords()` 函数返回 (Molecule, MolMetadata, Vec<(f64, f64)>)
  - 从 MOL V3000 原子块恢复 2D 坐标（以前被丢弃）
  - 匹配 V2000 `parse_mol_with_coords()` API 模式
  - 使往返旅行 2D 坐标保存成为可能

#### `chematic-perception` — 芳香族性模型细化

- **NEW**: `RingAromaticity` enum 区分 Aromatic/Antiaromatic/NonAromatic
- **NEW**: Hückel 4n+2 规则及反芳香族性（4n）检测
  - C、N（H 依赖）、O、S 原子的 π 电子计数
  - `ring_classifications()`、`antiaromatic_rings()`、`has_antiaromaticity()` 方法
  - 检测奇异体系：环丁二烯、环八四烯、杜瓦苯

### Fixed — WASM 构建可靠性及坐标系清晰度

#### `chematic-3d` — WASM RNG 种子设定（MD 必需）

- **CRITICAL WASM**: `fastrand = "2.4"` 现在为 wasm32-unknown-unknown 目标指定 `features = ["js"]`
  - 修复 MD 速度初始化以使用密码学随机性而不是固定种子
  - WASM 构建现在产生非确定性（物理上有意义的）轨迹
  - 原生构建不受影响；feature 仅对 WASM 目标激活

#### `chematic-depict`、`chematic-mol` — 坐标系文档

- **澄清**: `compute_layout()` 生成 SVG Y-down 像素坐标（不是化学 Y-up）
- **澄清**: `parse_cml()` 返回化学 Y-up 惯例；调用者必须为 SVG 渲染否定 Y
- **澄清**: `parse_cdxml()` 返回 ChemDraw Y-down（SVG 兼容，无需转换）
- 添加了详细的文档字符串注释以防止坐标系错误

### Changed — 错误处理完整性

#### 13 个错误类型现在实现 `std::error::Error` 特征

- **高优先级**: `SmartsError`、`ValenceError`、`StereoError` — 添加 Display + Error
- **其他**: `CmlError`、`CdxmlError`、`Mol2Error`、`RxnParseError`、`MolError`、`IupacError`、`ConformerError`、`RxnError`、`TransformError` — 添加 Error 特征
- 启用标准错误处理模式（`.source()`、`Box<dyn Error>` 等）

### 测试覆盖

- **+ 44 个新测试**：约束满足（12）、芳香族性（16）、V3000 坐标（2）、错误类型（14）
- **总计**: 全部 crate 171 个测试通过
- 已验证 WASM 构建（js feature 无回归）

---

## [0.1.30] — 2026-06-07

### Fixed — 关键物理学和静电学更正

#### `chematic-ewald` — PME 网格索引错误

- **CRITICAL**: 修复了 `interpolate_charges_to_mesh()` 中的 3D 到 1D 网格索引计算，该计算对于非立方网格丢失了 97% 的电荷数据
- 将错误的 `isqrt()` 近似替换为正确的 3D 索引公式：`linear_idx = ix + iy*M0 + iz*M0*M1`
- 现在在整个倒易空间网格中正确分布电荷

#### `chematic-3d` — 分子动力学和力场

- **CRITICAL**: 使用正确的单位转换因子（0.01038 kcal/mol → amu·Ų/fs²）修复了 Maxwell-Boltzmann 速度初始化
  - 之前的代码生成的速度小 347 倍，导致初始动能不正确
  - 现在生成与目标温度相匹配的物理正确热速度
  
- **CRITICAL**: 修复了 VDW 能量计算以使用 DREIDING 参数而不是硬编码值
  - 将 `r_eq = 2.0 Å`（所有对）替换为特定元素的 DREIDING VDW 半径
  - 实现了完整的 Lennard-Jones 12-6 势，包括吸引力分散项
  - 实现了原子对相互作用的 Lorentz-Berthelot 结合规则
  - 实现了 1-2 和 1-3 键合对排斥

---
