# chematic

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![docs.rs](https://docs.rs/chematic/badge.svg)](https://docs.rs/chematic)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)

[Open in Colab](https://colab.research.google.com/github/kent-tokyo/chematic/blob/main/notebooks/quickstart.ipynb)

[日本語](README_ja.md) | [中文](README_zh.md)

A cheminformatics library for Python, Rust, and the browser.

**Cheminformatics that's fast by default, safe by design.**  
Pure Rust · Zero C/C++ · Python · WebAssembly · [Website](https://chematic.io/) · [Live Demo](https://kent-tokyo.github.io/chematic/playground/)

| | chematic | RDKit (Python) | RDKit.js (WASM) |
|---|---|---|---|
| **Get started** | `pip install chematic` | `pip install rdkit` (official prebuilt wheels) or conda | `npm install @rdkit/rdkit`, no Python bindings |
| **Browser bundle** | **2.94 MB raw / 1.10 MB gzip** | not applicable (Python/C++ library) | 6.91 MB raw* |
| **Batch fingerprints** | **~78 µs/mol** (2–3× faster) | ~160–235 µs/mol | — |
| **Memory safety** | compiler-enforced (Rust) | C++ | C++ |
| **Build from source** | `cargo build` only | cmake + clang + Boost | Emscripten SDK |

\* RDKit.js gzip-over-the-wire size was not independently measured; raw figures are compared
on a like-for-like basis. RDKit.js is currently in a maintainer transition (see its repo for
current status).

All numbers are reproducible — see [benchmark details](https://kent-tokyo.github.io/chematic/benchmark/).  
WASM sizes (raw, measured 2026-08-21 from a clean `wasm-pack build --target web --release`
+ `wasm-opt -O3`, commit `ef7dc25`): chematic **2.94 MB** (**1.10 MB gzip**) · RDKit.js **6.91 MB**
(`@rdkit/rdkit@2025.3.4-1.0.0`'s `RDKit_minimal.wasm`, via unpkg.com) · Indigo (Ketcher build)
**11.24 MB** (`indigo-ketcher@1.45.1`'s main `.wasm`, via jsDelivr) — chematic's raw WASM binary
is currently about 2.3× smaller than RDKit.js's and about 3.8× smaller than Indigo's Ketcher-oriented
build, on a raw-to-raw basis.

The separate 2026-08-23 benchmark rebuild reports 2.98 MB raw / 1.11 MB gzip;
both figures are retained with their measurement dates because build outputs
can vary slightly by toolchain and build environment.

**Feature maturity at a glance:**

| Feature | Status |
|---|---|
| SMILES / SMARTS / fingerprints / descriptors | Stable |
| 3D conformer generation (DG + MMFF94) | Experimental |
| pKa / ADMET | Rule-based screening (not for clinical use) |
| IUPAC name generation | Partial (25+ classes) |
| Pure-Rust InChI | Approximate (enable `native-inchi` feature for exact) |

---

## What you get

```
$ python -c "import chematic; print(chematic.from_smiles('CC(=O)Oc1ccccc1C(=O)O').describe())"
Molecular weight 180.2 Da, formula C9H8O4.
LogP 1.31 (mildly lipophilic), TPSA 63.6 Å².
HBD 1, HBA 3, 3 rotatable bond(s), 1 aromatic ring(s).
Drug-likeness: no Lipinski rule-of-5 violations. likely orally bioavailable (passes Veber criteria).
QED 0.56 (0 = non-drug-like, 1 = ideal).
Structural alerts: Brenk alert.
```

One `pip install`. No RDKit, no conda, no C compiler. Works in Python, Rust, the browser, and AI agents.

```python
# HTML report — self-contained, opens in any browser and renders in Jupyter
mols = [chematic.from_smiles(s) for s in smiles_list]
report = chematic.report(mols, names=compound_names)
report.save("report.html")   # or: display(report) in Jupyter

# Side-by-side comparison
cmp = chematic.compare(aspirin, ibuprofen, names=("Aspirin", "Ibuprofen"))
cmp.save("compare.html")
```

---

## Common Use Cases

| Scenario | How chematic helps |
|---|---|
| **HTML report** | `chematic.report(mols, output="report.html")` — self-contained compound grid, no server needed |
| **Drug screening** | 190+ descriptors, ADMET, PAINS/Brenk, QED — batch over thousands of compounds |
| **Molecule search** | ECFP4/MACCS fingerprints, opt-in RDKit-compatible chiral Morgan fingerprints, Tanimoto, LSH approximate nearest-neighbour |
| **AI agent / MCP** | Built-in MCP server — Claude Desktop can call chemistry tools directly |
| **Browser app** | 1.10 MB gzip WASM bundle, zero backend required, React/Vue/Svelte ready |
| **Jupyter notebook** | `mol` renders SVG inline; `descriptors_df()` returns a pandas DataFrame |
| **Batch analysis** | Rayon-parallel descriptor/fingerprint/3D pipelines; SDF/CSV in, CSV out |
| **Rust server** | Pure-Rust crates with no C/C++ toolchain; Axum/Actix compatible |

Full worked examples → [Use cases](https://kent-tokyo.github.io/chematic/use-cases/)

---

## When to use chematic

**Use chematic if:**

- You want chemistry in the browser (WASM, 1.10 MB gzip, no server required)
- You need a pure Rust stack with no C++ toolchain dependencies
- You deploy to environments where installing RDKit is impractical or unsupported (Cloudflare Workers, Lambda, embedded — RDKit itself ships official `pip install rdkit` wheels, but those still assume a standard CPython environment)
- You build AI agents and want native MCP tool integration
- You process molecules in batch at high throughput (ECFP4: 2–3× faster than RDKit, Rayon-parallel)
- You want `pip install chematic` to just work — anywhere, no compiler needed

**Use RDKit if:**

- You need maximum ecosystem compatibility and 20+ years of production validation
- You need publication-quality 3D structures with ML-assisted torsion corrections (RDKit's ETKDGv3)
- You need bit-exact standard InChI without enabling the `native-inchi` feature
- You depend on community plugins written against the RDKit Python API

---

## Choose your interface

- [Rust](#quick-start)
- [Python](#quick-start)
- [WebAssembly / Node.js](#javascript--typescript-webassembly)
- [Materials and simulation formats](docs/format-capabilities.md) — mmCIF, PQR, QCSchema, ORCA, Gaussian Cube, OpenDX, LAMMPS
- [Migrating from RDKit](docs/rdkit-migration.md) — feature-by-feature Supported/Partial/Not-supported breakdown
- [Compatibility scope](docs/compatibility-scope.md) — final v0.89 boundary for RWMol, CDXML, polymer expansion, and RDKit/Morgan compatibility

---

## Quick Start

### Installation

```bash
# Python — no C/C++ compiler required
pip install chematic

# Rust
cargo add chematic --features "smiles,perception,chem,3d,fp"

# JavaScript/TypeScript
npm install @kent-tokyo/chematic
```

### Python

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin

# In Jupyter, type `mol` in a cell — 2D structure renders automatically
mol

# Access 190+ descriptors as properties
print(mol.mw, mol.logp, mol.tpsa)           # 180.16  1.31  63.6
print(mol.lipinski_passes, mol.pains_passes) # True   True

# Substructure search
mol.has_substructure("[OH]")   # True
mol.find_matches("[CX3](=O)O") # → [[1, 2, 3], [7, 8, 9]]

# Natural-language summary (one paragraph)
print(mol.describe())

# Structured Markdown report — paste into LLM, Jupyter, or save as .md
print(mol.review())
# → # Molecular Review\n## Structure\n## Physical Properties\n## Drug-likeness\n## ADMET...

# Structural diff between two molecules
ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")
d = mol.diff(ibuprofen)  # {"summary": "+C7, -O2. ΔLogP +2.75 ...", "delta_mw": 66.1, ...}

# Batch processing — parallel, numpy-ready
fps = chematic.bulk.ecfp4(["CCO", "c1ccccc1", "CC(=O)O"])  # (3, 2048) uint8

# One-liner DataFrame
df = chematic.descriptors_df(["CCO", "c1ccccc1", "CC(=O)O"])
df[["mw", "logp", "tpsa", "qed"]]
```

For Rust and JavaScript/TypeScript examples, see the [documentation](https://kent-tokyo.github.io/chematic/).

### Migrating from RDKit

`chematic.rdkit_compat` provides a lightweight RDKit-compatible subset so existing scripts port with minimal changes:

```python
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import Descriptors, rdMolDescriptors, DataStructs

mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
Descriptors.MolWt(mol)                       # 180.16
fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
DataStructs.TanimotoSimilarity(fp, fp)       # 1.0
```

It is **not a full RDKit clone**, and unsupported options fail loudly. See the
[RDKit migration guide](docs/rdkit-migration.md) for the compatibility matrix,
differential-validation results vs RDKit, and runnable examples.

### Diagnostics

```python
import chematic
chematic.doctor()
# chematic v0.89.0 candidate
# Python 3.12.x  |  darwin arm64
#
# Descriptor accuracy (benchmark 2026-07-17, v0.4.29 vs RDKit 2026.03.3 --
# descriptor calculation paths unchanged through v0.8.0, not re-measured since):
#   MW / HBA / HBD / ARC  100%   (4,999-mol ChEMBL subset)
#   TPSA                  100%   within ±0.1 Å²
#   LogP (Crippen)        100%*  (max Δ = 1.1×10⁻¹³)
#   Stereocenter count    99.96% (legacy) / 98.6% (new CIP FindPotentialStereo)
#   CIP R/S label         96.30% vs modern rdCIPLabeler (96.83% vs legacy)
# ...
```

---

## For AI / LLM Developers

chematic ships a native **MCP (Model Context Protocol) server** — the first cheminformatics library with built-in AI agent integration.

```json
// Claude Desktop (~/.config/claude/claude_desktop_config.json)
{
  "mcpServers": {
    "chematic": { "command": "chematic-mcp" }
  }
}
```

20 chemistry tools are callable from any MCP-compatible agent (full list in the
[`chematic-mcp` README](crates/chematic-mcp/README.md)):

| Tool | What it does |
|---|---|
| `name_to_smiles` | Resolve "aspirin", "caffeine", … to SMILES via PubChem (the only tool that makes a network call) |
| `calc_properties` | MW, exact mass, Crippen LogP, TPSA, HBD, HBA, rotatable bonds, QED |
| `smarts_match` | Substructure search |
| `pains_check` / `brenk_check` | Flag assay interference or reactive groups |
| `generate_3d` | 3D coordinates via rule-based placement + DREIDING force-field minimization |
| `find_mcs` | Maximum common substructure |
| + 13 more | `ecfp4`, `tanimoto`, `canonical_smiles`, `admet_profile`, `boiled_egg`, `sa_score`, `lipinski_check`, `retrosynthesis`, `smiles_to_moljson`, `moljson_to_smiles`, `representation_router`, `molecule_context_pack`, `parse_smiles` |

**Transport**: stdio (JSON-RPC 2.0 over stdin/stdout) only. Runs as a local
process; there is no hosted Remote MCP endpoint, no authentication, and no
public service SLA — a remote-ready refactor is under consideration but not
implemented.

**Protocol**: speaks both the legacy (`2024-11-05`-style `initialize`
handshake) and the modern MCP `2026-07-28` stateless dialect
(`server/discover`, per-request `_meta`, cacheable `tools/list`,
`structuredContent`) on the same stdio connection — see the
[`chematic-mcp` README](crates/chematic-mcp/README.md#protocol-eras)
for the protocol details.
Remote HTTP, OAuth, the Tasks extension, and MCP Apps remain unsupported.

---

## Why Pure Rust?

### Fast

Rust's zero-cost abstractions and ownership model eliminate overhead at the source.
chematic's ECFP4 fingerprint batch pipeline runs at **~78 µs/mol** on a diverse
molecule corpus — 2–3× faster than RDKit's Python API on the same hardware, via
Rayon parallelism across all CPU cores. No GIL, no interpreter overhead, no FFI
call overhead hidden inside a `_sys` crate.

### Safe

chematic's own ~180,700 lines of Rust (tokei-measured code lines, all 20 crates,
2026-08-21) contain **zero `unsafe` blocks outside one file**: 9 `unsafe {}` blocks
plus 1 `unsafe extern "C"` FFI declaration, all in the optional `native-inchi` layer
(below). No C++ heap corruptions. No segfaults from malformed SMILES input. No
platform-specific build failures from `-sys` crates. The compiler enforces memory
safety at every call site chematic itself wrote.

> The `native-inchi` feature is the single opt-in exception — it vendors the IUPAC InChI
> C library (v1.07.5) for bit-exact standard InChI. All other chematic crates stay
> FFI-free and unsafe-free. This count is chematic's own source only, not its
> dependency tree — the optional `depict` feature (SVG/PDF/EPS rendering) pulls in a
> font/image-rendering stack (resvg/usvg/rustybuzz/tiny-skia/zune-jpeg) that is **not**
> unsafe-free; see the comparison table footnote below for a measured count.

### Anywhere

Pure Rust compiles to `wasm32-unknown-unknown` natively — no Emscripten, no `cmake`,
no `clang`. The npm package `@kent-tokyo/chematic` is **1.10 MB gzip** (2.94 MB raw) —
roughly 2.3× smaller than RDKit.js's `RDKit_minimal.wasm` (6.91 MB raw) on a like-for-like
raw-size basis. One codebase runs on Linux, macOS, Windows, and in every browser.

---

## Benchmarks & Validation

| Metric | Result | Corpus |
|--------|--------|--------|
| ECFP4 throughput | **~78 µs/mol** (2–3× vs RDKit, diverse corpus) | 5,000-mol ChEMBL subset |
| HBA / HBD / aromatic ring count | **100% RDKit agreement** | 4,999-mol ChEMBL subset |
| TPSA | **100% RDKit agreement** within ±0.1 Å² | 4,999-mol ChEMBL subset |
| LogP (Crippen) | **100% RDKit agreement**\* | 4,999-mol ChEMBL subset |
| Stereocenter count | **99.96%** vs legacy†; 98.6% vs new CIP | 4,999-mol ChEMBL subset |
| CIP R/S label agreement | **96.30%** vs modern `rdCIPLabeler`‡; 96.83% vs legacy | 5,000-mol ChEMBL subset |
| WASM bundle | **1.10 MB** gzip (2.94 MB raw) | measured 2026-08-21, commit `ef7dc25` |

\*LogP max Δ = 1.1×10⁻¹³ across 4,999 molecules — within float64 rounding error.  
†Stereocenter count: ~99.96% vs legacy `CalcNumAtomStereoCenters` (a handful of molecules where chematic matches `FindPotentialStereo` and legacy under-counts); ~98.6% vs new-CIP `FindPotentialStereo` (cage/bridgehead molecules where both chematic and legacy correctly return fewer than the new oracle). chematic is calibrated between both extremes. This measures whether an atom is *flagged* as a stereocenter, not whether its R/S label is correct — see the next row.  
‡CIP R/S label agreement measures, for atoms both oracles agree are stereocenters, whether the assigned R/S descriptor matches — a stricter, separate check from stereocenter count agreement above. This row is chematic's *default* `assign_cip()` path. The separate `chematic-cip` engine now reaches 99.38% raw / 99.64% oracle-stable (Milestone 4 gate closed) and is reachable opt-in via `assign_cip_with_mode(mol, CipMode::Accurate)` (Rust), `Mol.cip_stereo(mode="accurate")` (Python), or `cip_assignments_accurate_json` (WASM). No default path changed; this row's 96.30% is unaffected.

All numbers are reproducible with the scripts in this repo.  
Full history → [benchmarks/](benchmarks/) · Methodology → [validation/](validation/)

---

## Comparison with Other Cheminformatics Libraries

| Feature                 | **chematic**                              | RDKit (rdkit-sys)  | OpenBabel FFI  | RDKit.js (WASM)    |
|-------------------------|-------------------------------------------|--------------------|----------------|--------------------|
| **C/C++ dependencies**  | **None (default)**†                       | Extensive C++      | Extensive C++  | C++ via Emscripten |
| **WASM binary size**    | **2.94 MB raw** (1.10 MB gzip)             | N/A (no WASM)      | N/A (no WASM)  | 6.91 MB raw        |
| **Build requirement**   | `cargo build` only                        | cmake + clang      | cmake + clang  | Emscripten SDK     |
| **WASM target support** | **Full (native)**                         | No                 | No             | Yes (Emscripten)   |
| **Python bindings**     | **Yes** (`pip install chematic`, PyO3)    | Yes (rdkit-sys)    | Yes            | No                 |
| **Unsafe Rust**         | **None in own crates**‡                   | Extensive          | Extensive      | N/A                |

See the [format capability matrix](docs/format-capabilities.md) and the
[RDKit migration guide](docs/rdkit-migration.md) for detailed support
differences. The table above is intentionally limited to deployment-level
differences; detailed feature claims belong in those maintained pages.

---

## JavaScript / TypeScript (WebAssembly)

**1.10 MB gzip — roughly 2.3× smaller than RDKit.js's raw WASM.** No Emscripten, no cmake. Drop-in for browser or Node.js.

```sh
npm install @kent-tokyo/chematic
```

```js
import init, { parse_smiles, get_descriptors_json, tanimoto_ecfp4,
               generate_3d_minimized_pdb, enumerate_stereo_isomers_json,
               maxmin_picks_ecfp4_json } from '@kent-tokyo/chematic';

await init();

const mol = parse_smiles('CC(=O)Oc1ccccc1C(=O)O'); // aspirin
console.log(mol.molecular_weight(), mol.qed(), mol.lipinski_passes());

// All descriptors as a JSON object
const desc = JSON.parse(get_descriptors_json(mol));

// Fingerprint similarity
const caffeine = parse_smiles('Cn1cnc2c1c(=O)n(c(=O)n2C)C');
console.log(tanimoto_ecfp4(mol, caffeine));  // 0.26

// 3D coordinates, stereoisomers, diversity picking
const pdb = generate_3d_minimized_pdb(mol);
const isomers = JSON.parse(enumerate_stereo_isomers_json(parse_smiles('C(F)(Cl)Br')));
const picks = JSON.parse(maxmin_picks_ecfp4_json('["CC","c1ccccc1","CCO","CCCC"]', 2));
```

218+ exported functions (plus `MolHandle`/`DepictOptions` class methods, measured 2026-08-21) cover descriptors, fingerprints, 3D geometry, reactions (incl. `retro_disconnect_json` — single-step retrosynthetic disconnection), diversity picking, and SDF round-trips.
See the [full WASM API reference](https://kent-tokyo.github.io/chematic/) for all exports.
---

## Crate Reference

| Crate                 | Description                                                                                              | Tests |
|-----------------------|----------------------------------------------------------------------------------------------------------|-------|
| `chematic-core`       | Atom, Bond, Molecule, Element, kekulization (no deps); mutable `add/remove_atom/bond`, `fragments()`, `is_connected()`, `formula_with_isotopes`, `validate_valence`; `StereoGroup`/`StereoGroupKind` | 132    |
| `chematic-smiles`     | OpenSMILES parser, writer, canonical SMILES; **stereo parity correction** (pre-solves RDKit #8775 — `@`/`@@` auto-flipped on odd permutations); **allene cumulated double bond stereo** (`C=C=C` `@`/`@@`, round-trip stable) | 202    |
| `chematic-perception` | SSSR, Hückel aromaticity + antiaromaticity (4n+2 rule), `apply_aromaticity`, `aromatize`/`kekulize_inplace`, `assign_stereo_from_2d`, `assign_ez_from_2d`, `cip_ez_descriptor`; **zero-order/dative bonds excluded from ring perception** | 194    |
| `chematic-mol`        | MOL/SDF V2000+V3000 (R/W with 2D coords, +partial charge writing), CML (R/W), CDXML (R); `SdfRecord` with coords+props; MDL RXN R/W; V3000 stereo-group COLLECTION R/W; **AutoDock PDBQT** (parse + write); **ChemicalJSON** (`parse_cjson`/`write_cjson`, Avogadro/MolSSI format); **2D wedge/hash tetrahedral parity + E/Z double-bond direction now perceived automatically on read** (`read_mol_with_diagnostics`/`read_mol_v3000_with_diagnostics`, typed opt-in diagnostics); **PDBx/mmCIF** (R/W, chain/altloc/insertion-code/model/occupancy/B-factor preserved — Open Babel's own mmCIF support is read-only); **PQR** (R/W); **QCSchema JSON** (`Molecule`/`AtomicInput`/`AtomicResult`, MolSSI schema, Bohr↔Å conversion); **ORCA** (input R/W with lossless unknown-block preservation, output R — final geometry/trajectory/energy/frequencies/termination/convergence as independent typed fields); new shared `VolumetricGrid` type + **Gaussian Cube** (R/W, streaming-*input* `CubeFileReader` for large grids — the parsed voxel array is still fully in-memory, non-orthogonal axes, explicit Bohr/Ångström unit tag) + **OpenDX/APBS scalar field** (R/W) — single-dataset only, multi-dataset Cube typed-rejected rather than silently truncated | 476   |
| `chematic-depict`     | 2D SVG (CPK colors, highlighting, grid), DepictData, `detect_crossings`, `render_svg_with_metadata`, reaction SVG; **PDF output** (`depict_pdf`/`depict_pdf_opts` via svg2pdf); **EPS output** (`depict_eps`/`depict_eps_opts`, pure Rust); `tiny_skia` PNG is optional `png` feature (default on, disabled for WASM) | 75    |
| `chematic-chem`       | 190+ descriptor values (71 functions), tautomers, scaffold, BRICS, QED, standardize, CIP; **pKa prediction** (23 SMARTS rules); **ADMET profile** (BBB/Caco-2/hERG/CYP3A4); **HBA 100% RDKit agreement** (4 999 / 4 999 mol benchmark); **TPSA 100% ±0.1 Å² / LogP 100%\* / HBD 100% / stereocenter count 99.96% (legacy) / 98.6% (new CIP)** vs RDKit (4,999-mol ChEMBL); **CIP R/S label agreement 96.30% (default), 99.64% oracle-stable via opt-in `CipMode::Accurate`** (5,000-mol ChEMBL); **topological descriptors** (`petitjean_index`, `graph_diameter`, `graph_radius`, `graph_eccentricities`, `eccentric_connectivity_index`, `hosoya_index`, `moran_autocorr`, `geary_autocorr`); **`schultz_mti`, `gutman_mti`, `vabc` (Bondi radii vdW volume), `gravitational_index`**; `clean_stereo_groups()` in standardize | 724   |
| `chematic-fp`         | ECFP2/4/6, FCFP4/6, MACCS, TopoPF, AtomPair, Torsion, Layered, Pattern, Pharmacophore, Reaction, **MAP4** (Minervini 2020, not in RDKit) — Tanimoto/Dice; bulk similarity | 266    |
| `chematic-ff`         | **Experimental MMFF94 implementation** with all 7 term families (Halgren 1996): Bond/Angle/Torsion/vdW/Elec + **OOP** + **Stretch-Bend**; steepest-descent + L-BFGS optimizer, torsion scan, energy breakdown; incomplete typing/parameter coverage remains observable as failure; DREIDING typing; **UFF** (metals/organometallics: Zn, Fe, Cu, …) | 198    |
| `chematic-smarts`     | SMARTS, VF2, MCS with chirality matching; **SmartsCache** (LRU compilation cache, 5–20×); **named_pattern()** library (20 functional group patterns); **atom map `:N` in SMARTS** (`[O;D1;H0:3]` — stored as metadata, not a match criterion); **`[kN]` ring-size primitive**; **VF2 early-exit** when query > target atom count; **`find_matches_with_rings`** — share SSSR across multi-pattern batches | 169   |
| `chematic-3d`         | 3D coordinate generation, distance geometry constraints, ETKDG KB (40 torsion patterns, adaptive noise), force-field minimization, shape descriptors, ConformerEnsemble with RMSD pruning, PDB/XYZ; **GETAWAY HATS-matrix** (full 19-dim implementation); **`whim_getaway_combined()`** now 29-dim | 540    |
| `chematic-rxn`        | Reaction SMILES/SMIRKS, `run_reactants`/`run_reactants_strict`; **`retro_disconnect()`** — 60 retro-SMIRKS templates (AmideBond/Ester/Ether/CNBond/CCBond/CSBond) + SA Score ranking; **parity-aware `@`/`@@` SMIRKS stereo filtering**; **E/Z double-bond stereo filtering** in `run_reactants` (`ez_stereo_outward`, `smirks_ez_stereo_ok`) | 180    |
| `chematic-inchi`      | InChI/InChIKey: pure-Rust approximation (WASM) **+ IUPAC-standard** via `native-inchi` feature (vendored C lib 1.07.5, bit-exact); **parse_inchi** reader; **verified canonical-SMILES dedup** (`dedup::{group_candidates, deduplicate_verified}`, fail-closed on legacy-CIP-unresolved specified tetrahedral stereo); **accurate-CIP dedup preflight** (issue #161) recovering verified-comparison capability on legacy-CIP-unresolved stereocentres; **indexed graph relation API** (`compare_indexed_graph_relation`, orthogonal `GraphStrictness`/`AtomMapPolicy` axes) | 108 (+16*)    |
| `chematic-cip`        | Opt-in accurate CIP engine (`assign_cip_accurate_experimental`, hierarchical digraph, Rules 1a/1b/2/4b/5, RDKit-compatible MANCUDE fractional atomic numbers) — the default `assign_cip()`/`CipMode::LegacyFast` is unchanged | —     |
| `chematic-wasm`       | **218+ WASM exports** (plus class methods; measured 2026-08-21) — npm: `@kent-tokyo/chematic` (published in lockstep with crates.io/PyPI); pKa/ADMET/BBB/Caco-2/hERG/CYP3A4; `smiles_to_pdbqt`, `minimize_uff_json`, **`retro_disconnect_json`** (issue #91) | 276   |
| `chematic-iupac`      | Local IUPAC name generation — **25+ compound classes**: alkanes, cycloalkanes, alkenes/alkynes, alcohols, amines, halides, aldehydes, ketones, acids, esters, amides, **piperidine, morpholine, piperazine, naphthalene, sulfides** | 56    |
| `chematic-mcp`        | **MCP (Model Context Protocol) server** — AI agent integration; **20 tools**: parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, name_to_smiles, retrosynthesis, smiles_to_moljson, moljson_to_smiles, representation_router, **molecule_context_pack**; dual-era protocol (legacy `2024-11-05` + modern `2026-07-28` stateless dialect), `structuredContent`/`outputSchema` on all 20 tools | 82    |
| `chematic-py`         | PyO3 Python bindings (`pip install chematic`); 300+ API endpoints: `from_smiles()`, `Mol.descriptors()`, `Mol.minimize_dreiding()`, `from_cxsmiles()`, `from_rxn_file()`/`to_rxn_file()`, `parse_sdf_with_coords()`, `Mol.ring_families()`, `tanimoto_matrix()`, `iter_sdf()`, `SimilarityIndex`; **`mol.to_pdf()`/`mol.to_eps()`** (depict); **`from_cjson()`/`mol.to_cjson()`** (ChemicalJSON); **`mol.schultz_mti`, `mol.gutman_mti`, `mol.vabc`, `mol.gravitational_index`**; **`bulk.substructure_match(smarts, mols)`** (parallel VF2 on pre-parsed Mol objects); **`mol.describe()`** (LLM/MCP-ready natural-language summary); **`mol.diff(other)`** (element + descriptor diff); **`PeriodicStructure.from_cif()`/`.from_poscar()`, `Lattice`, `Site`** (periodic/crystal structures — `chematic-crystal`'s first host-language binding); **`from_cif(text, expand_symmetry=True)`** expands a CIF's own literal symmetry-operation list into a full unit cell by default (`expand_symmetry=False` for the asymmetric unit only — no space-group database, no name/number-to-operations generation); Sprint 18–27 coverage | 300+  |
| `chematic-ewald`      | PME Ewald summation, B-spline interpolation (cubic, phase-corrected)                                     | 16    |
| `chematic`            | Umbrella crate with feature flags (all sub-crates, incl. `iupac`, `inchi`)                              | 1     |

```
cargo test --workspace --lib --quiet                                          # 3,912 tests, all passing (2026-08-21)
cargo test -p chematic-inchi --features native-inchi --test standard_inchi  # +16 IUPAC-exact InChI tests
```

---

## Recent Development

**v0.23.0** (2026-08-30): **MCS accuracy fix (behavior change), two more RDKit fingerprint ports, full MCS bindings**
- `chematic-smarts`: **behavior change** — `find_mcs`'s default `AtomCompare::Elements` no longer requires matching aromaticity, matching RDKit's identically-named `rdFMCS.AtomCompare.CompareElements` exactly (confirmed via live oracle: RDKit never encodes aromaticity as a per-atom constraint, only via bond-type queries). Agreement vs. a live RDKit oracle rose from 74.6%/68.2%/70.4% to 88.4%/88.5%/97.0% across three established corpora. There is currently no `AtomCompare` mode that restores the old strict element+aromaticity match
- `chematic-py`/`chematic-wasm`/`chematic-mcp`/`chematic-chem`: fixed `find_mcs` result reconstruction silently losing heteroatoms and/or aromaticity across all 4 binding surfaces (`QueryMolecule` → concrete `Molecule` conversion never unwrapped the compound atom query correctly) — surfaced while measuring the fix above
- `chematic-fp`/`chematic-py`: `rdkit_rdk_fp`/`rdkit_layered_fp` — RDKit-compatible `Chem.RDKFingerprint`/`Chem.LayeredFingerprint` ports, completing a 6-fingerprint parity series (100%/100%/99.44% and 100%/100%/99.46% bit-exact across 3 corpora vs. a live RDKit oracle)
- `chematic-py`/`chematic-wasm`: full `McsConfig`/`McsOutcome` exposed to `find_mcs` bindings (`match_charge`/`match_isotope`/`atom_compare`/`bond_compare`/`timeout_ms`/etc., previously Rust-only)
- Full details in `CHANGELOG.md`'s `[0.23.0]` section

Current development is tracked in [`CHANGELOG.md`](CHANGELOG.md). The latest
`v0.31.0` work adds Parent identity bindings for WASM, following the Python
bindings in `v0.30.0`; both expose bounded, status-aware operations.

The older release notes below are retained as a short historical summary.

**v0.22.0** (2026-08-29): **New WASM ensemble binding, canonicalization-hang fix, 3-membered-ring embedding fix**
- `chematic-wasm`: new `embed_ensemble_v2_json` binding for `chematic_3d::embed_ensemble_v2`, mirroring the Python binding (`Mol.conformer_ensemble_v2()`) via the existing `pipeline_v2.rs` conventions (camelCase JSON keys, `schemaVersion: 1` envelope) — purely additive
- `chematic-smiles`: `canonical_smiles`/`canonical_atom_order` could hang on molecules with several simultaneously-unresolved symmetric regions (issue #421) — the automorphism backtracking search had no internal step bound; fixed with an always-on step ceiling that safely falls back to "not proven automorphic" rather than searching unbounded
- `chematic-3d`: 3-membered rings (cyclopropane/epoxide/aziridine/thiirane) failed closed at the distance-geometry embedding stage — a generic angle bound was overwriting the correct, tighter bond-length bound for a ring-closing pair that's simultaneously "1-3" and directly bonded; fixed by skipping the angle bound for any already-bonded neighbor pair; strict-MMFF94 3D corpus 252/265 → 263/265
- Full details in `CHANGELOG.md`'s `[0.22.0]` section

**v0.21.0** (2026-08-27): **`McsConfig` charge/isotope matching + typed timeout outcome, four correctness fixes**
- `chematic-smarts`: `McsConfig` gains `match_charge`/`match_isotope` fields (mirroring the existing `match_chiral_tag`, default `false`) and a new `McsOutcome` enum (`Exhaustive`/`TimedOut`) via `find_mcs_with_config_checked`, reporting whether a timeout cut the search short rather than silently returning a possibly-non-optimal result — purely additive, `find_mcs`/`find_mcs_with_config` unchanged
- `chematic-smarts`: `find_mcs`'s branch-and-bound search was incomplete — `grow()` only ever tried the first frontier atom with no way to exclude it and try another, silently missing a strictly larger common substructure in some cases (minimal repro: `OC(N)N` vs `NC(N)` returned 2 atoms instead of the true 3); fixed via standard include/exclude branch-and-bound
- `chematic-chem`: `disconnect_metals` left a dative-bond-derived formal charge unneutralized after severing the metal bond (issue #403) — 34/4999 molecules in RDKit's own bundled NCI Diversity Set holdout were non-idempotent; fixed by recomputing the affected atom's H count via valence inference immediately after disconnection; NCI holdout 34 → 0 failures, new 11-fixture metal-complex holdout added
- `chematic-chem`: `normalize_zwitterion` invented a proton on the negative atom of a permanently charge-separated group with no transferable proton anywhere (e.g. a diazo-N,N'-dioxide), silently changing the molecule's formula (issue #407) — fixed by gating the transfer on both atoms actually having a proton to move; dev-corpus residual 4 → 1 (the remaining one is unrelated, see issue #402/#415)
- `chematic-chem`: `canonical_tautomer` could produce a chemically invalid, over-valent nitrogen in a fused/bridged ring system (issue #415) — both aromatic H-shift mechanisms now validate their output via kekulization before accepting it
- Full details in `CHANGELOG.md`'s `[0.21.0]` section

**v0.20.1** (2026-08-26): **Three canonical-SMILES/standardization correctness fixes (patch release, no breaking changes)**
- `chematic-smiles`: coupled E/Z canonicalization could silently change geometry on re-canonicalization (issue #390) — two independent defects in the canonical writer's E/Z marker machinery, both required to reproduce the filed witness; fixed together, verified against a real 290-compound corpus (290/290 idempotent, 290/290 matching independent RDKit InChIKeys, up from 289/290)
- `chematic-chem`: `standardize()` silently dropped stereo tables on several rebuild paths (issue #399) — 8 functions in `standardize.rs` rebuilt the molecule via a bare `MoleculeBuilder` without carrying `stereo_neighbor_order`/`bond_directions`/`stereo_groups` forward, flipping `@`/`@@` depending on ring-open/close role; dev-corpus standardize-path idempotency 615/519 → 68/60 (the exact pre-#392 baseline), NCI holdout (4,999 unused real molecules, run once) 0 stereo-related failures
- `chematic-smiles`: canonical-writer ring-closure bond markers ignored the closure partner's aromaticity (issue #395) — a genuinely non-aromatic ring-closure "fusion" bond between two individually-aromatic atoms (e.g. `c1-2`) silently became aromatic on re-parse; dev-corpus bare-parse idempotency 73/57 → 0/0, a complete fix, independent RDKit InChI oracle 0 mismatches across all 10,000 corpus lines
- Combined, the two standardize-path fixes bring that corpus to 0/4 residual — the 4 remaining failures are traced to two newly-filed, not-yet-fixed issues (#407, #402-class), not folded into this release
- Full details in `CHANGELOG.md`'s `[0.20.1]` section

**v0.20.0** (2026-08-25): **Stereo-safe 3D generation, a connectivity-ordered coordinate engine, and identity-correctness fixes for `remove_hydrogens`**
- `chematic-3d`/`chematic-py`/`chematic-wasm`: `PipelineV2Config::stereo_safe(force_field_policy)` — a single-call configuration that resolves a real gap for ring-fused declared stereocenters (testosterone, cholesterol, and similar), where `repair_tetrahedral_center` previously had no coordinate to reflect an implicit H against. Measured on a 29-molecule × 5-seed corpus: 144/145 (99.3%) correct_and_ok, 0 silently wrong, testosterone/cholesterol both 5/5 seeds on every declared stereocenter
- `chematic-3d`: `generate_coords_connectivity_ordered`, a new public alternative 3D placement engine (issues #256/#255) — places rings and chain atoms in true connectivity order rather than "all rings, then all chains." Measured: raw-geometry soundness 10/33 → 33/33 on a differential corpus with zero regressions, post-UFF bond-violation rate reaching 0.0000 (better than the legacy engine's own baseline). `generate_coords` itself is completely unchanged — ships as an available alternative, not a default-behavior switch; no existing caller is routed to it
- `chematic-3d`: `rescue_with_distance_geometry_v2` (the UFF-catastrophic-blowup rescue bridge) now enforces declared chirality on retry (issue #210, partial fix) — zero regressions on its 58-molecule corpus, one of 5 named residual molecules newly succeeds; 4 residuals remain unfixed via this specific bridge (though resolved via `stereo_safe` above)
- `chematic-chem`: `remove_hydrogens` no longer destroys isotope-labeled hydrogen (`[2H]`, `[3H]`) or silently drops declared stereo/E-Z information on every call — both were real, shipped correctness bugs found via a downstream consumer's 9.47M-compound real-world corpus scan, unaffected molecules confirmed unchanged, 289/290 of the originating investigation's identity mismatches resolved
- `chematic-py`: `Mol.conformer_ensemble_v2(config)` exposes `embed_ensemble_v2` (deterministic multi-conformer generation, energy ranking scoped within force field, full per-attempt provenance) — added alongside the existing `conformer_ensemble()`, not replacing it. New best-of-10 benchmark arm confirms it works robustly at scale (~250/265 molecules, median RMSD 2.147 Å / TFD 0.344 vs. RDKit) — RDKit conformer-*selection* parity is a separate, unestablished claim
- Known limitations: `generate_coords` not yet routed to the new engine; issue #210's 4 residuals remain open via that specific bridge; issue #390 (a single-molecule E/Z correctness residual, unrelated) remains open; no fresh full-corpus re-measurement was run solely for this release
- Full details in `CHANGELOG.md`'s `[0.20.0]` section

**v0.19.0** (2026-08-23): **Round 2C aromatic lactam/lactim tautomer fix, plus a benchmark/validation refresh**
- `chematic-chem`: Tautomer & Parent Identity round 2C (ROADMAP.md Phase 2) — `canonical_tautomer`/`tautomer_parent` now cover the aromatic lactam/lactim class for 2-pyridone, 4-pyridone, uracil, cytosine, guanine, methylpyrimidinone, and the primary/N9-methyl hypoxanthine cases. Aromatic exocyclic tautomer edges are traversed in both directions by bounded enumeration, while canonical selection retains the lactam preference. The former tp2-39 and tp2-holdout-06 residuals were corrected after RDKit InChIKey review exposed positional-isomer fixture errors; the unrelated nitroso/oxime defect is fixed. Python/WASM Parent-API bindings are implemented.
- Benchmark/validation refresh: every number in `docs/benchmark.md`/`docs/validation.md` was pinned to chematic v0.4.29/RDKit 2026.03.3 (~14 releases stale) — re-measured fresh against RDKit 2026.03.4. The 4,999-mol accuracy corpus is now committed (`scripts/chembl_accuracy_corpus_4999.smi`, previously an uncommitted personal path); molecular weight has a real corpus-wide check for the first time (99.82%, not the previously-unmeasured "175-mol"/100% placeholder); CIP R/S/E/Z label agreement re-measured at 99.74–99.78% (up from a stale 96.30–96.83%); WASM bundle size rebuilt clean; the ECFP4 "diverse corpus" figure now has a reproducible source (`benchmark_vs_rdkit.py --corpus`, previously none existed); 3D conformer generation's "Good (ETKDG rules)" framing corrected to "Experimental," matching the migration guide's own honest characterization
- Full details in `CHANGELOG.md`'s `[0.19.0]` section

**v0.18.0** (2026-08-20): **Python/WASM bindings for the 7 v0.17.0 formats, plus an MMFF94 atom-typing fix and a binding-quality/cross-language-consistency pass**
- `chematic-ff`: fixed the aryl-isothiocyanate cumulated-double-bond CSP carbon mistyping from issue #337 (`getTotalDegree() == 2` replacing a `triple_bonds > 0`-only check — a strict superset, RDKit's real rule); the other 6/8 molecules behind that issue were re-diagnosed as a genuine RDKit Kekulization/MMFF-aromaticity-perception artifact (confirmed via direct negative-control fragments) rather than a locally-fixable typing rule, and left as an honestly-disclosed residual
- `chematic-py`: Python bindings for all 7 v0.17.0 formats (mmCIF, PQR, ORCA, QCSchema, Gaussian Cube, OpenDX, LAMMPS data/dump) — previously Rust-only; `VolumetricGrid`/`LammpsDumpFrame` pyclasses with numpy-array properties, `to_opendx`/`to_opendx_lossy` fail-closed split preserved faithfully; a `py.typed` marker verified to actually ship in the built wheel (`mypy --strict` passes against a fresh-venv wheel install, not just the source tree)
- `chematic-wasm`: WASM (wasm-bindgen) bindings for the same 7 formats, plus 5 additive `js_sys::Float64Array`/`Uint32Array`-returning functions alongside the existing JSON-string API (large numeric grid/row data without a full JSON round trip) — this crate's first typed-array precedent
- Cross-language parity: the same 4 small fixtures (Cube, OpenDX, mmCIF, LAMMPS triclinic dump) independently verified to produce identical results from Rust, Python, and WASM entry points — no discrepancy found
- Full details in `CHANGELOG.md`'s `[0.18.0]` section

**v0.17.0** (2026-08-17): **Format/Python/materials-interop breadth, plus two MMFF94 charge/bond-order accuracy fixes**
- `chematic-mol`: square-planar (`@SP1`/`@SP2`/`@SP3`-equivalent) stereo read/write for MOL/SDF via 3D-coordinate-derived reperception; PDBx/mmCIF, PQR, QCSchema JSON, and ORCA input/output; CIF explicit symmetry-operation expansion into a full unit cell (Rust + Python); a shared `VolumetricGrid` type plus Gaussian Cube and OpenDX (APBS-scoped) I/O; LAMMPS data-file (`read_data` format) and dump/trajectory-file I/O as standalone document types, not integrated with `Molecule`
- `chematic-py`: new Python bindings for `chematic-crystal`'s `Lattice`/`PeriodicStructure`/`Site` — found and fixed a real pre-existing bug along the way (`to_cif()` silently re-declared an unexpanded-symmetry CIF as false P1)
- `chematic-ff`: MMFF94 bond-order-classification fix (`assign_mmff94_numeric_types_with_view` — production energy/gradient entry points and the coverage gate now agree on classification; `torsions_missing` 257→0 on the 265-molecule Wave 1 corpus) and an MMFF94 BCI partial-charge fix (own wrong `bond_type_for`, then root-caused a second bug to RDKit's atom-type-*derived* formal charge never being computed at all) with one post-minimization stereo-repair addition the first fix surfaced. Production `pipeline_v2_mmff94_strict`: 240/265 → 241/265
- Release-hygiene: corrected an MSRV declaration that didn't match reality (`rust-version` raised to 1.88, now continuously verified by a dedicated CI job)
- Full details in `CHANGELOG.md`'s `[0.17.0]` section

**v0.16.0** (2026-08-15): **Periodic-structure interoperability (CIF/POSCAR/FPS) and generalized stereochemistry foundation**
- `chematic-mol`: new optional `crystal` feature bridges the existing CIF reader/writer to `chematic_crystal::PeriodicStructure` (`parse_cif_periodic_structure`/`write_cif_periodic_structure`) — cell parameters to `Lattice`, `_atom_site_occupancy` to `Occupancy`, disorder-sharing atom-site rows merged into one `PeriodicSite`'s multi-species list. New `CifSymmetryStatus` enum distinguishes genuinely-P1 CIFs from CIFs that declared symmetry this parser doesn't expand, rather than silently treating the latter as P1. `chematic-crystal` itself remains independent of `chematic-mol`/`Molecule` (dependency direction is one-way: `chematic-mol` → `chematic-crystal`, optional)
- `chematic-crystal`: native POSCAR/CONTCAR (VASP structure format) read/write — `parse_poscar`/`parse_contcar`/`write_poscar`, VASP 5 only, both scale-factor conventions, Direct/Cartesian coordinates, selective dynamics, ion velocities, and CONTCAR's predictor-corrector MD-restart section preserved verbatim (VASP's own docs don't specify its numeric layout)
- `chematic-fp`: new `fps` module — streaming read/write for the FPS ("Fingerprint file format") text-based interchange format popularized by chemfp/OpenBabel, hex bit-ordering verified against the chemfp spec, reuses `BitVec2048`/`BitVecN` as the sole bit-vector representation
- `chematic-core`: new `stereo_geometry` module — stereo configuration modeled as a coordination geometry (`Tetrahedral`/`SquarePlanar`, `#[non_exhaustive]` for future TBP/octahedral) plus the equivalence class of ligand-slot permutations under that geometry's proper rotation group (A4, order 12, for tetrahedral; the order-8 S4-stabilizer of a trans-pair partition, not the naive order-4 in-plane-only group, for square-planar). Replaces two independent hand-written stereo-remapping algorithms in `chematic-smiles`; `@`/`@@`/`@SP1`/`@SP2`/`@SP3` semantics fully preserved (88-fixture byte-identical canonical-SMILES regression). Fixed a real bug found along the way: a square-planar-tagged atom in `chematic-3d` could be silently coerced into a tetrahedral chiral-volume check decided by floating-point noise; also fixed a transient allene-end-carbon parity regression surfaced during development, pinned by an exact golden-value test.
- Release-grade re-measurement of the `pipeline_v2` vs RDKit 2026.03.4 benchmark (superseding stale 2026-08-06 numbers): `mmff94_strict` 149/265 → 239/265. New finding: torsion parameter coverage, not bond/angle, is now the dominant remaining MMFF94 gap (71% of `complete_bonded_term_gated` failures cite missing torsion parameters, 0% OOP, 0% bonds) — direct evidence for the project's next MMFF94 roadmap item
- Full details in `CHANGELOG.md`'s `[0.16.0]` section

**v0.15.0** (2026-08-14): **`chematic-crystal` — periodic (crystal) structure foundation crate, MMFF94 Bond/Angle empirical-rule fallback (issue #227)**
- New crate `chematic-crystal`: periodic (crystal) structure representation and geometry — `Lattice` (triclinic-capable, validated matrix/inverse/reciprocal vectors), `FractionalCoord`/`CartesianCoord`, `PeriodicSite`/`SiteSpecies`/`Occupancy` (multi-species disorder-ready), and `PeriodicStructure` with exact (not `round()`-approximate) periodic minimum-image distance — equidistant periodic images resolve deterministically to the lexicographically smallest image — cutoff neighbor enumeration, and diagonal supercells. Deliberately **not** an extension of `chematic_core::Molecule` (a bond graph). Optional `serde` feature; optional `crystal` feature on the `chematic` facade, included in `full` (does not change `default`, which stays empty). No symmetry, no CIF parser changes, no Python/WASM/MCP bindings yet
- `chematic-ff`: ported Halgren's MMFF.V eq. 18-20 empirical Bond-stretch/Angle-bend rule (`mmff94_bond_energy_resolved`/`mmff94_angle_energy_resolved`, new additive functions — the existing `mmff94_bond_energy`/`mmff94_angle_energy` keep their original signatures), tried strictly after the existing exact-table/`eqLevel`-ladder lookup so it never overrides a real table hit. Along the way, found and fixed a real data gap: 97 rows present in RDKit's real Angle table (generic central-atom-type-only `theta0` defaults) were missing from chematic's port. One triple is deliberately left unresolved (fails closed) rather than guessed — the outer atom type has no equivalence-class entry and RDKit's own real code dereferences that unchecked (undefined behavior), so its live-oracle answer couldn't be attributed to any well-defined mechanism. Also fixed 5 pre-existing MMFF94 atom-typing gaps and ported RDKit's `eqLevel` atom-type-equivalence ladder for Angle lookup. Net effect on the 265-molecule Wave 1 corpus (production minimization path), reported as two separately-verified numbers (both via a full per-molecule join, zero regressions either way): full v0.14.1→v0.15.0 change 158/265 → 248/265 (107 → 17 failing); the empirical-rule work specifically (isolated from the atom-typing/`eqLevel` prerequisites merged earlier in this same release) 178/265 → 248/265 (87 → 17 failing). The 3 molecules still `MinimizationFailed` in the final state were already non-`Ok` in v0.14.1 — a pre-existing geometry issue newly exposed once real parameters became available, not a regression
- Full details in `CHANGELOG.md`'s `[0.15.0]` section

**v0.14.1** (2026-08-12): **Anticancer platinum coordination-chemistry compatibility fixes, Extended XYZ (extxyz) read/write**
- `chematic-core`: `valence_inferred_hcount` treated a `BondOrder::Dative` bond's donor side exactly like a covalent single bond when computing implicit hydrogen count — an un-bracketed dative donor like `N->[Pt]Cl` computed as `NH2` instead of the chemically correct `NH3`. Donor-side dative bonds now contribute 0 to the valence sum; found via a platinum coordination-chemistry benchmark but general (verified against Fe/Co/Pd/Ru acceptors too), not platinum-specific
- `chematic-mol`: MDL bond type 9 (dative/coordinate — RDKit's own V3000 convention for `Bond::BondType.DATIVE`) silently mapped to `BondOrder::Single` in both V2000 and V3000 readers, quietly discarding coordination-bond semantics on read. Both readers now map code 9 to `BondOrder::Dative`; V3000's writer now emits code 9 instead of collapsing to plain single
- `chematic-chem`: `avg_mass`/`mono_mass` covered only ~24 light main-group elements and silently fell back to `atomic_number as f64` for every other element — every transition metal, lanthanide, actinide, and heavy post-transition element (platinum: atomic number 78, real mass ~195 Da, previously returned "78.0 Da") got a wildly wrong mass with no error. Extended to all 118 `Element` values, sourced from RDKit's periodic table data, with the ~24 previously-covered values kept as-is where they differ (selenium: this project's value is the current IUPAC standard, RDKit ships the superseded pre-2013 value)
- `chematic-mol`: new Extended XYZ (extxyz) format support — `parse_extxyz`/`write_extxyz`, `ExtxyzReader`/`ExtxyzWriter`, `parse_extxyz_all`, built as an extension of the existing multi-frame `XyzFrame` type (ASE's `Lattice=` cell matrix, typed per-atom `Properties=` columns, arbitrary `key=value` frame metadata); a plain XYZ file round-trips through the extxyz reader/writer unchanged. Python: `from_extxyz`/`from_extxyz_all`/`to_extxyz`. WASM: `mol_from_extxyz`/`extxyz_frame_json`/`to_extxyz_json`. **Breaking (Rust API only)**: `XyzFrame` gained three public fields, `XyzError` gained seven variants, `write_extxyz` now returns `Result<String, XyzError>` — a real break to the `chematic-mol` v0.14.0 Rust API already published to crates.io, not merely an unreleased-API change
- Platinum coordination-chemistry stereochemistry (square-planar cis/trans identity, e.g. cisplatin vs. transplatin) remains unrepresented.
- Full details in `CHANGELOG.md`'s `[0.14.1]` section

**v0.14.0** (2026-08-11): **Stereo-aware distance geometry — declared E/Z enforced as a bound-matrix constraint, `enforce_chirality` composable with post-minimization stereo verification, Python/WASM exposure**
- `chematic-3d`: root-caused and fixed the issue #285 release-gate waiver from v0.13.0 — `apply_vdw_bounds`'s generic non-bonded Van der Waals lower bound was being applied to a declared-E/Z alkene's own 1-4 substituent pair regardless of declared stereochemistry, structurally excluding the correct cis geometry from ever being sampled. New `apply_declared_ez_bounds` (`enforce_chirality`-only) intersects an analytic same-side/opposite-side 1-4 distance bound into the bond matrix *before* the generic Van der Waals floor applies, so the correct geometry is reachable by construction, not by post-hoc repair/retry/reflection. Unlike tetrahedral chirality (which a pairwise distance matrix can never encode — a molecule and its mirror image have identical pairwise distances), declared E/Z is genuinely distance-representable, since cis/trans are two different scalar separations, not mirror images. Measured on the 265-molecule corpus's declared-E/Z subset (39 molecules): stereo-satisfied 22 → 42, violated 23 → 3, pipeline success/soundness unchanged
- `chematic-3d`: `embed_pipeline_v2`'s config-validation gate previously rejected `enforce_chirality: true` for any `stereo_policy` other than `Ignore`. Corpus measurement found this wrong — `enforce_chirality` protects embedding-time correctness only, and force-field minimization (which has no notion of declared stereo) can walk a correctly-embedded E/Z bond back across its boundary afterward (found on 2 real molecules, confirmed by re-running with no force field). `enforce_chirality: true` is now also allowed with `StereoPolicy::VerifyOnly`, whose existing post-minimization gate catches exactly this failure mode as a typed error instead of a silent wrong-stereo `success`
- `chematic-py`, `chematic-wasm`: `enforce_chirality` (default `false`) is now a real, settable parameter/field on `PipelineV2Config`/`PipelineV2Config.safe()` (Python) and the `enforceChirality` JSON field (WASM) — neither binding had ever threaded the field through before, so the fix above was previously unreachable from Python or WASM callers
- `chematic-rxn`: fixed `suzuki_biaryl`'s retro-template (issue #294) — `[c:1][c:2]` never matched a real biaryl bond, only intra-ring aromatic bonds, since two adjacent aromatic atoms with no explicit bond token default to aromatic in this crate's SMILES convention. Fixed to `[c:1]-[c:2]`. Found along the way: 14 of 59 `DEFAULT_TEMPLATES` entries silently never parse at all — filed as issue #296, not fixed here
- Opt-in only — `enforce_chirality: false` remains the default everywhere; the default conformer path (`generate_coords_etkdg`/`Mol.conformer_ensemble()`) is untouched
- Full details in `CHANGELOG.md`'s `[0.14.0]` section

**v0.13.0** (2026-08-10): **MMFF94 stretch-bend + torsion parameter-selection parity (both breaking), per-atom stereocenter API, E/Z completeness, macrocycle detection, notation-invariant atropisomer detection/assignment, XYZ I/O**
- `chematic-ff`: `mmff94_stbn`/`mmff94_stbn_type_only` now key the `MMFF94_STBN` table lookup on RDKit's real, finer-grained "stretch-bend type" (`getMMFFStretchBendType`, 0-11) instead of the coarser angle type (0-8) previously used as a stand-in (issue #227) — 220 of 427 stretch-bend routing candidates on the 265-molecule Wave 1 corpus move from RDKit's generic Dfsb periodic-row default to the correct, specific parameter; `angle_type_for`'s ring-offset formula also corrected to match RDKit's real `getMMFFAngleType`. **Breaking**: `mmff94_stbn`/`mmff94_stbn_type_only`'s leading `u8` parameter is now `stretch_bend_type`, not `angle_type` — same shape, different required value; use the new `pub stretch_bend_type_for` to compute it
- `chematic-ff`: `torsion_type_for` now classifies from the real j-k bond's MMFF bond type (reusing `bond_type_for`) plus RDKit's real local-bond-adjacency ring-4/5 override, instead of atom-type-membership alone — corrects 76.9% of the 1,107 previously-missing torsion instances via classification alone, and, corpus-wide, corrects 1,792 of 13,530 torsion instances that resolved to a *silently wrong* parameter value before (not just missing coverage) — 99.1% of the corrected values independently confirmed against a live RDKit oracle, 0 newly lost. **Breaking**: `torsion_type_for`'s signature changed from `(rings, i, j, k, l, tj, tk)` to `(mol, i, j, k, l, ti, tj, tk, tl)`
- `chematic-mol`: XYZ / multi-frame XYZ read-write (`parse_xyz`/`write_xyz`, `XyzReader`/`XyzWriter`) — explicit hydrogens kept as real atoms, no connectivity/bond-order inference, fails closed on atom-count mismatch or non-finite coordinates
- `chematic-perception`: `stereo_centers(&Molecule) -> Vec<(AtomIdx, bool)>` exposes per-atom tetrahedral-stereocenter classification (issue #263), previously only available as an aggregate count; fixed two bugs found while adding it — a `u64` overflow for negatively-charged atoms in the shared Morgan-rank helper (issue #267), and an implicit-hydrogen rank-0 sentinel colliding with a real atom's normalized rank 0 that silently dropped genuine specified stereocenters
- `chematic-chem`: `ez_completeness(&Molecule) -> EzCompleteness` (issue #264) reports specified/unspecified/total declared E/Z double bonds, matching RDKit's own stereo-bond-eligibility rules (terminal/symmetric bonds excluded, ring bonds <8 atoms excluded via BFS shortest-cycle, not SSSR alone — correctly handles bridged-bicyclic cases like norbornene)
- `chematic-chem`: `detect_atropisomers`/`assign_atropisomer_chirality` are now fully notation-invariant (issues #262, #276) — detection is SSSR-based (two aromatic carbons in separate rings, both with ortho substitution) rather than keyed off whether the SMILES wrote the inter-ring bond explicitly or left it implicit; chirality assignment's own redundant bond-order gate is now aligned with detection's own classification instead of re-deriving a separate, notation-sensitive check
- `chematic-perception`: `is_macrocycle(ring: &[AtomIdx]) -> bool` (issue #266) — a single shared ≥9-atom-ring predicate, replacing duplicated hardcoded thresholds in `chematic-3d`
- **v0.13.0 release-gate note**: 2 of 265 Wave 1 corpus molecules (`chembl_tier_b_0126`/`0168`) show a 1-stereocenter-satisfaction regression, root-caused to a pre-existing distance-geometry embedding defect (present since v0.12.0) that used to be accidentally masked by the now-fixed torsion classification bug — confirmed via RDKit's own MMFF94, given identical starting coordinates, exhibiting the same behavior. Shipped under an explicit waiver (issue #285); not a new defect introduced by this release's MMFF94 fixes
- Full details in `CHANGELOG.md`'s `[0.13.0]` section

**v0.12.0** (2026-08-09): **MMFF94 stretch-bend production fix (breaking), 3D starting-geometry fix for fused/multi-ring molecules**
- `chematic-ff`: `mmff94_stbn` now falls back to RDKit's real 29-row periodic-table-row stretch-bend defaults when the specific/generic MMFF-type table has no row — unconditional production behavior for every MMFF94 policy, not behind an opt-in flag. Missing stretch-bend instances on the 265-molecule Wave 1 corpus: 2,107 → 0. **Breaking**: `mmff94_stbn` gained 3 required `atomic_num_{i,j,k}: u8` parameters (prior type-only behavior kept as new `mmff94_stbn_type_only`); Python's raw `PipelineV2Config(...)` constructor gained a new required `gate_mmff94_stretch_bend` argument (`.safe(...)` unaffected)
- `chematic-3d`: `dg::generate_coords` no longer produces atom-coincident or wildly-stretched starting geometry for several multi-ring topologies (issue #185/#252) — root/ring-vertex collision, ring-fusion-order mismatch, and fixed-offset ring-island anchoring were all independent bugs, not a UFF minimizer defect as issue #185 originally suspected. All 28 `MinimizationFailed` cases on the 265-molecule corpus now resolve to `Ok`, 0 regressions. Two known, separately-tracked residual limitations remain unfixed: fused-ring seam orientation (issue #255) and chain-bridged ring islands (issue #256)
- Full details in `CHANGELOG.md`'s `[0.12.0]` section

**v0.11.0** (2026-08-04): **MMFF94 O2CM typing coverage, SMIRKS/CDXML stereo correctness, 2D/3D layout fixes**
- `chematic-ff`: closed the O2CM terminal-oxygen typing gap (issue #227 Priority 1A-3) — atom-type parity 98.82% → 99.37% on the 265-molecule Wave 1 corpus, oxygen-element parity 95.88% → 100%, strict-gate minimization success 123 → 130/265, 0 cross-element mismatches (unchanged). Issue #227 stays open
- `chematic-rxn`: SMIRKS product chirality assignment made parity-aware — a reordered mapped template neighbor order now correctly inverts/validates the product's `@`/`@@` flag instead of copying it verbatim; inherited (non-template) chirality now fails closed to `Chirality::None` when its neighbor order or mapped topology can't be validated
- `chematic-mol`: CDXML reader now perceives tetrahedral stereo from directional wedges (RDKit issue #9359), wired into the same shared mechanism MOL/MRV already use; non-directional `Bold`/`Hash`/`Dash` displays are opt-in via `CdxmlParseOptions` and, when enabled, invariant to CDXML's B/E bond-atom ordering
- `chematic-depict`: independent (non-fused) ring systems no longer collide at identical/near-identical 2D coordinates
- `chematic-3d`: ETKDG macrocyclic amide 1-4 distance bounds now split by true cis/trans ring-continuation role instead of blanket-pinning all four combinatorial pairs to cis; abstains to a relaxed band when a central amide bond is shared by multiple eligible macrocycles at once
- Full details in `CHANGELOG.md`'s `[0.11.0]` section

**v0.10.1** (2026-08-02): **MMFF94 numeric-typing correctness hotfix**
- `chematic-ff`: fixed a class of bug where MMFF94 could silently resolve an atom against a parameter row belonging to a different element and report the resulting physically-wrong energy as success (issue #227's "furan collision") — the aromatic atom typer never implemented RDKit's real 5-/6-ring alpha/beta-heteroatom classification. Ported from a pinned RDKit source with a new provenance-cited numeric-type registry, plus a construction-time semantic-compatibility invariant that makes this bug class fail closed (`NumericTypeError`) instead of silently wrong, going forward. That invariant caught two more instances of the identical bug: a protonated amine N and an anionic O were each being typed as the *other* element's parameter row. Measured on the 265-molecule Wave 1 corpus (production API): 44 → 102 successful MMFF94 minimizations, 0 cross-element type mismatches across 6693 comparable atoms vs. a pinned RDKit oracle (91.83% exact match)
- This is a correctness hotfix, not a coverage-completion release — issue #227 stays open (140 `MissingParameters` and 22 `MinimizationFailed` cases remain, stretch-bend is not yet gated, no full-corpus energy/gradient parity harness exists yet). See Migration notes in `CHANGELOG.md`'s `[0.10.1]` section if you cache MMFF94 results
- Full details in `CHANGELOG.md`'s `[0.10.1]` section

**v0.10.0** (2026-08-01): **Match-level SMIRKS reaction application, MRV 2D stereo, shared E/Z carrier bond fix**
- `chematic-rxn`: `find_reaction_matches`/`apply_reaction_match` (issue #225) — a public seam between enumerating a SMIRKS's matches against reactant molecules and applying one of them, for callers that need to accept some matches and reject others (e.g. based on whether the matched bond is a ring bond) without discarding the whole `run_reactants` call. `run_reactants`/`run_reactants_strict` are now implemented in terms of these two functions, unchanged in cost (still one SMIRKS parse + one VF2 match pass per call)
- `chematic-mol`: MRV reader now perceives 2D wedge/hash tetrahedral and E/Z stereo (issue #202) — `parse_mrv` previously read wedge/dash bonds and 2D coordinates into `coords_2d` but never converted them into `Atom.chirality`/bond E/Z direction, silently dropping stereochemistry present in the file
- `chematic-smiles`: shared E/Z carrier bonds now resolved via a joint component solver (issue #149) — 10 of 18 previously-abstained fixtures become fully permutation-invariant; the remaining 8 are a documented, RDKit-verified semantically-safe residual (endocyclic double bonds in 5-/6-membered rings, where marker choice has no free degree). Issue #149 stays open pending a scoped fix for the ring-constrained residual
- Full details in `CHANGELOG.md`'s `[0.10.0]` section

**v0.9.0** (2026-08-01): **Opt-in 3D embedding pipeline v2 in Python + WASM, WASM-portable monotonic clock**
- `chematic-py`: `Mol.embed_pipeline_v2(config)` — Python binding for the Rust-only `pipeline_v2::embed_pipeline_v2` (torsion-knowledge-aware distance geometry + stereo verification/repair + policy-gated force field), returning full per-stage evidence (never just final coordinates) and a typed `PipelineV2Error` with structured, diagnostic-only partial evidence on failure. Applies directly to the caller's own atom order — no canonicalize/reparse. Additive; no existing default 3D API changed
- `chematic-wasm`: `embed_pipeline_v2_json(mol, configJson)` — WASM mirror of the above, same config/evidence shape as a tagged-union JSON envelope. New CI job builds both `wasm-pack` targets (`nodejs`/`web`) and runs the Node integration suite on every push/PR — this repo had zero WASM-runtime CI coverage before this
- `chematic-3d`/`chematic-smarts`: fixed `std::time::Instant::now()` panicking unconditionally under real `wasm32-unknown-unknown` (issues #219, #221) — the pipeline v2 binding's own first real-runtime run is what surfaced this pre-existing gap. Fixed via a small crate-internal `clock` module (`web_time::Instant` on wasm32, `std::time::Instant` elsewhere) in each affected crate; no chemistry/geometry/torsion/force-field or timeout-contract change
- `chematic-3d`: `generate_and_minimize_uff()` deprecated (issue #204) — it never ran chematic-ff's real UFF despite its name; kept, not removed, not behavior-changed
- Full details in `CHANGELOG.md`'s `[0.9.0]` section

**v0.8.1** (2026-07-30): **`canonical_smiles()` explicit/implicit hydrogen-count correctness fix**
- `chematic-smiles`/`chematic-core`: two representations of the same molecule that differ only in whether an atom's H count came from bracket notation (`[Cl]`) or organic-subset notation (`Cl`) — when the explicit value merely repeats what valence inference gives anyway — now canonicalize identically. Some canonical SMILES output strings for existing inputs will change as a direct, intended consequence — see `CHANGELOG.md`'s `[0.8.1]` Migration notes if you depend on exact string stability across versions
- Full details in `CHANGELOG.md`'s `[0.8.1]` section

**v0.8.0** (2026-07-29): **Opt-in fail-closed 3D embedding pipeline, canonical-SMILES automorphism-orbit pruning**
- `chematic-3d`: new opt-in `pipeline_v2::embed_pipeline_v2` — stochastic distance geometry + torsion knowledge + stereo verification/repair + typed force-field minimization in one 12-stage pipeline, with a fail-closed stereo re-check *after* minimization. Existing default behavior is unchanged
- `chematic-smiles`: fixed the `canonical_smiles()` performance regression reported by RENKIN (~5x geomean speedup on high-symmetry molecules); also fixes a Dative-bond round-trip bug (issue #194)
- Full details, benchmark numbers, and known limitations in `CHANGELOG.md`

Full version history back to v0.1 — every release, corpus-level before/after
numbers, root causes, and migration notes — is in [CHANGELOG.md](CHANGELOG.md).

---

## Built with chematic

Using chematic in a project? [Share it in Discussions](https://github.com/kent-tokyo/chematic/discussions) or open a PR to add it here.

---

## Reliability by Feature

Not all features have the same validation depth. This table tells you what to trust.

| Feature | Status | Validation |
|---|---|---|
| SMILES parse / write | **Stable** | 4,999-mol ChEMBL comparison; OpenSMILES corpus (parse *correctness*, not canonical-form self-stability — see Canonical SMILES row) |
| Canonical SMILES (structural correctness) | **Stable** | `canonical_smiles(parse(x))` always represents the same molecule as `x`: **100%** across 5,000-mol ChEMBL worst-of-10 *and* a 33-compound acyclic-polyene corpus (retinoids/carotenoids/prostaglandins/leukotrienes/macrolides), each with a verified positive control — was 4.28% corrupting to a different stereoisomer. **Not yet a dedup/cache key** — see Known Limitations below |
| MW / HBA / HBD | **Stable** | 100% RDKit agreement on 4,999 mol |
| TPSA | **Stable** | **100%** on 4,999-mol ChEMBL subset (±0.1 Å²) — see [`docs/validation.md`](https://kent-tokyo.github.io/chematic/validation/) |
| LogP (Crippen) | **Stable** | **100%** on 4,999-mol corpus (max Δ = 1.1×10⁻¹³, within float64 rounding error) |
| ECFP4 / MACCS fingerprints | **Stable** | RDKit comparison + benchmark |
| Tanimoto similarity | **Stable** | RDKit comparison |
| SDF / MOL V2000/V3000 I/O | **Stable** | round-trip tests |
| Substructure search (SMARTS / VF2) | **Stable** | internal test suite |
| PAINS / Brenk filters | **Stable** | rule matching stable; ring-size SMARTS (`[r5]`/`[r6]`) now 0% instability across 5,000-mol worst-of-10 (was ~29–55% before the SSSR fix) |
| Ring perception (SSSR) | **Stable** | Horton algorithm, minimal + deterministic; 0% self-instability across 5,000-mol worst-of-10 (was 50.6%) — see Known Limitations below |
| Murcko scaffold | **Stable** (normalized) | normalized string output **100%** stable across 5,000-mol worst-of-10 (was 0.8% unstable, same root cause as the canonical-SMILES corruption above, now fixed); raw `.smiles` inherits the still-partially-open direction-normalization gap — normalize before comparing (see Known Limitations) |
| 2D SVG depiction | **Stable** | visual spot-checks; not publication-quality |
| 3D conformer (DG + MMFF94) | **Experimental** | reasonable geometry; not equivalent to RDKit ETKDGv3 quality |
| pKa prediction | **Rule-based screening** | 23 SMARTS rules; early triage only, not clinical |
| ADMET (BBB / Caco-2 / hERG / CYP3A4) | **Rule-based screening** | empirical models; directional, not validated on clinical endpoints |
| IUPAC name generation | **Partial** | common compound classes; complex structures may fail |
| Pure-Rust InChI | **Approximate** | enable `native-inchi` feature for bit-exact IUPAC InChI |

Full benchmark methodology → [validation/](validation/) · History → [benchmarks/](benchmarks/)

---

## Known Limitations

- **`canonical_smiles()` is now partially normalized for E/Z stereochemistry — still not safe as a dedup or cache key.** Isolated/simple E/Z double bonds have two equally correct `/`/`\` spellings (e.g. `/N=N/` vs `\N=N\`); the writer previously never normalized between them. Fixed for the general case: every connected E/Z system (a double bond plus every directional bond geometrically tied to it, including whole conjugated chains) is now normalized so its first directional bond in canonical write order is always `/`, regardless of input spelling. Measured on the 5,000-mol ChEMBL corpus, worst-of-10: E/Z-only self-instability (tetrahedral stereo stripped) improved **9.76% → 5.50%** (275/5000 still unstable); structural correctness unaffected by this change (re-verified **0/5000** ChEMBL and **0/33** acyclic-polyene corpus). The residual 275 are confirmed **100% cosmetic** — every unstable case's variants represent the same molecule per RDKit, zero corruption — but the cause is a **mixed pool, not fully root-caused**: about half match a specific motif (a small ring bearing two or more exocyclic double bonds, e.g. cross-conjugated cyclic diimines) where which physical bonds count as "one system" is not yet input-spelling-invariant; the other half is uncharacterized. Until this closes fully, ~1 in 18 stereo-bearing molecules (down from ~1 in 10) can still produce two different, individually valid `canonical_smiles()` strings for the same molecule — do not use it as a dedup or cache key today; document your own dedup key as `apply_aromaticity()`-normalized in the meantime if this matters for your use case.
- **Fail-closed canonical identity key:** `canonical_smiles_stable_key()` reparses and re-canonicalizes the candidate, returning `None` if the result is not idempotent or if multiple independent E/Z systems remain coupled. This is the recommended API for deduplication/cache identity; callers must handle `None` rather than fall back to raw `canonical_smiles()`. It does not claim the historical 275/5000 E/Z-only residual is eliminated. The recovered Issue #11 corpus is pinned in `validation/canonical_original_corpus_manifest.json`; the broader current diagnostic found 31 permutation-sensitive cases (22 semantic-only, 9 with structural differences) under its documented projection.
- **Canonical SMILES structural corruption — fixed.** Before this fix, `canonical_smiles(parse(x))` could silently emit a *different stereoisomer* (not just a differently-spelled but equivalent string) depending on `x`'s input traversal order. Measured on a 5,000-mol ChEMBL subset, worst-of-10 independently-traversed representations per molecule, RDKit-verified structural correctness: **4.28% (214/5000)** of molecules had at least one variant round-trip to the wrong molecule. Root-caused to two independent parser bugs (not the originally-suspected "conjugated double-bond markers are geometrically coupled across bonds" — that diagnosis was disproven, see below), each confirmed via a real found molecule and a minimal regression test: (1) a ring-closure directional-bond (`/`/`\`) marker read at the *closing* occurrence of a ring digit was stored raw instead of flipped to the opening→closing sense, corrupting a conjugated E/Z chain whenever its connecting bond happened to be routed through a ring closure; (2) a stereocenter that opens a ring whose partner closes *inside its own branch* had its neighbor-order resolution keyed by the reusable ring *digit* rather than a unique per-occurrence id, so a later, unrelated reuse of the same digit elsewhere in the SMILES could silently steal and corrupt the stereocenter's neighbor order. **After both fixes: structural correctness is 100% (0/5000) on ChEMBL**, confirmed three times over via independently-ordered reconstructions of the fix (with and without an unrelated third ranking fix, to rule out a hidden dependency). Because both root causes are ring-closure-specific — and retinoids, carotenoids, prostaglandins, leukotrienes, and polyene macrolides carry their long conjugated systems in *acyclic* chains, essentially absent from ChEMBL-random sampling — this was independently re-verified on a dedicated 33-compound corpus of exactly those classes (tretinoin, β-carotene, lycopene, amphotericin B, leukotriene B4, and 28 others; `scripts/polyene_corpus.csv`): **0/33 (0.00%) at worst-of-30**, with a positive control confirming 12/33 (36.36%) corruption on the pre-fix code for this same corpus (all 12 failures were ring-closure-heavy structures; zero purely-acyclic examples — including fully acyclic lycopene — ever failed, even unpatched). This directly disproves the original "any conjugated chain" diagnosis and closes the investigation with no remaining corruption class identified. Skeleton-only and tetrahedral-only self-*stability* also reached 0% (were 0.16% and 4.36%); raw combined self-stability (all stereo intact) improved 86.02% → 90.28% (13.98% → 9.72% unstable) — the entire remainder is the separate, non-corrupting direction-normalization gap described above, not residual corruption. Round-trip invariance (`canonical(parse(canonical(m))) == canonical(m)`) improved slightly, 98.26% → 98.32%, since it was never measuring the corruption class directly.
- **Ring perception (SSSR) was non-deterministic and non-minimal — fixed.** The old `find_sssr` built a single spanning tree and took one fundamental cycle per non-tree edge, with no redundancy to recover a smaller ring when the tree's shape made one unnecessarily large (naphthalene, `c1ccc2ccccc2c1`, deterministically returned ring sizes `[6, 10]` instead of `[6, 6]`). `find_sssr` now uses Horton's algorithm (candidate cycles from every vertex × every edge via shortest-path trees, O(V·E) candidates, canonical-rank tie-break for determinism), giving a genuinely minimum-weight, deterministic basis. Measured on a 5,000-mol ChEMBL subset, worst-of-10 independently-traversed representations per molecule: self-stability **100%** (was 50.6%); single-parse ring-size agreement with RDKit **98.9%** (was 72.4%) — the residual ~1.1% gap is RDKit's own `GetSymmSSSR` legitimately returning *more* rings than the topological minimum for symmetric fused systems (e.g. cubane: μ=5, RDKit=6), not a chematic bug; full symmetrization (Vismara relevant cycles) is future work, not required for correctness. Downstream wins, same corpus: ring-size SMARTS `[r5]`/`[r6]` **0%** instability (was 29–55%), `NumAromaticRings` **0%** (was ~4%), `RingCount`/MW/TPSA/HBA/HBD/LogP/MR unaffected (were already 0%). Two known-narrow exceptions where the *old* SSSR bug had accidentally compensated for a separate, still-open aromaticity bug — see the Aromaticity model bullet below. Full methodology: `scripts/ringinfo_parity.py`.
- **Murcko scaffold: ring topology and normalized string output are now fully stable.** The previously-reported "100% traversal-order instability" was itself a measurement-harness bug (comparing `Mol` objects by Python identity instead of value — always reported "unstable" regardless of the real result); that script bug is fixed (`scripts/ring_collateral_damage.py`). Re-measured on a 5,000-mol worst-of-10 run after the canonical-SMILES corruption fixes above: after normalizing (`apply_aromaticity().canonical_smiles_mode("nostereo")`), self-stability is **100% (0/5000 unstable)**, down from a 0.8% residual — confirming that residual was the same canonical-SMILES structural corruption, not a Murcko ring-selection bug, and it is now fully resolved. Raw isomeric `scaffold().smiles` string comparison (no normalization) is **79.30%** stable (20.70% unstable, was ~45%, essentially unchanged by the partial E/Z-normalization fix above — scaffolds strip most of the side-chain motifs that fix improves) — the remainder is the still-partially-open, non-corrupting `/`/`\` direction-normalization gap described above, not a scaffold-specific issue. `scaffold()` extracts the correct ring system reliably; compare via `mol.apply_aromaticity().canonical_smiles_mode("nostereo")` rather than raw `.smiles` if you need string equality across differently-ordered input.
- **Aromaticity model**: chematic applies Hückel 4n+2 per SSSR ring independently; RDKit uses fused-ring electron delocalization. Visible differences in N-heterocycles (pyridone, quinolone, indolizine). Current benchmark on 4,999-mol ChEMBL subset: HBA/HBD/aromatic ring count **100%**; TPSA **100%** (±0.1 Å²); LogP **100%** (±0.01). Aromaticity-flag parity on Kekulized input measured worst-of-10-representations: **96.3%** (`scripts/aromaticity_atom_parity.py`) — the default-path gap is root-caused to an `aromatic_context` bypass mechanism. The opt-in `AromaticityAlgorithm::RdkitLike` gate covers purine and azulene (9 and 10 aromatic atoms respectively), but this is a model-specific regression boundary, not universal RDKit parity; bridgehead-N and other fused/non-alternant residuals remain explicit.
- **Explicit fused/non-alternant model:** `AromaticityAlgorithm::RdkitLike` is the opt-in whole-graph model for fused and non-alternant systems. Its public regression gate covers purine (9 aromatic atoms) and azulene (10); the compatibility-preserving per-SSSR Hückel default remains intentionally distinct. Applications needing RDKit-like parity must select the model explicitly and record that choice.

---

## Repository Structure

```
chematic/
├── Cargo.toml                    workspace root (v0.89.0 candidate)
├── CHANGELOG.md
├── crates/
│   ├── chematic-core/            Atom, Bond, Molecule, Element, kekulization (4-pass + blossom)
│   ├── chematic-smiles/          OpenSMILES parser/writer, canonical SMILES
│   ├── chematic-perception/      SSSR, 2-pass Hückel aromaticity, CIP stereo
│   ├── chematic-smarts/          SMARTS parser, VF2 subgraph isomorphism, MCS, LRU cache
│   ├── chematic-chem/            190+ descriptors, pKa, ADMET, BOILED-Egg, QED, SA Score,
│   │                             PAINS/Brenk filters, scaffold, standardization, BRICS/RECAP
│   ├── chematic-fp/              ECFP/FCFP, MACCS, MAP4, AtomPair, Torsion, MHFP, ERG
│   ├── chematic-ff/              MMFF94 full stack (7 terms), DREIDING, L-BFGS minimizer
│   ├── chematic-3d/              ETKDG, MD, SASA, USR shape screen, WHIM, GETAWAY, XYZ/PDB I/O
│   ├── chematic-depict/          2D SVG rendering, grid layout, CPK colors, highlighting
│   ├── chematic-rxn/             Reaction SMILES/SMIRKS, RunReactants, RECAP/BRICS
│   ├── chematic-mol/             SDF/MOL V2000+V3000, CML, CDXML parser/writer
│   ├── chematic-inchi/           InChI/InChIKey (pure-Rust approx + IUPAC-exact via native-inchi)
│   ├── chematic-iupac/           IUPAC name generation (25+ compound classes)
│   ├── chematic-mcp/             MCP server — 20 AI-callable tools (JSON-RPC 2.0 over stdio)
│   ├── chematic-wasm/            218+ WASM exports → npm @kent-tokyo/chematic
│   ├── chematic-py/              PyO3 Python bindings → pip install chematic
│   ├── chematic-ewald/           PME Ewald summation, B-spline interpolation
│   ├── chematic-crystal/         Periodic crystal structures: lattice, PBC, neighbors, supercells, POSCAR/CONTCAR I/O (not Molecule)
│   └── chematic/                 Umbrella crate with feature flags
├── demo/                         Interactive WASM playground (→ /playground/ on GitHub Pages)
│   ├── index.html
│   └── pkg/                      Pre-built WASM bundle (rebuilt on each release)
└── docs/                         MkDocs documentation site source
    ├── cookbook.md
    ├── getting_started/
    └── api/
```

---

## Development Commands

```bash
cargo build --workspace                                                   # build all crates
cargo test --workspace --lib --quiet                                      # 3,912 lib tests
cargo test -p chematic-inchi --features native-inchi --test standard_inchi  # +16 InChI tests
cargo clippy --workspace -- -D warnings                                   # lints (zero warnings)
```

---

## Citation

If you use chematic in academic or research work, please cite:

```bibtex
@software{chematic,
  author    = {Kentaro Tanabe (kent-tokyo)},
  title     = {chematic: A pure-Rust cheminformatics toolkit},
  url       = {https://github.com/kent-tokyo/chematic},
  version   = {0.31.0},
  year      = {2026},
}
```

---

## License

Licensed under either of Apache License 2.0 or MIT License, at your option.
Copyright attribution: Kentaro Tanabe (kent-tokyo). See [`NOTICE`](NOTICE) for the
redistribution attribution notice.

---

If chematic saves you time, a [GitHub star](https://github.com/kent-tokyo/chematic) helps others discover it.
