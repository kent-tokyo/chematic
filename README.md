# chematic

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![docs.rs](https://docs.rs/chematic/badge.svg)](https://docs.rs/chematic)

![Pure Rust](https://img.shields.io/badge/Pure%20Rust-zero%20C%2B%2B-orange?logo=rust)
![WASM](https://img.shields.io/badge/WASM-504%20KB-blueviolet?logo=webassembly)
![MCP](https://img.shields.io/badge/MCP-agent%20ready-purple)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)
[![Demo](https://img.shields.io/badge/demo-live-brightgreen)](https://kent-tokyo.github.io/chematic/playground/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/kent-tokyo/chematic/blob/main/notebooks/quickstart.ipynb)

[日本語](README_ja.md) | [中文](README_zh.md)

A cheminformatics library for Python, Rust, and the browser.

**Cheminformatics that's fast by default, safe by design.**  
Pure Rust · Zero C/C++ · Python · WebAssembly · [Live Demo](https://kent-tokyo.github.io/chematic/playground/)

| | chematic | RDKit (Python) | RDKit.js (WASM) |
|---|---|---|---|
| **Get started** | `pip install chematic` | conda / cmake required | no Python bindings |
| **Browser bundle** | **504 KB** | not available | ~30 MB (60× larger) |
| **Batch fingerprints** | **3.6 µs/mol** (5–14× faster) | 20–50 µs/mol | — |
| **Memory safety** | compiler-enforced (Rust) | C++ | C++ |
| **Build from source** | `cargo build` only | cmake + clang + Boost | Emscripten SDK |

All numbers are reproducible — see [benchmark details](https://kent-tokyo.github.io/chematic/benchmark/).  
WASM sizes: chematic **504 KB** · RDKit.js ~30 MB · Indigo WASM ~40 MB

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
| **Molecule search** | ECFP4/MACCS fingerprints, Tanimoto, LSH approximate nearest-neighbour |
| **AI agent / MCP** | Built-in MCP server — Claude Desktop can call chemistry tools directly |
| **Browser app** | 504 KB WASM bundle, zero backend required, React/Vue/Svelte ready |
| **Jupyter notebook** | `mol` renders SVG inline; `descriptors_df()` returns a pandas DataFrame |
| **Batch analysis** | Rayon-parallel descriptor/fingerprint/3D pipelines; SDF/CSV in, CSV out |
| **Rust server** | Pure-Rust crates with no C/C++ toolchain; Axum/Actix compatible |

Full worked examples → [Use cases](https://kent-tokyo.github.io/chematic/use-cases/)

---

## When to use chematic

**Use chematic if:**

- You want chemistry in the browser (WASM, 504 KB, no server required)
- You need a pure Rust stack with no C++ toolchain dependencies
- You deploy to environments where `pip install rdkit` is impractical (Cloudflare Workers, Lambda, embedded)
- You build AI agents and want native MCP tool integration
- You process molecules in batch at high throughput (ECFP4: 5–14× faster than RDKit)
- You want `pip install chematic` to just work — anywhere, no compiler needed

**Use RDKit if:**

- You need maximum ecosystem compatibility and 20+ years of production validation
- You need publication-quality 3D structures with ML-assisted torsion corrections (RDKit's ETKDGv3)
- You need bit-exact standard InChI without enabling the `native-inchi` feature
- You depend on community plugins written against the RDKit Python API

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
[RDKit compatibility guide](docs/rdkit_compat.md) for the compatibility matrix,
differential-validation results vs RDKit, and runnable examples.

### Diagnostics

```python
import chematic
chematic.doctor()
# chematic v0.4.29
# Python 3.12.x  |  darwin arm64
#
# Descriptor accuracy (benchmark 2026-06, v0.4.29 vs RDKit 2026.03.3):
#   MW / HBA / HBD / ARC  100%   (4,999-mol ChEMBL subset)
#   TPSA                  100%   within ±0.1 Å²
#   LogP (Crippen)        100%*  (max Δ = 1.1×10⁻¹³)
#   Num stereocenters     99.98% (legacy) / 98.7% (new CIP FindPotentialStereo)
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

15 chemistry tools are callable from any MCP-compatible agent:

| Tool | What it does |
|---|---|
| `name_to_smiles` | Resolve "aspirin", "caffeine", … to SMILES via PubChem |
| `calc_properties` | MW, LogP, TPSA, HBA/HBD, QED, SA Score, pKa, ADMET |
| `smarts_match` | Substructure search |
| `pains_check` / `brenk_check` | Flag assay interference or reactive groups |
| `generate_3d` | 3D coordinates (ETKDG + MMFF94) |
| `find_mcs` | Maximum common substructure |
| + 9 more | `ecfp4`, `tanimoto`, `canonical_smiles`, `admet_profile`, `boiled_egg`, `sa_score`, `lipinski_check` … |

---

## Why Pure Rust?

### Fast

Rust's zero-cost abstractions and ownership model eliminate overhead at the source.
chematic's ECFP4 fingerprint batch pipeline runs at **3.6 µs/mol** — 5–14× faster
than RDKit's Python API on the same hardware. No GIL, no interpreter overhead, no
FFI call overhead hidden inside a `_sys` crate.

### Safe

The entire default dependency tree contains **~6 `unsafe` blocks** across 15,000+ lines
of Rust. No C++ heap corruptions. No segfaults from malformed SMILES input. No
platform-specific build failures from `-sys` crates. The compiler enforces memory
safety at every call site.

> The `native-inchi` feature is the single opt-in exception — it vendors the IUPAC InChI
> C library (v1.07.5) for bit-exact standard InChI. All other crates stay FFI-free.

### Anywhere

Pure Rust compiles to `wasm32-unknown-unknown` natively — no Emscripten, no `cmake`,
no `clang`. The npm package `@kent-tokyo/chematic` is **504 KB gzip** — 60× smaller
than RDKit.js. One codebase runs on Linux, macOS, Windows, and in every browser.

---

## Benchmarks & Validation

| Metric | Result | Corpus |
|--------|--------|--------|
| ECFP4 throughput | **3.6 µs/mol** (5–14× vs RDKit) | 4,999-mol ChEMBL subset |
| HBA / HBD / aromatic ring count | **100% RDKit agreement** | 4,999-mol ChEMBL subset |
| TPSA | **100% RDKit agreement** within ±0.1 Å² | 4,999-mol ChEMBL subset |
| LogP (Crippen) | **100% RDKit agreement**\* | 4,999-mol ChEMBL subset |
| Num stereocenters | **99.98%** vs legacy†; 98.7% vs new CIP | 4,999-mol ChEMBL subset |
| WASM bundle | **504 KB** gzip | — |

\*LogP max Δ = 1.1×10⁻¹³ across 4,999 molecules — within float64 rounding error.  
†Stereocenters: 99.98% vs legacy `CalcNumAtomStereoCenters` (1 molecule where chematic matches `FindPotentialStereo`=4 and legacy under-counts at 2); 98.7% vs new-CIP `FindPotentialStereo` (67 cage/bridgehead molecules where both chematic and legacy correctly return fewer than the new oracle). chematic is calibrated between both extremes.

All numbers are reproducible with the scripts in this repo.  
Full history → [benchmarks/](benchmarks/) · Methodology → [validation/](validation/)

---

## Comparison with Other Cheminformatics Libraries

| Feature                 | **chematic**                              | RDKit (rdkit-sys)  | OpenBabel FFI  | RDKit.js (WASM)    |
|-------------------------|-------------------------------------------|--------------------|----------------|--------------------|
| **C/C++ dependencies**  | **None (default)**†                       | Extensive C++      | Extensive C++  | C++ via Emscripten |
| **WASM binary size**    | **~500 KB** (504 KB gzip)                 | N/A (no WASM)      | N/A (no WASM)  | ~30 MB             |
| **Build requirement**   | `cargo build` only                        | cmake + clang      | cmake + clang  | Emscripten SDK     |
| **WASM target support** | **Full (native)**                         | No                 | No             | Yes (Emscripten)   |
| **Python bindings**     | **Yes** (`pip install chematic`, PyO3)    | Yes (rdkit-sys)    | Yes            | No                 |
| **Unsafe Rust**         | **None (default)**†                       | Extensive          | Extensive      | N/A                |

<details>
<summary>Full feature comparison (30+ capabilities)</summary>

| Feature                                      | **chematic**                                     | RDKit (rdkit-sys)   | OpenBabel FFI  | RDKit.js (WASM)   |
|----------------------------------------------|--------------------------------------------------|---------------------|----------------|-------------------|
| OpenSMILES parser                            | Full                                             | Full                | Full           | Full              |
| SMILES writer / canonical                    | Yes                                              | Yes                 | Yes            | Yes               |
| Kekulization                                 | **4-pass (incl. Edmonds' blossom)**              | Yes                 | Yes            | Yes               |
| Ring perception (SSSR)                       | Yes + iterative augmentation                     | Yes                 | Yes            | Yes               |
| SDF/MOL V2000+V3000 + SD fields              | Yes                                              | Yes                 | Yes            | Yes               |
| Tripos MOL2 format                           | **Yes** (parser + writer)                        | Yes                 | Yes            | No                |
| 2D depiction (SVG, CPK colors, **PDF, EPS**) | Yes                                              | Yes                 | Yes            | Yes               |
| ECFP/FCFP fingerprints (2/4/6)               | **All variants + bitvec**                        | Yes                 | Yes            | Yes               |
| AtomPair / Torsion / MACCS FP                | Yes                                              | Yes                 | Yes            | Yes               |
| **MAP4 fingerprint**                         | **Yes** (Minervini 2020)                         | No (external pkg)   | No             | No                |
| Molecular descriptors                        | **190+ descriptor values** (71 functions; MQN×42, BCUT2D, autocorr2d return multi-value arrays) | ~30  | ~20            | ~30               |
| **Topological descriptors**                  | **Yes** (Petitjean, Hosoya Z, ECI, Moran, Geary) | Partial            | Partial        | No                |
| BRICS / RECAP fragmentation                  | Yes                                              | Yes                 | No             | Yes               |
| Murcko scaffold                              | Yes                                              | Yes                 | No             | Yes               |
| Tautomer normalisation                       | Yes                                              | Yes                 | No             | Yes               |
| MCS                                          | Yes                                              | Yes                 | No             | Yes               |
| Stereoisomer enumeration                     | **Yes**                                          | Yes                 | No             | Yes               |
| CIP stereo (R/S, E/Z) detail                 | **Yes (per-atom JSON)**                          | Yes                 | Yes            | Yes               |
| Allene cumulated stereo (`C=C=C`)            | **Yes** (`@`/`@@`, round-trip stable)            | Yes                 | Partial        | No                |
| 3D coordinate generation                     | Yes (DG + MMFF94/DREIDING + L-BFGS)             | Yes (ETKDG)         | Yes            | Yes               |
| 3D shape descriptors (PMI/NPR/USR/…)         | **Yes**                                          | Yes                 | No             | Yes               |
| **3D GETAWAY descriptors (HATS-matrix)**     | **Yes** (19-dim; `whim_getaway_combined` 29-dim) | Yes                | No             | No                |
| MMFF94 force field (all 7 energy terms)      | **Yes**                                          | Yes                 | Yes            | No                |
| **UFF force field** (metals, organometallics)| **Yes**                                          | No                  | Yes            | No                |
| AutoDock PDBQT format (parse + write)        | **Yes** (docking pipeline ready)                 | Via Python API      | Yes            | No                |
| SDF with partial charges                     | **Yes** (`write_sdf_with_charges`)               | Yes                 | Yes            | No                |
| MaxMin / Butina diversity picking            | **Yes**                                          | Yes                 | No             | No                |
| Reaction SMILES/SMIRKS                       | Yes                                              | Yes                 | Yes            | Yes               |
| InChI / InChIKey                             | **Yes** — pure-Rust + **IUPAC-exact** via `native-inchi` | C lib required | C lib required | C lib required |
| **pKa prediction**                           | **Yes (15 SMARTS rules)**                        | No                  | No             | No                |
| **ADMET profile** (BBB/Caco-2/hERG/CYP3A4)  | **Yes + BOILED-Egg**                             | Partial             | No             | Partial           |
| **MCP server (AI agent API)**                | **Yes — 15 tools incl. Name→SMILES**            | No                  | No             | No                |
| IUPAC name generation                        | **Yes (25+ classes)**                            | No                  | No             | Partial           |
| Name → SMILES (PubChem proxy)                | **Yes** (`name_to_smiles` MCP tool)              | No                  | No             | No                |
| Maintenance (2026)                           | Active                                           | Active              | Minimal        | Active            |

</details>

† Default build only. The optional `native-inchi` feature adds a C-compiler dependency (and ~6 `unsafe` blocks, all confined to its FFI layer — see "Safe" above) for the vendored IUPAC InChI C library (v1.07.5). All other crates remain FFI-free and unsafe-free.

---

## JavaScript / TypeScript (WebAssembly)

**504 KB gzip — 60× smaller than RDKit.js.** No Emscripten, no cmake. Drop-in for browser or Node.js.

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

130+ exported functions cover descriptors, fingerprints, 3D geometry, reactions, diversity picking, and SDF round-trips.
See the [full WASM API reference](https://kent-tokyo.github.io/chematic/) for all exports.
---

## Crate Reference

| Crate                 | Description                                                                                              | Tests |
|-----------------------|----------------------------------------------------------------------------------------------------------|-------|
| `chematic-core`       | Atom, Bond, Molecule, Element, kekulization (no deps); mutable `add/remove_atom/bond`, `fragments()`, `is_connected()`, `formula_with_isotopes`, `validate_valence`; `StereoGroup`/`StereoGroupKind` | 71    |
| `chematic-smiles`     | OpenSMILES parser, writer, canonical SMILES; **stereo parity correction** (pre-solves RDKit #8775 — `@`/`@@` auto-flipped on odd permutations); **allene cumulated double bond stereo** (`C=C=C` `@`/`@@`, round-trip stable) | 109    |
| `chematic-perception` | SSSR, Hückel aromaticity + antiaromaticity (4n+2 rule), `apply_aromaticity`, `aromatize`/`kekulize_inplace`, `assign_stereo_from_2d`, `assign_ez_from_2d`, `cip_ez_descriptor`; **zero-order/dative bonds excluded from ring perception** | 101    |
| `chematic-mol`        | MOL/SDF V2000+V3000 (R/W with 2D coords, +partial charge writing), CML (R/W), CDXML (R); `SdfRecord` with coords+props; MDL RXN R/W; V3000 stereo-group COLLECTION R/W; **AutoDock PDBQT** (parse + write); **ChemicalJSON** (`parse_cjson`/`write_cjson`, Avogadro/MolSSI format) | 130    |
| `chematic-depict`     | 2D SVG (CPK colors, highlighting, grid), DepictData, `detect_crossings`, `render_svg_with_metadata`, reaction SVG; **PDF output** (`depict_pdf`/`depict_pdf_opts` via svg2pdf); **EPS output** (`depict_eps`/`depict_eps_opts`, pure Rust); `tiny_skia` PNG is optional `png` feature (default on, disabled for WASM) | 64    |
| `chematic-chem`       | 190+ descriptor values (71 functions), tautomers, scaffold, BRICS, QED, standardize, CIP; **pKa prediction** (15 SMARTS rules); **ADMET profile** (BBB/Caco-2/hERG/CYP3A4); **HBA 100% RDKit agreement** (4 999 / 4 999 mol benchmark); **TPSA 100% ±0.1 Å² / LogP 100%\* / HBD 100% / stereocenters 99.98% (legacy) / 98.7% (new CIP)** vs RDKit (4,999-mol ChEMBL); **topological descriptors** (`petitjean_index`, `graph_diameter`, `graph_radius`, `graph_eccentricities`, `eccentric_connectivity_index`, `hosoya_index`, `moran_autocorr`, `geary_autocorr`); **`schultz_mti`, `gutman_mti`, `vabc` (Bondi radii vdW volume), `gravitational_index`**; `clean_stereo_groups()` in standardize | 662   |
| `chematic-fp`         | ECFP2/4/6, FCFP4/6, MACCS, TopoPF, AtomPair, Torsion, Layered, Pattern, Pharmacophore, Reaction, **MAP4** (Minervini 2020, not in RDKit) — Tanimoto/Dice; bulk similarity | 185    |
| `chematic-ff`         | **MMFF94 all 7 terms** (Halgren 1996): Bond/Angle/Torsion/vdW/Elec + **OOP** (117 entries) + **Stretch-Bend** (282 entries); steepest-descent + L-BFGS optimizer, torsion scan, energy breakdown; DREIDING typing; **UFF** (metals/organometallics: Zn, Fe, Cu, …) | 98    |
| `chematic-smarts`     | SMARTS, VF2, MCS with chirality matching; **SmartsCache** (LRU compilation cache, 5–20×); **named_pattern()** library (20 functional group patterns); **atom map `:N` in SMARTS** (`[O;D1;H0:3]` — stored as metadata, not a match criterion); **`[kN]` ring-size primitive**; **VF2 early-exit** when query > target atom count; **`find_matches_with_rings`** — share SSSR across multi-pattern batches | 142   |
| `chematic-3d`         | 3D coordinate generation, distance geometry constraints, ETKDG KB (40 torsion patterns, adaptive noise), force-field minimization, shape descriptors, ConformerEnsemble with RMSD pruning, PDB/XYZ; **GETAWAY HATS-matrix** (full 19-dim implementation); **`whim_getaway_combined()`** now 29-dim | 265    |
| `chematic-rxn`        | Reaction SMILES/SMIRKS, `run_reactants`/`run_reactants_strict`; **`retro_disconnect()`** — 60 retro-SMIRKS templates (AmideBond/Ester/Ether/CNBond/CCBond/CSBond) + SA Score ranking; **parity-aware `@`/`@@` SMIRKS stereo filtering**; **E/Z double-bond stereo filtering** in `run_reactants` (`ez_stereo_outward`, `smirks_ez_stereo_ok`) | 137    |
| `chematic-inchi`      | InChI/InChIKey: pure-Rust approximation (WASM) **+ IUPAC-standard** via `native-inchi` feature (vendored C lib 1.07.5, bit-exact); **parse_inchi** reader | 96 (+16*)    |
| `chematic-wasm`       | **130+ WASM exports** — npm: `@kent-tokyo/chematic` v0.4.19 (~500 KB, 504 KB gzip); pKa/ADMET/BBB/Caco-2/hERG/CYP3A4; `smiles_to_pdbqt`, `minimize_uff_json` | 211   |
| `chematic-iupac`      | Local IUPAC name generation — **25+ compound classes**: alkanes, cycloalkanes, alkenes/alkynes, alcohols, amines, halides, aldehydes, ketones, acids, esters, amides, **piperidine, morpholine, piperazine, naphthalene, sulfides** | 47    |
| `chematic-mcp`        | **MCP (Model Context Protocol) server** — AI agent integration; **20 tools**: parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, name_to_smiles, retrosynthesis, smiles_to_moljson, moljson_to_smiles, representation_router, **molecule_context_pack** | 31    |
| `chematic-py`         | PyO3 Python bindings (`pip install chematic`); 300+ API endpoints: `from_smiles()`, `Mol.descriptors()`, `Mol.minimize_dreiding()`, `from_cxsmiles()`, `from_rxn_file()`/`to_rxn_file()`, `parse_sdf_with_coords()`, `Mol.ring_families()`, `tanimoto_matrix()`, `iter_sdf()`, `SimilarityIndex`; **`mol.to_pdf()`/`mol.to_eps()`** (depict); **`from_cjson()`/`mol.to_cjson()`** (ChemicalJSON); **`mol.schultz_mti`, `mol.gutman_mti`, `mol.vabc`, `mol.gravitational_index`**; **`bulk.substructure_match(smarts, mols)`** (parallel VF2 on pre-parsed Mol objects); **`mol.describe()`** (LLM/MCP-ready natural-language summary); **`mol.diff(other)`** (element + descriptor diff); Sprint 18–27 coverage | 300+  |
| `chematic-ewald`      | PME Ewald summation, B-spline interpolation (cubic, phase-corrected)                                     | 16    |
| `chematic`            | Umbrella crate with feature flags (all sub-crates, incl. `iupac`, `inchi`)                              | 1     |

```
cargo test --workspace --lib --quiet                                          # 2,366 tests, all passing
cargo test -p chematic-inchi --features native-inchi --test standard_inchi  # +16 IUPAC-exact InChI tests
```

---

## Recent Development (v0.4.x Era)

**v0.4.19** (2026-06-23): **PDF/EPS output, ChemicalJSON, new descriptors, WASM −38.5%**
- `chematic-depict`: `depict_pdf()` / `depict_eps()` — PDF and EPS output; pure Rust, no external tools
- `chematic-mol`: **ChemicalJSON** — `parse_cjson()` / `write_cjson()` for Avogadro2 / MolSSI interop
- `chematic-chem`: 4 new descriptors — `schultz_mti()`, `gutman_mti()`, `vabc()` (Bondi vdW volume), `gravitational_index()`
- `chematic-3d`: **Spectrophores** 3D fingerprints (pharmacophore shell encoding)
- `chematic-py`: `mol.to_pdf()`, `mol.to_eps()`, `mol.to_cjson()`, `from_cjson()`; `bulk.substructure_match(smarts, mols)` parallel VF2; `estate_all()` and `ring_bundle` in bulk
- **WASM bundle: 819 → 504 KB gzip (−38.5%)** — `tiny_skia` made optional, inline SHA-256, `opt-level="z" lto=true codegen-units=1`

**v0.4.18** (2026-06-23): **Python API expansion + benchmark docs**
- `chematic-py`: **Jupyter auto-display** — writing `mol` in a cell renders 2D structure via `_repr_svg_()`; `mol.has_substructure(smarts)`, `mol.find_matches(smarts)`; `from_smiles_list()`, `descriptors_df()`
- `chematic-chem`: `chi_all()` — all 10 Hall-Kier connectivity indices in a single pass; `cns_mpo_from_parts()`; `pains_passes_and_matches()` / `brenk_passes_and_matches()` — combined pass/match in one scan
- Docs: benchmark page added (ECFP4 5–14× vs RDKit, 100% descriptor accuracy on 4,999-mol ChEMBL corpus)

**v0.4.16–v0.4.17** (2026-06-22–23): **SSSR sharing performance sprint**
- `chematic-smarts`: `find_matches_with_rings()` — share a pre-computed `RingSet` across all patterns in a batch
- `chematic-chem`: Crippen 117 SSSR → 1 per `logp_crippen` call; PAINS ~480 → 1; QED 113 → 1; pKa 42 → 1; new `logp_and_mr()`, `logd_from_logp()`, `pka_both()` to avoid redundant passes
- `chematic-fp`: MHFP incremental BFS — 3N → N BFS operations per molecule at radius=2

**v0.4.15** (2026-06-21): **TPSA calibration + E/Z stereo in reactions**
- `chematic-chem`: TPSA ±0.1 Å² calibration sprint — **HBA 100%, HBD 100%, aromatic ring count 100%** on 4,999-mol ChEMBL subset; TPSA 86.7% → 93.3% (4,999-mol), 100% on 175-mol drug-like set
- `chematic-rxn`: E/Z double-bond stereo filtering in `run_reactants` — SMIRKS `/`/`\` geometry matching via `smirks_ez_stereo_ok()` / `ez_stereo_outward()`

**v0.4.14** (2026-06-21): **Topological descriptors + stereo correctness**
- `chematic-chem`: 8 topological descriptors — `petitjean_index()`, `graph_eccentricities()`, `graph_diameter()`, `graph_radius()`, `eccentric_connectivity_index()`, `hosoya_index()`, `moran_autocorr()`, `geary_autocorr()`
- `chematic-3d`: GETAWAY HATS-matrix (19-dim); `whim_getaway_combined()` now 29-dim
- `chematic-smiles`: allene cumulated stereo `C=C=C` `@`/`@@` — round-trip stable
- `chematic-smarts`: `[kN]` ring-size primitive; VF2 early-exit when query > target atom count
- `chematic-rxn`: parity-aware SMIRKS chirality matching; product bracket cleanup (`[O:1]` → `O`)
- `chematic-perception`: zero-order/dative bonds excluded from SSSR; `count_aromatic_rings()` handles Kekulé input

**v0.4.13** (2026-06-21): **Template retrosynthesis + descriptor fixes**
- `chematic-rxn`: `retro_disconnect()` — 60 retro-SMIRKS templates (AmideBond / Ester / Ether / CNBond / CCBond / CSBond) with SA Score ranking; Python `mol.retro_disconnect(reaction_class=...)`
- `chematic-3d`: ETKDG torsion KB 28 → 40 patterns; adaptive bond-flexibility noise scaling
- `chematic-chem`: `hbd_count()` now includes S-H (thiol); TPSA nitro-N / aromatic oxide bridge / Kekulé-N corrections

**v0.4.9–v0.4.12** (2026-06-19–21): **AutoDock, UFF, SMARTS atom-map, ring augmentation**
- `chematic-mol`: AutoDock PDBQT parse/write; `write_sdf_with_charges`
- `chematic-ff`: UFF force field for metals/organometallics (Zn, Fe, Cu, …)
- `chematic-smarts`: atom map `:N` in SMARTS (`[O;D1;H0:3]` — stored as metadata)
- `chematic-perception`: iterative `augmented_ring_set` for fused polycyclic aromatic ring counting (222/222 bench5k fixes)
- MCP: 15th tool `name_to_smiles` via PubChem REST proxy

**v0.4.5–v0.4.7** (2026-06-19): **Kekulization blossom + BOILED-Egg + InChI E/Z**
- Edmonds' blossom algorithm for non-bipartite aromatic graphs (128→2 failures)
- InChI `/b` E/Z layer, 6 new MCP tools, BOILED-Egg descriptor + Python/WASM bindings

**v0.4.0–v0.4.4** (2026-06-17–18): **PyO3 Python bindings + native-inchi**
- `chematic-py`: PyO3/maturin bindings — `from_smiles()`, `Mol.aromatic_ring_count`, `Mol.descriptors()`
- `native-inchi` feature: IUPAC-exact InChI via vendored C lib v1.07.5
- HBA rewrite: 99.98% agreement with RDKit (4,999-mol ChEMBL benchmark)


Full changelog: [CHANGELOG.md](CHANGELOG.md)

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
| TPSA | **Stable** | 100% on 175-mol drug-like set; **99.7%** on 4,999-mol ChEMBL subset (±0.1 Å²) |
| LogP (Crippen) | **Stable** | **100%** on 4,999-mol corpus (±0.01); ~99% on 175-mol drug-like set (±0.3) |
| ECFP4 / MACCS fingerprints | **Stable** | RDKit comparison + benchmark |
| Tanimoto similarity | **Stable** | RDKit comparison |
| SDF / MOL V2000/V3000 I/O | **Stable** | round-trip tests |
| Substructure search (SMARTS / VF2) | **Stable** | internal test suite |
| PAINS / Brenk filters | **Stable** | rule matching stable; ring-size SMARTS (`[r5]`/`[r6]`) now 0% instability across 5,000-mol worst-of-10 (was ~29–55% before the SSSR fix) |
| Ring perception (SSSR) | **Stable** | Horton algorithm, minimal + deterministic; 0% self-instability across 5,000-mol worst-of-10 (was 50.6%) — see Known Limitations below |
| Murcko scaffold | **Stable** (normalized) | normalized string output **100%** stable across 5,000-mol worst-of-10 (was 0.8% unstable, same root cause as the canonical-SMILES corruption above, now fixed); raw `.smiles` inherits the still-open direction-normalization gap — normalize before comparing (see Known Limitations) |
| 2D SVG depiction | **Stable** | visual spot-checks; not publication-quality |
| 3D conformer (DG + MMFF94) | **Experimental** | reasonable geometry; not equivalent to RDKit ETKDGv3 quality |
| pKa prediction | **Rule-based screening** | 15 SMARTS rules; early triage only, not clinical |
| ADMET (BBB / Caco-2 / hERG / CYP3A4) | **Rule-based screening** | empirical models; directional, not validated on clinical endpoints |
| IUPAC name generation | **Partial** | common compound classes; complex structures may fail |
| Pure-Rust InChI | **Approximate** | enable `native-inchi` feature for bit-exact IUPAC InChI |

Full benchmark methodology → [validation/](validation/) · History → [benchmarks/](benchmarks/)

---

## Known Limitations

- **`canonical_smiles()` is not yet a canonical form for molecules with E/Z stereochemistry — do not use it as a dedup or cache key today.** Structural correctness is unaffected (see below), but ~9.7% of a 5,000-mol ChEMBL corpus produce two *different, individually valid* strings for the same molecule depending on input atom order (isolated E/Z double bonds have two equally correct `/`/`\` spellings — `/N=N/` and `\N=N\` — and the writer does not yet normalize which one it emits). Any code path relying on `canonical_smiles()` output for deduplication, hashing, or a cache/DB key will silently treat ~1 in 10 stereo-bearing molecules as distinct when they are the same molecule. Not yet fixed — tracked for a future normalization pass; document your own dedup key as `apply_aromaticity()`-normalized in the meantime if this matters for your use case.
- **Canonical SMILES structural corruption — fixed.** Before this fix, `canonical_smiles(parse(x))` could silently emit a *different stereoisomer* (not just a differently-spelled but equivalent string) depending on `x`'s input traversal order. Measured on a 5,000-mol ChEMBL subset, worst-of-10 independently-traversed representations per molecule, RDKit-verified structural correctness: **4.28% (214/5000)** of molecules had at least one variant round-trip to the wrong molecule. Root-caused to two independent parser bugs (not the originally-suspected "conjugated double-bond markers are geometrically coupled across bonds" — that diagnosis was disproven, see below), each confirmed via a real found molecule and a minimal regression test: (1) a ring-closure directional-bond (`/`/`\`) marker read at the *closing* occurrence of a ring digit was stored raw instead of flipped to the opening→closing sense, corrupting a conjugated E/Z chain whenever its connecting bond happened to be routed through a ring closure; (2) a stereocenter that opens a ring whose partner closes *inside its own branch* had its neighbor-order resolution keyed by the reusable ring *digit* rather than a unique per-occurrence id, so a later, unrelated reuse of the same digit elsewhere in the SMILES could silently steal and corrupt the stereocenter's neighbor order. **After both fixes: structural correctness is 100% (0/5000) on ChEMBL**, confirmed three times over via independently-ordered reconstructions of the fix (with and without an unrelated third ranking fix, to rule out a hidden dependency). Because both root causes are ring-closure-specific — and retinoids, carotenoids, prostaglandins, leukotrienes, and polyene macrolides carry their long conjugated systems in *acyclic* chains, essentially absent from ChEMBL-random sampling — this was independently re-verified on a dedicated 33-compound corpus of exactly those classes (tretinoin, β-carotene, lycopene, amphotericin B, leukotriene B4, and 28 others; `scripts/polyene_corpus.csv`): **0/33 (0.00%) at worst-of-30**, with a positive control confirming 12/33 (36.36%) corruption on the pre-fix code for this same corpus (all 12 failures were ring-closure-heavy structures; zero purely-acyclic examples — including fully acyclic lycopene — ever failed, even unpatched). This directly disproves the original "any conjugated chain" diagnosis and closes the investigation with no remaining corruption class identified. Skeleton-only and tetrahedral-only self-*stability* also reached 0% (were 0.16% and 4.36%); raw combined self-stability (all stereo intact) improved 86.02% → 90.28% (13.98% → 9.72% unstable) — the entire remainder is the separate, non-corrupting direction-normalization gap described above, not residual corruption. Round-trip invariance (`canonical(parse(canonical(m))) == canonical(m)`) improved slightly, 98.26% → 98.32%, since it was never measuring the corruption class directly.
- **Ring perception (SSSR) was non-deterministic and non-minimal — fixed.** The old `find_sssr` built a single spanning tree and took one fundamental cycle per non-tree edge, with no redundancy to recover a smaller ring when the tree's shape made one unnecessarily large (naphthalene, `c1ccc2ccccc2c1`, deterministically returned ring sizes `[6, 10]` instead of `[6, 6]`). `find_sssr` now uses Horton's algorithm (candidate cycles from every vertex × every edge via shortest-path trees, O(V·E) candidates, canonical-rank tie-break for determinism), giving a genuinely minimum-weight, deterministic basis. Measured on a 5,000-mol ChEMBL subset, worst-of-10 independently-traversed representations per molecule: self-stability **100%** (was 50.6%); single-parse ring-size agreement with RDKit **98.9%** (was 72.4%) — the residual ~1.1% gap is RDKit's own `GetSymmSSSR` legitimately returning *more* rings than the topological minimum for symmetric fused systems (e.g. cubane: μ=5, RDKit=6), not a chematic bug; full symmetrization (Vismara relevant cycles) is future work, not required for correctness. Downstream wins, same corpus: ring-size SMARTS `[r5]`/`[r6]` **0%** instability (was 29–55%), `NumAromaticRings` **0%** (was ~4%), `RingCount`/MW/TPSA/HBA/HBD/LogP/MR unaffected (were already 0%). Two known-narrow exceptions where the *old* SSSR bug had accidentally compensated for a separate, still-open aromaticity bug — see the Aromaticity model bullet below. Full methodology: `scripts/ringinfo_parity.py`.
- **Murcko scaffold: ring topology and normalized string output are now fully stable.** The previously-reported "100% traversal-order instability" was itself a measurement-harness bug (comparing `Mol` objects by Python identity instead of value — always reported "unstable" regardless of the real result); that script bug is fixed (`scripts/ring_collateral_damage.py`). Re-measured on a 5,000-mol worst-of-10 run after the canonical-SMILES corruption fixes above: after normalizing (`apply_aromaticity().canonical_smiles_mode("nostereo")`), self-stability is **100% (0/5000 unstable)**, down from a 0.8% residual — confirming that residual was the same canonical-SMILES structural corruption, not a Murcko ring-selection bug, and it is now fully resolved. Raw isomeric `scaffold().smiles` string comparison (no normalization) is **79.24%** stable (20.76% unstable, was ~45%) — the remainder is the still-open, non-corrupting `/`/`\` direction-normalization gap described above, not a scaffold-specific issue. `scaffold()` extracts the correct ring system reliably; compare via `mol.apply_aromaticity().canonical_smiles_mode("nostereo")` rather than raw `.smiles` if you need string equality across differently-ordered input.
- **Aromaticity model**: chematic applies Hückel 4n+2 per SSSR ring independently; RDKit uses fused-ring electron delocalization. Visible differences in N-heterocycles (pyridone, quinolone, indolizine). Current benchmark on 4,999-mol ChEMBL subset: HBA/HBD/aromatic ring count **100%**; TPSA **99.7%** (±0.1 Å²); LogP **100%** (±0.01). Aromaticity-flag parity on Kekulized input measured worst-of-10-representations: **96.3%** (`scripts/aromaticity_atom_parity.py`) — bit-for-bit unchanged by the SSSR fix above, confirming the SSSR bug and the aromaticity gap are independent; the aromaticity gap is root-caused separately to an `aromatic_context` bypass mechanism, not yet fixed. Two molecules (azulene, purine) are known to have regressed by the SSSR fix specifically — the old, broken SSSR had accidentally been masking the `aromatic_context` bug for these non-alternant/bridgehead-heavy structures. They are **not present in the 5,000-mol measured corpus at all** (confirmed by direct search — ChEMBL-derived drug-like corpora don't contain bare azulene/purine), so the 96.3% figure is unchanged because it cannot see them, not because they have zero impact; both are pinned as `#[ignore]`d regressions in `chematic-perception`'s test suite with the root cause documented in-code, pending the `aromatic_context` fix.
- **TPSA edge cases**: remaining 0.3% discrepancy (16 of 4,999 molecules) concentrated in exotic phosphazene ring-N calibration and cyclic sulfurimide/S=N=P chemistry — not relevant for drug-like molecules.

---

## Repository Structure

```
chematic/
├── Cargo.toml                    workspace root (v0.4.26)
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
│   ├── chematic-wasm/            130+ WASM exports → npm @kent-tokyo/chematic
│   ├── chematic-py/              PyO3 Python bindings → pip install chematic
│   ├── chematic-ewald/           PME Ewald summation, B-spline interpolation
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
cargo test --workspace --lib --quiet                                      # 2,366 lib tests
cargo test -p chematic-inchi --features native-inchi --test standard_inchi  # +16 InChI tests
cargo clippy --workspace -- -D warnings                                   # lints (zero warnings)
```

---

## Citation

If you use chematic in academic or research work, please cite:

```bibtex
@software{chematic,
  author    = {kent-tokyo},
  title     = {chematic: A pure-Rust cheminformatics toolkit},
  url       = {https://github.com/kent-tokyo/chematic},
  version   = {0.4.26},
  year      = {2026},
}
```

---

## License

Licensed under either of Apache License 2.0 or MIT License, at your option.

---

If chematic saves you time, a [GitHub star](https://github.com/kent-tokyo/chematic) helps others discover it.
