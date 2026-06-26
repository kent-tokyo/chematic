# chematic

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/chematic.svg)](https://crates.io/crates/chematic)
[![PyPI](https://img.shields.io/pypi/v/chematic.svg)](https://pypi.org/project/chematic/)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic.svg)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](LICENSE-MIT)
[![Docs](https://img.shields.io/badge/docs-site-blue)](https://kent-tokyo.github.io/chematic/)
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

### Diagnostics

```python
import chematic
chematic.doctor()
# chematic v0.4.22
# Python 3.12.x  |  darwin arm64
#
# Descriptor accuracy (benchmark 2026-06, v0.4.22 vs RDKit 2026.03.3):
#   MW / HBA / HBD / ARC  100%   (4,999-mol ChEMBL subset)
#   TPSA                  93.3%
#   LogP (Crippen)        ~99%
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
| ECFP4 throughput | **3.6 µs/mol** (5–14× vs RDKit) | 5,000 mol |
| HBA / HBD / aromatic ring count | **100% RDKit agreement** | 4,999 mol |
| TPSA | **100%** within ±0.1 Å² | 175-mol drug-like set |
| TPSA | 93.3% within ±0.1 Å² | 4,999-mol ChEMBL subset |
| WASM bundle | **504 KB** gzip | — |

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
| **Unsafe Rust**         | **None**                                  | Extensive          | Extensive      | N/A                |

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
| Molecular descriptors                        | **190+ (incl. MQN×42, BCUT2D, ESOL, LogD, XLogP3, BOILED-Egg)** | ~30  | ~20            | ~30               |
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

† Default build only. The optional `native-inchi` feature adds a C-compiler dependency for the vendored IUPAC InChI C library (v1.07.5). All other crates remain FFI-free.

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
| `chematic-core`       | Atom, Bond, Molecule, Element, kekulization (no deps); mutable `add/remove_atom/bond`, `fragments()`, `is_connected()`, `formula_with_isotopes`, `validate_valence`; `StereoGroup`/`StereoGroupKind` | 69    |
| `chematic-smiles`     | OpenSMILES parser, writer, canonical SMILES; **stereo parity correction** (pre-solves RDKit #8775 — `@`/`@@` auto-flipped on odd permutations); **allene cumulated double bond stereo** (`C=C=C` `@`/`@@`, round-trip stable) | 48    |
| `chematic-perception` | SSSR, Hückel aromaticity + antiaromaticity (4n+2 rule), `apply_aromaticity`, `aromatize`/`kekulize_inplace`, `assign_stereo_from_2d`, `assign_ez_from_2d`, `cip_ez_descriptor`; **zero-order/dative bonds excluded from ring perception** | 34    |
| `chematic-mol`        | MOL/SDF V2000+V3000 (R/W with 2D coords, +partial charge writing), CML (R/W), CDXML (R); `SdfRecord` with coords+props; MDL RXN R/W; V3000 stereo-group COLLECTION R/W; **AutoDock PDBQT** (parse + write); **ChemicalJSON** (`parse_cjson`/`write_cjson`, Avogadro/MolSSI format) | 31    |
| `chematic-depict`     | 2D SVG (CPK colors, highlighting, grid), DepictData, `detect_crossings`, `render_svg_with_metadata`, reaction SVG; **PDF output** (`depict_pdf`/`depict_pdf_opts` via svg2pdf); **EPS output** (`depict_eps`/`depict_eps_opts`, pure Rust); `tiny_skia` PNG is optional `png` feature (default on, disabled for WASM) | 28    |
| `chematic-chem`       | 190+ descriptors, tautomers, scaffold, BRICS, QED, standardize, CIP; **pKa prediction** (15 SMARTS rules); **ADMET profile** (BBB/Caco-2/hERG/CYP3A4); **HBA 100% RDKit agreement** (4 999 / 4 999 mol benchmark); **TPSA ±0.1 Å² / LogP ±0.3 / HBD 100%** vs RDKit (175-mol bulk regression); **topological descriptors** (`petitjean_index`, `graph_diameter`, `graph_radius`, `graph_eccentricities`, `eccentric_connectivity_index`, `hosoya_index`, `moran_autocorr`, `geary_autocorr`); **`schultz_mti`, `gutman_mti`, `vabc` (Bondi radii vdW volume), `gravitational_index`**; `clean_stereo_groups()` in standardize | 211   |
| `chematic-fp`         | ECFP2/4/6, FCFP4/6, MACCS, TopoPF, AtomPair, Torsion, Layered, Pattern, Pharmacophore, Reaction, **MAP4** (Minervini 2020, not in RDKit) — Tanimoto/Dice; bulk similarity | 87    |
| `chematic-ff`         | **MMFF94 all 7 terms** (Halgren 1996): Bond/Angle/Torsion/vdW/Elec + **OOP** (117 entries) + **Stretch-Bend** (282 entries); steepest-descent + L-BFGS optimizer, torsion scan, energy breakdown; DREIDING typing; **UFF** (metals/organometallics: Zn, Fe, Cu, …) | 51    |
| `chematic-smarts`     | SMARTS, VF2, MCS with chirality matching; **SmartsCache** (LRU compilation cache, 5–20×); **named_pattern()** library (20 functional group patterns); **atom map `:N` in SMARTS** (`[O;D1;H0:3]` — stored as metadata, not a match criterion); **`[kN]` ring-size primitive**; **VF2 early-exit** when query > target atom count; **`find_matches_with_rings`** — share SSSR across multi-pattern batches | 142   |
| `chematic-3d`         | 3D coordinate generation, distance geometry constraints, ETKDG KB (40 torsion patterns, adaptive noise), force-field minimization, shape descriptors, ConformerEnsemble with RMSD pruning, PDB/XYZ; **GETAWAY HATS-matrix** (full 19-dim implementation); **`whim_getaway_combined()`** now 29-dim | 45    |
| `chematic-rxn`        | Reaction SMILES/SMIRKS, `run_reactants`/`run_reactants_strict`; **`retro_disconnect()`** — 60 retro-SMIRKS templates (AmideBond/Ester/Ether/CNBond/CCBond/CSBond) + SA Score ranking; **parity-aware `@`/`@@` SMIRKS stereo filtering**; **E/Z double-bond stereo filtering** in `run_reactants` (`ez_stereo_outward`, `smirks_ez_stereo_ok`) | 25    |
| `chematic-inchi`      | InChI/InChIKey: pure-Rust approximation (WASM) **+ IUPAC-standard** via `native-inchi` feature (vendored C lib 1.07.5, bit-exact); **parse_inchi** reader | 28 (+16*)    |
| `chematic-wasm`       | **130+ WASM exports** — npm: `@kent-tokyo/chematic` v0.4.18 (~500 KB, 504 KB gzip); pKa/ADMET/BBB/Caco-2/hERG/CYP3A4; `smiles_to_pdbqt`, `minimize_uff_json` | 209   |
| `chematic-iupac`      | Local IUPAC name generation — **25+ compound classes**: alkanes, cycloalkanes, alkenes/alkynes, alcohols, amines, halides, aldehydes, ketones, acids, esters, amides, **piperidine, morpholine, piperazine, naphthalene, sulfides** | 45    |
| `chematic-mcp`        | **MCP (Model Context Protocol) server** — AI agent integration; **15 tools**: parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, **name_to_smiles** | 28    |
| `chematic-py`         | PyO3 Python bindings (`pip install chematic`); 300+ API endpoints: `from_smiles()`, `Mol.descriptors()`, `Mol.minimize_dreiding()`, `from_cxsmiles()`, `from_rxn_file()`/`to_rxn_file()`, `parse_sdf_with_coords()`, `Mol.ring_families()`, `tanimoto_matrix()`, `iter_sdf()`, `SimilarityIndex`; **`mol.to_pdf()`/`mol.to_eps()`** (depict); **`from_cjson()`/`mol.to_cjson()`** (ChemicalJSON); **`mol.schultz_mti`, `mol.gutman_mti`, `mol.vabc`, `mol.gravitational_index`**; **`bulk.substructure_match(smarts, mols)`** (parallel VF2 on pre-parsed Mol objects); **`mol.describe()`** (LLM/MCP-ready natural-language summary); **`mol.diff(other)`** (element + descriptor diff); Sprint 18–27 coverage | 300+  |
| `chematic-ewald`      | PME Ewald summation, B-spline interpolation (cubic, phase-corrected)                                     | 12    |
| `chematic`            | Umbrella crate with feature flags (all sub-crates, incl. `iupac`, `inchi`)                              | 1     |

```
cargo test --workspace --lib --quiet                                          # 211 tests, all passing
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
- Docs: benchmark page added (ECFP4 5–14× vs RDKit, 100% descriptor accuracy on 5 000-mol corpus)

**v0.4.16–v0.4.17** (2026-06-22–23): **SSSR sharing performance sprint**
- `chematic-smarts`: `find_matches_with_rings()` — share a pre-computed `RingSet` across all patterns in a batch
- `chematic-chem`: Crippen 117 SSSR → 1 per `logp_crippen` call; PAINS ~480 → 1; QED 113 → 1; pKa 42 → 1; new `logp_and_mr()`, `logd_from_logp()`, `pka_both()` to avoid redundant passes
- `chematic-fp`: MHFP incremental BFS — 3N → N BFS operations per molecule at radius=2

**v0.4.15** (2026-06-21): **TPSA calibration + E/Z stereo in reactions**
- `chematic-chem`: TPSA ±0.1 Å² calibration sprint — **HBA 100%, HBD 100%, aromatic ring count 100%** on 5 000-mol corpus; TPSA 86.7% → 93.3% (5 000-mol), 100% on 175-mol drug-like set
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
- HBA rewrite: 99.98% agreement with RDKit (5,000 molecule benchmark)


Full changelog: [CHANGELOG.md](CHANGELOG.md)

---

## Built with chematic

Using chematic in a project? [Share it in Discussions](https://github.com/kent-tokyo/chematic/discussions) or open a PR to add it here.

---

## Reliability by Feature

Not all features have the same validation depth. This table tells you what to trust.

| Feature | Status | Validation |
|---|---|---|
| SMILES parse / write | **Stable** | 5,000-mol RDKit comparison; OpenSMILES corpus |
| MW / HBA / HBD | **Stable** | 100% RDKit agreement on 4,999 mol |
| TPSA | **Stable** | 100% on 175-mol drug-like set; 93.3% on 4,999-mol ChEMBL subset |
| LogP (Crippen) | **Stable** | ~99% on 175-mol drug-like set (±0.3) |
| ECFP4 / MACCS fingerprints | **Stable** | RDKit comparison + benchmark |
| Tanimoto similarity | **Stable** | RDKit comparison |
| SDF / MOL V2000/V3000 I/O | **Stable** | round-trip tests |
| Substructure search (SMARTS / VF2) | **Stable** | internal test suite |
| PAINS / Brenk filters | **Stable** | rule-based; matches public SMARTS databases |
| 2D SVG depiction | **Stable** | visual spot-checks; not publication-quality |
| 3D conformer (DG + MMFF94) | **Experimental** | reasonable geometry; not equivalent to RDKit ETKDGv3 quality |
| pKa prediction | **Rule-based screening** | 15 SMARTS rules; early triage only, not clinical |
| ADMET (BBB / Caco-2 / hERG / CYP3A4) | **Rule-based screening** | empirical models; directional, not validated on clinical endpoints |
| IUPAC name generation | **Partial** | common compound classes; complex structures may fail |
| Pure-Rust InChI | **Approximate** | enable `native-inchi` feature for bit-exact IUPAC InChI |

Full benchmark methodology → [validation/](validation/) · History → [benchmarks/](benchmarks/)

---

## Known Limitations

- **Aromaticity model**: chematic applies Hückel 4n+2 per SSSR ring independently; RDKit uses fused-ring electron delocalization. Visible differences in N-heterocycles (pyridone, quinolone, indolizine). Current benchmark on 5,000-molecule corpus: HBA/HBD/aromatic ring count **100%**; TPSA **98.1%** (±0.1 Å²).
- **TPSA edge cases**: remaining 1.9% discrepancy concentrated in azide groups, macrolide lactones, and exotic phosphazene chemistry.

---

## Repository Structure

```
chematic/
├── Cargo.toml                    workspace root (v0.4.22)
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
│   ├── chematic-mcp/             MCP server — 15 AI-callable tools (JSON-RPC 2.0 over stdio)
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
cargo test --workspace --lib --quiet                                      # 211 lib tests
cargo test -p chematic-inchi --features native-inchi --test standard_inchi  # +16 InChI tests
cargo clippy --workspace -- -D warnings                                   # lints (zero warnings)
```

---

## Citation

If you use chematic in academic or research work, please cite:

```bibtex
@software{chematic,
  author    = {Tanabe, Kent},
  title     = {chematic: A pure-Rust cheminformatics toolkit},
  url       = {https://github.com/kent-tokyo/chematic},
  version   = {0.4.21},
  year      = {2026},
}
```

---

## License

Licensed under either of Apache License 2.0 or MIT License, at your option.

---

If chematic saves you time, a [GitHub star](https://github.com/kent-tokyo/chematic) helps others discover it.
