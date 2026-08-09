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
Pure Rust · Zero C/C++ · Python · WebAssembly · [Live Demo](https://kent-tokyo.github.io/chematic/playground/)

| | chematic | RDKit (Python) | RDKit.js (WASM) |
|---|---|---|---|
| **Get started** | `pip install chematic` | conda / cmake required | no Python bindings |
| **Browser bundle** | **719 KB** | not available | ~30 MB (~42× larger) |
| **Batch fingerprints** | **~78 µs/mol** (2–3× faster) | ~160–235 µs/mol | — |
| **Memory safety** | compiler-enforced (Rust) | C++ | C++ |
| **Build from source** | `cargo build` only | cmake + clang + Boost | Emscripten SDK |

All numbers are reproducible — see [benchmark details](https://kent-tokyo.github.io/chematic/benchmark/).  
WASM sizes: chematic **719 KB** · RDKit.js ~30 MB · Indigo WASM ~40 MB

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
| **Browser app** | 719 KB WASM bundle, zero backend required, React/Vue/Svelte ready |
| **Jupyter notebook** | `mol` renders SVG inline; `descriptors_df()` returns a pandas DataFrame |
| **Batch analysis** | Rayon-parallel descriptor/fingerprint/3D pipelines; SDF/CSV in, CSV out |
| **Rust server** | Pure-Rust crates with no C/C++ toolchain; Axum/Actix compatible |

Full worked examples → [Use cases](https://kent-tokyo.github.io/chematic/use-cases/)

---

## When to use chematic

**Use chematic if:**

- You want chemistry in the browser (WASM, 719 KB, no server required)
- You need a pure Rust stack with no C++ toolchain dependencies
- You deploy to environments where `pip install rdkit` is impractical (Cloudflare Workers, Lambda, embedded)
- You build AI agents and want native MCP tool integration
- You process molecules in batch at high throughput (ECFP4: 2–3× faster than RDKit, Rayon-parallel)
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
# chematic v0.12.0
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
[`chematic-mcp` README](crates/chematic-mcp/README.md#protocol-eras) and
[`docs/mcp/2026-07-28-implementation-rfc.md`](docs/mcp/2026-07-28-implementation-rfc.md).
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

chematic's own ~149,000 lines of Rust (tokei-measured code lines, all 18 crates,
2026-08-02) contain **zero `unsafe` blocks outside one file**: 9 `unsafe {}` blocks
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
no `clang`. The npm package `@kent-tokyo/chematic` is **719 KB gzip** — ~42× smaller
than RDKit.js. One codebase runs on Linux, macOS, Windows, and in every browser.

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
| WASM bundle | **719 KB** gzip | — |

\*LogP max Δ = 1.1×10⁻¹³ across 4,999 molecules — within float64 rounding error.  
†Stereocenter count: ~99.96% vs legacy `CalcNumAtomStereoCenters` (a handful of molecules where chematic matches `FindPotentialStereo` and legacy under-counts); ~98.6% vs new-CIP `FindPotentialStereo` (cage/bridgehead molecules where both chematic and legacy correctly return fewer than the new oracle). chematic is calibrated between both extremes. This measures whether an atom is *flagged* as a stereocenter, not whether its R/S label is correct — see the next row.  
‡CIP R/S label agreement measures, for atoms both oracles agree are stereocenters, whether the assigned R/S descriptor matches — a stricter, separate check from stereocenter count agreement above. This row is chematic's *default* `assign_cip()` path. The separate `chematic-cip` engine now reaches 99.38% raw / 99.64% oracle-stable (Milestone 4 gate closed) and is reachable opt-in via `assign_cip_with_mode(mol, CipMode::Accurate)` (Rust), `Mol.cip_stereo(mode="accurate")` (Python), or `cip_assignments_accurate_json` (WASM) — see [`docs/cip_accurate_rfc.md`](docs/cip_accurate_rfc.md). No default path changed; this row's 96.30% is unaffected.

All numbers are reproducible with the scripts in this repo.  
Full history → [benchmarks/](benchmarks/) · Methodology → [validation/](validation/)

---

## Comparison with Other Cheminformatics Libraries

| Feature                 | **chematic**                              | RDKit (rdkit-sys)  | OpenBabel FFI  | RDKit.js (WASM)    |
|-------------------------|-------------------------------------------|--------------------|----------------|--------------------|
| **C/C++ dependencies**  | **None (default)**†                       | Extensive C++      | Extensive C++  | C++ via Emscripten |
| **WASM binary size**    | **~1.9 MB** (719 KB gzip)                 | N/A (no WASM)      | N/A (no WASM)  | ~30 MB             |
| **Build requirement**   | `cargo build` only                        | cmake + clang      | cmake + clang  | Emscripten SDK     |
| **WASM target support** | **Full (native)**                         | No                 | No             | Yes (Emscripten)   |
| **Python bindings**     | **Yes** (`pip install chematic`, PyO3)    | Yes (rdkit-sys)    | Yes            | No                 |
| **Unsafe Rust**         | **None in own crates**‡                   | Extensive          | Extensive      | N/A                |

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
| **MCP server (AI agent API)**                | **Yes — 20 tools incl. Name→SMILES (stdio only)** | No                  | No             | No                |
| IUPAC name generation                        | **Yes (25+ classes)**                            | No                  | No             | Partial           |
| Name → SMILES (PubChem proxy)                | **Yes** (`name_to_smiles` MCP tool)              | No                  | No             | No                |
| Maintenance (2026)                           | Active                                           | Active              | Minimal        | Active            |

</details>

† Default build only. The optional `native-inchi` feature adds a C-compiler dependency for the vendored IUPAC InChI C library (v1.07.5). This is about C/C++ FFI specifically — the `depict` feature below pulls in pure-Rust rendering crates, so it doesn't add a C compiler dependency even though it isn't unsafe-free (see ‡).

‡ chematic's own ~149,000 lines of Rust (tokei-measured): unsafe-free outside `native-inchi`'s 9 FFI blocks (see "Safe" above) — a real, verifiable claim about code chematic wrote, and categorically different from RDKit/OpenBabel's *C++ FFI* unsafe (uncheckable by any compiler at that boundary) even where the raw count is comparable. It is **not** true of the full dependency tree: the optional `depict` feature (SVG/PDF/EPS rendering) pulls in resvg/usvg/rustybuzz/tiny-skia/zune-jpeg, pure-Rust crates that are themselves not unsafe-free — measured directly (`unsafe fn`/`impl`/`trait`/`{` openings): tiny-skia 151, zune-jpeg 79, rustybuzz 14, image 8, fontdb 3, tiny-skia-path 3 (258 total in this set alone). `chematic-py` (`pip install chematic`) and the npm package both depend on `chematic-depict` directly, so this applies to both real-world install paths, not just an edge case.

---

## JavaScript / TypeScript (WebAssembly)

**719 KB gzip — ~42× smaller than RDKit.js.** No Emscripten, no cmake. Drop-in for browser or Node.js.

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

130+ exported functions cover descriptors, fingerprints, 3D geometry, reactions (incl. `retro_disconnect_json` — single-step retrosynthetic disconnection), diversity picking, and SDF round-trips.
See the [full WASM API reference](https://kent-tokyo.github.io/chematic/) for all exports.
---

## Crate Reference

| Crate                 | Description                                                                                              | Tests |
|-----------------------|----------------------------------------------------------------------------------------------------------|-------|
| `chematic-core`       | Atom, Bond, Molecule, Element, kekulization (no deps); mutable `add/remove_atom/bond`, `fragments()`, `is_connected()`, `formula_with_isotopes`, `validate_valence`; `StereoGroup`/`StereoGroupKind` | 71    |
| `chematic-smiles`     | OpenSMILES parser, writer, canonical SMILES; **stereo parity correction** (pre-solves RDKit #8775 — `@`/`@@` auto-flipped on odd permutations); **allene cumulated double bond stereo** (`C=C=C` `@`/`@@`, round-trip stable) | 109    |
| `chematic-perception` | SSSR, Hückel aromaticity + antiaromaticity (4n+2 rule), `apply_aromaticity`, `aromatize`/`kekulize_inplace`, `assign_stereo_from_2d`, `assign_ez_from_2d`, `cip_ez_descriptor`; **zero-order/dative bonds excluded from ring perception** | 101    |
| `chematic-mol`        | MOL/SDF V2000+V3000 (R/W with 2D coords, +partial charge writing), CML (R/W), CDXML (R); `SdfRecord` with coords+props; MDL RXN R/W; V3000 stereo-group COLLECTION R/W; **AutoDock PDBQT** (parse + write); **ChemicalJSON** (`parse_cjson`/`write_cjson`, Avogadro/MolSSI format); **2D wedge/hash tetrahedral parity + E/Z double-bond direction now perceived automatically on read** (`read_mol_with_diagnostics`/`read_mol_v3000_with_diagnostics`, typed opt-in diagnostics) | 130+   |
| `chematic-depict`     | 2D SVG (CPK colors, highlighting, grid), DepictData, `detect_crossings`, `render_svg_with_metadata`, reaction SVG; **PDF output** (`depict_pdf`/`depict_pdf_opts` via svg2pdf); **EPS output** (`depict_eps`/`depict_eps_opts`, pure Rust); `tiny_skia` PNG is optional `png` feature (default on, disabled for WASM) | 64    |
| `chematic-chem`       | 190+ descriptor values (71 functions), tautomers, scaffold, BRICS, QED, standardize, CIP; **pKa prediction** (15 SMARTS rules); **ADMET profile** (BBB/Caco-2/hERG/CYP3A4); **HBA 100% RDKit agreement** (4 999 / 4 999 mol benchmark); **TPSA 100% ±0.1 Å² / LogP 100%\* / HBD 100% / stereocenter count 99.96% (legacy) / 98.6% (new CIP)** vs RDKit (4,999-mol ChEMBL); **CIP R/S label agreement 96.30% (default), 99.64% oracle-stable via opt-in `CipMode::Accurate`** (5,000-mol ChEMBL, see `docs/cip_accurate_rfc.md`); **topological descriptors** (`petitjean_index`, `graph_diameter`, `graph_radius`, `graph_eccentricities`, `eccentric_connectivity_index`, `hosoya_index`, `moran_autocorr`, `geary_autocorr`); **`schultz_mti`, `gutman_mti`, `vabc` (Bondi radii vdW volume), `gravitational_index`**; `clean_stereo_groups()` in standardize | 662   |
| `chematic-fp`         | ECFP2/4/6, FCFP4/6, MACCS, TopoPF, AtomPair, Torsion, Layered, Pattern, Pharmacophore, Reaction, **MAP4** (Minervini 2020, not in RDKit) — Tanimoto/Dice; bulk similarity | 185    |
| `chematic-ff`         | **MMFF94 all 7 terms** (Halgren 1996): Bond/Angle/Torsion/vdW/Elec + **OOP** (117 entries) + **Stretch-Bend** (282 entries); steepest-descent + L-BFGS optimizer, torsion scan, energy breakdown; DREIDING typing; **UFF** (metals/organometallics: Zn, Fe, Cu, …) | 98    |
| `chematic-smarts`     | SMARTS, VF2, MCS with chirality matching; **SmartsCache** (LRU compilation cache, 5–20×); **named_pattern()** library (20 functional group patterns); **atom map `:N` in SMARTS** (`[O;D1;H0:3]` — stored as metadata, not a match criterion); **`[kN]` ring-size primitive**; **VF2 early-exit** when query > target atom count; **`find_matches_with_rings`** — share SSSR across multi-pattern batches | 142   |
| `chematic-3d`         | 3D coordinate generation, distance geometry constraints, ETKDG KB (40 torsion patterns, adaptive noise), force-field minimization, shape descriptors, ConformerEnsemble with RMSD pruning, PDB/XYZ; **GETAWAY HATS-matrix** (full 19-dim implementation); **`whim_getaway_combined()`** now 29-dim | 265    |
| `chematic-rxn`        | Reaction SMILES/SMIRKS, `run_reactants`/`run_reactants_strict`; **`retro_disconnect()`** — 60 retro-SMIRKS templates (AmideBond/Ester/Ether/CNBond/CCBond/CSBond) + SA Score ranking; **parity-aware `@`/`@@` SMIRKS stereo filtering**; **E/Z double-bond stereo filtering** in `run_reactants` (`ez_stereo_outward`, `smirks_ez_stereo_ok`) | 137    |
| `chematic-inchi`      | InChI/InChIKey: pure-Rust approximation (WASM) **+ IUPAC-standard** via `native-inchi` feature (vendored C lib 1.07.5, bit-exact); **parse_inchi** reader; **verified canonical-SMILES dedup** (`dedup::{group_candidates, deduplicate_verified}`, fail-closed on legacy-CIP-unresolved specified tetrahedral stereo); **accurate-CIP dedup preflight** (issue #161) recovering verified-comparison capability on legacy-CIP-unresolved stereocentres; **indexed graph relation API** (`compare_indexed_graph_relation`, orthogonal `GraphStrictness`/`AtomMapPolicy` axes) | 108 (+16*)    |
| `chematic-cip`        | Opt-in accurate CIP engine (`assign_cip_accurate_experimental`, hierarchical digraph, Rules 1a/1b/2/4b/5, RDKit-compatible MANCUDE fractional atomic numbers) — the default `assign_cip()`/`CipMode::LegacyFast` is unchanged | —     |
| `chematic-wasm`       | **131+ WASM exports** — npm: `@kent-tokyo/chematic` (published in lockstep with crates.io/PyPI); pKa/ADMET/BBB/Caco-2/hERG/CYP3A4; `smiles_to_pdbqt`, `minimize_uff_json`, **`retro_disconnect_json`** (issue #91) | 223   |
| `chematic-iupac`      | Local IUPAC name generation — **25+ compound classes**: alkanes, cycloalkanes, alkenes/alkynes, alcohols, amines, halides, aldehydes, ketones, acids, esters, amides, **piperidine, morpholine, piperazine, naphthalene, sulfides** | 47    |
| `chematic-mcp`        | **MCP (Model Context Protocol) server** — AI agent integration; **20 tools**: parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, name_to_smiles, retrosynthesis, smiles_to_moljson, moljson_to_smiles, representation_router, **molecule_context_pack**; dual-era protocol (legacy `2024-11-05` + modern `2026-07-28` stateless dialect), `structuredContent`/`outputSchema` on all 20 tools | 82    |
| `chematic-py`         | PyO3 Python bindings (`pip install chematic`); 300+ API endpoints: `from_smiles()`, `Mol.descriptors()`, `Mol.minimize_dreiding()`, `from_cxsmiles()`, `from_rxn_file()`/`to_rxn_file()`, `parse_sdf_with_coords()`, `Mol.ring_families()`, `tanimoto_matrix()`, `iter_sdf()`, `SimilarityIndex`; **`mol.to_pdf()`/`mol.to_eps()`** (depict); **`from_cjson()`/`mol.to_cjson()`** (ChemicalJSON); **`mol.schultz_mti`, `mol.gutman_mti`, `mol.vabc`, `mol.gravitational_index`**; **`bulk.substructure_match(smarts, mols)`** (parallel VF2 on pre-parsed Mol objects); **`mol.describe()`** (LLM/MCP-ready natural-language summary); **`mol.diff(other)`** (element + descriptor diff); Sprint 18–27 coverage | 300+  |
| `chematic-ewald`      | PME Ewald summation, B-spline interpolation (cubic, phase-corrected)                                     | 16    |
| `chematic`            | Umbrella crate with feature flags (all sub-crates, incl. `iupac`, `inchi`)                              | 1     |

```
cargo test --workspace --lib --quiet                                          # 3,235 tests, all passing (2026-08-02)
cargo test -p chematic-inchi --features native-inchi --test standard_inchi  # +16 IUPAC-exact InChI tests
```

---

## Recent Development

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
| TPSA | **Stable** | 100% on 175-mol drug-like set; **99.7%** on 4,999-mol ChEMBL subset (±0.1 Å²) |
| LogP (Crippen) | **Stable** | **100%** on 4,999-mol corpus (±0.01); ~99% on 175-mol drug-like set (±0.3) |
| ECFP4 / MACCS fingerprints | **Stable** | RDKit comparison + benchmark |
| Tanimoto similarity | **Stable** | RDKit comparison |
| SDF / MOL V2000/V3000 I/O | **Stable** | round-trip tests |
| Substructure search (SMARTS / VF2) | **Stable** | internal test suite |
| PAINS / Brenk filters | **Stable** | rule matching stable; ring-size SMARTS (`[r5]`/`[r6]`) now 0% instability across 5,000-mol worst-of-10 (was ~29–55% before the SSSR fix) |
| Ring perception (SSSR) | **Stable** | Horton algorithm, minimal + deterministic; 0% self-instability across 5,000-mol worst-of-10 (was 50.6%) — see Known Limitations below |
| Murcko scaffold | **Stable** (normalized) | normalized string output **100%** stable across 5,000-mol worst-of-10 (was 0.8% unstable, same root cause as the canonical-SMILES corruption above, now fixed); raw `.smiles` inherits the still-partially-open direction-normalization gap — normalize before comparing (see Known Limitations) |
| 2D SVG depiction | **Stable** | visual spot-checks; not publication-quality |
| 3D conformer (DG + MMFF94) | **Experimental** | reasonable geometry; not equivalent to RDKit ETKDGv3 quality |
| pKa prediction | **Rule-based screening** | 15 SMARTS rules; early triage only, not clinical |
| ADMET (BBB / Caco-2 / hERG / CYP3A4) | **Rule-based screening** | empirical models; directional, not validated on clinical endpoints |
| IUPAC name generation | **Partial** | common compound classes; complex structures may fail |
| Pure-Rust InChI | **Approximate** | enable `native-inchi` feature for bit-exact IUPAC InChI |

Full benchmark methodology → [validation/](validation/) · History → [benchmarks/](benchmarks/)

---

## Known Limitations

- **`canonical_smiles()` is now partially normalized for E/Z stereochemistry — still not safe as a dedup or cache key.** Isolated/simple E/Z double bonds have two equally correct `/`/`\` spellings (e.g. `/N=N/` vs `\N=N\`); the writer previously never normalized between them. Fixed for the general case: every connected E/Z system (a double bond plus every directional bond geometrically tied to it, including whole conjugated chains) is now normalized so its first directional bond in canonical write order is always `/`, regardless of input spelling. Measured on the 5,000-mol ChEMBL corpus, worst-of-10: E/Z-only self-instability (tetrahedral stereo stripped) improved **9.76% → 5.50%** (275/5000 still unstable); structural correctness unaffected by this change (re-verified **0/5000** ChEMBL and **0/33** acyclic-polyene corpus). The residual 275 are confirmed **100% cosmetic** — every unstable case's variants represent the same molecule per RDKit, zero corruption — but the cause is a **mixed pool, not fully root-caused**: about half match a specific motif (a small ring bearing two or more exocyclic double bonds, e.g. cross-conjugated cyclic diimines) where which physical bonds count as "one system" is not yet input-spelling-invariant; the other half is uncharacterized. Until this closes fully, ~1 in 18 stereo-bearing molecules (down from ~1 in 10) can still produce two different, individually valid `canonical_smiles()` strings for the same molecule — do not use it as a dedup or cache key today; document your own dedup key as `apply_aromaticity()`-normalized in the meantime if this matters for your use case.
- **Canonical SMILES structural corruption — fixed.** Before this fix, `canonical_smiles(parse(x))` could silently emit a *different stereoisomer* (not just a differently-spelled but equivalent string) depending on `x`'s input traversal order. Measured on a 5,000-mol ChEMBL subset, worst-of-10 independently-traversed representations per molecule, RDKit-verified structural correctness: **4.28% (214/5000)** of molecules had at least one variant round-trip to the wrong molecule. Root-caused to two independent parser bugs (not the originally-suspected "conjugated double-bond markers are geometrically coupled across bonds" — that diagnosis was disproven, see below), each confirmed via a real found molecule and a minimal regression test: (1) a ring-closure directional-bond (`/`/`\`) marker read at the *closing* occurrence of a ring digit was stored raw instead of flipped to the opening→closing sense, corrupting a conjugated E/Z chain whenever its connecting bond happened to be routed through a ring closure; (2) a stereocenter that opens a ring whose partner closes *inside its own branch* had its neighbor-order resolution keyed by the reusable ring *digit* rather than a unique per-occurrence id, so a later, unrelated reuse of the same digit elsewhere in the SMILES could silently steal and corrupt the stereocenter's neighbor order. **After both fixes: structural correctness is 100% (0/5000) on ChEMBL**, confirmed three times over via independently-ordered reconstructions of the fix (with and without an unrelated third ranking fix, to rule out a hidden dependency). Because both root causes are ring-closure-specific — and retinoids, carotenoids, prostaglandins, leukotrienes, and polyene macrolides carry their long conjugated systems in *acyclic* chains, essentially absent from ChEMBL-random sampling — this was independently re-verified on a dedicated 33-compound corpus of exactly those classes (tretinoin, β-carotene, lycopene, amphotericin B, leukotriene B4, and 28 others; `scripts/polyene_corpus.csv`): **0/33 (0.00%) at worst-of-30**, with a positive control confirming 12/33 (36.36%) corruption on the pre-fix code for this same corpus (all 12 failures were ring-closure-heavy structures; zero purely-acyclic examples — including fully acyclic lycopene — ever failed, even unpatched). This directly disproves the original "any conjugated chain" diagnosis and closes the investigation with no remaining corruption class identified. Skeleton-only and tetrahedral-only self-*stability* also reached 0% (were 0.16% and 4.36%); raw combined self-stability (all stereo intact) improved 86.02% → 90.28% (13.98% → 9.72% unstable) — the entire remainder is the separate, non-corrupting direction-normalization gap described above, not residual corruption. Round-trip invariance (`canonical(parse(canonical(m))) == canonical(m)`) improved slightly, 98.26% → 98.32%, since it was never measuring the corruption class directly.
- **Ring perception (SSSR) was non-deterministic and non-minimal — fixed.** The old `find_sssr` built a single spanning tree and took one fundamental cycle per non-tree edge, with no redundancy to recover a smaller ring when the tree's shape made one unnecessarily large (naphthalene, `c1ccc2ccccc2c1`, deterministically returned ring sizes `[6, 10]` instead of `[6, 6]`). `find_sssr` now uses Horton's algorithm (candidate cycles from every vertex × every edge via shortest-path trees, O(V·E) candidates, canonical-rank tie-break for determinism), giving a genuinely minimum-weight, deterministic basis. Measured on a 5,000-mol ChEMBL subset, worst-of-10 independently-traversed representations per molecule: self-stability **100%** (was 50.6%); single-parse ring-size agreement with RDKit **98.9%** (was 72.4%) — the residual ~1.1% gap is RDKit's own `GetSymmSSSR` legitimately returning *more* rings than the topological minimum for symmetric fused systems (e.g. cubane: μ=5, RDKit=6), not a chematic bug; full symmetrization (Vismara relevant cycles) is future work, not required for correctness. Downstream wins, same corpus: ring-size SMARTS `[r5]`/`[r6]` **0%** instability (was 29–55%), `NumAromaticRings` **0%** (was ~4%), `RingCount`/MW/TPSA/HBA/HBD/LogP/MR unaffected (were already 0%). Two known-narrow exceptions where the *old* SSSR bug had accidentally compensated for a separate, still-open aromaticity bug — see the Aromaticity model bullet below. Full methodology: `scripts/ringinfo_parity.py`.
- **Murcko scaffold: ring topology and normalized string output are now fully stable.** The previously-reported "100% traversal-order instability" was itself a measurement-harness bug (comparing `Mol` objects by Python identity instead of value — always reported "unstable" regardless of the real result); that script bug is fixed (`scripts/ring_collateral_damage.py`). Re-measured on a 5,000-mol worst-of-10 run after the canonical-SMILES corruption fixes above: after normalizing (`apply_aromaticity().canonical_smiles_mode("nostereo")`), self-stability is **100% (0/5000 unstable)**, down from a 0.8% residual — confirming that residual was the same canonical-SMILES structural corruption, not a Murcko ring-selection bug, and it is now fully resolved. Raw isomeric `scaffold().smiles` string comparison (no normalization) is **79.30%** stable (20.70% unstable, was ~45%, essentially unchanged by the partial E/Z-normalization fix above — scaffolds strip most of the side-chain motifs that fix improves) — the remainder is the still-partially-open, non-corrupting `/`/`\` direction-normalization gap described above, not a scaffold-specific issue. `scaffold()` extracts the correct ring system reliably; compare via `mol.apply_aromaticity().canonical_smiles_mode("nostereo")` rather than raw `.smiles` if you need string equality across differently-ordered input.
- **Aromaticity model**: chematic applies Hückel 4n+2 per SSSR ring independently; RDKit uses fused-ring electron delocalization. Visible differences in N-heterocycles (pyridone, quinolone, indolizine). Current benchmark on 4,999-mol ChEMBL subset: HBA/HBD/aromatic ring count **100%**; TPSA **99.7%** (±0.1 Å²); LogP **100%** (±0.01). Aromaticity-flag parity on Kekulized input measured worst-of-10-representations: **96.3%** (`scripts/aromaticity_atom_parity.py`) — bit-for-bit unchanged by the SSSR fix above, confirming the SSSR bug and the aromaticity gap are independent; the aromaticity gap is root-caused separately to an `aromatic_context` bypass mechanism, not yet fixed. Two molecules (azulene, purine) are known to have regressed by the SSSR fix specifically — the old, broken SSSR had accidentally been masking the `aromatic_context` bug for these non-alternant/bridgehead-heavy structures. They are **not present in the 5,000-mol measured corpus at all** (confirmed by direct search — ChEMBL-derived drug-like corpora don't contain bare azulene/purine), so the 96.3% figure is unchanged because it cannot see them, not because they have zero impact; both are pinned as `#[ignore]`d regressions in `chematic-perception`'s test suite with the root cause documented in-code, pending the `aromatic_context` fix.
- **TPSA edge cases**: remaining 0.3% discrepancy (16 of 4,999 molecules) concentrated in exotic phosphazene ring-N calibration and cyclic sulfurimide/S=N=P chemistry — not relevant for drug-like molecules.

---

## Repository Structure

```
chematic/
├── Cargo.toml                    workspace root (v0.12.0)
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
cargo test --workspace --lib --quiet                                      # 3,235 lib tests
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
  version   = {0.12.0},
  year      = {2026},
}
```

---

## License

Licensed under either of Apache License 2.0 or MIT License, at your option.

---

If chematic saves you time, a [GitHub star](https://github.com/kent-tokyo/chematic) helps others discover it.
