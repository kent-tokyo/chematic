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
# chematic v0.10.0
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

chematic's own ~15,000 lines of Rust contain **~6 `unsafe` blocks**, all confined to the
optional `native-inchi` FFI layer (below). No C++ heap corruptions. No segfaults from
malformed SMILES input. No platform-specific build failures from `-sys` crates. The
compiler enforces memory safety at every call site chematic itself wrote.

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

‡ chematic's own ~15,000 lines of Rust: unsafe-free outside `native-inchi`'s ~6 FFI blocks (see "Safe" above) — a real, verifiable claim about code chematic wrote, and categorically different from RDKit/OpenBabel's *C++ FFI* unsafe (uncheckable by any compiler at that boundary) even where the raw count is comparable. It is **not** true of the full dependency tree: the optional `depict` feature (SVG/PDF/EPS rendering) pulls in resvg/usvg/rustybuzz/tiny-skia/zune-jpeg, pure-Rust crates that are themselves not unsafe-free — measured directly (`unsafe fn`/`impl`/`trait`/`{` openings): tiny-skia 151, zune-jpeg 79, rustybuzz 14, image 8, fontdb 3, tiny-skia-path 3 (258 total in this set alone). `chematic-py` (`pip install chematic`) and the npm package both depend on `chematic-depict` directly, so this applies to both real-world install paths, not just an edge case.

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

131+ exported functions cover descriptors, fingerprints, 3D geometry, reactions (incl. `retro_disconnect_json` — single-step retrosynthetic disconnection), diversity picking, and SDF round-trips.
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
cargo test --workspace --lib --quiet                                          # 3,223 tests, all passing (2026-08-01)
cargo test -p chematic-inchi --features native-inchi --test standard_inchi  # +16 IUPAC-exact InChI tests
```

---

## Recent Development

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
- `chematic-smiles`/`chematic-core`: two representations of the same molecule that differ only in whether an atom's H count came from bracket notation (`[Cl]`) or organic-subset notation (`Cl`) — when the explicit value merely repeats what valence inference gives anyway — now canonicalize identically; `initial_invariant`'s Morgan-rank seed and the canonical writer's bracket-necessity check previously read the raw `hydrogen_count` field directly instead of the crate's own `implicit_hcount()` unifier. Isotope, charge, real stereochemistry, and any genuinely disambiguating explicit H count remain fully distinguishing — only the redundant case is unified. Found via [kent-tokyo/renkin PR #65](https://github.com/kent-tokyo/renkin/pull/65) (a `run_reactants` product from a bracket-notation SMIRKS template diverging from a directly user-typed SMILES of the identical compound). Some canonical SMILES output strings for existing inputs will change as a direct, intended consequence — see `CHANGELOG.md`'s `[0.8.1]` Migration notes if you depend on exact string stability across versions
- Full details in `CHANGELOG.md`'s `[0.8.1]` section

**v0.8.0** (2026-07-29): **Opt-in fail-closed 3D embedding pipeline, canonical-SMILES automorphism-orbit pruning**
- `chematic-3d`: new opt-in `pipeline_v2::embed_pipeline_v2` integrates stochastic distance geometry, macrocycle 1-4 distance-bound adjustments, validated ETKDG-style torsion knowledge, chiral-volume/E-Z stereo verification and repair, and typed force-field minimization into one 12-stage pipeline, with a fail-closed stereo re-check *after* force-field minimization (measured: UFF/DREIDING can reintroduce a stereo violation repair had just fixed). `RingTorsionApplicationPolicy::FailClosed` (default) rejects mechanically applying ring/macrocycle torsion potentials the current optimizer can only score, not apply — `DiagnosticOnly` is an explicit, non-default opt-in. Existing default behavior (Rust/Python/WASM/MCP) is unchanged
- `chematic-smiles`: fixed the `canonical_smiles()` performance regression reported by RENKIN — canonicalization search now prunes branches via exact, independently-verified colored-graph automorphism checks, never rank/hash equality alone. Measured (audited, conservative figures): approximately 5x geomean speedup on high-symmetry molecules (search-leaf count 6625 → 13), approximately 1.1–1.2x geomean speedup on low-symmetry negative controls with no regression, approximately 1.09–1.10x per-molecule geomean speedup on a 5,000-molecule external corpus (approximately 5x total-elapsed improvement there, tail-dominated by a few high-symmetry molecules, not representative of typical per-molecule cost), 0 correctness mismatches. `canonical_smiles(&Molecule) -> String`'s public signature is unchanged; new fallible `canonical_smiles_with_limits` added alongside it. Also fixes a Dative-bond round-trip bug (issue #194): the writer was silently truncating `->` to a plain single bond, and — a deeper, separate root cause — the parser could not distinguish `->` from `<-` at all, silently collapsing both to the same donor/acceptor assignment; both are now correctly direction-aware
- Full details, benchmark numbers, and known limitations in `CHANGELOG.md`

**v0.7.1** (2026-07-27): **`run_reactants`/`canonical_smiles()` performance fix**
- `chematic-smiles`: fixed `canonical_smiles()` calling `CanonicalWriter::write_all()` a second, fully redundant time on every call (the winner was already written once while resolving individualize-refine ties, then thrown away and rewritten) — a real regression introduced by `c219ee7`, first shipped in 0.4.30 (`be5dbb1`, previously misattributed as the origin, is a correctness fix for combinatorial enumeration and is not the commit that introduced the double write — corrected via independent audit), most visible on highly symmetric molecules (plain rings, cages, `CF3`/`tBu`-style substituents: 45-48x slower `canonical_smiles()` in isolation on chematic 0.4.30 vs 0.4.25, reported via an external consumer's `run_reactants`/`apply_retro` regression). Pure refactor — output is byte-identical (verified against all existing tests, including the issue #50 E/Z regression suite). New `chematic-rxn` `perf-instrumentation` feature + `reaction_transform_perf_report` example benchmark added alongside. Full bisect and remaining known gap (genuine automorphism-orbit pruning, not attempted) in `docs/reaction_transform_perf.md`
- Full details in `CHANGELOG.md`'s `[0.7.1]` section

**v0.7.0** (2026-07-26): **2D wedge/hash + E/Z stereo now perceived automatically from MOL/SDF, verified canonical-SMILES dedup, CIP Rule-5 phosphorus fix, native InChI explicit-H/isotope fix**
- `chematic-mol`/`chematic-perception`: MOL V2000/V3000/SDF readers now perceive both tetrahedral wedge/hash parity (PR #154) and E/Z double-bond direction (PR #162) automatically on read, CIP-independent, via `read_mol_with_diagnostics`/`read_mol_v3000_with_diagnostics` — typed opt-in diagnostics (`StereoDiagnostic`/`EzDirectionDiagnostic`) never guess on malformed/ambiguous input. Broad-corpus validation (4,999 molecules vs. RDKit 2026.03.3): E/Z — 622 RDKit-resolved double bonds, 0 semantic inversions, 0 false positives. Also fixed a V2000 MDL-code-4 bug, two V3000 `CFG` bugs, and a V2000 writer bug found while building this. PR #162's own corpus reconciliation surfaced a **new, not-yet-fixed** gap: `canonical_smiles()` can drop an already-correctly-`write()`-encoded E/Z marker for some molecules (aromaticity-unaware carrier grouping) — not yet filed as an issue
- `chematic-inchi`: new verified canonical-SMILES dedup (`dedup` module) — fast canonical-SMILES candidate bucketing reconciled against native-InChI verified identity; fails closed (`VerificationUnavailable`) rather than risk a false merge whenever a specified stereocenter's legacy-CIP ranking is unresolved (closes a real false-`VerifiedDuplicate` found via live 5,000-molecule corpus verification). Follow-up tracked as [#161](https://github.com/kent-tokyo/chematic/issues/161) (an accurate-CIP preflight could recover most of the conservative cases)
- `chematic-cip`: CIP Rule-5 pseudoasymmetric r/s phosphorus fix
- `chematic-inchi`: native InChI (`native-inchi` feature) explicit-hydrogen/isotope conversion fix
- Full details in `CHANGELOG.md`'s `[Unreleased]` section

**v0.6.0** (2026-07-25): **RDKit-bit-exact ECFP4 cross-language stable API, canonical SMILES E/Z marker consistency, opt-in authoritative aromaticity demotion**
- `chematic-fp`/Python/WASM: promoted the RDKit-bit-exact Morgan/ECFP4 path to a documented, cross-language opt-in API (`Mol.rdkit_ecfp4()` in Python, `rdkit_ecfp4_bitvec()` in WASM), generalized to an independently oracle-verified `(radius, fpSize)` matrix (4×5 = 20 cells, closed enums so an unsupported combination can't be constructed). Rust/Python/WASM verified byte-identical by actually building and running all three, not "identical by construction." `ecfp4()`'s default behavior is unchanged
- `chematic-smiles`: fixed canonical output placing an E/Z direction marker on a different substituent bond depending on input atom order — permutation invariance for marker placement **93.0% → 98.1%** (264/282 known-divergent molecules now converge; 18 residual cases deliberately left unresolved to avoid corrupting a shared-carrier-bond edge case, tracked as [#149](https://github.com/kent-tokyo/chematic/issues/149)). Also investigated and **ruled out** a previously-suspected bridged-bicyclic canonicalization bug — the two SMILES in question turned out to be genuinely different molecules per independent RDKit InChI verification, not a chematic defect
- `chematic-perception`: new opt-in `apply_aromaticity_authoritative_experimental` — makes aromatic-flag promotion/demotion fully authoritative to the Hückel model's actual verdict in both directions (the default `apply_aromaticity`/`apply_aromaticity_ex` are unchanged, verified byte-identical). Fixed a ring-fusion-bond misclassification affecting fused diazines (quinazoline/quinoxaline/purine-shaped rings) as part of this work, with a beneficial side effect of also resolving 32 pre-existing false-positive regression pins
- Full details and known limitations in `CHANGELOG.md`

**v0.5.0** (2026-07-23): **CIP-independent 2D wedge parity, charge-aware kekulization (6 previously-hard-failing molecule classes), non-silent PAINS/Brenk budget matching**
- `chematic-perception`: new `local_parity_from_wedges`/`apply_local_parity_from_wedges` — computes `Atom.chirality`/`stereo_neighbor_order` directly from wedge/hash bonds and 2D coordinates without any CIP ranking, so a CIP tie can't erase a known local parity; sign convention measured against RDKit's raw chiral tag, not derived by analogy. Not yet wired into any reader's default parse path
- `chematic-core`: `kekulize()`'s atom-matching rules were charge-blind and missing Tellurium — tropylium, imidazolium, pyridinium, pyrylium, tellurophene, and phosphole now kekulize successfully, bond-for-bond identical to RDKit, with zero regressions; `Element::normal_valences()` gained a source-verified Tellurium entry, fixing a resulting ECFP4 aromaticity mismatch
- `chematic-perception`: charge-aware Hückel π-electron counting — tropylium, imidazolium, pyridinium, and pyrylium now match RDKit's aromatic atom/bond flags exactly (tellurophene/phosphole and a broader authoritative-flag-demotion fix remain open, tracked separately)
- `chematic-smiles`: fixed two writer bugs — bracket-forced atoms (isotope/charge/atom-map) were silently dropping implicit hydrogens (`[NH4+]` → `[N+]`), and standalone wedge bonds with no adjacent double bond were written as meaningless SMILES `/`/`\` tokens
- `chematic-smarts`/`chematic-chem`: PAINS/Brenk substructure matching could hang for minutes on symmetric targets; VF2 now takes an explicit visit budget with a genuine three-way outcome (`Found`/`NotFound`/`BudgetExhausted`) that's never silently folded into a false negative — see [#139](https://github.com/kent-tokyo/chematic/issues/139) for the follow-up (symmetry-aware candidate ordering) that would let the shipped budget resolve some now-conservatively-flagged cases correctly instead of just safely
- Full details, corpus-level before/after numbers, and known limitations in `CHANGELOG.md`

**v0.4.30** (2026-07-17): **`chematic-cip` opt-in wired to every surface, SMARTS `[rN]` fix, new RDKit-parity aromaticity engine, 5 stereo-metadata bug fixes**
- `chematic-smarts`: fixed `[rN]` (ring-size SMARTS, e.g. `[r5]`/`[r6]`) being wrongly aliased to `[kN]`'s any-ring semantics — RDKit's real `[rN]` means "this atom's *smallest* ring is exactly size N", a materially different predicate (confirmed empirically: on a fusion atom shared between a 5-ring and 6-ring, RDKit's `[k6]` matches but `[r6]` doesn't). No ring-model change — `[rN]` now has its own `MinRingSize` primitive computed from chematic's existing SSSR. SMARTS match-set agreement vs RDKit **96.9% → 99.93%** on a 5,000-molecule corpus, 0 regressions; see `docs/rdkit_compat.md`'s "SMARTS-R0"/"SMARTS-R1" entries
- Milestone 5A: opt-in access to the accurate engine from every public surface — `chematic_chem::assign_cip_with_mode(mol, CipMode::Accurate)` (Rust), `Mol.cip_stereo(mode="accurate")` + `Mol.cip_stereo_unresolved()` (Python), `cip_assignments_accurate_json`/`cip_unresolved_json` (WASM). Every default (`assign_cip()`, `cip_stereo()`, the un-suffixed WASM functions) is unchanged — this is additive only, not a default switch; see `docs/cip_accurate_rfc.md`'s Milestone 5A entry for the merge semantics (accurate tetrahedral R/S + legacy E/Z, since the accurate engine doesn't compute bond stereo) and the "never guess" contract (ties/budget-outs surface explicitly, never silently backfilled)
- Milestone 4 gate closed: 99.64% oracle-stable agreement (raw 99.38%, 4160/4186) on a full-corpus, representation-stability-stratified score — the last residual (11 phosphorus cyclophosphazene rows) turned out to be an oracle instability (RDKit's own labels change under a chemically-neutral Kekulé respelling of the identical molecule), not a chematic defect; the 15-row Rule 5 cage family remains deferred, unaffected by this gate
- The accurate engine (a provenance-carrying, sphere-by-sphere digraph comparator — Rules 1a/1b/2 — plus RDKit-compatible MANCUDE fractional atomic numbers for aromatic ring stereocenters) is available through the opt-in APIs above but is not yet the default implementation behind `assign_cip()`
- Found and fixed a real ~10-14x perf regression (SSSR misused for a boolean ring-bond check, replaced with an O(V+E) bridge-edge DFS); CI Criterion-gate bootstrap fix; a Criterion-gate reliability finding (pseudo-replication, [#70](https://github.com/kent-tokyo/chematic/issues/70)) — process-level redesign (independent process-run observations, two-stage screening, same-binary null control) landed. Follow-up ([#117](https://github.com/kent-tokyo/chematic/pull/117)) fixed two further, more specific bugs found via real CI runs: Stage 1's 3-block sign-test could structurally never fail (replaced with a pure magnitude-routing screen, threshold chosen from an offline evaluation against 28 historical no-op runs) and Stage 2's sign-test was independently magnitude-blind to small build/codegen variance (fixed with a combined sign-test + practical-effect-threshold gate). The gate stays non-required; the same-binary null control's blind spot to cross-binary codegen variance is a known, reported, not-yet-fixed gap
- Milestone 4A: `CipCode::LowerR`/`LowerS` — Rule 5 (pseudoasymmetry), scoped to 2 verified-independent rows; a three-armed symmetric-cage family (15 rows) was found to be provably unreachable by this pairwise architecture and deferred as Milestone 4A-2 (needs symmetry/automorphism detection)
- Milestone 4A-0: re-froze the residual fresh at 34 rows and mechanically classified 100% of it (0 unexplained) — 15 Rule 5/pseudoasymmetry (the 4A-2 cage family), 8 Rule 4 candidate (positively confirmed via a structural-identity check, not inferred), 11 phosphorus (9 comparator-bug "wrong" + 2 genuinely-tied)
- `chematic-perception`: new opt-in `assign_aromaticity_rdkit_parity_experimental`/`apply_aromaticity_rdkit_parity_experimental` — a source-verified port of RDKit's actual aromaticity algorithm, **100.0000% atom/bond agreement** with real RDKit on 4,999/5,000 comparable molecules. Not wired into the default path (`RdkitLike`/`Huckel` unchanged); default-promotion is blocked on a pre-existing, unrelated canonical-SMILES-writer sensitivity, not this engine
- Fixed 5 instances of the same missing-metadata-copy bug (a `MoleculeBuilder` rebuild not calling `copy_stereo_groups_from`/`copy_stereo_from`/`copy_bond_directions_from`), each silently dropping `stereo_neighbor_order` or worse: `apply_kekule` (P0), `enumerate_stereoisomers` (could silently flip a newly-assigned stereocenter's CIP code), `transfer_hydrogen_aromatic`/`clone_mol` (now deleted, replaced by `Molecule::clone`), `transfer_hydrogen`, and `invert_stereocenter` (which turned out to be a functional no-op on plain `@`/`@@` SMILES input, a separate and more severe bug)
- `chematic-smiles`: consolidated aromatic bond-direction stashing across all 3 parser bond-creation paths (chain-edge/ring-closure/branch-attachment) into one shared helper — fixes a canonical-round-trip representation instability (4,994/5,000 → 5,000/5,000 stable); does not fix `assign_ez`'s pre-existing blindness to this side channel (tracked as follow-up)
- Full-corpus accuracy on the experimental CIP engine 96.68% → 99.38% raw / 99.64% oracle-stable vs modern RDKit `rdCIPLabeler` (0 regressions) — see `docs/cip_accurate_rfc.md` for the full milestone history
- Benchmarks refreshed (`benchmarks/2026-07-17.md`, Apple M4): the previous ECFP4 throughput headline (3.6 µs/mol, 5–14× vs RDKit) does not reproduce on a clean remeasurement — updated to today's measured numbers (~78 µs/mol / 2–3× on a diverse corpus) throughout this README and `docs/`; descriptor accuracy numbers reproduce cleanly

**v0.4.29** (2026-07-10): **Kabsch rotation bug fix + SDF V3000/CDXML write, Avalon FP, O3A**
- `chematic-3d`: fixed `align_coords`'s Kabsch rotation computed in the wrong direction — was giving grossly inflated RMSD for any non-pure-translation alignment (live on v0.4.28 across crates.io/PyPI/npm before this patch); `correspondence_search` for O3A atom correspondence
- `chematic-mol`: SDF V3000 write wiring; CDXML write
- `chematic-fp`: Avalon fingerprint

**v0.4.28** (2026-07-09): **SMARTS perf, registry re-sync**
- `chematic-smarts`: existence-check short-circuit — `bulk.substructure_search` 2.2× faster than RDKit
- No git tag had been pushed for v0.4.23–v0.4.27 (crates.io stayed current via manual `cargo publish`, but PyPI/npm/GitHub Releases fell behind) — this release re-syncs all three registries

**v0.4.27** (2026-07-04): **Descriptor fixes, RWMol/FCFP, veridict CI gates**
- `chematic-chem`: `kappa1-3`, `balaban_j`, `labute_asa`, `bcut2d`, `hall_kier_alpha` descriptor fixes
- `chematic-fp`: `useFeatures=True` FCFP
- `chematic-mol`: RWMol in-place editing
- CI: veridict-based performance/Criterion/accuracy-drift regression gates; integration-test CI coverage gap fix

**v0.4.26** (2026-06-29): **E/Z stereo transfer in reactions + validation Sprint 6/7**
- `chematic-rxn`: reaction products now preserve `/`/`\` double-bond geometry from reactants in `run_reactants()` (previously lost on transformation)
- Validation: canonical SMILES differential validation vs RDKit (Sprint 6); SMARTS/aromaticity differential tests + I/O compatibility (rdkit_compat Sprint 7); root-caused remaining RDKit canonical divergence to aromaticity round-trip, not Morgan ranks

**v0.4.25** (2026-06-29): **`chematic.rdkit_compat` layer**
- `chematic-py`: RDKit API compatibility surface (Sprints 1–5) — Morgan `bitInfo`, Fingerprint/Mol/Atom/Bond/RingInfo compatibility, differential tests against RDKit; streaming `SDMolSupplier`/`SDWriter`/`Mol.GetProp`
- `chematic-perception`: `AromaticityAlgorithm::RdkitLike` — Se/Te chalcogen aromaticity matching RDKit's model

**v0.4.24** (2026-06-29): **CIP Rule 5, bridgehead/rotatable-bonds/TPSA/MR to 100%, HDF fingerprints**
- `chematic-chem`: CIP Rule 5 stereo tie-breaking (stereocenters 99.8% → 99.98% vs RDKit); bridgehead detection 98.5% → 100%; rotatable bonds 99.1% → 100%; TPSA 100%; molar refractivity 97.5% → 100% (3-ring XOR augmentation) — all on the 5,000-mol ChEMBL corpus
- `chematic-py`: `bulk.descriptors_array()` columnar numpy output; true-streaming SDF (`SdfFileReader`/`iter_sdf_batched`); `screen()` compound-filter workflow
- LLM/RAG: representation router (`to_llm_text`, `best_representation`), molecule context pack, **Hyper-Dimensional Fingerprints (HDF)** — training-free dense molecular vectors

**v0.4.23** (2026-06-26): **LogP 96.5% → 99.7%**
- `chematic-chem`: `crippen_anchor_sets` fixed to use `uniquify: false`, so symmetric triple bonds (internal alkynes) yield both VF2 match orientations instead of one falling back to the generic `[#6]` value

**v0.4.22** (2026-06-26): **CITATION.cff + `chematic.doctor()`**
- `chematic-py`: `doctor()` self-diagnostic; Reliability-by-Feature matrix added to README

**v0.4.21** (2026-06-25): **HTML/Markdown reporting for LLM/Jupyter**
- `chematic-py`: `chematic.report()` self-contained HTML compound grid, `chematic.compare()`, `mol.review()` Markdown analysis
- Docs: `benchmarks/`/`validation/` reproducible accuracy history

**v0.4.20** (2026-06-25): **ETKDG torsion KB 44 → 80 rules, `mol.describe()`/`diff()`**
- `chematic-3d`: chair/envelope ring conformations for 6/5-membered aliphatic rings; SMARTS-based torsion rules as a high-precision pre-check layer
- `chematic-py`: `mol.describe()`/`mol.diff(other)` for LLM/MCP agents; `bulk.generate_3d`/`tanimoto_matrix`/`standardize`

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

Earlier v0.4.x work (template retrosynthesis, AutoDock/UFF, Kekulization blossom
algorithm, PyO3 bindings, native-inchi) and the full v0.1–v0.3 history:
[CHANGELOG.md](CHANGELOG.md)

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
├── Cargo.toml                    workspace root (v0.10.0)
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
cargo test --workspace --lib --quiet                                      # 2,812 lib tests
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
  version   = {0.10.0},
  year      = {2026},
}
```

---

## License

Licensed under either of Apache License 2.0 or MIT License, at your option.

---

If chematic saves you time, a [GitHub star](https://github.com/kent-tokyo/chematic) helps others discover it.
