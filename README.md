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
Pure Rust by default · optional native InChI C FFI · Python · WebAssembly · [Website](https://chematic.io/) · [Live Demo](https://kent-tokyo.github.io/chematic/playground/)

| | chematic | RDKit (Python) | RDKit.js (WASM) |
|---|---|---|---|
| **Get started** | `pip install chematic` | `pip install rdkit` (official prebuilt wheels) or conda | `npm install @rdkit/rdkit`, no Python bindings |
| **Browser bundle** | **3.30 MB raw / 1.21 MB gzip** | not applicable (Python/C++ library) | 6.91 MB raw* |
| **ECFP4 batch** | **54.7 µs/mol** | 94.3 µs/mol | — |
| **Canonical SMILES** | **24.95 / 18.27 µs/mol** | 25.58 / 26.82 µs/mol | — |
| **SDF graph read / serialization-only write** | **9.48 / 7.62 µs/mol** | 99.96 / 79.54 µs/mol | — |
| **Memory safety** | compiler-enforced (Rust) | C++ | C++ |
| **Build from source** | `cargo build` only | cmake + clang + Boost | Emscripten SDK |

\* RDKit.js gzip-over-the-wire size was not independently measured; raw figures are compared
on a like-for-like basis. RDKit.js is currently in a maintainer transition (see its repo for
current status).

The canonical and SDF rows are scoped 2026-09-04 macOS arm64 medians, not
cross-platform claims; see the exact corpora and operation boundaries in the
[benchmark details](docs/benchmark.md).
The chematic WASM size was measured 2026-09-04 from the v1.0.2 release candidate with
`wasm-pack 0.13.1` + `wasm-opt 130 -O3`: **3.30 MB raw** (**1.21 MB gzip**). The pinned
historical comparators are RDKit.js **6.91 MB**
(`@rdkit/rdkit@2025.3.4-1.0.0`'s `RDKit_minimal.wasm`, via unpkg.com) · Indigo (Ketcher build)
**11.24 MB** (`indigo-ketcher@1.45.1`'s main `.wasm`, via jsDelivr) — chematic's raw WASM binary
is currently about 2.1× smaller than RDKit.js's and about 3.8× smaller than Indigo's Ketcher-oriented
build, on a raw-to-raw basis. See the [artifact record](benchmarks/2026-09-04-wasm-size.md).

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

### v1.0.7 release boundary

The v1.0.7 release retains the v1.0.0 bounded compatibility contract while
adding typed reaction documents, document-level CDXML edits, explicit bounded
Markush/polymer expansion, crystal composition summaries, safer UFF rescue,
and canonical/SDF hot-path improvements. Spectrophores is intentionally
removed from the Rust and Python APIs while its patent/FTO status remains
independently uncleared. The
complete compatibility contract and reproducible local release gate are in
[`docs/compatibility-scope.md`](docs/compatibility-scope.md) and
[`docs/v1.0-local-release-gate.md`](docs/v1.0-local-release-gate.md). The
algorithm and third-party provenance boundary is recorded in
[`docs/implementation-provenance.md`](docs/implementation-provenance.md).

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
| **Browser app** | 1.21 MB gzip WASM bundle, zero backend required, React/Vue/Svelte ready |
| **Jupyter notebook** | `mol` renders SVG inline; `descriptors_df()` returns a pandas DataFrame |
| **Batch analysis** | Rayon-parallel descriptor/fingerprint/3D pipelines; SDF/CSV in, CSV out |
| **Rust server** | Pure-Rust crates with no C/C++ toolchain; Axum/Actix compatible |

Full worked examples → [Use cases](https://kent-tokyo.github.io/chematic/use-cases/)

---

## When to use chematic

**Use chematic if:**

- You want chemistry in the browser (WASM, 1.21 MB gzip, no server required)
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
- [Compatibility scope](docs/compatibility-scope.md) — v1.0 boundary for RWMol, CDXML, polymer expansion, and RDKit/Morgan compatibility

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
# chematic v1.0.7
# Python 3.12.x  |  darwin arm64
#
# Descriptor accuracy (benchmark 2026-08-23, v0.18.0 vs RDKit 2026.03.4):
#   MW                    99.82% within ±0.01 Da
#   HBA / HBD / ARC       100%   (4,999-mol ChEMBL subset)
#   TPSA                  100%   within ±0.1 Å²
#   LogP (Crippen)        100%*  (max Δ = 1.1×10⁻¹³)
#   Stereocenter count    99.96% (legacy) / 98.6% (new CIP FindPotentialStereo)
#   CIP R/S/E/Z labels    99.64% stable-oracle agreement (15 P rows fail closed)
# ...
```

---

## For AI / LLM Developers

chematic ships a native **MCP (Model Context Protocol) server** for local AI
agent integration.

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
The recorded v0.18.0 ECFP4 batch median is **54.7 µs/mol** versus RDKit's
94.3 µs/mol on the same 5,000-molecule corpus and Apple M4 environment. This
is a dated, corpus-specific result; current claims and reproduction details are
kept in [the benchmark guide](docs/benchmark.md).

### Safe

The common chemistry core is safe Rust and public untrusted-input paths use
finite defaults and typed failures. The optional `native-inchi` feature vendors
the IUPAC InChI C library and is the documented FFI exception. Dependencies may
contain their own unsafe code; see the [security policy](SECURITY.md) and unsafe-
surface gate for the exact boundary.

### Anywhere

Pure Rust compiles to `wasm32-unknown-unknown` natively — no Emscripten, no `cmake`,
no `clang`. The npm package `@kent-tokyo/chematic` is **1.21 MB gzip** (3.30 MB raw) —
roughly 2.1× smaller than RDKit.js's `RDKit_minimal.wasm` (6.91 MB raw) on a like-for-like
raw-size basis. One codebase targets Linux, macOS, Windows, and browser WASM;
Chromium, Firefox, and WebKit are covered by the browser CI lane.

---

## Benchmarks & Validation

| Metric | Recorded result | Scope |
|---|---|---|
| Canonical SMILES | 24.95 vs 25.58 µs/mol; 18.27 vs 26.82 µs/mol | chematic/RDKit, two 5,000-entry corpora, macOS arm64 |
| SDF graph read | 9.48 vs 99.96 µs/mol | chematic/RDKit, 365 records, graph-only |
| SDF serialization-only write | 7.62 vs 79.54 µs/mol | chematic/RDKit, same corpus, layout disabled |
| Molecular weight | 99.82% within ±0.01 Da | 4,999-molecule ChEMBL-derived corpus |
| HBA/HBD/TPSA/LogP | 100% at documented tolerances | same corpus |
| CIP R/S/E/Z | 99.64% | opt-in accurate engine; 15 representation-unstable P rows fail closed |
| WASM artifact | 3.30 MB raw / 1.21 MB gzip | v1.0.2 candidate, dated build |

These are dated, operation-specific measurements rather than universal performance or parity claims. See the [benchmark guide](docs/benchmark.md), [validation report](docs/validation.md), and [dated records](benchmarks/) for versions, corpus hashes, hardware, tolerances, and commands.

---

## Comparison with Other Cheminformatics Libraries

| Feature                 | **chematic**                              | RDKit (rdkit-sys)  | OpenBabel FFI  | RDKit.js (WASM)    |
|-------------------------|-------------------------------------------|--------------------|----------------|--------------------|
| **C/C++ dependencies**  | **None (default)**†                       | Extensive C++      | Extensive C++  | C++ via Emscripten |
| **WASM binary size**    | **3.30 MB raw** (1.21 MB gzip)             | N/A (no WASM)      | N/A (no WASM)  | 6.91 MB raw        |
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

**1.21 MB gzip — roughly 2.1× smaller than RDKit.js's raw WASM.** No Emscripten, no cmake. Drop-in for browser or Node.js.

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

The WASM binding exposes selected descriptors, fingerprints, 2D/3D operations,
reactions, diversity picking, and molecular-format conversions. See the
[WASM README](crates/chematic-wasm/README.md) and generated documentation for
the current export surface.
---

## Crate Reference

| Area | Crates |
|---|---|
| Molecular graph and identity | `chematic-core`, `chematic-smiles`, `chematic-perception`, `chematic-cip` |
| Queries, descriptors, and fingerprints | `chematic-smarts`, `chematic-chem`, `chematic-fp` |
| File and reaction models | `chematic-mol`, `chematic-rxn`, `chematic-inchi`, `chematic-iupac` |
| 2D, 3D, and materials | `chematic-depict`, `chematic-3d`, `chematic-ff`, `chematic-crystal`, `chematic-ewald` |
| User interfaces | `chematic`, `chematic-py`, `chematic-wasm`, `chematic-cli`, `chematic-mcp` |

See [format capabilities](docs/format-capabilities.md), [language bindings](docs/language-bindings.md), and the individual crate READMEs for supported operations and limitations.

---

## Recent Development

**Unreleased:** closed the remaining #210 legacy UFF stereo-rescue cases and continued canonical SMILES and SDF hot-path work. The fixed-version measurements are recorded in [benchmarks](benchmarks/).

**v1.0.7 (2026-09-05):** carries forward the v1.0.6 descriptor provenance,
shared cross-binding contracts, fused/non-alternant aromaticity and held-out
CIP boundaries, and records the #149/#337 residuals as fail-closed or
diagnostic-only contracts. It also adds the measured canonical/SDF hot-path
improvements documented in [the v1.0.7 benchmark record](benchmarks/2026-09-05-hotpath-110.md).
Spectrophores remains excluded pending independent patent/FTO review.

For public release summaries, see the [changelog](CHANGELOG.md); detailed
development notes are retained in its linked archive.

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
| Molecular weight | **Stable** | 99.82% within ±0.01 Da on 4,999 mol |
| HBA / HBD | **Stable** | 100% RDKit agreement on 4,999 mol |
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

- `canonical_smiles()` is a representation, not a cache or deduplication key. Use the fail-closed `canonical_smiles_stable_key()` and handle `None`; coupled E/Z systems using aromatic direction stashes are intentionally rejected until their spelling stability is proven.
- Aromaticity and CIP have explicit default and opt-in models; the default
  Hückel path has a bounded all-carbon odd/odd fused-envelope fallback, while
  other fused/non-alternant rings and symmetric cages are not claimed as
  universal RDKit parity. Accurate phosphorus CIP is fail-closed when the
  oracle is representation-unstable.
- 3D generation and MMFF94 are Experimental. Successful output is sanity-checked but does not promise ETKDGv3 quality or complete force-field coverage.
- Python `RWMol`, CDXML editing, and Markush/polymer expansion intentionally expose bounded subsets, not complete RDKit or ChemDraw compatibility.
- Pure-Rust InChI is approximate; enable `native-inchi` for standard IUPAC InChI.

See [compatibility scope](docs/compatibility-scope.md), [validation](docs/validation.md), and [error and resource limits](docs/error-and-limits.md) for the precise contracts.

---

## Repository Structure

```
chematic/
├── Cargo.toml                    workspace root (v1.0.7)
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
│   ├── chematic-mcp/             MCP JSON-RPC server over stdio
│   ├── chematic-wasm/            WASM/Node bindings → npm @kent-tokyo/chematic
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
cargo test --workspace --all-targets --locked
cargo test -p chematic-inchi --features native-inchi --test standard_inchi
cargo clippy --workspace --all-targets --locked -- -D warnings
```

---

## Citation

If you use chematic in academic or research work, please cite:

```bibtex
@software{chematic,
  author    = {Kentaro Tanabe (kent-tokyo)},
  title     = {chematic: A pure-Rust cheminformatics toolkit},
  url       = {https://github.com/kent-tokyo/chematic},
  version   = {1.0.7},
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
