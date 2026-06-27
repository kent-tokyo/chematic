# chematic — Agent Guide

Pure-Rust cheminformatics library. Zero C/C++ dependencies; **504 KB WASM bundle** (gzip).

@tasks/todo.md
@tasks/lessons.md
@README.md
@CHANGELOG.md

---

## Current State (v0.4.23)

| Crate | Purpose | Tests |
|-------|---------|-------|
| `chematic-core` | Atom, Bond, Molecule, Element, kekulization (blossom) | 69 |
| `chematic-smiles` | OpenSMILES parser/writer, canonical SMILES | 48 |
| `chematic-perception` | SSSR, 2-pass Hückel aromaticity, CIP stereo, `count_aromatic_rings` | 34 |
| `chematic-smarts` | SMARTS, VF2 subgraph, MCS (McGregor), LRU SMARTS cache | 142 |
| `chematic-chem` | 190+ descriptors, ADMET, BOILED-Egg, QED, SA Score, PAINS/Brenk, pKa | 659 |
| `chematic-fp` | ECFP/FCFP, MACCS, MAP4, AtomPair, Torsion, MHFP, ERG, Tanimoto | 87 |
| `chematic-ff` | MMFF94 full stack (7 terms), DREIDING, L-BFGS minimizer | 51 |
| `chematic-3d` | ETKDG (80 torsion rules), WHIM (22-dim), GETAWAY (19-dim), RDF (20-dim), SASA | 45 |
| `chematic-mol` | SDF/MOL V2000+V3000, CML, CDXML, **MolJSON**, PDBQT, CIF, RXN, KET, etc. | 126 |
| `chematic-depict` | 2D SVG, PDF, EPS output | 34 |
| `chematic-rxn` | Reaction SMILES/SMIRKS, RECAP/BRICS, retrosynthesis (60 retro-SMIRKS) | 25 |
| `chematic-inchi` | InChI/InChIKey (pure-Rust + native-inchi feature) | 28 |
| `chematic-iupac` | IUPAC naming (25+ compound classes) | 45 |
| `chematic-mcp` | MCP server — **18 tools** for AI agents | 31 |
| `chematic-ewald` | PME Ewald summation, B-spline interpolation | 12 |
| `chematic-wasm` | **160 WASM exports** → npm `@kent-tokyo/chematic` (504 KB gzip) | 211 |
| `chematic-py` | PyO3 Python bindings (`pip install chematic`) | — |
| `chematic` | Umbrella crate (feature-gated re-exports) | — |

**Total: ~2,319 tests (lib only), all passing.**

---

## Essential Commands

```bash
# Always run before committing
bash scripts/check.sh        # fmt + clippy + test + version consistency

# Individual checks
cargo test --workspace --lib --quiet
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check

# Python bindings (dev rebuild)
.venv/bin/maturin develop --release -m crates/chematic-py/Cargo.toml

# Benchmark vs RDKit (requires RDKit in .venv)
.venv/bin/python scripts/bench5k.py ~/Downloads/SMILES.csv
```

### Version bump (release process)

```bash
# 1. Edit Cargo.toml workspace.package.version
# 2. Sync all docs in one shot:
python3 scripts/bump_version.py
# 3. Verify:
bash scripts/check.sh
```

### Branch strategy

`feat/*` new feature · `fix/*` bug fix · `docs/*` docs only · `release/*` version bump

PRs are required (main is protected: Test, Clippy, Format Check, Build Check must pass).

---

## Code Conventions

**Hard rules:**
- No `unsafe` — `#![forbid(unsafe_code)]` in every crate
- No C/C++ deps — no `*-sys`, no `cc` build scripts (except `native-inchi` feature)
- WASM-compatible — no `std::fs`, no threads, no platform I/O in core/wasm crates
- No petgraph — implement graph algorithms directly

**Style:**
- Index newtypes: `AtomIdx(u32)`, `BondIdx(u32)` — never raw `usize` in public API
- WASM functions return `String` (JSON) or primitives — no `Vec<T>` across the boundary
- `atom.charge` is `i8`; shift `(charge.wrapping_add(8)) as u8` when hashing

**Testing:**
- Unit tests in the same file under `#[cfg(test)] mod tests`
- `bash scripts/check.sh` must pass before every commit (runs fmt + clippy + test)

---

## Codebase Map

```
crates/
├── chematic-core/src/
│   ├── atom.rs          Atom, CipCode, Chirality, Element
│   ├── molecule.rs      Molecule, MoleculeBuilder
│   └── kekulize.rs      Kekulization (Edmonds blossom + Pass 3 bridgehead-N)
├── chematic-smiles/src/
│   ├── parser.rs        parse() → Molecule
│   └── canonical.rs     canonical_smiles() with stereo parity correction
├── chematic-perception/src/
│   ├── sssr.rs          find_sssr() → RingSet, augmented_ring_set()
│   └── aromaticity.rs   apply_aromaticity(), count_aromatic_rings()
├── chematic-chem/src/
│   ├── descriptors.rs   190+ descriptors: mw, logp, tpsa, hba, hbd, …
│   ├── pka.rs           rule-based pKa (15+ rules)
│   └── admet.rs         BBB, Caco-2, hERG, CYP3A4, AMES prediction
├── chematic-mol/src/
│   ├── moljson.rs       parse_moljson / write_moljson (LLM-friendly JSON)
│   ├── cjson.rs         ChemicalJSON (Avogadro format)
│   ├── cml.rs           Chemical Markup Language
│   ├── mol2000.rs       MOL V2000 / SDF
│   └── mol3000.rs       MOL V3000
├── chematic-fp/src/
│   ├── ecfp.rs          ecfp4(), ecfp6(); FCFP
│   └── map4.rs          MAP4 fingerprint
├── chematic-mcp/src/
│   └── tools.rs         18 MCP tools (smiles_to_moljson, moljson_to_smiles, …)
└── chematic-wasm/src/
    └── lib.rs           MolHandle + 160 #[wasm_bindgen] exports
```

---

## Adding a New Format (Pattern: follow moljson.rs)

1. `crates/chematic-mol/src/myformat.rs` — `parse_myformat` / `write_myformat` / `MyFormatError`
2. `crates/chematic-mol/src/lib.rs` — `pub mod myformat; pub use myformat::…;`
3. `crates/chematic-py/src/lib.rs` — `#[pyfunction] fn from_myformat` + `Mol::to_myformat`
4. `crates/chematic-wasm/src/lib.rs` — `#[wasm_bindgen] pub fn mol_from_myformat`
5. `crates/chematic-mcp/src/tools.rs` — add schema to `list_tools()` + arm in `call_tool()`

---

## Adding a New Descriptor (Pattern: follow descriptors.rs)

```rust
// crates/chematic-chem/src/descriptors.rs
pub fn my_descriptor(mol: &Molecule) -> f64 { … }

// crates/chematic-py/src/lib.rs — in #[pymethods] impl Mol:
#[getter]
fn my_descriptor(&self) -> f64 { chematic_chem::my_descriptor(&self.inner) }
```

---

## RDKit Parity (v0.4.23, 4999-mol ChEMBL corpus)

| Descriptor | Agreement | Tolerance |
|---|---|---|
| HBA / HBD / ARC / MW | **100%** | exact |
| [nH] SMARTS | **100%** | exact |
| TPSA | **99.4%** | ±0.1 Å² |
| LogP (Crippen) | **99.7%** | ±0.01 |

Benchmark: `scripts/bench5k.py` vs RDKit `CalcTPSA(includeSandP=True)`.

---

## MolJSON (LLM-friendly molecular representation)

```json
{"atoms":[{"id":"a1","element":"C","charge":0,"isotope":null,"hydrogens":3,"aromatic":false}],
 "bonds":[{"id":"b1","source_id":"a1","target_id":"a2","order":1.0,"aromatic":false}]}
```

Python: `mol.to_moljson()` / `chematic.from_moljson(s)`  
WASM:  `to_moljson(mol)` / `mol_from_moljson(json)`  
MCP:   `smiles_to_moljson` / `moljson_to_smiles`  
Note: `hydrogens` is informational; bracket-H notation ([nH]) may differ after roundtrip.

---

## MCP Tools (18 total)

`parse_smiles` · `canonical_smiles` · `calc_properties` · `lipinski_check`  
`sa_score` · `pains_check` · `brenk_check` · `admet_profile` · `boiled_egg`  
`ecfp4` · `tanimoto` · `smarts_match` · `find_mcs` · `generate_3d`  
`name_to_smiles` · `retrosynthesis`  
**`smiles_to_moljson`** · **`moljson_to_smiles`**

---

## Known Lessons (see tasks/lessons.md for detail)

- **Pre-commit**: `bash scripts/check.sh` before every push. CI hardened (no continue-on-error).
- **cargo fmt before push**: CI fmt check fails if code is unformatted.
- **CodeQL**: GitHub Default setup conflicts with custom codeql.yml — use Default setup only.
- **cargo-deny**: internal tool crates need `license.workspace = true`; CDLA-Permissive-2.0 is allowed.
- **TPSA calibration**: P=N bridging N → 12.36, ring-N-external-P → 3.01, trivalent P (no H) → 13.59.
- **Version sync**: run `python3 scripts/bump_version.py` after editing `Cargo.toml` version.
- **Kekulization**: always kekulize before `implicit_hcount` or H-sensitive descriptors.
- **BFS double-buffer**: in ECFP/Morgan, snapshot previous-round IDs before updating.
- **CIP**: requires iterative sphere expansion for tie-breaking.
- **canonical stereo parity**: adjacency list order is bond-addition order, not SMILES text order.
