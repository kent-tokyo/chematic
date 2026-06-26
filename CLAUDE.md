# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Pre-commit (必須)

コミット前に必ず実行すること。CI と同等のチェックを1コマンドで走らせる:

```bash
bash scripts/check.sh
```

CI に新しいチェックを追加した場合は、まずこのスクリプトがローカルでパスすることを確認してからコミットすること。

### Build & Check
```bash
cargo build --workspace
cargo check --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all
```

### Test
```bash
# All tests
cargo test --workspace --lib --quiet

# Single crate
cargo test -p chematic-chem --lib --quiet

# Single test by name pattern
cargo test -p chematic-chem --lib -- aromatic_ring_count

# native-inchi feature (requires C compiler)
cargo test -p chematic-inchi --features native-inchi --lib --quiet
cargo test -p chematic-inchi --features native-inchi --test standard_inchi --quiet
```

### Python Bindings (maturin + PyO3)
```bash
# Development build (installs into active venv)
.venv/bin/maturin develop --release -m crates/chematic-py/Cargo.toml

# Run benchmark against RDKit reference
.venv/bin/python scripts/bench5k.py ~/Downloads/SMILES.csv
.venv/bin/python scripts/bench5k.py ~/Downloads/SMILES.csv --detail
```

---

## Architecture

### Crate Layers

```
chematic-core          (no deps)                   Atom, Bond, Molecule, Element, kekulization
chematic-smiles        → core                      SMILES parser/writer, canonical SMILES
chematic-perception    → core                      SSSR, Hückel aromaticity, CIP stereo
chematic-smarts        → core, perception          SMARTS, VF2 isomorphism, MCS
chematic-chem          → core, perception,         70+ descriptors, ADMET
                          smiles, smarts, fp, iupac
chematic-fp            → core, smarts              ECFP/FCFP, MACCS, MAP4, Tanimoto
chematic-ff            → core, perception          MMFF94 / DREIDING atom typing
chematic-3d            → core, perception,         ETKDG, MD, SASA, WHIM
                          ff, chem, fp, smarts
chematic-depict        → core, perception,         2D SVG
                          rxn, smiles
chematic-rxn           → core, smiles, smarts      Reaction SMILES/SMIRKS
chematic-inchi         → core, smiles, chem        InChI/InChIKey
chematic-iupac         → core, perception          IUPAC naming
chematic-py            → nearly all                PyO3 Python bindings
chematic-wasm          → nearly all                wasm-bindgen JS bindings
chematic               (umbrella, feature-gated re-exports)
```

### Core Data Model

`Molecule` (`crates/chematic-core/src/molecule.rs`) stores only heavy atoms; implicit H is computed on demand via `implicit_hcount()`. Key access:

```rust
mol.atom_count()              // heavy atoms only
mol.atom(idx: AtomIdx)        // → &Atom  (element, charge, aromatic flag, chirality)
mol.bond(idx: BondIdx)        // → &Bond  (atom1, atom2, BondOrder)
mol.neighbors(idx: AtomIdx)   // → Iterator<(AtomIdx, BondIdx)>
mol.bond_between(a, b)        // → Option<(BondIdx, &Bond)>
```

`MoleculeBuilder` is used for programmatic construction.

### Aromaticity

Aromatic SMILES atoms (`c`, `n`, `o`, …) set `atom.aromatic = true` directly during parsing. Kekulé input requires explicit `apply_aromaticity(&mol)` (`chematic-perception`) to set aromatic flags via 2-pass Hückel (4n+2 π-electron rule).

For **aromatic ring counting**, use `chematic_perception::count_aromatic_rings(mol)` rather than directly filtering `find_sssr().rings()`. The SSSR algorithm can return a large fundamental cycle (e.g. a 9-ring for indolizine) instead of the two smaller component rings; `count_aromatic_rings` applies `augmented_ring_set` to recover missing small rings and then strips "envelope" rings (rings that equal the bond-XOR of two smaller aromatic rings) to prevent double-counting.

### Ring Perception

`find_sssr(mol)` (`chematic-perception`) returns a `RingSet`; call `.rings()` to get `&[Vec<AtomIdx>]`. `augmented_ring_set(mol, sssr_rings)` is used internally by aromaticity and ring-count code to correct SSSR decomposition artefacts in fused systems.

### Descriptor Pipeline

`chematic-chem/src/descriptors.rs` contains all physicochemical functions as free functions taking `&Molecule`. SMARTS-based descriptors (TPSA, HBA/HBD, PAINS, pKa) compile patterns via `parse_smarts` and match with `find_matches`. The LRU cache in `chematic-smarts` gives 5–20× speedup on repeated patterns.

### Python API Shape

```python
import chematic
mol = chematic.from_smiles("c1ccccc1")  # returns Mol (not Molecule)
mol.mw                          # property, not method
mol.aromatic_ring_count         # property
mol.descriptors()               # → dict of 70+ keys
chematic.smarts_match("[nH]", mol)  # module-level function
```

`Mol` wraps `Arc<chematic_core::Molecule>`. Properties are computed on each access (no caching in the Python layer).

### Benchmark / RDKit Reference

`scripts/bench5k.py` compares chematic vs RDKit on 5,000 molecules for HBA, aromatic ring count, and `[nH]` SMARTS. Current agreement: HBA 100%, aromatic ring count ~95.6%. `scripts/rdkit_benchmark.py` and `scripts/gen_rdkit_reference.py` generate the reference TSV files in `scripts/`.
