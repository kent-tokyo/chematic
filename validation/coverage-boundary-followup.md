# Coverage boundary follow-up

The coverage workflow, production algorithms, and measured file scope stay unchanged. These tests extend the 94.00% baseline on main (`2530290ea62654b7fb76b85e697f82b559240f3c`). The processed Linux report remains the authority for the next percentage; this change alone does not claim 96%.

- NMR interchange: accumulate all diagnostic categories with field paths, validate serialized size limits, and round-trip normalization and vendor metadata. These tests use the existing `serde` feature.
- Fingerprints: check bit provenance against the actual set bits for both invariant modes, radii including the documented cap, different output sizes, chirality, and double folding; verify public shortcuts.
- Graphs and descriptors: check CIP diagnostic node topology, pharmacophore feature encoding, normalized explicit-H spellings, analytic path matching counts, resource sentinels, isolated-atom isotope envelopes, and branched alcohol output compatibility.
- Force fields and geometry: check finite positive UFF/DREIDING parameter scales, periodic torsion scores, and introsort output against Rust's integer sorting oracle over adversarial and duplicate-heavy sequences.
- Formats and transport: check CIF and POSCAR diagnostics, CLI molecular format round trips, MCP argument errors and stdio output channels, raw InChI parity and isotope reconstruction, and reaction primitive expansion with shared maps and charge/H annotations.

External fixtures remain pinned to RDKit 2026.03.1. Bounds matrices expand from 216 to 728 cases using 64 molecules selected at a fixed stride from the existing benchmark corpus; every matrix entry is compared at the existing tolerance. PDB cases add complete protein, RNA, and DNA topologies with connectivity-only CONECT records, so the reader must infer bond orders from standard residue names. The oracle records each reader option separately.

Regenerate with an environment containing exactly RDKit 2026.03.1:

```sh
python scripts/generate_bounds_boundary_fixtures.py
python scripts/generate_pdb_boundary_fixtures.py
```

Validate with the unchanged command in `.github/workflows/coverage.yml`, `cargo clippy --workspace --all-targets -- -D warnings`, and `cargo fmt --all -- --check`. MCP HTTP integration needs a permitted loopback socket. Coverage demonstrates execution; it does not establish complete chemical correctness or full RDKit compatibility.
