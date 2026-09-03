<div class="chm-hero" markdown>

# Cheminformatics that runs entirely in your browser

<p class="chm-subhead">A compact Rust and WebAssembly chemistry engine for interactive tools,
local analysis, and serverless applications — plus native Rust and Python bindings from the
same codebase. No backend required for supported browser workflows.</p>

<div class="chm-cta-row">
  <a class="chm-btn chm-btn-primary" href="https://kent-tokyo.github.io/chematic/explorer/">Open Local Compound Explorer</a>
  <a class="chm-btn chm-btn-secondary" href="https://kent-tokyo.github.io/chematic/playground/">Try the Playground</a>
  <a class="chm-btn chm-btn-secondary" href="getting_started/installation/">Install chematic</a>
</div>

<p class="chm-links-row">
  <a href="https://github.com/kent-tokyo/chematic">View on GitHub</a>
  <a href="use-cases/browser-app/">Read the browser integration guide</a>
  <a href="rdkit-comparison/#wasm-deployment">View benchmark methodology</a>
</p>

</div>

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions)
[![PyPI](https://img.shields.io/pypi/v/chematic)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic)](https://www.npmjs.com/package/@kent-tokyo/chematic)

---

## Runs locally, ships light, one core everywhere

**Runs locally.** Supported analysis (parsing, descriptors, fingerprints, similarity search,
2D depiction) executes inside the browser's own WASM sandbox — the molecule data you type or
upload is never sent to a chematic server. The [Local Compound Explorer](https://kent-tokyo.github.io/chematic/explorer/)
and [Playground](https://kent-tokyo.github.io/chematic/playground/) are both static pages with
no backend of their own. (This describes chematic's own browser tools; if you build a product
on top of chematic-wasm that calls other network APIs, that's your own code's choice, not
something chematic does on your behalf.)

**Lightweight deployment.** The WASM bundle is **3.30 MB raw / 1.21 MB gzip**, measured
2026-09-04 from the v1.0.2 release candidate (`wasm-pack 0.13.1` +
`wasm-opt 130 -O3`) — see the [artifact record](../benchmarks/2026-09-04-wasm-size.md)
for the digest and reproduction commands.

**One Rust core, multiple interfaces.** The same `chematic-*` Rust crates back the native Rust
API, the Python bindings (`pip install chematic`), and the WASM/JavaScript bindings
(`npm install @kent-tokyo/chematic`) — one implementation, not three ports to keep in sync.

---

## 30 seconds of chematic

=== "JavaScript / WASM"

    ```js
    import init, { parse_smiles } from "@kent-tokyo/chematic";
    await init();

    const mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O");  // aspirin
    console.log(mol.molecular_weight(), mol.tpsa(), mol.lipinski_passes());
    // 180.16  63.6  true
    mol.free();
    ```

=== "Rust"

    ```rust
    use chematic::{smiles, chem};

    let mol = smiles::parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();  // aspirin
    println!("{:.2} {:.1}", chem::molecular_weight(&mol), chem::tpsa(&mol));
    // 180.16 63.6
    ```

=== "Python"

    ```python
    import chematic

    mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
    print(mol.mw, mol.tpsa, mol.lipinski_passes)
    # 180.16  63.6  True
    ```

---

## Pick your entry point

<div class="chm-card-grid" markdown>

<div class="chm-card" markdown>
### Browser developers
- JavaScript / TypeScript, native WASM (no Emscripten)
- No backend for supported local workflows
- SVG 2D depiction, descriptors, fingerprints, similarity search
- [Browser integration guide →](use-cases/browser-app.md)
</div>

<div class="chm-card" markdown>
### Rust developers
- Native Rust API, `cargo add chematic`
- Zero C/C++ toolchain in the standard pure-Rust path
- Embeds in servers, CLIs, and embedded targets
- [Rust server guide →](use-cases/rust-server.md)
</div>

<div class="chm-card" markdown>
### AI developers
- Built-in MCP server, 20 structured chemistry tools
- Runs locally over stdio — no hosted service
- Does not implement remote/HTTP MCP transports
- [AI-assisted analysis guide →](use-cases/ai-drug-discovery.md)
</div>

<div class="chm-card" markdown>
### Python users
- `pip install chematic` — prebuilt wheels, no C/C++ compiler needed
- Jupyter-friendly inline SVG rendering, pandas DataFrame export
- RDKit-familiar API subset — not a full drop-in replacement
- [Python notebook guide →](use-cases/python-notebook.md)
</div>

</div>

---

## Common use cases

| Scenario | How chematic helps |
|---|---|
| **Local compound triage** | [Local Compound Explorer](https://kent-tokyo.github.io/chematic/explorer/) — load a CSV/SDF, filter, sort, and export, entirely client-side |
| **Browser app** | 1.21 MB gzip WASM bundle, zero backend required, React/Vue/Svelte ready |
| **Drug screening** | 190+ descriptor values, ADMET, PAINS/Brenk, QED — batch over thousands of compounds |
| **AI agent / MCP** | Built-in MCP server — Claude Desktop can call chemistry tools directly |
| **Batch analysis** | Rayon-parallel descriptor/fingerprint/3D pipelines; SDF/CSV in, CSV out |
| **Rust server** | Pure-Rust crates with no C/C++ toolchain; Axum/Actix compatible |

Full worked examples → [Use cases](use-cases/)

---

## Honest comparison

| | chematic | RDKit (Python) | RDKit.js (WASM) |
|---|---|---|---|
| Install | `pip install chematic` | `pip install rdkit` (official prebuilt wheels) or conda | `npm install @rdkit/rdkit`, no Python bindings |
| C/C++ toolchain | Not required, even building from source | Not required for the prebuilt wheel; required building from source | Not required by consumers of the published package |
| Browser / WASM | Yes — 3.30 MB raw / 1.21 MB gzip | Not applicable (Python/C++ library) | Yes — 6.91 MB raw (`RDKit_minimal.wasm`; a separate community project, currently in a maintainer transition) |
| pKa / ADMET prediction | Built-in, rule-based screening — not for clinical use | External tool required | External tool required |
| AI agent / MCP integration | Built-in, 20 tools (stdio only) | — | — |
| Ecosystem maturity | Growing (2024–) | Established (2006–) | Established, but the WASM distribution specifically is community-maintained |

The chematic bundle was measured 2026-09-04 from the v1.0.2 candidate; RDKit.js is a
pinned historical raw-size comparator because its gzip-over-the-wire size was not independently
measured. Full detail, including where chematic is weaker: [Detailed RDKit comparison](rdkit-comparison.md).

---

## Validation

Descriptor accuracy is measured against RDKit on a 4,999-molecule ChEMBL-derived corpus:
MW, HBA, HBD, TPSA, LogP (Crippen), molar refractivity, Fsp3, and ring/stereocenter counts all
reach 100% agreement (LogP within float64 rounding error). Full breakdown, known residuals, and
reproduction commands: [Validation report](validation.md).

---

## When to use chematic

- You want chemistry in the browser (WASM, 1.21 MB gzip, no server required)
- You need a pure Rust stack with no C++ toolchain dependencies
- You deploy to environments where installing RDKit is impractical (Cloudflare Workers, Lambda, embedded)
- You build AI agents and want native MCP tool integration
- You want `pip install chematic` to just work, anywhere, no compiler needed

## When to use RDKit

- You need maximum ecosystem compatibility and 20+ years of production validation
- You need publication-quality 3D structures with ML-assisted torsion corrections (ETKDGv3)
- You depend on community plugins written against the RDKit Python API
- You need bit-exact standard InChI without enabling an opt-in feature

---

## Quick links

- [Local Compound Explorer](https://kent-tokyo.github.io/chematic/explorer/) — analyze a batch of compounds entirely in your browser
- [Playground](https://kent-tokyo.github.io/chematic/playground/) — interactive single-molecule WASM demo
- [Cookbook](cookbook.md) — 20 copy-paste-ready tasks
- [Use cases](use-cases/) — AI agent workflows, notebooks, browser apps, Rust servers, batch analysis
- [Benchmark](benchmark.md) — performance vs RDKit, descriptor accuracy
- [RDKit migration guide](rdkit_cheatsheet.md) — side-by-side API comparison
- [API Reference](api/chematic.md) — full Python API
- [GitHub](https://github.com/kent-tokyo/chematic)
- [crates.io](https://crates.io/crates/chematic)
- [PyPI](https://pypi.org/project/chematic/)
