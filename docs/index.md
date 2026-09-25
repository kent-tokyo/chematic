<div class="chm-hero" markdown>

# Cheminformatics that runs locally

<p class="chm-subhead">A typed Rust chemistry core for browser, Python,
Node, and native applications. Supported browser workflows run entirely in
WebAssembly, without a chematic backend.</p>

<div class="chm-cta-row">
  <a class="chm-btn chm-btn-primary" href="https://kent-tokyo.github.io/chematic/explorer/">Open Local Compound Explorer</a>
  <a class="chm-btn chm-btn-secondary" href="https://kent-tokyo.github.io/chematic/playground/">Try the Playground</a>
  <a class="chm-btn chm-btn-secondary" href="getting_started/installation/">Install chematic</a>
</div>

<p class="chm-links-row">Current release: <strong>v1.0.25</strong> · <a href="changelog/">release notes</a></p>

</div>

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions)
[![PyPI](https://img.shields.io/pypi/v/chematic)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic)](https://www.npmjs.com/package/@kent-tokyo/chematic)

## Why chematic

- **Local-first:** parsing, descriptors, fingerprints, similarity search, and
  2D depiction can run in the browser. The Explorer and Playground are static
  sites and do not upload molecule data to a chematic service.
- **One core:** Rust, Python, Node, and WebAssembly use the same Rust
  implementation instead of separate ports.
- **Typed boundaries:** unsupported chemistry, malformed input, cancellation,
  and resource limits are explicit outcomes.
- **Reproducible claims:** compatibility and performance results identify the
  version, corpus, operation, options, and failure policy.

## Quick start

=== "JavaScript / WASM"

    ```js
    import init, { parse_smiles } from "@kent-tokyo/chematic";
    await init();

    const mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O");
    console.log(mol.molecular_weight(), mol.tpsa());
    mol.free();
    ```

=== "Python"

    ```python
    import chematic

    mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    print(mol.mw, mol.tpsa)
    ```

=== "Rust"

    ```rust
    use chematic::{chem, smiles};

    let mol = smiles::parse("CC(=O)Oc1ccccc1C(=O)O").expect("valid SMILES");
    println!("{:.2} {:.1}", chem::molecular_weight(&mol), chem::tpsa(&mol));
    ```

Installation: `pip install chematic`, `cargo add chematic`, or
`npm install @kent-tokyo/chematic`.

## Choose an entry point

| You are building | Start here |
|---|---|
| A local browser tool | [Browser integration](use-cases/browser-app.md) or the [Explorer](https://kent-tokyo.github.io/chematic/explorer/) |
| A Python notebook or report | [Python notebook guide](use-cases/python-notebook.md) |
| A Rust service or CLI | [Rust server guide](use-cases/rust-server.md) |
| An MCP-enabled agent | [AI-assisted analysis](use-cases/ai-drug-discovery.md) |
| A reproducible research workflow | [Researcher guide](researchers.md) |
| A migration from RDKit | [RDKit migration guide](rdkit-migration.md) |

## Current evidence boundary

The current published release is v1.0.25.

- It adds explicit atom-output order and source-atom provenance for fragments
  and Rust reaction products, and fixes the Python RDKit-compatible HBA profile.
- Its release-source differential covers six existing-output operations across
  46,736 molecules with zero differences. The aromatic/non-aromatic ring-closure
  writer fix is intentionally outside that unchanged-output claim.

- Its registry-installed browser package records 1.398x parse-inclusive and
  3.511x prepared compatible-Morgan speedups against official
  `@rdkit/rdkit@2026.03.6`, with 9,999/9,999 supported rows bit-exact.
- Its fixed MMFF94 stereo-safe lane produces 265/265 sound, stereo-clean,
  clash-free outputs, but measures 0.944x RDKit speed.
- The published v1.0.15 WASM asset is 4,005,280 bytes raw / 1,460,499 bytes
  gzip; official RDKit.js 2026.03.6 is 7,333,095 / 2,379,975 bytes under the
  same local compression method. Feature surfaces differ.

These are scoped measurements, not general claims that chematic is always
faster, smaller, or more accurate. The dated v1.0.24 differential remains
shared-VM source evidence, not a WASM or cross-platform performance claim. It
does not relabel the v1.0.20 package performance record.

See [validation](validation.md), [benchmark methodology](benchmark.md), and
[compatibility scope](compatibility-scope.md) for exact conditions.

## Stable and bounded areas

Stable selected paths include SMILES/SMARTS, descriptors, fingerprints,
similarity and substructure search, common MOL/SDF I/O, and the documented
Rust/Python/Node/WASM bindings.

3D, pKa/ADMET screening, IUPAC naming, Markush/polymer expansion, rich CDXML
editing, and the RDKit-style API are experimental or intentionally bounded.
`canonical_smiles()` is not always a safe identity key; use the fail-closed
`canonical_smiles_stable_key()` where the documented domain is sufficient.

RDKit remains the better choice for maximum ecosystem coverage, mature ETKDG
and force-field workflows, or APIs that chematic marks unsupported.

## Reference

- [Cookbook](cookbook.md)
- [Format support](format-capabilities.md)
- [API reference](api/chematic.md)
- [Errors and resource limits](error-and-limits.md)
- [Benchmark records](https://github.com/kent-tokyo/chematic/tree/main/benchmarks)
- [GitHub](https://github.com/kent-tokyo/chematic)
