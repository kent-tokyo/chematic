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

<p class="chm-links-row">Current release: <strong>v1.0.41</strong> · <a href="changelog/">release notes</a></p>

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

The v1.0.41 source comparison covers 161 RDKit 2026.03.1-compatible operations
on ChEMBL 5k and RDKit.js 10k, plus a 4,072-row unusual-SMILES stress corpus.
The published 83-row reaction gate is 80 exact graph/origin/map rows and three
inputs invalid in both engines. SMARTS/SMIRKS dialect work remains open.

Each result applies only to its named artifact, corpus, options, and failure
policy. Release availability is not accuracy evidence, source results are not
package results, and typed refusals are not matches.

See [validation](validation.md) for denominators, [benchmarks](benchmark.md)
for timed boundaries, and [compatibility scope](compatibility-scope.md) for
unsupported chemistry. Historical results remain in the
[record index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks).

## Stable and bounded areas

Stable selected paths include SMILES/SMARTS, descriptors, fingerprints,
similarity and substructure search, common MOL/SDF I/O, and the documented
Rust/Python/Node/WASM bindings.

3D, pKa/ADMET screening, IUPAC naming, Markush/polymer expansion, rich CDXML
editing, and the RDKit-style API are experimental or intentionally bounded.
`canonical_smiles()` is not always a safe identity key; use the fail-closed
`canonical_smiles_stable_key()` where the documented domain is sufficient.

RDKit remains the better choice for maximum ecosystem coverage, broader 3D
workflows, or APIs that chematic marks unsupported.

## Reference

- [Cookbook](cookbook.md)
- [Format support](format-capabilities.md)
- [API reference](api/chematic.md)
- [Errors and resource limits](error-and-limits.md)
- [MCP protocol and transport](mcp-protocol-reference.md)
- [Benchmark records](https://github.com/kent-tokyo/chematic/tree/main/benchmarks)
- [GitHub](https://github.com/kent-tokyo/chematic)
