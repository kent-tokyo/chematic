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

<p class="chm-links-row">Current release: <strong>v1.2.1</strong> · <a href="changelog/">release notes</a></p>

</div>

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions)
[![PyPI](https://img.shields.io/pypi/v/chematic)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic)](https://www.npmjs.com/package/@kent-tokyo/chematic)

## Choose chematic when local execution and explicit boundaries matter

- **Local-first:** parsing, descriptors, fingerprints, similarity search, and
  2D depiction can run in the browser. The Explorer and Playground are static
  sites and do not upload molecule data to a chematic service.
- **One core:** Rust, Python, Node, and WebAssembly use the same Rust
  implementation instead of separate ports.
- **Typed boundaries:** unsupported chemistry, malformed input, cancellation,
  and resource limits are explicit outcomes.
- **Reproducible claims:** compatibility and performance results identify the
  version, corpus, operation, options, and failure policy.

RDKit remains the better default for maximum chemistry and ecosystem breadth.
Open Babel remains a strong choice for broad command-line format conversion.
Use the [selection guide](reference/choosing-chematic.md) to compare deployment,
failure handling, supported chemistry, and evidence before choosing a toolkit.

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
| A reproducible research workflow | [Researcher guide](reference/researchers.md) |
| A migration from RDKit | [RDKit migration guide](reference/rdkit-migration.md) |

## Current evidence boundary

The v1.2.0 wheels passed the hosted three-platform runtime smoke and 11-job
chemistry workflow. That establishes package pass/fail evidence, not a general
speed or accuracy claim. Source experiments and older package measurements keep
their original versions, denominators, and host scope; typed refusals are not
counted as matches.

Use [validation](reference/validation.md) for current denominators,
[benchmark method](reference/benchmark.md) for timing rules, and
[compatibility scope](reference/compatibility-scope.md) for unsupported or
bounded chemistry.

## Stable and bounded areas

Stable selected paths include SMILES/SMARTS, descriptors, fingerprints,
similarity and substructure search, common MOL/SDF I/O, and the documented
Rust/Python/Node/WASM bindings.

3D, pKa/ADMET screening, IUPAC naming, Markush/polymer expansion, rich CDXML
editing, and the RDKit-style API are experimental or intentionally bounded.
`canonical_smiles()` is not always a safe identity key; use the fail-closed
`canonical_smiles_stable_key()` where the documented domain is sufficient.

Rust line coverage is measured separately from chemistry validation. CI
enforces at least 96.7% for its declared Rust scope; this shows which code ran,
not whether every chemical result is correct.

## Reference

- [Cookbook](reference/cookbook.md)
- [Format support](reference/format-capabilities.md)
- [API reference](api/chematic.md)
- [Errors and resource limits](reference/error-and-limits.md)
- [Product naming](reference/product-naming.md)
- [MCP protocol and transport](reference/mcp-protocol-reference.md)
- [Benchmark method and evidence](reference/benchmark.md)
- [GitHub](https://github.com/kent-tokyo/chematic)
