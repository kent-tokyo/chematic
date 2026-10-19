# When to choose chematic

Choose chematic when a chemistry feature must run locally in a browser, fit
inside a Rust application, or expose the same core behavior through Python,
Node, and WebAssembly. It is designed for bounded operations whose failures,
unsupported inputs, and information loss must remain visible.

Do not choose it only because a microbenchmark is fast. Start from the
operation, required chemistry, deployment target, and acceptable failure
policy, then verify the exact package on your own data.

## Good fits

| Need | Why chematic fits | Start here |
|---|---|---|
| A browser tool that keeps structures on the user's device | Supported parsing, descriptors, fingerprints, search, and depiction run in WebAssembly without a chematic backend | [Browser integration](use-cases/browser-app.md) · [Explorer](https://kent-tokyo.github.io/chematic/explorer/) |
| A Rust service without a Python or C++ runtime | The Rust crates can be embedded directly and expose typed errors and resource limits | [Rust server guide](use-cases/rust-server.md) |
| A Python batch pipeline with explicit rejected rows | Wheels share the Rust core; batch and streaming APIs can retain source indices and failures | [Batch analysis](use-cases/batch-analysis.md) |
| One chemistry core across native and browser products | Rust, Python, Node, and WebAssembly bindings are built around the same implementation | [Language bindings](language-bindings.md) |
| An agent that must receive structured chemistry errors | The local MCP server returns named tools and machine-readable failures | [MCP guide](use-cases/ai-drug-discovery.md) |
| A measured migration of selected RDKit calls | Named compatibility profiles, validation records, and typed unsupported cases make the boundary reviewable | [RDKit migration](rdkit-migration.md) |

## When another toolkit is the better default

- Choose **RDKit** when maximum chemistry breadth, its Python ecosystem, or an
  API outside chematic's declared compatibility profiles is required.
- Choose **Open Babel** when broad command-line file conversion or one of its
  long-tail formats matters more than chematic's typed semantic-loss reports.
- Compare **COSMolKit** and chematic on your own corpus when a pure-Rust
  implementation is the main requirement; do not transfer results between
  different operations, versions, or hosts.
- Use a specialized library when you need production pKa/ADMET prediction,
  broad Markush/polymer chemistry, rich CDXML editing, or an established 3D
  workflow that chematic still marks experimental.

This is a selection guide, not a claim that one toolkit is universally
better. A project can also use chematic for a local browser or Rust path while
retaining RDKit or Open Babel for unsupported work.

## Verify one operation before adopting

1. Pin the chematic and comparator versions.
2. Run the exact operation on representative inputs, including invalid and
   unsupported cases.
3. Count successes, refusals, unsupported rows, and errors separately.
4. Compare outputs before comparing speed.
5. Record the host, options, corpus hash, and package hash.

The [validation report](validation.md) lists current denominators. The
[benchmark method](benchmark.md) explains equal-work timing, and the
[compatibility scope](compatibility-scope.md) identifies bounded or
experimental chemistry. Dated records are retained in the
[benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks).

## Minimal Python check

```bash
pip install chematic==1.1.8
```

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
print({"mw": mol.mw, "tpsa": mol.tpsa, "qed": mol.qed})
print(mol.has_substructure("[OH]"))
```

For production evaluation, keep the input identifier beside every result and
do not silently replace a rejected structure with a modified one.
