# Use case: Browser-first chemistry app with WASM

## Problem

You want to ship a chemistry tool to users who won't install anything — a web app for medicinal chemists, a public screening tool, or an internal dashboard. Server-side chemistry APIs add latency, infrastructure cost, and data-privacy concerns. RDKit's official JavaScript/WASM package (`@rdkit/rdkit`) is the pinned browser comparator.

## Solution

The npm package runs supported parsing, descriptors, fingerprints, search, and
depiction locally in WebAssembly. No chematic server is required, and the app
can work offline after its assets are cached.

## Output / What you get

The public [Playground](https://kent-tokyo.github.io/chematic/playground/) provides a zero-setup
reference app: start from a guided workflow, inspect 2D/3D structures and descriptors, compare
similarity, run reaction examples, or parse materials-science formats. The same WASM boundary
can be embedded in your own UI as a React component that renders a 2D structure + property card
from a SMILES string, entirely client-side.

## Why browser-first matters

- Zero installation for end users
- Works offline after first load
- No data leaves the browser (suitable for proprietary structures)
- Embeds in any web app with a single script tag

## Setup

```bash
npm install @kent-tokyo/chematic
```

## SMILES to descriptors in the browser

```js
import init, { parse_smiles } from "@kent-tokyo/chematic";
await init();

const mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O");
console.log(mol.molecular_weight());   // 180.16
console.log(mol.tpsa());               // 63.6
console.log(mol.lipinski_passes());    // true
console.log(mol.qed());               // 0.55
mol.free();   // release the WASM-side handle once its data has been extracted
```

## Similarity search in the browser

```js
import init, { MhfpLshHandle } from "@kent-tokyo/chematic";
await init();

// num_hashes must be a multiple of 16
const idx = new MhfpLshHandle(128);
for (const smiles of ["CCO", "c1ccccc1", "CC(=O)O", "CCCCCC", "c1cccnc1"]) {
  idx.add_smiles(smiles);
}
const hits = JSON.parse(idx.query_json("CC(=O)Oc1ccccc1C(=O)O", 0.3));
// [{index: 2, similarity: 0.38}, ...]
```

## React component example

```jsx
import { useState, useEffect } from "react";
import init, { parse_smiles } from "@kent-tokyo/chematic";

let wasmReady = false;

function svgToDataUrl(svgString) {
  return "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgString)));
}

export function MoleculeCard({ smiles }) {
  const [info, setInfo] = useState(null);

  useEffect(() => {
    (async () => {
      if (!wasmReady) { await init(); wasmReady = true; }
      let mol;
      try { mol = parse_smiles(smiles); }   // throws a string on invalid SMILES, never returns null
      catch (e) { return; }
      setInfo({
        mw:     mol.molecular_weight().toFixed(2),
        logp:   mol.logp_crippen().toFixed(2),
        tpsa:   mol.tpsa().toFixed(1),
        passes: mol.lipinski_passes(),
        svgUrl: svgToDataUrl(mol.depict_svg()),
      });
      mol.free();   // free the WASM-side handle once its data has been extracted
    })();
  }, [smiles]);

  if (!info) return <div>Loading...</div>;

  return (
    <div>
      <img src={info.svgUrl} alt="2D structure" width="200" />
      <dl>
        <dt>MW</dt>    <dd>{info.mw} Da</dd>
        <dt>LogP</dt>  <dd>{info.logp}</dd>
        <dt>TPSA</dt>  <dd>{info.tpsa} A2</dd>
        <dt>Lipinski</dt><dd>{info.passes ? "Pass" : "Fail"}</dd>
      </dl>
    </div>
  );
}
```

## SDF upload and analysis in the browser

```js
import init, { sdf_to_records_json, parse_smiles } from "@kent-tokyo/chematic";
await init();

document.getElementById("file-input").addEventListener("change", async (e) => {
  const text = await e.target.files[0].text();
  // sdf_to_records_json returns each record's name + canonical SMILES + SD properties
  // (sdf_to_smiles_json, by contrast, returns a bare array of SMILES strings with no names)
  const records = JSON.parse(sdf_to_records_json(text));

  const results = records.map((record) => {
    // a record that failed to parse is the bare JSON value `null`, not an object
    if (!record) return { name: null, mw: null, passes: null };
    const mol = parse_smiles(record.smiles);
    const result = { name: record.name, mw: mol.molecular_weight(), passes: mol.lipinski_passes() };
    mol.free();
    return result;
  });

  renderResultsTable(results);
});
```

## Performance evidence

Browser results depend on the package version, browser, host, operation, and
JavaScript boundary. Use the
[three-browser Morgan packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-linux-three-browser-paired20.md)
for the retained browser comparison. Exact package size is release-specific;
measure the installed npm artifact rather than reusing an older release's size.
Native Python/Rust timings do not transfer to browser workflows.

## Related APIs

- `parse_smiles(smiles)` — returns a `MolHandle` with all descriptor methods (throws a string on invalid input; call `.free()` once done with it)
- `new MhfpLshHandle(numHashes)` / `.add_smiles(smiles)` / `.query_json(smiles, threshold)` — MinHash LSH approximate nearest-neighbour index
- `sdf_to_records_json(text)` — parse SDF file contents to an array of `{smiles, name, properties, stereo_diagnostics}` (or `null` for a record that failed to parse)
- `mol.depict_svg()` / `depict_svg_grid_highlighted(smilesBlock, cols, matchSmarts)` — 2D structure rendering
- [Live demo](https://kent-tokyo.github.io/chematic/playground/) — try WASM in the browser now
- [Local Compound Explorer](https://kent-tokyo.github.io/chematic/explorer/) — load, filter, and compare a batch of compounds entirely in the browser
