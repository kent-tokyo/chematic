# Use case: Browser-first chemistry app with WASM

## Problem

You want to ship a chemistry tool to users who won't install anything — a web app for medicinal chemists, a public screening tool, or an internal dashboard. Server-side chemistry APIs add latency, infrastructure cost, and data-privacy concerns. RDKit.js at ~30 MB is too heavy for a smooth page load.

## Solution

chematic compiles to WebAssembly at **504 KB gzip** — roughly 60× smaller than RDKit.js. No server required: descriptor calculation, fingerprint generation, and similarity search run entirely in the browser, offline-capable after first load.

## Output / What you get

A React component that renders a 2D structure + property card from a SMILES string, entirely client-side. Load time for the WASM module: ~150 ms on a 4G connection.

## Why browser-first matters

- Zero installation for end users
- Works offline after first load
- No data leaves the browser (suitable for proprietary structures)
- Embeds in any web app with a single script tag

## Setup

```bash
npm install @kent-tokyo/chematic-wasm
```

## SMILES to descriptors in the browser

```js
import init, { parse_smiles } from "@kent-tokyo/chematic-wasm";
await init();

const mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O");
console.log(mol.molecular_weight());   // 180.16
console.log(mol.tpsa());               // 63.6
console.log(mol.lipinski_passes());    // true
console.log(mol.qed());               // 0.55
```

## Similarity search in the browser

```js
import init, { SimilarityIndex } from "@kent-tokyo/chematic-wasm";
await init();

const library = ["CCO", "c1ccccc1", "CC(=O)O", "CCCCCC", "c1cccnc1"];
const idx = SimilarityIndex.from_smiles(library);
const hits = idx.search("CC(=O)Oc1ccccc1C(=O)O", 0.3, 5);
// [{index: 2, score: 0.38}, ...]
```

## React component example

```jsx
import { useState, useEffect } from "react";
import init, { parse_smiles } from "@kent-tokyo/chematic-wasm";

let wasmReady = false;

function svgToDataUrl(svgString) {
  return "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgString)));
}

export function MoleculeCard({ smiles }) {
  const [info, setInfo] = useState(null);

  useEffect(() => {
    (async () => {
      if (!wasmReady) { await init(); wasmReady = true; }
      const mol = parse_smiles(smiles);
      if (!mol) return;
      setInfo({
        mw:     mol.molecular_weight().toFixed(2),
        logp:   mol.logp().toFixed(2),
        tpsa:   mol.tpsa().toFixed(1),
        passes: mol.lipinski_passes(),
        svgUrl: svgToDataUrl(mol.svg()),
      });
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
import init, { sdf_to_smiles_json, parse_smiles } from "@kent-tokyo/chematic-wasm";
await init();

document.getElementById("file-input").addEventListener("change", async (e) => {
  const text = await e.target.files[0].text();
  const parsed = JSON.parse(sdf_to_smiles_json(text));

  const results = parsed.map(({ smiles, name }) => {
    const mol = parse_smiles(smiles);
    return {
      name,
      mw:     mol?.molecular_weight(),
      passes: mol?.lipinski_passes(),
    };
  });

  renderResultsTable(results);
});
```

## Performance

| Task | chematic WASM | RDKit.js |
|------|--------------|----------|
| Bundle size (gzip) | 504 KB | ~30 MB |
| Parse SMILES | ~0.5 us | ~2 us |
| ECFP4 fingerprint | ~2 us | ~6 us |
| Tanimoto (pair) | ~0.1 us | ~0.3 us |

All benchmarks run in Chrome 124 on M2 MacBook Pro.

## Related APIs

- `parse_smiles(smiles)` — returns a `Mol` with all descriptor methods
- `SimilarityIndex.from_smiles(library)` / `.search(query, threshold, k)` — LSH nearest-neighbour
- `sdf_to_smiles_json(text)` — parse SDF file contents to `[{smiles, name}]`
- `mol.svg()` / `mol.svg_highlighted(atoms, color)` — 2D structure rendering
- [Live demo](https://kent-tokyo.github.io/chematic/playground/) — try WASM in the browser now
