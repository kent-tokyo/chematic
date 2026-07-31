// Runtime smoke test for the --target web wasm-pack artifact specifically
// (as opposed to pipeline_v2.test.mjs, which uses --target nodejs).
//
// .github/workflows/publish-npm.yml and pages.yml both build with
// `--target web` -- that's what actually ships to npm and the Pages demo, so
// it's exercised here directly rather than assumed to behave like the
// nodejs target. Loads the built artifact's raw .wasm bytes via `initSync`
// (Node has no `fetch`-a-relative-file story the `web` target's default init
// path expects), which is sufficient to exercise the wasm32-unknown-unknown
// binary itself -- the same one a browser would load.
//
// Run manually after building the web target:
//
//   wasm-pack build crates/chematic-wasm --target web --out-dir pkg-web --release
//   node crates/chematic-wasm/tests/pipeline_v2_web_target.test.mjs

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const pkgDir = path.join(__dirname, "..", "pkg-web");

const wasm = await import(path.join(pkgDir, "chematic_wasm.js"));
wasm.initSync({ module: readFileSync(path.join(pkgDir, "chematic_wasm_bg.wasm")) });

const { parse_smiles, embed_pipeline_v2_json } = wasm;

const SAFE_CONFIG = {
  embedSeed: 7,
  maxAttempts: 8,
  embedTimeoutMs: null,
  useExpTorsions: false,
  useSmallRingTorsions: false,
  useMacrocycleTorsions: false,
  useMacrocycle14Bounds: false,
  includeLegacyTorsionHeuristic: false,
  stereoPolicy: "ignore",
  failOnUnevaluableStereo: false,
  forceFieldPolicy: "none",
  forceFieldMaxIterations: 200,
  gateMmff94TorsionOop: false,
  ringTorsionPolicy: "fail_closed",
  totalTimeoutMs: null,
};

{
  const mol = parse_smiles("CCCCCCCCCC"); // decane
  const response = JSON.parse(embed_pipeline_v2_json(mol, JSON.stringify(SAFE_CONFIG)));
  assert.equal(response.ok, true, "web-target success");
  assert.equal(response.result.coords.length, 10);
  for (const [x, y, z] of response.result.coords) {
    assert.ok(Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z));
  }
  console.log("web-target success: ok (decane, 10 finite coords)");
}

{
  const mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
  const response = JSON.parse(
    embed_pipeline_v2_json(mol, JSON.stringify({ ...SAFE_CONFIG, totalTimeoutMs: 0 })),
  );
  assert.equal(response.ok, false, "web-target typed timeout");
  assert.equal(response.error.cause.kind, "timeout");
  console.log("web-target typed timeout: ok (totalTimeoutMs: 0)");
}

console.log("pipeline_v2_web_target.test.mjs: all assertions passed");
