// Tests for embed_pipeline_v2_json (crates/chematic-wasm/src/pipeline_v2.rs),
// the WASM mirror of the Python binding (crates/chematic-py/src/pipeline_v2.rs).
//
// Reads the same validation/pipeline_v2_wasm_parity_fixtures.json the Rust test
// suite's wasm_binding_matches_python_reference_fixtures /
// raw_rust_api_matches_python_reference_fixtures read -- generated once via
// scripts/gen_pipeline_v2_wasm_parity_fixtures.py against the Python binding.
//
// KNOWN BLOCKER (github.com/kent-tokyo/chematic/issues/219): any call that
// reaches chematic_3d::pipeline_v2::embed_pipeline_v2 traps under real
// wasm32-unknown-unknown, because that function (and distance_geometry_v2)
// call std::time::Instant::now() unconditionally, which panics with "time
// not implemented on this platform" outside a native target. This is a
// pre-existing gap in chematic-3d, not in this binding. Sections below that
// would otherwise exercise the real pipeline instead assert that the call
// traps (documenting today's reality); they should be converted back to real
// assertions once #219 is fixed. Sections that are rejected at the WASM
// input-validation boundary (before any Instant::now() call) are unaffected
// and assert real behavior today.
//
// Not wired into any CI workflow (this repo has no WASM-test CI job). Run
// manually after building the Node target:
//
//   wasm-pack build crates/chematic-wasm --target nodejs --out-dir pkg-node --release
//   node crates/chematic-wasm/tests/pipeline_v2.test.mjs

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");

const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
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

// Any call whose config passes WASM-level validation reaches the real
// pipeline and traps per issue #219. Assert exactly that, rather than
// crashing the whole test run on an uncaught WebAssembly.RuntimeError.
function assertTrapsPendingIssue219(fn, label) {
  assert.throws(fn, /unreachable/, `${label}: expected issue #219 trap`);
}

// ---------------------------------------------------------------------------
// Basic envelope shape (success) -- BLOCKED by issue #219
// ---------------------------------------------------------------------------

{
  const mol = parse_smiles("CCCCCCCCCC"); // decane
  assertTrapsPendingIssue219(
    () => embed_pipeline_v2_json(mol, JSON.stringify(SAFE_CONFIG)),
    "success envelope shape",
  );
  console.log("success envelope shape: traps per issue #219 (as expected today)");
}

// ---------------------------------------------------------------------------
// Failure envelope shape (ring-torsion fail-closed) -- BLOCKED by issue #219
// (the real embedding stage runs, and traps, before ever reaching the
// ring-torsion-specific fail-closed branch)
// ---------------------------------------------------------------------------

{
  const mol = parse_smiles("C1CCCCC1CCCCCCCCCCCC"); // cyclohexane + acyclic tail
  const config = { ...SAFE_CONFIG, useSmallRingTorsions: true, forceFieldPolicy: "dreiding" };
  assertTrapsPendingIssue219(
    () => embed_pipeline_v2_json(mol, JSON.stringify(config)),
    "failure envelope shape",
  );
  console.log("failure envelope shape: traps per issue #219 (as expected today)");
}

// ---------------------------------------------------------------------------
// Malformed / incomplete config JSON -- rejected at the WASM input-validation
// boundary, before any Instant::now() call -- unaffected by issue #219, still
// real assertions: never throws, never silently defaults.
// ---------------------------------------------------------------------------

{
  const mol = parse_smiles("CC");

  const cases = [
    ["malformed JSON", "{not valid json"],
    ["unknown field", JSON.stringify({ ...SAFE_CONFIG, notARealField: true })],
    [
      "missing required field",
      (() => {
        const c = { ...SAFE_CONFIG };
        delete c.maxAttempts;
        return JSON.stringify(c);
      })(),
    ],
    [
      "missing nullable timeout field (key must still be present)",
      (() => {
        const c = { ...SAFE_CONFIG };
        delete c.embedTimeoutMs;
        return JSON.stringify(c);
      })(),
    ],
    ["unknown enum value", JSON.stringify({ ...SAFE_CONFIG, stereoPolicy: "not_a_real_policy" })],
    ["out-of-range integer", JSON.stringify({ ...SAFE_CONFIG, maxAttempts: -1 })],
  ];

  for (const [label, configJson] of cases) {
    const json = embed_pipeline_v2_json(mol, configJson);
    const response = JSON.parse(json); // must always be valid JSON, never throw
    assert.equal(response.ok, false, label);
    assert.equal(response.error.stage, "wasm_input_validation", label);
    assert.equal(response.error.cause.kind, "invalid_config", label);
  }
  console.log("malformed/incomplete config JSON: ok (never throws, never silently defaults)");
}

// ---------------------------------------------------------------------------
// null timeout / a present numeric timeout both pass WASM-level validation
// (i.e. are not rejected as invalid config) -- BLOCKED by issue #219 past
// that point, since both then reach the real pipeline and trap.
// ---------------------------------------------------------------------------

{
  const mol = parse_smiles("CC");
  assertTrapsPendingIssue219(
    () => embed_pipeline_v2_json(mol, JSON.stringify(SAFE_CONFIG)),
    "null timeout accepted past validation",
  );
  assertTrapsPendingIssue219(
    () => embed_pipeline_v2_json(mol, JSON.stringify({ ...SAFE_CONFIG, totalTimeoutMs: 60000 })),
    "present timeout accepted past validation",
  );
  console.log(
    "nullable timeout fields (null and present): accepted past validation, then trap per issue #219",
  );
}

// ---------------------------------------------------------------------------
// Safety limits: oversized config JSON -- rejected before any Instant::now()
// call -- unaffected by issue #219, still a real assertion.
// ---------------------------------------------------------------------------

{
  const mol = parse_smiles("CC");
  const padded = " ".repeat(1_000_001); // WASM_MAX_INPUT_BYTES + 1
  const json = embed_pipeline_v2_json(mol, padded);
  const response = JSON.parse(json);
  assert.equal(response.ok, false);
  assert.equal(response.error.cause.kind, "input_too_large");
  console.log("oversized config JSON: fail-closed, ok");
}

// ---------------------------------------------------------------------------
// Finite/count/order checks across 4 fixtures -- BLOCKED by issue #219
// (each of these is a full real run)
// ---------------------------------------------------------------------------

{
  const fixtures = ["CCCCCCCCCC", "c1ccc2ccccc2c1", "CC(=O)Oc1ccccc1C(=O)O", "CCC(C)C"];
  for (const smiles of fixtures) {
    const mol = parse_smiles(smiles);
    assertTrapsPendingIssue219(
      () => embed_pipeline_v2_json(mol, JSON.stringify(SAFE_CONFIG)),
      `finite/count/order: ${smiles}`,
    );
  }
  console.log("finite/count/order checks across 4 fixtures: trap per issue #219 (as expected today)");
}

// ---------------------------------------------------------------------------
// Cross-binding parity vs. the frozen Python-derived reference fixtures --
// BLOCKED by issue #219 for the actual pipeline run; atom-count parsing
// itself (mol.atom_count() vs. fixture.atomCount) does not touch the
// pipeline and is still checked for real.
// ---------------------------------------------------------------------------

{
  const fixturesPath = path.join(repoRoot, "validation", "pipeline_v2_wasm_parity_fixtures.json");
  const { fixtures } = JSON.parse(readFileSync(fixturesPath, "utf8"));

  for (const fixture of fixtures) {
    const mol = parse_smiles(fixture.smiles);
    assert.equal(mol.atom_count(), fixture.atomCount, `${fixture.name}: atom count`);
    assertTrapsPendingIssue219(
      () => embed_pipeline_v2_json(mol, JSON.stringify(fixture.config)),
      `parity fixture: ${fixture.name}`,
    );
  }
  console.log(
    `cross-binding parity vs. Python reference (${fixtures.length} fixtures): atom counts ok, ` +
      "pipeline run traps per issue #219 (as expected today)",
  );
}

console.log("pipeline_v2.test.mjs: all assertions passed (see issue #219 for the pending blocker)");
