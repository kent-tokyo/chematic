// Tests for embed_pipeline_v2_json (crates/chematic-wasm/src/pipeline_v2.rs),
// the WASM mirror of the Python binding (crates/chematic-py/src/pipeline_v2.rs).
//
// Reads the same validation/pipeline_v2_wasm_parity_fixtures.json the Rust test
// suite's wasm_binding_matches_python_reference_fixtures /
// raw_rust_api_matches_python_reference_fixtures read -- generated once via
// scripts/gen_pipeline_v2_wasm_parity_fixtures.py against the Python binding.
//
// Previously blocked by issue #219 (std::time::Instant::now() panicking under
// real wasm32-unknown-unknown); fixed in #222 via chematic-3d's crate::clock
// abstraction. Every case below now exercises the real pipeline end to end.
//
// Not wired into the general CI job (this is a targeted job, see
// .github/workflows/ci.yml's `test-wasm` job). Run manually after building
// the Node target:
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

// ---------------------------------------------------------------------------
// Basic envelope shape (success)
// ---------------------------------------------------------------------------

{
  const mol = parse_smiles("CCCCCCCCCC"); // decane
  const response = JSON.parse(embed_pipeline_v2_json(mol, JSON.stringify(SAFE_CONFIG)));
  assert.equal(response.ok, true, "success envelope shape");
  assert.equal(response.schemaVersion, 1);
  assert.equal(response.result.coords.length, 10);
  for (const [x, y, z] of response.result.coords) {
    assert.ok(Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z));
  }
  console.log("success envelope shape: ok (decane, 10 finite coords)");
}

// ---------------------------------------------------------------------------
// Failure envelope shape (ring-torsion fail-closed)
// ---------------------------------------------------------------------------

{
  const mol = parse_smiles("C1CCCCC1CCCCCCCCCCCC"); // cyclohexane + acyclic tail
  const config = { ...SAFE_CONFIG, useSmallRingTorsions: true, forceFieldPolicy: "dreiding" };
  const response = JSON.parse(embed_pipeline_v2_json(mol, JSON.stringify(config)));
  assert.equal(response.ok, false, "failure envelope shape");
  assert.equal(response.error.stage, "torsion_optimization");
  assert.equal(response.error.cause.kind, "ring_torsion_application_unsupported");
  console.log("failure envelope shape: ok (ring-torsion fail-closed typed failure)");
}

// ---------------------------------------------------------------------------
// Malformed / incomplete config JSON -- rejected at the WASM input-validation
// boundary: never throws, never silently defaults.
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
// null timeout / a present numeric timeout both succeed for real.
// ---------------------------------------------------------------------------

{
  const mol = parse_smiles("CC");
  const r1 = JSON.parse(embed_pipeline_v2_json(mol, JSON.stringify(SAFE_CONFIG)));
  assert.equal(r1.ok, true, "null timeout succeeds");

  const r2 = JSON.parse(
    embed_pipeline_v2_json(mol, JSON.stringify({ ...SAFE_CONFIG, totalTimeoutMs: 60000 })),
  );
  assert.equal(r2.ok, true, "present (large) timeout succeeds");

  console.log("nullable timeout fields (null and present): both succeed for real");
}

// ---------------------------------------------------------------------------
// An actual (near-zero) timeout fails closed with a typed Timeout cause --
// never a trap, never partial coords reported as success.
// ---------------------------------------------------------------------------

{
  const mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
  const response = JSON.parse(
    embed_pipeline_v2_json(mol, JSON.stringify({ ...SAFE_CONFIG, totalTimeoutMs: 0 })),
  );
  assert.equal(response.ok, false, "zero timeout must fail closed");
  assert.equal(response.error.cause.kind, "timeout");
  assert.equal("result" in response, false, "no partial coords reported as success on timeout");
  console.log("actual timeout (totalTimeoutMs: 0): ok (typed timeout failure, no partial success)");
}

// ---------------------------------------------------------------------------
// Safety limits: oversized config JSON and atom-count limit -- both rejected
// before the pipeline runs.
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

{
  // embed_pipeline_v2_json's own `mol.inner.atom_count() > WASM_MAX_ATOMS`
  // check (-> typed `atom_limit_exceeded`) is unreachable from JS today:
  // every MolHandle-producing entry point in this crate (parse_smiles,
  // mol_from_*, etc.) already enforces the same WASM_MAX_ATOMS limit at
  // construction time, so no oversized MolHandle can ever reach this
  // function. That branch is verified directly at the Rust unit-test layer
  // instead (pipeline_v2.rs's `atom_limit_exceeded` test, built via an
  // internal-only oversized-mol constructor bypassing parse_smiles). What's
  // testable from JS is that the size limit is enforced somewhere on the
  // path -- confirmed here via parse_smiles itself.
  assert.throws(
    () => parse_smiles("C".repeat(10_001)), // WASM_MAX_ATOMS + 1
    /[Ee]xceeds maximum atom count|too large/,
    "oversized molecule must be rejected before any MolHandle exists",
  );
  console.log("atom-count limit: enforced at MolHandle construction (parse_smiles), ok");
}

// ---------------------------------------------------------------------------
// Finite/count/order checks across 4 fixtures
// ---------------------------------------------------------------------------

{
  const fixtures = ["CCCCCCCCCC", "c1ccc2ccccc2c1", "CC(=O)Oc1ccccc1C(=O)O", "CCC(C)C"];
  for (const smiles of fixtures) {
    const mol = parse_smiles(smiles);
    const response = JSON.parse(embed_pipeline_v2_json(mol, JSON.stringify(SAFE_CONFIG)));
    assert.equal(response.ok, true, smiles);
    assert.equal(response.result.coords.length, mol.atom_count(), smiles);
    for (const [x, y, z] of response.result.coords) {
      assert.ok(Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z), smiles);
    }
  }
  console.log("finite/count/order checks across 4 fixtures: ok");
}

// ---------------------------------------------------------------------------
// Cross-binding parity vs. the frozen Python-derived reference fixtures.
// Covers decane/naphthalene/aspirin/branched success, declared tetrahedral
// stereo, E/Z, ring-torsion fail-closed typed failure, and force-field
// fallback -- structural fields only (atom count, coords length, stereo
// counts, force-field actual/fallback, final soundness); timing excluded.
// ---------------------------------------------------------------------------

{
  const fixturesPath = path.join(repoRoot, "validation", "pipeline_v2_wasm_parity_fixtures.json");
  const { fixtures } = JSON.parse(readFileSync(fixturesPath, "utf8"));

  for (const fixture of fixtures) {
    const mol = parse_smiles(fixture.smiles);
    assert.equal(mol.atom_count(), fixture.atomCount, `${fixture.name}: atom count`);

    const response = JSON.parse(embed_pipeline_v2_json(mol, JSON.stringify(fixture.config)));
    assert.equal(response.ok, fixture.ok, `${fixture.name}: ok`);

    if (fixture.ok) {
      const { result } = response;
      assert.equal(result.coords.length, fixture.coordsLength, `${fixture.name}: coords length`);
      assert.equal(
        result.stereoBefore.nDeclared,
        fixture.stereoBeforeDeclared,
        `${fixture.name}: stereoBefore.nDeclared`,
      );
      assert.equal(
        result.finalStereo.nDeclared,
        fixture.stereoAfterDeclared,
        `${fixture.name}: finalStereo.nDeclared`,
      );
      assert.equal(
        result.finalStereo.nViolations,
        fixture.stereoAfterViolations,
        `${fixture.name}: finalStereo.nViolations`,
      );
      assert.equal(
        result.forceField.requestedForceField,
        fixture.forceFieldRequested,
        `${fixture.name}: forceField.requestedForceField`,
      );
      assert.equal(
        result.forceField.actualForceFieldUsed,
        fixture.forceFieldActual,
        `${fixture.name}: forceField.actualForceFieldUsed`,
      );
      assert.equal(
        result.forceField.fallbackReason !== null,
        fixture.hasFallback,
        `${fixture.name}: forceField fallback presence`,
      );
      assert.equal(result.finalValidation.sound, fixture.sound, `${fixture.name}: finalValidation.sound`);
    } else {
      assert.equal(response.error.stage, fixture.stage, `${fixture.name}: error.stage`);
      assert.equal(response.error.cause.kind, fixture.causeKind, `${fixture.name}: error.cause.kind`);
    }
  }
  console.log(
    `cross-binding parity vs. Python reference (${fixtures.length} fixtures): all structural fields match`,
  );
}

console.log("pipeline_v2.test.mjs: all assertions passed");
