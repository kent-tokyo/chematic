// Compare the Explorer Worker boundary with native scalar and batch APIs.
//
// The browser side is intentionally exercised through the real module Worker;
// importing wasm-bindgen directly in Node would skip structured cloning,
// initialization and handle-release behavior that users depend on.

import assert from "node:assert/strict";
import { spawnSync } from "node:child_process";
import { chromium } from "playwright";

const cli = process.env.CHEMATIC_CLI ?? "target/release/chematic";
const records = [
  { name: "ethanol", smiles: "CCO" },
  { name: "invalid", smiles: "C1CC" },
  { name: "ethylamine", smiles: "CCN" },
];
const numericFields = new Map([
  ["mw", "molecular_weight"],
  ["exactMass", "exact_mass"],
  ["tpsa", "tpsa"],
  ["logP", "logp"],
  ["hbd", "hbd"],
  ["hba", "hba"],
  ["rotatableBonds", "rotatable_bonds"],
  ["heavyAtomCount", "heavy_atoms"],
]);

function invoke(args, input) {
  const completed = spawnSync(cli, args, { encoding: "utf8", input });
  if (completed.error) throw completed.error;
  return completed;
}

function jsonCommand(args, input) {
  const completed = invoke(args, input);
  assert.equal(completed.status, 0, completed.stderr);
  return JSON.parse(completed.stdout);
}

function approximatelyEqual(actual, expected, field) {
  assert.ok(Math.abs(actual - expected) <= 1e-9, `${field}: ${actual} !== ${expected}`);
}

const batch = jsonCommand(["batch-descriptors"], records.map((record) => record.smiles).join("\n") + "\n");
assert.equal(batch.status, "complete");
assert.equal(batch.records.length, records.length);

const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage();
  const errors = [];
  page.on("pageerror", (error) => errors.push(String(error)));
  page.on("console", (message) => {
    if (message.type() === "error") errors.push(message.text());
  });
  await page.goto("http://127.0.0.1:8765/explorer/index.html?native-worker-parity=1", {
    waitUntil: "networkidle",
  });
  await page.locator("#loading-overlay").waitFor({ state: "hidden" });
  const workerRecords = await page.evaluate(async (workerInput) => {
    const worker = new Worker("./worker.js", { type: "module" });
    let nextId = 1;
    const request = (type, payload = {}) => new Promise((resolve, reject) => {
      const requestId = nextId++;
      const timer = setTimeout(() => reject(new Error(`worker ${type} timed out`)), 15_000);
      worker.addEventListener("message", function listener({ data }) {
        if (data.requestId !== requestId) return;
        worker.removeEventListener("message", listener);
        clearTimeout(timer);
        if (data.type === "error") reject(new Error(data.message));
        else resolve(data);
      });
      worker.postMessage({ type, requestId, ...payload });
    });
    try {
      await request("init");
      return (await request("parse", { records: workerInput, startIndex: 0 })).records;
    } finally {
      worker.terminate();
    }
  }, records);

  assert.equal(workerRecords.length, batch.records.length);
  for (const [index, workerRecord] of workerRecords.entries()) {
    const nativeRecord = batch.records[index];
    assert.equal(workerRecord.inputSmiles, nativeRecord.input_smiles);
    if (nativeRecord.error) {
      assert.equal(workerRecord.status, "error", `record ${index} must preserve native rejection`);
      continue;
    }
    assert.equal(workerRecord.status, "ok");
    assert.equal(workerRecord.canonicalSmiles, nativeRecord.descriptors.smiles);
    for (const [workerField, nativeField] of numericFields) {
      approximatelyEqual(workerRecord.descriptors[workerField], nativeRecord.descriptors[nativeField], `${index}:${workerField}`);
    }
    const scalar = jsonCommand(["descriptors", workerRecord.inputSmiles]);
    assert.equal(scalar.smiles, workerRecord.canonicalSmiles);
    for (const [workerField, nativeField] of numericFields) {
      approximatelyEqual(workerRecord.descriptors[workerField], scalar[nativeField], `scalar ${index}:${workerField}`);
    }
  }
  assert.deepEqual(errors, []);
} finally {
  await browser.close();
}

console.log("Explorer native scalar/batch/Worker parity passed");
