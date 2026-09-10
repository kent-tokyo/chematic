#!/usr/bin/env node

/** Run the WASM comparison gate in a real Chromium page without Playwright. */

import { createHash } from "node:crypto";
import { createServer } from "node:http";
import { mkdtempSync, readFileSync, statSync, writeFileSync } from "node:fs";
import { basename, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { spawnSync } from "node:child_process";

const ROOT = resolve(join(resolve(fileURLToPath(import.meta.url), ".."), ".."));
const DEFAULT_CORPUS = join(ROOT, "scripts", "descriptor_census_corpus.smi");

function option(args, name, fallback = undefined) {
  const index = args.indexOf(name);
  return index >= 0 ? args[index + 1] : fallback;
}

function required(args, name) {
  const value = option(args, name);
  if (!value) throw new Error(`${name} is required`);
  return value;
}

function sha256(value) {
  return createHash("sha256").update(value).digest("hex");
}

function percentile(values, fraction) {
  const ordered = [...values].sort((a, b) => a - b);
  const position = (ordered.length - 1) * fraction;
  const low = Math.floor(position);
  const high = Math.min(low + 1, ordered.length - 1);
  return ordered[low] + (ordered[high] - ordered[low]) * (position - low);
}

function summary(values) {
  return {
    count: values.length,
    p50_ms: Number(percentile(values, 0.5).toFixed(6)),
    p95_ms: Number(percentile(values, 0.95).toFixed(6)),
    mean_ms: Number((values.reduce((a, b) => a + b, 0) / values.length).toFixed(6)),
  };
}

function html({ corpus, schematicDir, rdkitDir, warmup, initTimeoutMs }) {
  const corpusJson = JSON.stringify(corpus).replace(/</g, "\\u003c");
  return `<!doctype html>
<meta charset="utf-8">
<pre id="result">running</pre>
<script src="/rdkit/RDKit_minimal.js"></script>
<script type="module">
import initSchematic, { parse_smiles, rdkit_ecfp4_bitvec } from "/schematic/chematic_wasm.js";
const smiles = ${corpusJson};
const warmup = ${warmup};
const initTimeoutMs = ${initTimeoutMs};
const out = document.querySelector("#result");
const withTimeout = (promise, label) => Promise.race([
  promise,
  new Promise((_, reject) => setTimeout(() => reject(new Error(label + " timed out after " + initTimeoutMs + " ms")), initTimeoutMs)),
]);
const digest = async (value) => {
  const bytes = typeof value === "string" ? new TextEncoder().encode(value) : value;
  const hash = await crypto.subtle.digest("SHA-256", bytes);
  return [...new Uint8Array(hash)].map((v) => v.toString(16).padStart(2, "0")).join("");
};
const packedBits = (bytes) => {
  let result = "";
  for (const byte of bytes) for (let bit = 0; bit < 8; bit += 1) result += (byte >> bit) & 1;
  return result;
};
const summary = (values) => {
  const ordered = [...values].sort((a, b) => a - b);
  const percentile = (fraction) => {
    const position = (ordered.length - 1) * fraction;
    const low = Math.floor(position);
    const high = Math.min(low + 1, ordered.length - 1);
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low);
  };
  return { count: values.length, p50_ms: percentile(0.5), p95_ms: percentile(0.95), mean_ms: values.reduce((a, b) => a + b, 0) / values.length };
};
const time = async (values, fn) => {
  for (const value of values.slice(0, warmup)) await fn(value);
  const samples = [];
  for (const value of values) {
    const start = performance.now();
    await fn(value);
    samples.push(performance.now() - start);
  }
  return summary(samples);
};
const run = async () => {
  const schematicStart = performance.now();
  await withTimeout(initSchematic("/schematic/chematic_wasm_bg.wasm"), "chematic WASM initialization");
  const schematicInit = performance.now() - schematicStart;
  const schematicParse = await time(smiles, (value) => { const mol = parse_smiles(value); mol.free(); });
  const schematicWrite = await time(smiles, (value) => { const mol = parse_smiles(value); mol.canonical_smiles(); mol.free(); });
  const schematicFp = await time(smiles, (value) => { const mol = parse_smiles(value); rdkit_ecfp4_bitvec(mol); mol.free(); });
  const schematicHashes = [];
  for (const value of smiles) { const mol = parse_smiles(value); schematicHashes.push(await digest(packedBits(rdkit_ecfp4_bitvec(mol)))); mol.free(); }

  const rdkitStart = performance.now();
  const RDKit = await withTimeout(window.initRDKitModule({ locateFile: (name) => name.endsWith(".wasm") ? "/rdkit/RDKit_minimal.wasm" : name }), "RDKit WASM initialization");
  const rdkitInit = performance.now() - rdkitStart;
  const getMol = (value) => { const mol = RDKit.get_mol(value); if (!mol) throw new Error("RDKit rejected " + value); return mol; };
  const rdkitParse = await time(smiles, (value) => { getMol(value).delete(); });
  const rdkitWrite = await time(smiles, (value) => { const mol = getMol(value); mol.get_smiles(); mol.delete(); });
  const rdkitFp = await time(smiles, (value) => { const mol = getMol(value); mol.get_morgan_fp(JSON.stringify({ radius: 2, nBits: 2048 })); mol.delete(); });
  const rdkitHashes = [];
  for (const value of smiles) { const mol = getMol(value); rdkitHashes.push(await digest(mol.get_morgan_fp(JSON.stringify({ radius: 2, nBits: 2048 })))); mol.delete(); }
  out.textContent = JSON.stringify({
    environment: { user_agent: navigator.userAgent, platform: navigator.platform },
    schematic: { init_ms: schematicInit, parse: schematicParse, smiles_write: schematicWrite, rdkit_compatible_ecfp4: schematicFp },
    rdkit: { version: RDKit.version(), init_ms: rdkitInit, parse: rdkitParse, smiles_write: rdkitWrite, morgan_radius2_2048: rdkitFp },
    fingerprint_parity: { compared_rows: smiles.length, exact_matches: schematicHashes.filter((v, i) => v === rdkitHashes[i]).length },
  });
};
run().catch((error) => { out.textContent = JSON.stringify({ error: String(error), stack: error.stack }); });
</script>`;
}

function serve(files) {
  const server = createServer((request, response) => {
    const routes = {
      "/": { path: files.html, type: "text/html" },
      "/schematic/chematic_wasm.js": { path: join(files.schematic, "chematic_wasm.js"), type: "text/javascript" },
      "/schematic/chematic_wasm_bg.wasm": { path: join(files.schematic, "chematic_wasm_bg.wasm"), type: "application/wasm" },
      "/rdkit/RDKit_minimal.js": { path: join(files.rdkit, "RDKit_minimal.js"), type: "text/javascript" },
      "/rdkit/RDKit_minimal.wasm": { path: join(files.rdkit, "RDKit_minimal.wasm"), type: "application/wasm" },
    };
    const route = routes[request.url];
    if (!route) { response.writeHead(404); response.end(); return; }
    response.writeHead(200, { "Content-Type": route.type, "Cache-Control": "no-store" });
    response.end(readFileSync(route.path));
  });
  return server;
}

async function main() {
  const args = process.argv.slice(2);
  if (args.includes("--help") || args.includes("-h")) {
    console.log("Usage: node scripts/bench_browser_wasm_vs_rdkit.mjs --rdkit-package PATH --output PATH [--rows N] [--warmup N] [--chromium PATH] [--browser-timeout-ms N]");
    return;
  }
  const corpusPath = resolve(option(args, "--corpus", DEFAULT_CORPUS));
  const rows = Number(option(args, "--rows", "1000"));
  const warmup = Number(option(args, "--warmup", "20"));
  const initTimeoutMs = Number(option(args, "--init-timeout-ms", "10000"));
  const corpusBytes = readFileSync(corpusPath);
  const corpus = corpusBytes.toString("utf8").split(/\r?\n/).map((line) => line.trim()).filter(Boolean).slice(0, rows);
  if (corpus.length !== rows) throw new Error(`corpus contains ${corpus.length} rows, expected ${rows}`);
  const schematic = resolve(option(args, "--schematic-dir", join(ROOT, "demo", "pkg")));
  const rdkitPackage = resolve(required(args, "--rdkit-package"));
  const rdkit = join(rdkitPackage, "dist");
  if (!statSync(join(rdkit, "RDKit_minimal.js")).isFile()) throw new Error(`RDKit dist not found: ${rdkit}`);
  const output = resolve(required(args, "--output"));
  const htmlPath = "/private/tmp/chematic-browser-benchmark-runner.html";
  writeFileSync(htmlPath, html({ corpus, schematicDir: schematic, rdkitDir: rdkit, warmup, initTimeoutMs }));
  const server = serve({ html: htmlPath, schematic, rdkit });
  await new Promise((resolveServer) => server.listen(0, "127.0.0.1", resolveServer));
  const port = server.address().port;
  try {
    const chromium = option(args, "--chromium", "/opt/homebrew/bin/chromium");
    const browserTimeoutMs = Number(option(args, "--browser-timeout-ms", "60000"));
    const profile = mkdtempSync("/private/tmp/chematic-browser-gate-profile-");
    const browser = spawnSync(chromium, ["--headless=new", "--no-sandbox", "--disable-gpu", "--disable-dev-shm-usage", "--no-first-run", `--user-data-dir=${profile}`, "--virtual-time-budget=30000", "--dump-dom", `http://127.0.0.1:${port}/`], { encoding: "utf8", maxBuffer: 10 * 1024 * 1024, timeout: browserTimeoutMs });
    if (browser.error) throw browser.error;
    if (browser.status !== 0) throw new Error(browser.stderr || "Chromium exited unsuccessfully");
    const marker = browser.stdout.match(/<pre id="result">([\s\S]*?)<\/pre>/);
    if (!marker) throw new Error(`Chromium did not return a result: ${browser.stdout.slice(-1000)}`);
    const measurement = JSON.parse(marker[1].replaceAll("&quot;", '"').replaceAll("&amp;", "&"));
    if (measurement.error) throw new Error(`${measurement.error}\n${measurement.stack || ""}`);
    const result = {
      schema_version: 1,
      gate: "wasm-vs-rdkit-browser-chromium",
      corpus: { path: corpusPath.startsWith(`${ROOT}/`) ? corpusPath.slice(ROOT.length + 1) : corpusPath, rows, sha256: sha256(corpusBytes) },
      configuration: { warmup_rows: warmup, browser: "Chromium headless", timing: "same browser page, sequential lanes" },
      ...measurement,
      interpretation: "Browser measurements exclude peak RSS and are not equivalent to Node process measurements; they are a Chromium-specific WASM/API comparison.",
    };
    writeFileSync(output, `${JSON.stringify(result, null, 2)}\n`);
    console.log(JSON.stringify(result, null, 2));
  } finally {
    server.close();
  }
}

main().catch((error) => { console.error(error.stack || error.message || error); process.exitCode = 1; });
