#!/usr/bin/env node

/**
 * Compare the generated chematic WASM package with an installed
 * @rdkit/rdkit package under one Node.js process contract.
 *
 * This is a measurement harness, not a claim that the two libraries expose
 * identical chemistry. Fingerprint parity is checked only through the
 * explicitly named RDKit-compatible Morgan/ECFP4 APIs.
 */

import { gzipSync } from "node:zlib";
import { createHash } from "node:crypto";
import { readFileSync, statSync, writeFileSync } from "node:fs";
import { basename, dirname, extname, join, resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import { performance } from "node:perf_hooks";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const DEFAULT_CORPUS = join(ROOT, "scripts", "descriptor_census_corpus.smi");

function option(args, name, fallback = undefined) {
  const index = args.indexOf(name);
  return index >= 0 ? args[index + 1] : fallback;
}

function requiredOption(args, name) {
  const value = option(args, name);
  if (!value) throw new Error(`${name} is required`);
  return value;
}

function percentile(values, fraction) {
  const ordered = [...values].sort((a, b) => a - b);
  const position = (ordered.length - 1) * fraction;
  const low = Math.floor(position);
  const high = Math.min(low + 1, ordered.length - 1);
  return ordered[low] + (ordered[high] - ordered[low]) * (position - low);
}

function timingSummary(samples) {
  return {
    count: samples.length,
    p50_ms: Number((percentile(samples, 0.5)).toFixed(6)),
    p95_ms: Number((percentile(samples, 0.95)).toFixed(6)),
    mean_ms: Number((samples.reduce((a, b) => a + b, 0) / samples.length).toFixed(6)),
    min_ms: Number(Math.min(...samples).toFixed(6)),
    max_ms: Number(Math.max(...samples).toFixed(6)),
  };
}

function measure(operation, values, warmup, fn) {
  for (const value of values.slice(0, warmup)) fn(value);
  const samples = [];
  for (const value of values) {
    const start = performance.now();
    fn(value);
    samples.push(performance.now() - start);
  }
  return timingSummary(samples);
}

function fileSizes(path) {
  const bytes = readFileSync(path);
  return {
    file: basename(path),
    raw_bytes: bytes.length,
    gzip_bytes: gzipSync(bytes, { level: 9 }).length,
    sha256: digest(bytes),
  };
}

function findRdkitFiles(packagePath) {
  const candidate = resolve(packagePath);
  const directories = statSync(candidate).isDirectory() ? [candidate] : [dirname(candidate)];
  const jsCandidates = directories.flatMap((dir) => [
    join(dir, "dist", "RDKit_minimal.js"),
    join(dir, "Code", "MinimalLib", "dist", "RDKit_minimal.js"),
    join(dir, "RDKit_minimal.js"),
  ]);
  const wasmCandidates = directories.flatMap((dir) => [
    join(dir, "dist", "RDKit_minimal.wasm"),
    join(dir, "Code", "MinimalLib", "dist", "RDKit_minimal.wasm"),
    join(dir, "RDKit_minimal.wasm"),
  ]);
  const js = jsCandidates.find((path) => {
    try { return statSync(path).isFile(); } catch { return false; }
  });
  const wasm = wasmCandidates.find((path) => {
    try { return statSync(path).isFile(); } catch { return false; }
  });
  if (!js || !wasm) throw new Error(`could not locate RDKit_minimal.js/.wasm under ${candidate}`);
  const packageJsonCandidates = [join(candidate, "package.json"), join(dirname(js), "..", "package.json")];
  const packageJson = packageJsonCandidates.find((path) => {
    try { return statSync(path).isFile(); } catch { return false; }
  });
  let packageVersion = null;
  if (packageJson) {
    try { packageVersion = JSON.parse(readFileSync(packageJson, "utf8")).version ?? null; } catch { /* digest remains authoritative */ }
  }
  const declarationCandidates = [
    join(candidate, "dist", "index.d.ts"),
    join(candidate, "dist", "RDKit_minimal.d.ts"),
    join(candidate, "index.d.ts"),
  ];
  const declaration = declarationCandidates.find((path) => {
    try { return statSync(path).isFile(); } catch { return false; }
  });
  return { js, wasm, packageJson, packageVersion, declaration };
}

function packedBits(bytes) {
  let result = "";
  for (const byte of bytes) {
    for (let bit = 0; bit < 8; bit += 1) result += (byte >> bit) & 1;
  }
  return result;
}

function digest(value) {
  return createHash("sha256").update(value).digest("hex");
}

function fileDigest(path) {
  return path ? digest(readFileSync(path)) : null;
}

async function runSchematic(args, smiles, warmup) {
  const packageDir = resolve(option(args, "--schematic-dir", join(ROOT, "demo", "pkg")));
  const module = await import(pathToFileURL(join(packageDir, "chematic_wasm.js")));
  const initStart = performance.now();
  const wasmPath = join(packageDir, "chematic_wasm_bg.wasm");
  if (typeof module.initSync === "function") {
    module.initSync({ module: readFileSync(wasmPath) });
  } else if (typeof module.default === "function") {
    await module.default(pathToFileURL(wasmPath));
  } else if (typeof module.parse_smiles === "function") {
    // The wasm-pack `nodejs` target instantiates the module while it is
    // imported and intentionally exports no default initializer.  Keep this
    // path explicit so a missing or incompatible artifact cannot be mistaken
    // for a successful initialization.
  } else {
    throw new Error("schematic WASM package exports neither an initializer nor parse_smiles");
  }
  const initMs = performance.now() - initStart;
  const parse = (value) => module.parse_smiles(value);
  const parseTiming = measure("parse", smiles, warmup, (value) => {
    const mol = parse(value);
    mol.free();
  });
  const writeTiming = measure("write", smiles, warmup, (value) => {
    const mol = parse(value);
    mol.canonical_smiles();
    mol.free();
  });
  const fingerprintTiming = measure("fingerprint", smiles, warmup, (value) => {
    const mol = parse(value);
    module.rdkit_ecfp4_bitvec(mol);
    mol.free();
  });
  const fingerprintDigests = smiles.map((value) => {
    const mol = parse(value);
    const bits = packedBits(module.rdkit_ecfp4_bitvec(mol));
    mol.free();
    return digest(bits);
  });
  return {
    implementation: "chematic",
    package_dir: basename(packageDir),
    init_ms: Number(initMs.toFixed(6)),
    parse: parseTiming,
    canonical_smiles: writeTiming,
    rdkit_compatible_ecfp4: fingerprintTiming,
    fingerprint_digests: fingerprintDigests,
    peak_rss_bytes: process.resourceUsage().maxRSS * 1024,
  };
}

async function runRdkit(args, smiles, warmup) {
  const files = findRdkitFiles(requiredOption(args, "--rdkit-package"));
  const initStart = performance.now();
  const imported = await import(pathToFileURL(files.js));
  const init = imported.default ?? imported.initRDKitModule;
  if (typeof init !== "function") throw new Error("RDKit package does not export an initialization function");
  const RDKit = await init({ locateFile: (name) => name.endsWith(".wasm") ? files.wasm : name });
  const initMs = performance.now() - initStart;
  const parse = (value) => {
    const mol = RDKit.get_mol(value);
    if (!mol) throw new Error(`RDKit rejected SMILES: ${value}`);
    return mol;
  };
  const parseTiming = measure("parse", smiles, warmup, (value) => {
    const mol = parse(value);
    mol.delete();
  });
  const writeTiming = measure("write", smiles, warmup, (value) => {
    const mol = parse(value);
    mol.get_smiles();
    mol.delete();
  });
  const fingerprintTiming = measure("fingerprint", smiles, warmup, (value) => {
    const mol = parse(value);
    mol.get_morgan_fp('{"radius":2,"nBits":2048}');
    mol.delete();
  });
  const fingerprintDigests = smiles.map((value) => {
    const mol = parse(value);
    const bits = mol.get_morgan_fp('{"radius":2,"nBits":2048}');
    mol.delete();
    return digest(bits);
  });
  return {
    implementation: "rdkit",
    version: typeof RDKit.version === "function" ? RDKit.version() : "unknown",
    package_version: files.packageVersion,
    package_dir: basename(dirname(files.js)),
    wasm_file: basename(files.wasm),
    init_ms: Number(initMs.toFixed(6)),
    parse: parseTiming,
    smiles_write: writeTiming,
    morgan_radius2_2048: fingerprintTiming,
    fingerprint_digests: fingerprintDigests,
    peak_rss_bytes: process.resourceUsage().maxRSS * 1024,
  };
}

async function child(args) {
  const corpusPath = resolve(option(args, "--corpus", DEFAULT_CORPUS));
  const rows = Number(option(args, "--rows", "1000"));
  const warmup = Number(option(args, "--warmup", "20"));
  const smiles = readFileSync(corpusPath, "utf8").split(/\r?\n/).map((line) => line.trim()).filter(Boolean).slice(0, rows);
  if (smiles.length !== rows) throw new Error(`corpus contains ${smiles.length} rows, expected ${rows}`);
  const result = option(args, "--lane") === "rdkit"
    ? await runRdkit(args, smiles, warmup)
    : await runSchematic(args, smiles, warmup);
  result.corpus_rows = rows;
  result.corpus_sha256 = "see parent report";
  process.stdout.write(`${JSON.stringify(result)}\n`);
}

async function main() {
  const args = process.argv.slice(2);
  if (args.includes("--help") || args.includes("-h")) {
    console.log("Usage: node scripts/bench_wasm_vs_rdkit.mjs --rdkit-package PATH --output PATH [--rows N] [--warmup N] [--corpus PATH] [--schematic-dir PATH] [--same-process]");
    return;
  }
  if (args.includes("--child")) return child(args);
  const corpusPath = resolve(option(args, "--corpus", DEFAULT_CORPUS));
  const output = resolve(requiredOption(args, "--output"));
  const rows = Number(option(args, "--rows", "1000"));
  const warmup = Number(option(args, "--warmup", "20"));
  const corpus = readFileSync(corpusPath);
  const smiles = corpus.toString("utf8").split(/\r?\n/).map((line) => line.trim()).filter(Boolean).slice(0, rows);
  if (smiles.length !== rows) throw new Error(`corpus contains ${smiles.length} rows, expected ${rows}`);
  const base = ["--child", "--corpus", corpusPath, "--rows", String(rows), "--warmup", String(warmup)];
  const { spawnSync } = await import("node:child_process");
  const run = (lane, extra) => {
    const result = spawnSync(process.execPath, [fileURLToPath(import.meta.url), ...base, "--lane", lane, ...extra], { encoding: "utf8" });
    if (result.status !== 0) throw new Error(`${lane} lane failed:\n${result.stderr || result.stdout}`);
    return JSON.parse(result.stdout.trim().split(/\n/).at(-1));
  };
  const schematicDir = resolve(option(args, "--schematic-dir", join(ROOT, "demo", "pkg")));
  const rdkitPackage = requiredOption(args, "--rdkit-package");
  const sameProcess = args.includes("--same-process");
  const schematic = sameProcess
    ? await runSchematic(args, smiles, warmup)
    : run("chematic", ["--schematic-dir", schematicDir]);
  const rdkit = sameProcess
    ? await runRdkit(args, smiles, warmup)
    : run("rdkit", ["--rdkit-package", resolve(rdkitPackage)]);
  const schematicWasm = join(schematicDir, "chematic_wasm_bg.wasm");
  const rdkitFiles = findRdkitFiles(rdkitPackage);
  const result = {
    schema_version: 1,
    gate: sameProcess ? "wasm-vs-rdkit-node-paired" : "wasm-vs-rdkit-node",
    environment: {
      node: process.version,
      platform: process.platform,
      arch: process.arch,
    },
    corpus: { path: corpusPath.startsWith(`${ROOT}/`) ? corpusPath.slice(ROOT.length + 1) : corpusPath, rows, sha256: createHash("sha256").update(corpus).digest("hex") },
    configuration: { warmup_rows: warmup, morgan_radius: 2, fingerprint_bits: 2048, timing: sameProcess ? "both implementations loaded and measured in one Node process" : "one isolated Node process per implementation" },
    schematic,
    rdkit,
    artifacts: {
      schematic_wasm: fileSizes(schematicWasm),
      rdkit_minimal_wasm: fileSizes(rdkitFiles.wasm),
      rdkit_minimal_js: fileSizes(rdkitFiles.js),
      rdkit_package_json: rdkitFiles.packageJson ? { file: rdkitFiles.packageJson, sha256: fileDigest(rdkitFiles.packageJson) } : null,
      rdkit_typescript_declarations: rdkitFiles.declaration ? { file: rdkitFiles.declaration, sha256: fileDigest(rdkitFiles.declaration) } : null,
    },
    fingerprint_parity: {
      compared_rows: rows,
      exact_matches: schematic.fingerprint_digests.filter((hash, index) => hash === rdkit.fingerprint_digests[index]).length,
      definition: "chematic rdkit_ecfp4_bitvec vs RDKit.js Morgan radius=2, 2048-bit binary text",
    },
    interpretation: "Node/WASM measurements are implementation and artifact comparisons. They do not establish chemistry accuracy or universal browser performance.",
  };
  delete result.schematic.corpus_sha256;
  delete result.rdkit.corpus_sha256;
  delete result.schematic.fingerprint_digests;
  delete result.rdkit.fingerprint_digests;
  writeFileSync(output, `${JSON.stringify(result, null, 2)}\n`);
  process.stdout.write(`${JSON.stringify(result, null, 2)}\n`);
}

main().catch((error) => {
  console.error(error.stack || error.message || error);
  process.exitCode = 1;
});
