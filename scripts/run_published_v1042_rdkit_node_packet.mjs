#!/usr/bin/env node

/** Reproduce the v1.0.42 Node/WASM lanes that have an equivalent RDKit.js API. */

import { createHash } from "node:crypto";
import { readFileSync, readdirSync, writeFileSync } from "node:fs";
import { basename, join, resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { performance } from "node:perf_hooks";
import os from "node:os";

const OPS = ["labute_asa", "morgan2_chiral"];

function option(args, name, fallback = undefined) {
  const index = args.indexOf(name);
  return index >= 0 ? args[index + 1] : fallback;
}

function required(args, name) {
  const value = option(args, name);
  if (!value) throw new Error(`${name} is required`);
  return value;
}

function digest(value) {
  return createHash("sha256").update(value).digest("hex");
}

function artifactRecords(directory) {
  return readdirSync(directory, { withFileTypes: true })
    .filter((entry) => entry.isFile())
    .map((entry) => {
      const path = join(directory, entry.name);
      const bytes = readFileSync(path);
      return { file: entry.name, bytes: bytes.length, sha256: digest(bytes) };
    })
    .sort((left, right) => left.file.localeCompare(right.file));
}

function rowsFrom(path, limit = undefined) {
  const rows = readFileSync(path, "utf8")
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter(Boolean)
    .map((line) => line.split(/\s+/)[0]);
  return limit === undefined ? rows : rows.slice(0, limit);
}

function packedBits(bytes) {
  let result = "";
  for (const byte of bytes) {
    for (let bit = 0; bit < 8; bit += 1) result += (byte >> bit) & 1;
  }
  return result;
}

function quantile(sorted, fraction) {
  const index = Math.floor((sorted.length - 1) * fraction);
  return sorted[index];
}

function makeRng(seed) {
  let state = seed >>> 0;
  return () => {
    state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
    return state / 0x100000000;
  };
}

function pairedRatioCi(numerator, denominator, seed) {
  const ratios = numerator.map((value, index) => value / denominator[index]);
  const rng = makeRng(seed);
  const estimates = [];
  for (let sample = 0; sample < 10000; sample += 1) {
    const selected = [];
    for (let index = 0; index < ratios.length; index += 1) {
      selected.push(ratios[Math.floor(rng() * ratios.length)]);
    }
    selected.sort((a, b) => a - b);
    estimates.push(selected[Math.floor(selected.length / 2)]);
  }
  estimates.sort((a, b) => a - b);
  const ordered = [...ratios].sort((a, b) => a - b);
  return {
    median: ordered[Math.floor(ordered.length / 2)],
    ci95: [quantile(estimates, 0.025), quantile(estimates, 0.975)],
    winning_blocks: ratios.filter((value) => value > 1).length,
    blocks: ratios.length,
  };
}

async function loadEngines(schematicDir, rdkitDir) {
  const schematicPackage = JSON.parse(readFileSync(join(schematicDir, "package.json"), "utf8"));
  const schematic = await import(pathToFileURL(join(schematicDir, "chematic_wasm.js")));
  const schematicWasm = join(schematicDir, "chematic_wasm_bg.wasm");
  if (typeof schematic.initSync === "function") {
    schematic.initSync({ module: readFileSync(schematicWasm) });
  } else if (typeof schematic.default === "function") {
    await schematic.default(pathToFileURL(schematicWasm));
  } else {
    throw new Error("published chematic package has no WASM initializer");
  }
  const rdkitPackage = JSON.parse(readFileSync(join(rdkitDir, "package.json"), "utf8"));
  const imported = await import(pathToFileURL(join(rdkitDir, "dist", "RDKit_minimal.js")));
  const init = imported.default ?? imported.initRDKitModule;
  const rdkit = await init({
    locateFile: (name) => name.endsWith(".wasm")
      ? join(rdkitDir, "dist", "RDKit_minimal.wasm")
      : name,
  });
  if (schematicPackage.version !== "1.0.42" || rdkitPackage.version !== "2026.9.1") {
    throw new Error(`expected npm 1.0.42 / 2026.9.1, got ${schematicPackage.version} / ${rdkitPackage.version}`);
  }
  return { schematic, rdkit, schematicPackage, rdkitPackage };
}

function adapters(engines) {
  const ch = {
    parse: (smiles) => engines.schematic.parse_smiles(smiles),
    free: (mol) => mol.free(),
    labute_asa: (mol) => mol.labute_asa(),
    morgan2_chiral: (mol) => engines.schematic.rdkit_ecfp_config_chiral_bitvec(mol, 2, 2048),
  };
  const rd = {
    parse: (smiles) => {
      const mol = engines.rdkit.get_mol(smiles);
      if (!mol) throw new Error(`RDKit rejected ${smiles}`);
      return mol;
    },
    free: (mol) => mol.delete(),
    labute_asa: (mol) => JSON.parse(mol.get_descriptors()).labuteASA,
    morgan2_chiral: (mol) => mol.get_morgan_fp(JSON.stringify({
      radius: 2, nBits: 2048, useChirality: true,
    })),
  };
  return { chematic: ch, rdkit: rd };
}

function compareOutputs(rows, engines) {
  const api = adapters(engines);
  const report = Object.fromEntries(OPS.map((op) => [op, {
    compared: 0, matches: 0, mismatch_count: 0, error_count: 0,
    mismatches: [], errors: [], equivalent: false,
  }]));
  rows.forEach((smiles, row) => {
    let chMol;
    let rdMol;
    try {
      chMol = api.chematic.parse(smiles);
      rdMol = api.rdkit.parse(smiles);
      for (const op of OPS) {
        try {
          const actualRaw = api.chematic[op](chMol);
          const expected = api.rdkit[op](rdMol);
          const actual = op === "morgan2_chiral" ? packedBits(actualRaw) : actualRaw;
          const matches = op === "labute_asa"
            ? Math.abs(actual - expected) <= 1e-9
            : actual === expected;
          report[op].compared += 1;
          report[op].matches += Number(matches);
          if (!matches && report[op].mismatches.length < 50) {
            report[op].mismatches.push({ row, smiles, chematic: actual, rdkit: expected });
          }
        } catch (error) {
          report[op].errors.push({ row, error: String(error) });
        }
      }
    } catch (error) {
      for (const op of OPS) report[op].errors.push({ row, error: String(error) });
    } finally {
      if (chMol) api.chematic.free(chMol);
      if (rdMol) api.rdkit.free(rdMol);
    }
  });
  for (const op of OPS) {
    report[op].mismatch_count = report[op].compared - report[op].matches;
    report[op].error_count = report[op].errors.length;
    report[op].equivalent = report[op].compared === rows.length
      && report[op].matches === rows.length
      && report[op].error_count === 0;
  }
  return report;
}

function measureBlock(rows, api, block) {
  const result = {};
  const ordered = OPS.slice(block % OPS.length).concat(OPS.slice(0, block % OPS.length));
  for (const op of ordered) {
    const prepared = rows.map((smiles) => api.parse(smiles));
    let checksum = 0;
    let started = performance.now();
    for (const mol of prepared) {
      const value = api[op](mol);
      checksum ^= typeof value === "number" ? Math.trunc(value * 1e6) : value.length;
    }
    result[`${op}_prepared_seconds`] = (performance.now() - started) / 1000;
    result[`${op}_prepared_checksum`] = checksum;
    prepared.forEach((mol) => api.free(mol));

    checksum = 0;
    started = performance.now();
    for (const smiles of rows) {
      const mol = api.parse(smiles);
      const value = api[op](mol);
      checksum ^= typeof value === "number" ? Math.trunc(value * 1e6) : value.length;
      api.free(mol);
    }
    result[`${op}_pipeline_seconds`] = (performance.now() - started) / 1000;
    result[`${op}_pipeline_checksum`] = checksum;
  }
  return result;
}

async function main() {
  const args = process.argv.slice(2);
  const corpus = resolve(required(args, "--corpus"));
  const artifactDir = resolve(required(args, "--artifact-dir"));
  const schematicDir = resolve(required(args, "--schematic-dir"));
  const rdkitDir = resolve(required(args, "--rdkit-dir"));
  const output = resolve(required(args, "--output"));
  const blocks = Number(option(args, "--blocks", "21"));
  const limit = Number(option(args, "--limit", "1000"));
  const engines = await loadEngines(schematicDir, rdkitDir);
  const api = adapters(engines);
  const allRows = rowsFrom(corpus);
  const timingRows = allRows.slice(0, limit);
  const accuracy = compareOutputs(allRows, engines);
  const samples = { chematic: [], rdkit: [] };
  for (let block = 0; block < blocks; block += 1) {
    const order = block % 2 === 0 ? ["chematic", "rdkit"] : ["rdkit", "chematic"];
    for (const engine of order) samples[engine].push(measureBlock(timingRows, api[engine], block));
    console.error(`published Node block ${block + 1}/${blocks}`);
  }
  const speed = {};
  const equivalentWork = { labute_asa: false, morgan2_chiral: true };
  OPS.forEach((op, index) => {
    const eligible = accuracy[op].equivalent && equivalentWork[op];
    speed[op] = {
      eligible,
      ineligible_reason: !accuracy[op].equivalent
        ? "10000-row output gate did not pass"
        : (!equivalentWork[op]
          ? "RDKit.js exposes Labute ASA only through get_descriptors(), which computes a descriptor bundle"
          : null),
    };
    for (const lane of ["prepared", "pipeline"]) {
      const key = `${op}_${lane}_seconds`;
      speed[op][lane] = eligible
        ? pairedRatioCi(
          samples.rdkit.map((sample) => sample[key]),
          samples.chematic.map((sample) => sample[key]),
          0x1042 + index * 2 + Number(lane === "pipeline"),
        )
        : null;
    }
  });
  const corpusBytes = readFileSync(corpus);
  const report = {
    schema: "published-v1.0.42-rdkit-packet/v1",
    binding: "node-wasm",
    versions: {
      chematic: engines.schematicPackage.version,
      rdkit: engines.rdkitPackage.version,
      rdkit_runtime: engines.rdkit.version(),
    },
    artifacts: artifactRecords(artifactDir),
    corpus: {
      path: basename(corpus), rows: allRows.length, timing_rows: limit,
      sha256: digest(corpusBytes),
    },
    environment: {
      node: process.version, platform: process.platform, arch: process.arch,
      cpus: os.cpus().length, runner_os: process.env.RUNNER_OS ?? null,
      runner_arch: process.env.RUNNER_ARCH ?? null,
      runner_name: process.env.RUNNER_NAME ?? null,
    },
    method: {
      blocks,
      order: "alternating engines in one Node process for each block",
      prepared: "parse outside timing, then first operation call",
      pipeline: "parse plus operation inside timing",
      confidence_interval: "paired bootstrap of per-block RDKit/chematic ratios, 10000 resamples",
      speed_claim_policy: "withhold both intervals unless the 10000-row output gate passes and both APIs perform equivalent work",
      timed_result: "each operation returns its native published-package value; bit normalization is accuracy-only",
      unavailable_v1042_lanes: {
        tpsa: "npm v1.0.42 exposes native TPSA, not the Python rdkit_tpsa profile",
        num_rings: "npm v1.0.42 exposes ring_count, not the Python rdkit_num_rings profile",
      },
    },
    accuracy,
    raw_samples_seconds: samples,
    speedup_rdkit_over_chematic: speed,
  };
  writeFileSync(output, `${JSON.stringify(report, null, 2)}\n`);
}

main().catch((error) => {
  console.error(error.stack || error);
  process.exitCode = 1;
});
