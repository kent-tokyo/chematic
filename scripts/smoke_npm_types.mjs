#!/usr/bin/env node
/** Type-check a generated WASM package as an isolated npm consumer. */

import assert from "node:assert/strict";
import { cp, mkdtemp, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { resolve } from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const [packageDir, tscCommand] = process.argv.slice(2);
if (!packageDir || !tscCommand || process.argv.length !== 4) {
  console.error("usage: node scripts/smoke_npm_types.mjs <package-directory> <tsc-command>");
  process.exit(2);
}

const root = resolve(fileURLToPath(new URL("..", import.meta.url)));
const fixture = resolve(root, "scripts/tests/fixtures/npm_typing_smoke.ts");
const consumer = await mkdtemp(resolve(tmpdir(), "chematic-npm-types-"));
const installed = resolve(consumer, "node_modules/@kent-tokyo/chematic");

try {
  await cp(resolve(packageDir), installed, { recursive: true, errorOnExist: true });
  await cp(fixture, resolve(consumer, "npm_typing_smoke.ts"));
  await writeFile(
    resolve(consumer, "tsconfig.json"),
    JSON.stringify(
      {
        compilerOptions: {
          target: "ES2022",
          module: "NodeNext",
          moduleResolution: "NodeNext",
          strict: true,
          noEmit: true,
          skipLibCheck: false,
          lib: ["ES2022", "DOM", "ESNext.Disposable"],
        },
        files: ["npm_typing_smoke.ts"],
      },
      null,
      2,
    ) + "\n",
    "utf8",
  );
  const result = spawnSync(resolve(tscCommand), ["--project", "tsconfig.json"], {
    cwd: consumer,
    encoding: "utf8",
  });
  assert.equal(result.status, 0, `TypeScript consumer check failed:\n${result.stdout}${result.stderr}`);
  console.log("generated npm package TypeScript consumer check passed");
} finally {
  await rm(consumer, { recursive: true, force: true });
}
