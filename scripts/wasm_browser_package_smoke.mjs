#!/usr/bin/env node
/**
 * Exercise the browser-facing WASM package under a deployable CSP boundary.
 *
 * This deliberately does not visit the demo UI. It serves the generated
 * `demo/pkg` package as a consumer would, imports its ESM entrypoint, loads
 * the adjacent WASM asset, and releases a molecule handle. The page has no
 * network capability beyond same-origin asset loading.
 */

import assert from "node:assert/strict";
import { createServer } from "node:http";
import { readFile } from "node:fs/promises";
import { resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { chromium, firefox, webkit } from "playwright";

const browserName = process.argv[2] ?? "chromium";
const browsers = { chromium, firefox, webkit };
assert.ok(browsers[browserName], `unknown browser: ${browserName}`);

const root = resolve(fileURLToPath(new URL("..", import.meta.url)));
const packageDir = resolve(root, "demo", "pkg");
const csp = "default-src 'self'; script-src 'self' 'wasm-unsafe-eval'; connect-src 'self'; worker-src 'self'; style-src 'none'; img-src 'none'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'";

const page = "<!doctype html><meta charset=utf-8><title>WASM package smoke</title><script type=module src=/smoke.js></script>";
const smokeModule = `
import init, { parse_smiles } from "/pkg/chematic_wasm.js";
try {
  await init("/pkg/chematic_wasm_bg.wasm");
  const molecule = parse_smiles("c1ccccc1");
  const result = { atoms: molecule.atom_count(), formula: molecule.formula() };
  molecule.free();
  window.__chematicPackageSmoke = { ok: true, result };
} catch (error) {
  window.__chematicPackageSmoke = { ok: false, error: String(error) };
}
`;

function send(response, status, type, body) {
  response.writeHead(status, {
    "Content-Type": type,
    "Content-Security-Policy": csp,
    "Cache-Control": "no-store",
  });
  response.end(body);
}

const server = createServer(async (request, response) => {
  try {
    if (request.url === "/") return send(response, 200, "text/html; charset=utf-8", page);
    if (request.url === "/smoke.js") return send(response, 200, "text/javascript; charset=utf-8", smokeModule);
    if (request.url === "/pkg/chematic_wasm.js") {
      return send(response, 200, "text/javascript; charset=utf-8", await readFile(resolve(packageDir, "chematic_wasm.js")));
    }
    if (request.url === "/pkg/chematic_wasm_bg.wasm") {
      return send(response, 200, "application/wasm", await readFile(resolve(packageDir, "chematic_wasm_bg.wasm")));
    }
    return send(response, 404, "text/plain; charset=utf-8", "not found");
  } catch (error) {
    return send(response, 500, "text/plain; charset=utf-8", String(error));
  }
});

await new Promise((resolveListen, rejectListen) => {
  server.once("error", rejectListen);
  server.listen(0, "127.0.0.1", resolveListen);
});
const address = server.address();
assert.ok(address && typeof address === "object", "server address");
const url = `http://127.0.0.1:${address.port}/`;

const browser = await browsers[browserName].launch({ headless: true });
try {
  const browserPage = await browser.newPage();
  const errors = [];
  browserPage.on("pageerror", (error) => errors.push(String(error)));
  await browserPage.goto(url, { waitUntil: "networkidle" });
  await browserPage.waitForFunction(() => Boolean(window.__chematicPackageSmoke));
  const result = await browserPage.evaluate(() => window.__chematicPackageSmoke);
  assert.deepEqual(result, { ok: true, result: { atoms: 6, formula: "C6H6" } });
  assert.deepEqual(errors, [], `browser page errors: ${errors.join("; ")}`);
  console.log(`${browserName}: CSP package import, init, parse, and free passed`);
} finally {
  await browser.close();
  await new Promise((resolveClose, rejectClose) => server.close((error) => error ? rejectClose(error) : resolveClose()));
}
