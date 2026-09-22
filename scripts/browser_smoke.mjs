import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { chromium, firefox, webkit } from "playwright";

const browserName = process.argv[2] ?? "chromium";
const browsers = { chromium, firefox, webkit };
assert.ok(browsers[browserName], `unknown browser: ${browserName}`);
const workspaceManifest = readFileSync(new URL("../Cargo.toml", import.meta.url), "utf8");
const workspaceVersion = workspaceManifest.match(/^version\s*=\s*"([^"]+)"/m)?.[1];
assert.ok(workspaceVersion, "workspace Cargo.toml must declare a version");
const expectedVersionBadge = `v${workspaceVersion}`;

const ETHANE_MOL_BLOCK = `ethane
  chematic

  2  1  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END`;

  const browser = await browsers[browserName].launch({ headless: true });
try {
  const page = await browser.newPage();
  // Keep assertions deterministic across hosts whose navigator language differs.
  await page.addInitScript(() => localStorage.setItem("chematic-lang", "en"));
  const errors = [];
  const dialogs = [];
  const wasmAssetVersions = [];
  page.on("pageerror", (error) => errors.push(String(error)));
  page.on("console", (message) => {
    if (message.type() === "error") errors.push(message.text());
  });
  page.on("request", (request) => {
    const url = new URL(request.url());
    if (url.pathname.endsWith("/pkg/chematic_wasm.js") || url.pathname.endsWith("/pkg/chematic_wasm_bg.wasm")) {
      wasmAssetVersions.push(url.searchParams.get("v"));
    }
  });
  page.on("dialog", (dialog) => {
    dialogs.push(dialog.message());
    void dialog.dismiss();
  });
  await page.goto("http://127.0.0.1:8765/index.html?browser-smoke=0.89", {
    waitUntil: "networkidle",
  });
  await page.locator("#version-badge").waitFor({ state: "visible" });
  assert.equal(await page.locator("#version-badge").innerText(), expectedVersionBadge);
  assert.ok(wasmAssetVersions.length >= 2, "page must request both WASM assets");
  assert.ok(
    wasmAssetVersions.every((version) => version === workspaceVersion),
    `WASM asset cache version must match Cargo.toml: ${wasmAssetVersions}`,
  );
  await page.locator("#smiles-input").fill("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
  await page.locator("#btn-calc").click();
  const hba = page.locator("#desc-tbody tr").filter({ hasText: "HBA" }).locator("td").nth(1);
  await hba.waitFor({ state: "visible" });
  assert.equal(await hba.innerText(), "3");

  // Expose the existing breadth of the demo through task-oriented entry points.
  await page.locator('[data-workflow="molecule"]').click();
  await page.locator("#tb-2d").waitFor({ state: "visible" });
  await page.locator('[data-workflow="similarity"]').click();
  await page.locator("#tb-sim").waitFor({ state: "visible" });
  assert.equal(await page.locator("#sim-a").inputValue(), "CC(=O)Oc1ccccc1C(=O)O");
  await page.locator('[data-workflow="report"]').click();
  await page.locator("#tb-report").waitFor({ state: "visible" });
  await page.locator("#report-output").waitFor({ state: "visible" });
  await page.locator("#btn-report-json").click();
  await page.locator("#report-json-status").waitFor({ state: "visible" });
  await page.locator('[data-workflow="reaction"]').click();
  await page.locator("#tb-rxn").waitFor({ state: "visible" });
  await page.locator('[data-workflow="formats"]').click();
  await page.locator("#tb-formats").waitFor({ state: "visible" });
  await page.locator("#quickstart-guide-summary").click();
  await page.locator("#quickstart-guide-step3").waitFor({ state: "visible" });
  await page.locator("#tc-structure").click();
  await page.locator("#tb-2d").click();
  for (const [nextSmiles, expectedHba] of [
    ["CCO", "1"],
    ["CC(=O)O", "1"],
    ["c1ccccc1", "0"],
    ["Cn1cnc2c1c(=O)n(c(=O)n2C)C", "3"],
  ]) {
    await page.locator("#smiles-input").fill(nextSmiles);
    await page.locator("#btn-calc").click();
    await page.locator("#error-desc").waitFor({ state: "hidden" });
    assert.equal(await hba.innerText(), expectedHba);
  }
  await page.locator("#smiles-input").fill("C".repeat(1_000_001));
  await page.locator("#btn-calc").click();
  await page.locator("#error-desc").waitFor({ state: "visible" });
  assert.match(await page.locator("#error-desc").innerText(), /input|size|large|atom/i);
  await page.locator("#smiles-input").fill("C1CC");
  await page.locator("#btn-calc").click();
  await page.locator("#error-desc").waitFor({ state: "visible" });
  assert.match(await page.locator("#error-desc").innerText(), /parse|invalid|SMILES/i);
  await page.locator("#smiles-input").fill("CCO");
  await page.locator("#btn-calc").click();
  await page.locator("#error-desc").waitFor({ state: "hidden" });
  await hba.waitFor({ state: "visible" });
  assert.equal(await hba.innerText(), "1");
  await page.locator("#smiles-input").fill("  CCO  ");
  await page.locator("#btn-calc").click();
  await page.locator("#error-desc").waitFor({ state: "hidden" });
  assert.equal(await hba.innerText(), "1");
  await page.locator("#smiles-input").fill(" \t\n ");
  await page.locator("#btn-calc").click();
  await page.locator("#error-desc").waitFor({ state: "visible" });
  assert.match(await page.locator("#error-desc").innerText(), /enter|empty/i);
  await page.locator("#smiles-input").fill("CCO");
  await page.locator("#btn-calc").click();
  await page.locator("#error-desc").waitFor({ state: "hidden" });
  assert.equal(await hba.innerText(), "1");
  await page.locator("#smarts-input").fill("[");
  await page.getByRole("button", { name: "Highlight", exact: true }).click();
  await page.locator("#error-desc").waitFor({ state: "visible" });
  assert.match(await page.locator("#error-desc").innerText(), /SMARTS/i);
  await page.locator("#smarts-input").fill("c1ccccc1");
  await page.getByRole("button", { name: "Highlight", exact: true }).click();
  await page.locator("#error-desc").waitFor({ state: "hidden" });
  await page.locator("#tc-reactions").click();
  await page.locator("#tb-rxn").waitFor({ state: "visible" });
  await page.locator("#tb-rxn").click();
  await page.locator("#rxn-smirks").fill("[");
  await page.locator("#rxn-reactants").fill("CCO");
  await page.locator("#btn-rxn").click();
  await page.locator("#error-rxn").waitFor({ state: "visible" });
  assert.match(await page.locator("#error-rxn").innerText(), /parse|invalid|reaction|SMILES/i);
  await page.locator("#rxn-smirks").fill("[C:1]Br.[N:2]>>[C:1][N:2]");
  await page.locator("#rxn-reactants").fill("CCBr|CN");
  await page.locator("#btn-rxn").click();
  await page.locator("#rxn-svg-wrap svg").waitFor({ state: "visible" });
  await page.locator("#error-rxn").waitFor({ state: "hidden" });
  await page.locator("#rxn-eq-input").fill("not-a-reaction");
  await page.locator("#btn-rxn-eq").click();
  await page.locator("#error-rxn-eq").waitFor({ state: "visible" });
  assert.match(await page.locator("#error-rxn-eq").innerText(), /parse|invalid|reaction|SMILES/i);
  await page.locator("#rxn-eq-input").fill("CC(=O)O.CCO>>CC(=O)OCC.O");
  await page.locator("#btn-rxn-eq").click();
  await page.locator("#rxn-eq-svg-wrap svg").first().waitFor({ state: "visible" });
  await page.locator("#error-rxn-eq").waitFor({ state: "hidden" });
  await page.locator("#tc-structure").click();
  await page.locator("#tb-2d").click();
  await page.getByRole("button", { name: "日", exact: true }).click();
  await page.getByText("記述子計算機", { exact: true }).waitFor({ state: "visible" });
  await page.getByRole("button", { name: "EN", exact: true }).click();
  await page.getByText("Descriptor Calculator", { exact: true }).waitFor({ state: "visible" });
  await page.locator("#tc-analysis").click();
  await page.locator("#tb-sim").waitFor({ state: "visible" });
  await page.locator("#tb-sim").click();
  await page.locator("#sim-a").press("ControlOrMeta+A");
  await page.locator("#sim-a").press("Backspace");
  await page.locator("#sim-b").fill("CCO");
  await page.locator("#btn-sim").click();
  await page.locator("#error-sim").waitFor({ state: "visible" });
  assert.match(await page.locator("#error-sim").innerText(), /enter|empty/i);
  await page.locator("#sim-a").fill("CCO");
  await page.locator("#sim-b").fill("C1CC");
  await page.locator("#btn-sim").click();
  await page.locator("#error-sim").waitFor({ state: "visible" });
  assert.match(await page.locator("#error-sim").innerText(), /B:.*parse|invalid|SMILES/i);
  await page.locator("#sim-b").fill("CCN");
  await page.locator("#btn-sim").click();
  await page.locator("#error-sim").waitFor({ state: "hidden" });
  await page.locator("#sim-svgs").waitFor({ state: "visible" });
  await page.locator("#tc-structure").click();
  await page.locator("#tb-2d").click();
  await page.locator("#smiles-input").fill("CCCC");
  await page.locator("#btn-calc").click();
  await page.locator("#error-desc").waitFor({ state: "hidden" });
  await page.locator("#tc-labs").click();
  await page.locator("#tb-dynamics").waitFor({ state: "visible" });
  await page.locator("#tb-dynamics").click();
  await page.locator("#btn-torsion-scan").click();
  await page.locator("#torsion-scan-result").waitFor({ state: "visible" });
  assert.ok(Number.isFinite(Number(await page.locator("#scan-min").innerText())));
  assert.ok(Number.isFinite(Number(await page.locator("#scan-max").innerText())));
  assert.deepEqual(dialogs, []);
  await page.locator("#tc-structure").click();
  await page.locator("#tb-2d").waitFor({ state: "visible" });
  await page.locator("#tb-2d").click();
  const sdfInput = page.locator("#sdf-input");
  const sdfError = page.locator("#error-sdf");
  await sdfInput.fill("not-a-valid-mol");
  await page.locator("#btn-sdf-load").click();
  await sdfError.waitFor({ state: "visible" });
  assert.match(await sdfError.innerText(), /No valid|invalid|parse/i);
  await sdfInput.fill("");
  await page.locator("#btn-sdf-load").click();
  await sdfError.waitFor({ state: "visible" });
  assert.match(await sdfError.innerText(), /enter|empty/i);
  await sdfInput.fill(`${ETHANE_MOL_BLOCK}$$$$\n${ETHANE_MOL_BLOCK}$$$$`);
  await page.locator("#btn-sdf-load").click();
  await page.locator("#sdf-grid-output svg").waitFor({ state: "visible" });
  await sdfError.waitFor({ state: "hidden" });
  assert.equal(await page.locator("#sdf-grid-output svg").count(), 1);
  const sdfHba = page.locator("#desc-tbody tr").filter({ hasText: "HBA" }).locator("td").nth(1);
  await sdfHba.waitFor({ state: "visible" });
  assert.equal(await sdfHba.innerText(), "0");
  await sdfInput.fill("not-a-valid-mol");
  await page.locator("#btn-sdf-load").click();
  await sdfError.waitFor({ state: "visible" });
  assert.match(await sdfError.innerText(), /No valid|invalid|parse/i);
  await page.waitForFunction(
    () => typeof window.__browserSmoke?.sdfToSmilesJson === "function",
  );
  const sdfBoundaryResults = await page.evaluate(() =>
    [1_000_000, 1_000_001, 1_000_001].map((size) =>
      window.__browserSmoke.sdfToSmilesJson("x".repeat(size)),
    ),
  );
  assert.doesNotMatch(sdfBoundaryResults[0], /SDF input too large/i);
  assert.match(sdfBoundaryResults[1], /SDF input too large|input.*large|size/i);
  assert.equal(sdfBoundaryResults[1], sdfBoundaryResults[2]);
  await sdfInput.fill(ETHANE_MOL_BLOCK);
  await page.locator("#btn-sdf-load").click();
  await page.locator("#sdf-grid-output svg").waitFor({ state: "visible" });
  await sdfError.waitFor({ state: "hidden" });
  assert.equal(await sdfHba.innerText(), "0");

  // Sharing uses a bounded URL fragment and restores the molecule after a
  // fresh page load; the structure itself stays client-side.
  await page.locator("#tc-structure").click();
  await page.locator("#tb-2d").click();
  await page.locator("#smiles-input").fill("CCO");
  await page.locator("#btn-calc").click();
  await page.locator("#error-desc").waitFor({ state: "hidden" });
  await page.locator("#btn-share-molecule").click();
  await page.waitForFunction(() => location.hash.startsWith("#smiles="));
  assert.match(page.url(), /#smiles=CCO$/);
  await page.reload({ waitUntil: "networkidle" });
  await page.locator("#version-badge").waitFor({ state: "visible" });
  assert.equal(await page.locator("#smiles-input").inputValue(), "CCO");

  // Malformed or oversized fragments must fail closed to the normal example,
  // not prevent the WASM app from booting.
  await page.goto("http://127.0.0.1:8765/index.html#smiles=%E0%A4%A", {
    waitUntil: "networkidle",
  });
  await page.locator("#version-badge").waitFor({ state: "visible" });
  assert.doesNotMatch(await page.locator("#smiles-input").inputValue(), /%E0|%A4/);
  await page.goto(`http://127.0.0.1:8765/index.html#smiles=${"C".repeat(20001)}`, {
    waitUntil: "networkidle",
  });
  await page.locator("#version-badge").waitFor({ state: "visible" });
  assert.ok((await page.locator("#smiles-input").inputValue()).length <= 20000);

  await page.getByRole("tab", { name: "Data & Formats", exact: true }).click();
  const formatCases = [
    "Gaussian Cube",
    "OpenDX",
    "mmCIF",
    "PQR",
    "QCSchema",
    "ORCA Input",
    "ORCA Output",
    "LAMMPS Data",
    "LAMMPS Dump",
  ];
  const formatInput = page.locator("#formats-input");
  const formatError = page.locator("#error-formats");
  const formatOutput = page.locator("#formats-output");
  const loadExample = page.getByRole("button", { name: "Load Example", exact: true });
  const parseFormat = page.getByRole("button", { name: "Parse", exact: true });

  // Switching away from a built-in example must replace it with a valid example
  // for the selected parser. Previously the Cube text stayed in the textarea,
  // so ORCA, QCSchema, and LAMMPS parsing failed immediately after a tab click.
  await page.getByRole("button", { name: "Gaussian Cube", exact: true }).click();
  await loadExample.click();
  for (const formatName of formatCases.slice(1)) {
    await page.getByRole("button", { name: formatName, exact: true }).click();
    assert.notEqual(await formatInput.inputValue(), "");
    await parseFormat.click();
    await formatOutput.waitFor({ state: "visible" });
    await formatError.waitFor({ state: "hidden" });
  }

  // A hand-written document must not be replaced merely because the parser
  // selection changes; the user can intentionally choose another parser.
  await formatInput.fill("user supplied input");
  await page.getByRole("button", { name: "Gaussian Cube", exact: true }).click();
  assert.equal(await formatInput.inputValue(), "user supplied input");

  for (const formatName of formatCases) {
    await page.getByRole("button", { name: formatName, exact: true }).click();
    await loadExample.click();
    await parseFormat.click();
    await formatOutput.waitFor({ state: "visible" });
    await formatError.waitFor({ state: "hidden" });
    await formatInput.fill("malformed input");
    await parseFormat.click();
    if (formatName === "ORCA Output") {
      const orcaResult = JSON.parse(await page.locator("#formats-raw-json").textContent());
      assert.equal(orcaResult.termination.kind, "incomplete");
    } else {
      await formatError.waitFor({ state: "visible" });
      assert.match(
        await formatError.innerText(),
        /invalid|malformed|parse|unexpected|expected|found|no /i,
      );
    }
    await loadExample.click();
    await parseFormat.click();
    await formatOutput.waitFor({ state: "visible" });
    await formatError.waitFor({ state: "hidden" });
  }
  // The local explorer is a separate entry point but shares the generated
  // browser artifact. Cover its public workflow and error boundary so the
  // browser gate does not only exercise the demo tab shell.
  await page.goto("http://127.0.0.1:8765/explorer/index.html?browser-smoke=0.89", {
    waitUntil: "networkidle",
  });
  await page.locator("#loading-overlay").waitFor({ state: "hidden" });
  await page.locator("#explorer-btn-sample").click();
  await page.locator("#explorer-status").filter({ hasText: /loaded/i }).waitFor({ state: "visible" });
  assert.match(await page.locator("#explorer-result-count").innerText(), /^\d+ of \d+ shown$/);
  await page.locator("#explorer-paste-textarea").fill("CCO\nC1CC\nCCN");
  await page.locator("#explorer-btn-parse-paste").click();
  await page.locator("#explorer-status").filter({ hasText: /failed to parse/i }).waitFor({ state: "visible" });
  const explorerResultCount = page.locator("#explorer-result-count");
  await explorerResultCount.filter({ hasText: "3 of 3 shown" }).waitFor({ state: "visible" });
  assert.equal(await explorerResultCount.innerText(), "3 of 3 shown");
  assert.equal(await page.locator("#explorer-tbody tr").count(), 3);
  const cancellationInput = Array.from({ length: 2000 }, () => "CCO").join("\n");
  await page.locator("#explorer-paste-textarea").fill(cancellationInput);
  await page.locator("#explorer-btn-parse-paste").click();
  await page.locator("#explorer-cancel").waitFor({ state: "visible" });
  await page.locator("#explorer-cancel").click();
  await page.locator("#explorer-status").filter({ hasText: /cancelled; complete=false/i }).waitFor({ state: "visible" });
  assert.match(await page.locator("#explorer-result-count").innerText(), /^\d+ of \d+ shown$/);
  await page.locator("#explorer-btn-sample").click();
  // The cancelled batch's final status also contains "loaded". Wait for the
  // deterministic 16-row sample result rather than treating that stale text
  // as completion of the new import.
  await explorerResultCount.filter({ hasText: "16 of 16 shown" }).waitFor({ state: "visible" });
  await page.locator("#explorer-filter-text").fill("Aspirin");
  assert.match(await page.locator("#explorer-result-count").innerText(), /^1 of \d+ shown$/);
  await page.locator("#explorer-reference-smiles").fill("C1CC");
  await page.locator("#explorer-btn-similarity").click();
  await page.locator("#explorer-error").waitFor({ state: "visible" });
  assert.match(await page.locator("#explorer-error").innerText(), /Invalid|SMILES|parse/i);
  await page.locator("#explorer-reference-smiles").fill("CCO");
  await page.locator("#explorer-btn-similarity").click();
  await page.locator("#explorer-status").filter({ hasText: /complete/i }).waitFor({ state: "visible" });
  await page.locator("#explorer-filter-text").fill("");
  const exportDownload = page.waitForEvent("download");
  await page.locator("#explorer-btn-export").click();
  assert.equal((await exportDownload).suggestedFilename(), "chematic-explorer-export.csv");

  assert.deepEqual(errors, []);
} finally {
  await browser.close();
}
console.log(`${browserName}: browser smoke passed`);
