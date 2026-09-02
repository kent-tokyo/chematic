import assert from "node:assert/strict";
import { chromium, firefox, webkit } from "playwright";

const browserName = process.argv[2] ?? "chromium";
const browsers = { chromium, firefox, webkit };
assert.ok(browsers[browserName], `unknown browser: ${browserName}`);

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
  const errors = [];
  page.on("pageerror", (error) => errors.push(String(error)));
  page.on("console", (message) => {
    if (message.type() === "error") errors.push(message.text());
  });
  await page.goto("http://127.0.0.1:8765/index.html?browser-smoke=0.89", {
    waitUntil: "networkidle",
  });
  await page.locator("#version-badge").waitFor({ state: "visible" });
  assert.equal(await page.locator("#version-badge").innerText(), "v0.89.0");
  await page.locator("#smiles-input").fill("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
  await page.locator("#btn-calc").click();
  const hba = page.locator("#desc-tbody tr").filter({ hasText: "HBA" }).locator("td").nth(1);
  await hba.waitFor({ state: "visible" });
  assert.equal(await hba.innerText(), "6");
  for (const [nextSmiles, expectedHba] of [
    ["CCO", "1"],
    ["CC(=O)O", "1"],
    ["c1ccccc1", "0"],
    ["Cn1cnc2c1c(=O)n(c(=O)n2C)C", "6"],
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
  for (const formatName of formatCases) {
    await page.getByRole("button", { name: formatName, exact: true }).click();
    await page.getByRole("button", { name: "Load Example", exact: true }).click();
    await page.getByRole("button", { name: "Parse", exact: true }).click();
    await formatOutput.waitFor({ state: "visible" });
    await formatError.waitFor({ state: "hidden" });
    await formatInput.fill("malformed input");
    await page.getByRole("button", { name: "Parse", exact: true }).click();
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
    await page.getByRole("button", { name: "Load Example", exact: true }).click();
    await page.getByRole("button", { name: "Parse", exact: true }).click();
    await formatOutput.waitFor({ state: "visible" });
    await formatError.waitFor({ state: "hidden" });
  }
  assert.deepEqual(errors, []);
} finally {
  await browser.close();
}
console.log(`${browserName}: browser smoke passed`);
