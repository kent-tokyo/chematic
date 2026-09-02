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
  await page.locator("#smiles-input").fill("C".repeat(1_000_001));
  await page.locator("#btn-calc").click();
  await page.locator("#error-desc").waitFor({ state: "visible" });
  assert.match(await page.locator("#error-desc").innerText(), /input|size|large/i);
  await page.locator("#smarts-input").fill("[");
  await page.getByRole("button", { name: "Highlight", exact: true }).click();
  await page.locator("#error-desc").waitFor({ state: "visible" });
  assert.match(await page.locator("#error-desc").innerText(), /SMARTS/i);
  await page.locator("#smarts-input").fill("c1ccccc1");
  await page.getByRole("button", { name: "Highlight", exact: true }).click();
  await page.locator("#error-desc").waitFor({ state: "hidden" });
  await page.getByRole("button", { name: "日", exact: true }).click();
  await page.getByText("記述子計算機", { exact: true }).waitFor({ state: "visible" });
  await page.getByRole("button", { name: "EN", exact: true }).click();
  await page.getByText("Descriptor Calculator", { exact: true }).waitFor({ state: "visible" });
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
  assert.deepEqual(errors, []);
} finally {
  await browser.close();
}
console.log(`${browserName}: browser smoke passed`);
