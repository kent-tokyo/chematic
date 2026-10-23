import assert from "node:assert/strict";
import { chromium } from "playwright";

const recordArgument = process.argv.find((argument) => argument.startsWith("--records="));
const recordCount = recordArgument === undefined ? 10_000 : Number(recordArgument.slice("--records=".length));
if (!Number.isSafeInteger(recordCount) || recordCount < 1 || recordCount > 100_000) {
  throw new Error("--records must be an integer from 1 through 100000");
}
const capacityMode = recordArgument !== undefined;
const startedAt = performance.now();
if (capacityMode) console.log(`Explorer Worker capacity smoke starting (${recordCount.toLocaleString("en-US")} records)`);

const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage();
  const errors = [];
  page.on("pageerror", (error) => errors.push(String(error)));
  page.on("console", (message) => {
    if (message.type() === "error") errors.push(message.text());
  });
  await page.goto("http://127.0.0.1:8765/explorer/index.html?worker-10k=1", { waitUntil: "networkidle" });
  await page.locator("#loading-overlay").waitFor({ state: "hidden" });
  assert.equal(await page.locator("html").getAttribute("data-explorer-analysis"), "worker");
  if (capacityMode) console.log("Explorer Worker capacity smoke initialized");

  const records = Array.from({ length: recordCount }, () => "C").join("\n");
  await page.locator("#explorer-paste-textarea").evaluate((element, value) => {
    element.value = value;
    element.dispatchEvent(new Event("input", { bubbles: true }));
  }, records);
  await page.locator("#explorer-btn-parse-paste").click();
  if (capacityMode) console.log("Explorer Worker capacity smoke submitted");
  await page.locator("#explorer-status").waitFor({ hasText: `${recordCount} molecules loaded.` , timeout: 600_000 });
  const resultCount = page.locator("#explorer-result-count");
  const renderedCount = Math.min(250, recordCount);
  const expectedCount = renderedCount === recordCount
    ? `${recordCount} of ${recordCount} shown`
    : `${renderedCount} rendered of ${recordCount} matching (${recordCount} loaded)`;
  await page.waitForFunction(
    (expected) => document.querySelector("#explorer-result-count")?.textContent === expected,
    expectedCount,
  );
  assert.equal(await resultCount.textContent(), expectedCount);
  assert.deepEqual(errors, []);
} finally {
  await browser.close();
}
console.log(
  `Explorer Worker ${recordCount.toLocaleString("en-US")}-record smoke passed ` +
  `(${(performance.now() - startedAt).toFixed(0)} ms)`
);
