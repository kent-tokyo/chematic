import assert from "node:assert/strict";
import { chromium } from "playwright";

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

  const records = Array.from({ length: 10_000 }, () => "C").join("\n");
  await page.locator("#explorer-paste-textarea").evaluate((element, value) => {
    element.value = value;
    element.dispatchEvent(new Event("input", { bubbles: true }));
  }, records);
  await page.locator("#explorer-btn-parse-paste").click();
  await page.locator("#explorer-status").waitFor({ hasText: "10000 molecules loaded." , timeout: 120_000 });
  const resultCount = page.locator("#explorer-result-count");
  await resultCount.waitFor({ hasText: "250 rendered of 10000 matching (10000 loaded)" });
  assert.equal(await resultCount.textContent(), "250 rendered of 10000 matching (10000 loaded)");
  assert.deepEqual(errors, []);
} finally {
  await browser.close();
}
console.log("Explorer Worker 10,000-record smoke passed");
