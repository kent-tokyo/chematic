import assert from "node:assert/strict";
import { chromium, firefox, webkit } from "playwright";

const browserName = process.argv[2] ?? "chromium";
const browsers = { chromium, firefox, webkit };
assert.ok(browsers[browserName], `unknown browser: ${browserName}`);

const browser = await browsers[browserName].launch({ headless: true });
try {
  const page = await browser.newPage();
  const errors = [];
  page.on("pageerror", (error) => errors.push(String(error)));
  page.on("console", (message) => {
    if (message.type() === "error") errors.push(message.text());
  });

  await page.goto("http://127.0.0.1:8765/explorer/index.html?browser-smoke=1", {
    waitUntil: "networkidle",
  });
  await page.locator("#loading-overlay").waitFor({ state: "hidden" });

  const input = page.locator("#explorer-paste-textarea");
  const status = page.locator("#explorer-status");
  const error = page.locator("#explorer-error");
  await input.fill("CCO\nC1CC\nCCN");
  await page.locator("#explorer-btn-parse-paste").click();
  await status.waitFor({ state: "visible" });
  await page.waitForFunction(() => /loaded|failed/i.test(document.querySelector("#explorer-status")?.textContent ?? ""));
  assert.match(await status.innerText(), /2 loaded, 1 failed/);
  assert.equal(await page.locator("#explorer-result-count").innerText(), "3 of 3 shown");
  assert.equal(await error.isVisible(), false);

  const largeInput = Array.from({ length: 2000 }, () => "CCO").join("\n");
  await input.fill(largeInput);
  await page.locator("#explorer-btn-parse-paste").click();
  await status.waitFor({ hasText: /Parsing/ });
  await page.locator("#explorer-cancel").click();
  await status.waitFor({ hasText: /Cancelled after/ });
  assert.equal(await page.locator("#explorer-cancel").isVisible(), false);
  const cancelledCount = Number((await page.locator("#explorer-result-count").innerText()).match(/\d+/)?.[0]);
  assert.ok(cancelledCount < 2000, `cancel should stop before all records: ${cancelledCount}`);
  assert.deepEqual(errors, []);
} finally {
  await browser.close();
}
console.log(`${browserName}: explorer browser smoke passed`);
