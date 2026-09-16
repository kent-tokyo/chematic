#!/usr/bin/env python3
"""Verify that retrying failed SMILES preserves every Explorer row.

Start a static server for ``demo/`` first. The fixture has one valid and one
malformed SMILES. Retrying cannot repair malformed chemistry, but it must keep
the original error row, input index, and export contract intact.
"""

from __future__ import annotations

import argparse
import csv
import io
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default="http://127.0.0.1:8765/explorer/index.html")
    parser.add_argument("--engine", choices=("chromium", "firefox", "webkit"), default="chromium")
    parser.add_argument("--browser", type=Path, help="explicit Chromium executable; required only with --engine chromium")
    args = parser.parse_args()
    if args.engine == "chromium" and args.browser is None:
        parser.error("--browser is required with --engine chromium")

    try:
        from playwright.sync_api import sync_playwright
    except ImportError as exc:  # pragma: no cover - host-dependent browser test
        raise SystemExit(f"Python Playwright is required: {exc}") from exc

    with sync_playwright() as playwright:
        launch_options = {"headless": True}
        if args.browser is not None:
            launch_options["executable_path"] = str(args.browser)
        browser = getattr(playwright, args.engine).launch(**launch_options)
        page = browser.new_page()
        errors: list[str] = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        try:
            page.goto(args.url, wait_until="networkidle")
            page.locator("#loading-overlay").wait_for(state="hidden", timeout=30_000)
            page.locator("#explorer-smiles-input").fill("CCO\nC1")
            page.locator("#explorer-btn-parse").click()
            page.wait_for_function(
                'document.querySelector("#explorer-status").textContent.includes("failed to parse")',
                timeout=30_000,
            )
            assert page.locator("#explorer-result-count").text_content() == "2 of 2 shown"
            page.locator("#explorer-btn-retry-failed").click()
            page.wait_for_function(
                'document.querySelector("#explorer-status").textContent.startsWith("Retry complete:")',
                timeout=30_000,
            )
            assert page.locator("#explorer-status").text_content() == "Retry complete: 0 recovered; 1 still failed."
            assert page.locator("#explorer-result-count").text_content() == "2 of 2 shown"
            with page.expect_download() as download_info:
                page.locator("#explorer-btn-export").click()
            csv_text = Path(download_info.value.path()).read_text(encoding="utf-8")
            rows = list(csv.DictReader(io.StringIO(csv_text)))
            assert [row["input_index"] for row in rows] == ["0", "1"]
            assert [row["parse_status"] for row in rows] == ["ok", "error"]
            assert errors == [], errors
        finally:
            browser.close()
    print(f"Explorer failed-SMILES retry smoke passed ({args.engine})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
