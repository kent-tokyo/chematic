#!/usr/bin/env python3
"""Exercise the Explorer's resumable Worker SDF path in a real browser.

Start a static server for ``demo/`` first, then run this command against it.
The test uses a valid two-record fixture plus one malformed record and verifies
that the malformed input is preserved as an error row instead of disappearing.
"""

from __future__ import annotations

import argparse
import csv
import io
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


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
    except ImportError as exc:  # pragma: no cover - host-dependent integration test
        raise SystemExit(f"Python Playwright is required: {exc}") from exc

    valid = (ROOT / "benchmarks" / "fixtures" / "streaming.sdf").read_text(encoding="utf-8")
    mixed = f"{valid}malformed\n$$$$\n".encode()
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
            page.set_input_files(
                "#explorer-file-input",
                {
                    "name": "mixed.sdf",
                    "mimeType": "chemical/x-mdl-sdfile",
                    "buffer": mixed,
                },
            )
            page.wait_for_function(
                'document.querySelector("#explorer-status").textContent.includes("SDF records")',
                timeout=30_000,
            )
            assert page.locator("#explorer-status").text_content() == (
                "2 loaded, 1 rejected from 3 SDF records."
            )
            assert page.locator("#explorer-result-count").text_content() == "3 of 3 shown"
            with page.expect_download() as download_info:
                page.locator("#explorer-btn-export").click()
            exported = download_info.value
            csv_text = Path(exported.path()).read_text(encoding="utf-8")
            rows = list(csv.DictReader(io.StringIO(csv_text)))
            assert [row["input_index"] for row in rows] == ["0", "1", "2"]
            assert [row["parse_status"] for row in rows] == ["ok", "ok", "error"]
            assert errors == [], errors
        finally:
            browser.close()
    print(f"Explorer Worker SDF smoke passed ({args.engine})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
