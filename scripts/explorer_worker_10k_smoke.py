#!/usr/bin/env python3
"""Exercise the Explorer's 10,000-record Worker path in a real browser.

Start a static server for ``demo/`` first.  This test bypasses keystroke-based
textarea filling so that it measures the Explorer's parsing/rendering path,
not automation transport of a 20 KB string. It checks both successful loading
and a responsive cancellation request. By default the browser context goes
offline after the page and Worker have initialized, proving that subsequent
local file/paste analysis does not depend on a network request.
"""

from __future__ import annotations

import argparse
import csv
import io
from pathlib import Path


def set_records(page, records: str) -> None:
    page.locator("#explorer-paste-textarea").evaluate(
        """(node, value) => {
          node.value = value;
          node.dispatchEvent(new Event("input", { bubbles: true }));
        }""",
        records,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default="http://127.0.0.1:8765/explorer/index.html")
    parser.add_argument(
        "--engine",
        choices=("chromium", "firefox", "webkit"),
        default="chromium",
        help="Playwright engine to launch (default: chromium)",
    )
    parser.add_argument(
        "--browser",
        type=Path,
        help="explicit Chromium executable; required only with --engine chromium",
    )
    parser.add_argument(
        "--online-after-load",
        action="store_true",
        help="keep the browser online after initialization (diagnostic only)",
    )
    args = parser.parse_args()
    if args.engine == "chromium" and args.browser is None:
        parser.error("--browser is required with --engine chromium")

    try:
        from playwright.sync_api import sync_playwright
    except ImportError as exc:  # pragma: no cover - host-dependent integration test
        raise SystemExit(f"Python Playwright is required: {exc}") from exc

    records = ("C\n" * 10_000).rstrip()
    # Keep the completion lane inexpensive, but use a nontrivial valid molecule
    # for the cancellation lane so the control remains actionable long enough
    # to exercise an in-flight Worker request.
    cancel_records = ("CC(=O)Oc1ccccc1C(=O)O\n" * 10_000).rstrip()
    expected_count = "250 rendered of 10000 matching (10000 loaded)"
    with sync_playwright() as playwright:
        launcher = getattr(playwright, args.engine)
        launch_options = {"headless": True}
        if args.browser is not None:
            launch_options["executable_path"] = str(args.browser)
        browser = launcher.launch(**launch_options)
        context = browser.new_context()
        page = context.new_page()
        errors: list[str] = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        try:
            page.goto(args.url, wait_until="networkidle")
            page.locator("#loading-overlay").wait_for(state="hidden", timeout=30_000)
            if not args.online_after_load:
                # Static assets and both WASM instances are already ready. A fetch in
                # paste/file/cancel handling must fail this test instead of using a
                # network path that the product does not advertise.
                context.set_offline(True)
            set_records(page, records)
            page.locator("#explorer-btn-parse-paste").click()
            page.wait_for_function(
                'document.querySelector("#explorer-status").textContent.includes("10000 molecules loaded.")',
                timeout=120_000,
            )
            assert page.locator("#explorer-result-count").text_content() == expected_count

            csv_text = "smiles,name\n" + "\n".join(
                f"C,molecule-{index}" for index in range(10_000)
            )
            page.set_input_files(
                "#explorer-file-input",
                {
                    "name": "ten-thousand.csv",
                    "mimeType": "text/csv",
                    "buffer": csv_text.encode(),
                },
            )
            # The previous SMI run has the same terminal status text. Wait for a
            # CSV-specific rendered name before accepting that a new operation began;
            # otherwise a fast assertion could export the previous data set.
            page.wait_for_function(
                'document.querySelector("#explorer-tbody").textContent.includes("molecule-0")',
                timeout=120_000,
            )
            page.wait_for_function(
                'document.querySelector("#explorer-status").textContent.includes("10000 molecules loaded.")',
                timeout=120_000,
            )
            assert page.locator("#explorer-result-count").text_content() == expected_count
            # Export must retain every analysed row, not merely the 250-row DOM
            # window. The default input-order sort is part of the data contract.
            with page.expect_download() as download_info:
                page.locator("#explorer-btn-export").click()
            exported = download_info.value
            rows = list(csv.DictReader(io.StringIO(Path(exported.path()).read_text(encoding="utf-8"))))
            assert len(rows) == 10_000
            assert [row["input_index"] for row in rows] == [str(index) for index in range(10_000)]
            assert {row["parse_status"] for row in rows} == {"ok"}

            # CSV is an unknown-length stream. Cancellation must preserve the
            # observed-but-unprocessed prefix and state that the later suffix
            # was never read; it must not turn either into skipped/success rows.
            cancel_csv_text = "smiles,name\n" + "\n".join(
                f"CC(=O)Oc1ccccc1C(=O)O,cancel-{index}" for index in range(10_000)
            )
            page.set_input_files(
                "#explorer-file-input",
                {
                    "name": "cancel-stream.csv",
                    "mimeType": "text/csv",
                    "buffer": cancel_csv_text.encode(),
                },
            )
            page.locator("#explorer-cancel").wait_for(state="visible", timeout=30_000)
            page.locator("#explorer-cancel").dispatch_event("click")
            page.wait_for_function(
                'document.querySelector("#explorer-status").textContent.includes("unread input=unknown")',
                timeout=30_000,
            )
            cancelled_status = page.locator("#explorer-status").text_content()
            assert "complete=false" in cancelled_status
            assert "skipped" not in cancelled_status.lower()

            set_records(page, cancel_records)
            page.locator("#explorer-btn-parse-paste").click()
            page.locator("#explorer-cancel").wait_for(state="visible", timeout=30_000)
            # Dispatch immediately after visibility: an ordinary automation click waits
            # for layout stability, but even this deliberately heavier fixture can advance
            # a Worker batch during that wait.
            page.locator("#explorer-cancel").dispatch_event("click")
            page.wait_for_function(
                'document.querySelector("#explorer-status").textContent.startsWith("Cancelled after ")',
                timeout=30_000,
            )
            assert errors == [], errors
        finally:
            browser.close()
    print(f"Explorer Worker 10,000-record smoke passed ({args.engine})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
