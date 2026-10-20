#!/usr/bin/env python3
"""Run the generated WASM comparison page with Python Playwright.

This is the browser-engine fallback for hosts where Chrome's ``--dump-dom``
does not terminate after an asynchronous module page has finished.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--schematic-dir", type=Path, required=True)
    parser.add_argument("--rdkit-package", type=Path, required=True)
    parser.add_argument("--browser", type=Path, required=True)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--rows", type=int, required=True)
    parser.add_argument("--timeout-ms", type=int, default=60_000)
    args = parser.parse_args()

    try:
        from playwright.sync_api import sync_playwright
    except ImportError as exc:  # pragma: no cover - environment-dependent gate
        raise SystemExit(f"Python Playwright is required for this gate: {exc}") from exc

    html = args.html.read_text(encoding="utf-8")
    if '<base href="http://schematic.test/">' not in html:
        html = html.replace(
            '<meta charset="utf-8">',
            '<meta charset="utf-8"><base href="http://schematic.test/">',
            1,
        )
    schematic_dir = args.schematic_dir.resolve()
    rdkit_dist = (args.rdkit_package / "dist").resolve()
    files = {
        "/schematic/chematic_wasm.js": schematic_dir / "chematic_wasm.js",
        "/schematic/chematic_wasm_bg.wasm": schematic_dir / "chematic_wasm_bg.wasm",
        "/rdkit/RDKit_minimal.js": rdkit_dist / "RDKit_minimal.js",
        "/rdkit/RDKit_minimal.wasm": rdkit_dist / "RDKit_minimal.wasm",
    }

    def artifact(path: Path) -> dict:
        data = path.read_bytes()
        return {
            "file": path.name,
            "raw_bytes": len(data),
            "gzip_bytes": len(gzip.compress(data, compresslevel=9)),
            "sha256": hashlib.sha256(data).hexdigest(),
        }

    package_json = rdkit_dist.parent / "package.json"
    package_metadata = json.loads(package_json.read_text(encoding="utf-8"))
    corpus_bytes = args.corpus.read_bytes()
    corpus_rows = [line for line in corpus_bytes.decode().splitlines() if line.strip()]
    if len(corpus_rows) < args.rows:
        raise SystemExit(f"corpus contains {len(corpus_rows)} rows, expected at least {args.rows}")

    def fulfill(route) -> None:
        url = route.request.url
        path = "/" + url.split("/", 3)[-1] if "schematic.test" in url else "/"
        if path == "/":
            route.fulfill(content_type="text/html", body=html)
        elif path in files:
            content_type = "application/wasm" if path.endswith(".wasm") else "text/javascript"
            route.fulfill(content_type=content_type, body=files[path].read_bytes())
        else:
            route.abort()

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=str(args.browser),
            headless=True,
        )
        page = browser.new_page()
        page.route("**/*", fulfill)
        page.set_content(html, wait_until="domcontentloaded")
        page.wait_for_function(
            'document.querySelector("#result").textContent !== "running"',
            timeout=args.timeout_ms,
        )
        measurement = json.loads(page.locator("#result").text_content() or "")
        browser.close()

    if measurement.get("error"):
        raise SystemExit(f"browser page failed: {measurement['error']}\n{measurement.get('stack', '')}")
    result = {
        "schema_version": 1,
        "gate": "wasm-vs-rdkit-browser-playwright",
        "configuration": {
            "browser": "Playwright Chromium",
            "timing": "same browser page, sequential lanes",
            "timeout_ms": args.timeout_ms,
        },
        "corpus": {
            "path": str(args.corpus),
            "rows": args.rows,
            "sha256": hashlib.sha256(corpus_bytes).hexdigest(),
        },
        "artifacts": {
            "schematic_wasm": artifact(files["/schematic/chematic_wasm_bg.wasm"]),
            "rdkit_minimal_wasm": artifact(files["/rdkit/RDKit_minimal.wasm"]),
            "rdkit_minimal_js": artifact(files["/rdkit/RDKit_minimal.js"]),
            "rdkit_package_json": {
                "file": str(package_json),
                "version": package_metadata.get("version"),
                "sha256": hashlib.sha256(package_json.read_bytes()).hexdigest(),
            },
            "rdkit_typescript_declarations": {
                "file": str(rdkit_dist / "index.d.ts"),
                "sha256": hashlib.sha256((rdkit_dist / "index.d.ts").read_bytes()).hexdigest(),
            },
        },
        **measurement,
        "interpretation": "Browser measurements exclude peak RSS and are not equivalent to Node process measurements; they are a Playwright Chromium-specific WASM/API comparison.",
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
