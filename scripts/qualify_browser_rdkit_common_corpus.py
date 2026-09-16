#!/usr/bin/env python3
"""Freeze a 10k common-parse population using schematic WASM and RDKit.js.

The candidate input is already source-hash-pinned. This runner qualifies each
candidate in one browser page using both engines, retains the first requested
common-parse rows in deterministic candidate order, and records rejections.
It is a corpus-construction gate, not a timing or chemistry-parity claim.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Handler(BaseHTTPRequestHandler):
    files: dict[str, Path] = {}
    html: bytes = b""

    def do_GET(self) -> None:  # noqa: N802
        if self.path == "/":
            body, content_type = self.html, "text/html"
        elif self.path in self.files:
            path = self.files[self.path]
            body = path.read_bytes()
            content_type = "application/wasm" if path.suffix == ".wasm" else "text/javascript"
        else:
            self.send_error(404)
            return
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, _format: str, *_args: object) -> None:
        pass


def page_html(smiles: list[str]) -> str:
    data = json.dumps(smiles).replace("<", "\\u003c")
    return f"""<!doctype html><meta charset=\"utf-8\">
<script src=\"/rdkit/RDKit_minimal.js\"></script><pre id=\"result\">running</pre>
<script type=\"module\">
const smiles = {data};
const pause = () => new Promise((resolve) => setTimeout(resolve, 0));
const run = async () => {{
  const chematic = await import('/schematic/chematic_wasm.js');
  await chematic.default('/schematic/chematic_wasm_bg.wasm');
  const rdkit = await window.initRDKitModule({{ locateFile: (name) =>
    name.endsWith('.wasm') ? '/rdkit/RDKit_minimal.wasm' : name }});
  const accepted = [];
  const rejected = [];
  for (let index = 0; index < smiles.length; index += 1) {{
    const value = smiles[index];
    let left = null;
    let chematicOk = true;
    try {{ left = chematic.parse_smiles(value); }} catch (_) {{ chematicOk = false; }}
    if (left) left.free();
    let right = null;
    try {{ right = rdkit.get_mol(value); }} catch (_) {{ right = null; }}
    const rdkitOk = Boolean(right);
    if (right) right.delete();
    if (chematicOk && rdkitOk) accepted.push(index);
    else if (rejected.length < 64) rejected.push({{ index, chematic_ok: chematicOk, rdkit_ok: rdkitOk }});
    if ((index + 1) % 100 === 0) await pause();
  }}
  document.querySelector('#result').textContent = JSON.stringify({{
    rdkit_version: rdkit.version(), accepted, accepted_count: accepted.length,
    rejected_sample: rejected,
  }});
}};
run().catch((error) => {{ document.querySelector('#result').textContent = JSON.stringify({{ error: String(error) }}); }});
</script>"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rdkit-package", type=Path, required=True)
    parser.add_argument("--browser", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, default=ROOT / "validation/benchmark_corpora/rdkit-js-browser-candidates-v1.smi")
    parser.add_argument("--output", type=Path, default=ROOT / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")
    parser.add_argument("--rows", type=int, default=10_000)
    parser.add_argument("--schematic-dir", type=Path, default=ROOT / "demo/pkg")
    args = parser.parse_args()
    if args.rows <= 0:
        parser.error("--rows must be positive")
    candidate_bytes = args.candidate.read_bytes()
    smiles = [line.strip() for line in candidate_bytes.decode().splitlines() if line.strip()]
    if len(smiles) < args.rows:
        raise SystemExit(f"candidate has {len(smiles)} rows, expected at least {args.rows}")
    rdkit_dist = args.rdkit_package / "dist"
    Handler.files = {
        "/schematic/chematic_wasm.js": args.schematic_dir / "chematic_wasm.js",
        "/schematic/chematic_wasm_bg.wasm": args.schematic_dir / "chematic_wasm_bg.wasm",
        "/rdkit/RDKit_minimal.js": rdkit_dist / "RDKit_minimal.js",
        "/rdkit/RDKit_minimal.wasm": rdkit_dist / "RDKit_minimal.wasm",
    }
    missing = [str(path) for path in Handler.files.values() if not path.is_file()]
    if missing:
        raise SystemExit(f"missing browser artifacts: {missing}")
    Handler.html = page_html(smiles).encode()
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        from playwright.sync_api import sync_playwright
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(executable_path=str(args.browser), headless=True)
            try:
                page = browser.new_page()
                page.goto(f"http://127.0.0.1:{server.server_address[1]}/", wait_until="domcontentloaded")
                page.wait_for_function('document.querySelector("#result").textContent !== "running"', timeout=180_000)
                result = json.loads(page.locator("#result").text_content() or "{}")
            finally:
                browser.close()
    finally:
        server.shutdown()
        thread.join()
    if result.get("error"):
        raise SystemExit(f"common-parse qualification failed: {result['error']}")
    accepted = result.get("accepted", [])
    if not isinstance(accepted, list) or len(accepted) < args.rows:
        raise SystemExit(f"only {len(accepted) if isinstance(accepted, list) else 0} common-parse rows, expected {args.rows}")
    output_rows = [smiles[index] for index in accepted[: args.rows]]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(output_rows) + "\n", encoding="utf-8")
    package = json.loads((args.rdkit_package / "package.json").read_text(encoding="utf-8"))
    manifest = {
        "schema_version": 1,
        "id": "rdkit-js-browser-10k-v1",
        "purpose": "Fixed common-parse input for browser timing and configured-fingerprint checks.",
        "status": "exposed_not_sealed",
        "not_for": "accuracy tuning, model selection, or a claim of independent generalization",
        "candidate": {"path": str(args.candidate.relative_to(ROOT)), "sha256": hashlib.sha256(candidate_bytes).hexdigest(), "rows": len(smiles)},
        "qualification": {
            "browser": str(args.browser), "rdkit_package_version": package.get("version"),
            "rdkit_core_version": result.get("rdkit_version"), "accepted_count": len(accepted),
            "retained_rows": args.rows, "rejected_sample": result.get("rejected_sample", []),
            "criterion": "both browser-resident engines parse the original input without throwing or returning null",
        },
        "output": {"path": str(args.output.relative_to(ROOT)), "rows": len(output_rows), "sha256": digest(args.output)},
    }
    args.output.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
