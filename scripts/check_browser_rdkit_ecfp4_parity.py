#!/usr/bin/env python3
"""Check chematic RDKit-compatible ECFP4 against official RDKit.js in a browser.

This is a correctness gate, not a timing benchmark.  Both engines share one
browser page so their packed fingerprint bits can be compared directly.  The
page retains only counters and the first mismatch: unlike the older DOM-dump
runner, it never materializes a multi-megabyte array of per-row bit strings.
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


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def page_html(smiles: list[str], schematic_operation: str) -> str:
    data = json.dumps(smiles).replace("<", "\\u003c")
    operation = json.dumps(schematic_operation)
    return f"""<!doctype html><meta charset=\"utf-8\">
<script src=\"/rdkit/RDKit_minimal.js\"></script><pre id=\"result\">running</pre>
<script type=\"module\">
const smiles = {data};
const schematicOperation = {operation};
const yieldToBrowser = () => new Promise((resolve) => setTimeout(resolve, 0));
const run = async () => {{
  const chematic = await import('/schematic/chematic_wasm.js');
  await chematic.default('/schematic/chematic_wasm_bg.wasm');
  const rdkit = await window.initRDKitModule({{ locateFile: (name) =>
    name.endsWith('.wasm') ? '/rdkit/RDKit_minimal.wasm' : name }});
  let exact = 0;
  let unsupported = [];
  let firstMismatch = null;
  const mismatchSample = [];
  for (let index = 0; index < smiles.length; index += 1) {{
    const value = smiles[index];
    const left = chematic.parse_smiles(value);
    let leftBits;
    try {{
      if (schematicOperation === "prepared") {{
        const prepared = chematic.prepare_rdkit_ecfp4(left);
        try {{ leftBits = prepared.bitvec(); }} finally {{ prepared.free(); }}
      }} else {{
        leftBits = chematic.rdkit_ecfp4_bitvec(left);
      }}
    }} catch (error) {{
      left.free();
      const right = rdkit.get_mol(value);
      if (!right) throw new Error('RDKit rejected row ' + index);
      right.delete();
      const message = String(error);
      if (!message.includes('unsupported RDKit coordination sanitization')) throw error;
      unsupported.push({{ index, smiles: value, error: message }});
      if ((index + 1) % 100 === 0) await yieldToBrowser();
      continue;
    }}
    const right = rdkit.get_mol(value);
    if (!right) throw new Error('RDKit rejected row ' + index);
    const rightBits = right.get_morgan_fp(JSON.stringify({{ radius: 2, nBits: 2048 }}));
    let differingBits = Math.abs(leftBits.length * 8 - rightBits.length);
    const firstDifferingBits = [];
    if (leftBits.length * 8 === rightBits.length) {{
      for (let byte = 0; byte < leftBits.length; byte += 1) {{
        const packed = leftBits[byte];
        for (let bit = 0; bit < 8; bit += 1) {{
          if (((packed >> bit) & 1) !== (rightBits.charCodeAt(byte * 8 + bit) - 48)) {{
            differingBits += 1;
            if (firstDifferingBits.length < 16) firstDifferingBits.push(byte * 8 + bit);
          }}
        }}
      }}
    }}
    left.free();
    right.delete();
    if (differingBits === 0) exact += 1;
    else {{
      const mismatch = {{ index, smiles: value, differing_bits: differingBits, first_differing_bits: firstDifferingBits }};
      if (!firstMismatch) firstMismatch = mismatch;
      if (mismatchSample.length < 64) mismatchSample.push(mismatch);
    }}
    if ((index + 1) % 100 === 0) await yieldToBrowser();
  }}
  document.querySelector('#result').textContent = JSON.stringify({{
    rdkit_version: rdkit.version(), compared_rows: smiles.length,
    exact_matches: exact, supported_rows: smiles.length - unsupported.length,
    unsupported_rows: unsupported.length, unsupported_sample: unsupported.slice(0, 64),
    first_mismatch: firstMismatch, mismatch_sample: mismatchSample,
  }});
}};
run().catch((error) => {{ document.querySelector('#result').textContent = JSON.stringify({{ error: String(error) }}); }});
</script>"""


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


def artifact(path: Path) -> dict[str, object]:
    content = path.read_bytes()
    return {
        "path": str(path), "raw_bytes": len(content),
        "gzip9_bytes": len(gzip.compress(content, 9)), "sha256": hashlib.sha256(content).hexdigest(),
    }


def npm_package_provenance(package_dir: Path) -> dict[str, object]:
    """Return the npm tarball identity recorded beside an installed package."""
    node_modules = package_dir.parents[1]
    lock_path = node_modules / ".package-lock.json"
    if not lock_path.is_file():
        return {"lock_status": "not_found"}
    try:
        entry = json.loads(lock_path.read_text(encoding="utf-8")).get("packages", {}).get(
            "node_modules/@rdkit/rdkit"
        )
    except (json.JSONDecodeError, OSError):
        return {"lock_status": "unreadable"}
    if not isinstance(entry, dict):
        return {"lock_status": "entry_not_found"}
    return {
        "lock_status": "recorded",
        "lock_sha256": sha256(lock_path),
        "resolved": entry.get("resolved"),
        "integrity": entry.get("integrity"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rdkit-package", type=Path, required=True)
    parser.add_argument("--engine", choices=("chromium", "firefox", "webkit"), default="chromium")
    parser.add_argument("--browser", type=Path, help="explicit Chromium executable; required only with --engine chromium")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--corpus", type=Path, default=ROOT / "scripts" / "descriptor_census_corpus.smi")
    parser.add_argument("--schematic-dir", type=Path, default=ROOT / "demo" / "pkg")
    parser.add_argument(
        "--schematic-package-kind",
        choices=("published", "source_candidate"),
        default="published",
    )
    parser.add_argument("--schematic-source-revision")
    parser.add_argument("--schematic-source-diff-sha256")
    parser.add_argument("--rows", type=int, default=1000)
    parser.add_argument(
        "--schematic-operation", choices=("direct", "prepared"), default="direct"
    )
    parser.add_argument(
        "--allow-unsupported",
        action="store_true",
        help="permit only the declared RDKit coordination-sanitization typed refusal; other errors still fail",
    )
    args = parser.parse_args()
    if args.rows <= 0:
        raise SystemExit("--rows must be positive")
    if args.engine == "chromium" and args.browser is None:
        parser.error("--browser is required with --engine chromium")
    if args.schematic_package_kind == "source_candidate" and (
        not args.schematic_source_revision or not args.schematic_source_diff_sha256
    ):
        parser.error(
            "source candidates require --schematic-source-revision and "
            "--schematic-source-diff-sha256"
        )

    corpus_bytes = args.corpus.read_bytes()
    smiles = [line.strip() for line in corpus_bytes.decode().splitlines() if line.strip()][:args.rows]
    if len(smiles) != args.rows:
        raise SystemExit(f"corpus has {len(smiles)} usable rows, expected {args.rows}")
    rdkit_dist = args.rdkit_package / "dist"
    Handler.files = {
        "/schematic/chematic_wasm.js": args.schematic_dir / "chematic_wasm.js",
        "/schematic/chematic_wasm_bg.wasm": args.schematic_dir / "chematic_wasm_bg.wasm",
        "/rdkit/RDKit_minimal.js": rdkit_dist / "RDKit_minimal.js",
        "/rdkit/RDKit_minimal.wasm": rdkit_dist / "RDKit_minimal.wasm",
    }
    missing = [str(path) for path in Handler.files.values() if not path.is_file()]
    if missing:
        raise SystemExit(f"missing benchmark artifacts: {missing}")
    Handler.html = page_html(smiles, args.schematic_operation).encode()
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as exc:  # pragma: no cover - host dependent
        server.shutdown()
        raise SystemExit(f"Python Playwright is required: {exc}") from exc
    try:
        with sync_playwright() as playwright:
            launch_options = {"headless": True}
            if args.browser is not None:
                launch_options["executable_path"] = str(args.browser)
            browser = getattr(playwright, args.engine).launch(**launch_options)
            try:
                page = browser.new_page()
                page.goto(f"http://127.0.0.1:{server.server_address[1]}/", wait_until="domcontentloaded")
                page.wait_for_function('document.querySelector("#result").textContent !== "running"', timeout=120_000)
                result = json.loads(page.locator("#result").text_content() or "{}")
            finally:
                browser.close()
    finally:
        server.shutdown()
        thread.join()
    if result.get("error"):
        raise SystemExit(f"browser parity run failed: {result['error']}")
    package_json = json.loads((args.rdkit_package / "package.json").read_text())
    document = {
        "schema_version": 1,
        "gate": "browser-rdkit-ecfp4-bit-parity",
        "configuration": {
            "engine": args.engine, "browser": str(args.browser) if args.browser else None, "rows": args.rows,
            "allow_unsupported": args.allow_unsupported,
            "schematic_operation": args.schematic_operation,
            "operation": "chematic rdkit_ecfp4_bitvec vs RDKit.js Morgan radius=2, 2048 bits",
            "comparison": "same browser page; direct packed-bit comparison; yields every 100 rows",
        },
        "corpus": {"path": str(args.corpus), "sha256": hashlib.sha256(corpus_bytes).hexdigest(), "rows": args.rows},
        "artifacts": {
            "schematic_wasm": {
                **artifact(Handler.files["/schematic/chematic_wasm_bg.wasm"]),
                "kind": args.schematic_package_kind,
                "source_revision": args.schematic_source_revision,
                "source_diff_sha256": args.schematic_source_diff_sha256,
            },
            "rdkit_wasm": artifact(Handler.files["/rdkit/RDKit_minimal.wasm"]),
            "rdkit_package": {
                "version": package_json.get("version"),
                "package_json_sha256": sha256(args.rdkit_package / "package.json"),
                **npm_package_provenance(args.rdkit_package),
            },
        },
        "result": result,
        "boundary": (
            "Published package artifact. "
            if args.schematic_package_kind == "published"
            else "Locally built source candidate; not a registry release. "
        )
        + "This confirms an explicitly configured fingerprint bit-vector only for successfully evaluated rows. Typed refusals are separately counted and are not parity claims. It is not a general SMILES, stereo, canonicalization, search-ranking, memory, or performance claim.",
    }
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(document["result"], indent=2))
    unsupported = int(result.get("unsupported_rows", 0))
    expected_exact = args.rows - unsupported
    return 0 if result["exact_matches"] == expected_exact and (args.allow_unsupported or unsupported == 0) else 1


if __name__ == "__main__":
    raise SystemExit(main())
