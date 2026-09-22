#!/usr/bin/env python3
"""Measure chematic and RDKit.js in fresh, isolated browser processes.

This is deliberately separate from the older same-page gate.  Each arm receives
its own browser process and cache-free local server route, so initialization and
browser-exposed JS heap snapshots are not contaminated by loading the other arm.
It does not measure process RSS; that belongs to a separately instrumented lane.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import socket
import statistics
import subprocess
import tempfile
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse
from urllib.request import urlopen


ROOT = Path(__file__).resolve().parents[1]

# Measured by the direct browser bit-parity gate. RDKit sanitizes this exact
# Fe(II) degree-10, anionic-carbon coordination graph into directed coordinate
# bonds, which chematic intentionally rejects with a typed compatibility error.
# Keep the input in the corpus and record this exclusion in every performance
# result; do not silently remove it from the corpus or infer a wider metal rule.
FINGERPRINT_TYPED_UNSUPPORTED_SMILES = frozenset(
    {
        "CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10",
    }
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifact(path: Path) -> dict[str, object]:
    body = path.read_bytes()
    return {
        "path": str(path),
        "raw_bytes": len(body),
        "gzip9_bytes": len(gzip.compress(body, 9)),
        "sha256": hashlib.sha256(body).hexdigest(),
    }


def npm_package_provenance(
    package_dir: Path, package_json: Path, package_name: str
) -> dict[str, object]:
    """Capture npm's resolved tarball identity when the install kept its lockfile.

    A package.json hash proves only the installed metadata.  npm v7+ writes a
    hidden node_modules lock that also records the registry URL and SRI string;
    retain both when available without making an extracted/offline package
    unusable as a benchmark input.
    """
    record: dict[str, object] = {"sha256": sha256(package_json)}
    node_modules = next(
        (parent for parent in package_dir.parents if parent.name == "node_modules"),
        None,
    )
    if node_modules is None:
        record["lock_status"] = "not_found"
        return record
    lock_path = node_modules / ".package-lock.json"
    if not lock_path.is_file():
        record["lock_status"] = "not_found"
        return record
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
        entry = lock.get("packages", {}).get(f"node_modules/{package_name}")
    except (json.JSONDecodeError, OSError):
        record["lock_status"] = "unreadable"
        return record
    if not isinstance(entry, dict):
        record["lock_status"] = "entry_not_found"
        return record
    record.update(
        {
            "lock_status": "recorded",
            "lock_sha256": sha256(lock_path),
            "resolved": entry.get("resolved"),
            "integrity": entry.get("integrity"),
        }
    )
    return record


def summary(values: list[float]) -> dict[str, float | int]:
    ordered = sorted(values)

    def pct(fraction: float) -> float:
        position = (len(ordered) - 1) * fraction
        low = int(position)
        high = min(low + 1, len(ordered) - 1)
        return round(
            ordered[low] + (ordered[high] - ordered[low]) * (position - low), 6
        )

    return {
        "count": len(values),
        "p50_ms": pct(0.5),
        "p95_ms": pct(0.95),
        "mean_ms": round(statistics.fmean(values), 6),
    }


def free_local_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("127.0.0.1", 0))
        return int(probe.getsockname()[1])


def process_tree_rss_bytes(root_pid: int) -> tuple[int, int]:
    """Return summed RSS and descendant count from `ps`'s KiB RSS column.

    This is process-tree RSS, not de-duplicated physical memory: shared pages
    can appear in more than one process. Keeping that qualification makes the
    value useful for same-host arm comparison without misrepresenting it as a
    unique-resident-memory measurement.
    """
    output = subprocess.run(
        ["ps", "-axo", "pid=,ppid=,rss="],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    children: dict[int, list[int]] = {}
    rss_kib: dict[int, int] = {}
    for line in output.splitlines():
        fields = line.split()
        if len(fields) != 3:
            continue
        try:
            pid, parent, rss = map(int, fields)
        except ValueError:
            continue
        children.setdefault(parent, []).append(pid)
        rss_kib[pid] = rss
    pending, seen = [root_pid], set()
    while pending:
        pid = pending.pop()
        if pid in seen:
            continue
        seen.add(pid)
        pending.extend(children.get(pid, []))
    return sum(rss_kib.get(pid, 0) for pid in seen) * 1024, len(seen)


class ProcessTreeRssMonitor:
    def __init__(self, root_pid: int) -> None:
        self.root_pid = root_pid
        self._stop = threading.Event()
        self._lock = threading.Lock()
        self.peak_rss_bytes = 0
        self.peak_processes = 0
        self.samples = 0
        self._thread = threading.Thread(target=self._run, daemon=True)

    def start(self) -> None:
        self._thread.start()

    def _run(self) -> None:
        while not self._stop.is_set():
            try:
                rss_bytes, processes = process_tree_rss_bytes(self.root_pid)
            except (OSError, subprocess.SubprocessError):
                time.sleep(0.05)
                continue
            with self._lock:
                self.peak_rss_bytes = max(self.peak_rss_bytes, rss_bytes)
                self.peak_processes = max(self.peak_processes, processes)
                self.samples += 1
            time.sleep(0.05)

    def finish(self) -> dict[str, object]:
        self._stop.set()
        self._thread.join(timeout=2)
        with self._lock:
            if self.samples == 0:
                return {"status": "unavailable", "reason": "no successful ps samples"}
            return {
                "status": "measured",
                "metric": "summed process-tree RSS; shared pages are not de-duplicated",
                "ps_rss_unit": "KiB x 1024",
                "peak_rss_bytes": self.peak_rss_bytes,
                "peak_processes": self.peak_processes,
                "samples": self.samples,
            }


def launch_chromium_with_rss(
    playwright: object, executable: Path
) -> tuple[object, subprocess.Popen[bytes], tempfile.TemporaryDirectory[str]]:
    """Launch Chromium under our PID so its process tree can be measured."""
    temp_profile = tempfile.TemporaryDirectory(prefix="chematic-browser-rss-")
    port = free_local_port()
    process = subprocess.Popen(
        [
            str(executable),
            "--headless=new",
            "--no-first-run",
            "--no-default-browser-check",
            f"--remote-debugging-port={port}",
            f"--user-data-dir={temp_profile.name}",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    endpoint = f"http://127.0.0.1:{port}"
    deadline = time.monotonic() + 20
    while time.monotonic() < deadline:
        if process.poll() is not None:
            temp_profile.cleanup()
            raise RuntimeError(
                f"Chromium exited before CDP connection (exit {process.returncode})"
            )
        try:
            with urlopen(f"{endpoint}/json/version", timeout=0.5):
                browser = playwright.chromium.connect_over_cdp(endpoint)
                return browser, process, temp_profile
        except OSError:
            time.sleep(0.05)
    process.terminate()
    process.wait(timeout=5)
    temp_profile.cleanup()
    raise RuntimeError("Chromium CDP endpoint did not become ready")


def page_html(
    arm: str,
    smiles: list[str],
    warmup: int,
    fingerprint_excluded_smiles: frozenset[str],
) -> str:
    config = json.dumps(
        {
            "arm": arm,
            "smiles": smiles,
            "warmup": warmup,
            "fingerprint_excluded_smiles": sorted(fingerprint_excluded_smiles),
        }
    ).replace("<", "\\u003c")
    rdkit_script = (
        '<script src="/rdkit/RDKit_minimal.js"></script>' if arm == "rdkit" else ""
    )
    return f"""<!doctype html><meta charset="utf-8">{rdkit_script}<pre id="result">running</pre>
<script type="module">
const config = {config};
const memory = () => performance.memory ? {{status: "measured", used_js_heap_bytes: performance.memory.usedJSHeapSize, total_js_heap_bytes: performance.memory.totalJSHeapSize, js_heap_limit_bytes: performance.memory.jsHeapSizeLimit}} : {{status: "unavailable"}};
const summarize = (values) => {{
  const ordered = [...values].sort((a, b) => a - b);
  const percentile = (fraction) => {{ const pos = (ordered.length - 1) * fraction; const low = Math.floor(pos); const high = Math.min(low + 1, ordered.length - 1); return ordered[low] + (ordered[high] - ordered[low]) * (pos - low); }};
  return {{count: values.length, p50_ms: percentile(.5), p95_ms: percentile(.95), mean_ms: values.reduce((sum, x) => sum + x, 0) / values.length}};
}};
const time = (fn, values = config.smiles) => {{
  for (const value of values.slice(0, config.warmup)) fn(value);
  const samples = [];
  for (const value of values) {{ const start = performance.now(); fn(value); samples.push(performance.now() - start); }}
  return {{...summarize(samples), input_rows: values.length}};
}};
const MORGAN_OPTIONS = JSON.stringify({{radius: 2, nBits: 2048}});
const FP_BYTES = 256;
const fnv1a32 = (hash, bytes) => {{
  for (const byte of bytes) {{ hash ^= byte; hash = Math.imul(hash, 0x01000193); }}
  return hash >>> 0;
}};
const requirePackedFp = (bytes) => {{
  if (!(bytes instanceof Uint8Array) || bytes.length !== FP_BYTES) {{
    throw new Error(`expected packed ${{FP_BYTES}}-byte Morgan fingerprint`);
  }}
  return bytes;
}};
const rdkitBitsToPacked = (bits) => {{
  if (typeof bits !== "string" || bits.length !== FP_BYTES * 8) {{
    throw new Error(`RDKit returned ${{typeof bits}} length ${{bits?.length}}, expected ${{FP_BYTES * 8}} bits`);
  }}
  const packed = new Uint8Array(FP_BYTES);
  for (let bit = 0; bit < bits.length; bit += 1) {{
    const value = bits.charCodeAt(bit) - 48;
    if (value !== 0 && value !== 1) throw new Error(`RDKit fingerprint has non-bit at ${{bit}}`);
    packed[bit >> 3] |= value << (bit & 7);
  }}
  return packed;
}};
const timeFingerprint = (fn, values) => {{
  for (const value of values.slice(0, config.warmup)) requirePackedFp(fn(value));
  let digest = 0x811c9dc5;
  let setBits = 0;
  const samples = [];
  for (const value of values) {{
    const start = performance.now();
    const packed = requirePackedFp(fn(value));
    samples.push(performance.now() - start);
    digest = fnv1a32(digest, packed);
    for (const byte of packed) setBits += byte.toString(2).split("1").length - 1;
  }}
  return {{...summarize(samples), input_rows: values.length, output: {{format: "packed-lsb-first-2048-bit", bytes_per_row: FP_BYTES, fnv1a32: digest.toString(16).padStart(8, "0"), set_bits: setBits}}}};
}};
const withTimeout = (promise, label) => Promise.race([promise, new Promise((_, reject) => setTimeout(() => reject(new Error(label + " timed out")), 30000))]);
const run = async () => {{
  const before = memory();
  const fingerprintRows = config.smiles.filter((value) => !config.fingerprint_excluded_smiles.includes(value));
  if (config.arm === "chematic") {{
    const mod = await withTimeout(import("/schematic/chematic_wasm.js"), "chematic module import");
    const initialized = performance.now();
    await withTimeout(mod.default("/schematic/chematic_wasm_bg.wasm"), "chematic initialization");
    const init_ms = performance.now() - initialized;
    // This starts at navigationStart, before the HTML, JS module, and Wasm
    // assets are requested. The route is cache-disabled and each repetition
    // launches a new browser process, so the value is a cold local
    // download-to-ready measurement rather than init-only timing.
    const download_to_ready_ms = performance.now();
    const hasLinearMemoryMetric = typeof mod.wasm_linear_memory_bytes === "function";
    const linearMemoryAfterInit = hasLinearMemoryMetric ? mod.wasm_linear_memory_bytes() : null;
    const parse = time((value) => {{ const mol = mod.parse_smiles(value); mol.free(); }});
    const parse_write = time((value) => {{ const mol = mod.parse_smiles(value); mol.canonical_smiles(); mol.free(); }});
    const parse_fp = timeFingerprint((value) => {{
      const mol = mod.parse_smiles(value);
      let prepared;
      try {{ prepared = mod.prepare_rdkit_ecfp4(mol); return prepared.bitvec(); }}
      finally {{ if (prepared) prepared.free(); mol.free(); }}
    }}, fingerprintRows);
    const preparedMols = fingerprintRows.map((value) => {{
      const mol = mod.parse_smiles(value);
      try {{ return mod.prepare_rdkit_ecfp4(mol); }} finally {{ mol.free(); }}
    }});
    let prepared_fp;
    try {{
      prepared_fp = timeFingerprint((mol) => mol.bitvec(), preparedMols);
    }} finally {{
      for (const mol of preparedMols) mol.free();
    }}
    document.querySelector("#result").textContent = JSON.stringify({{arm: config.arm, init_ms, download_to_ready_ms, operations: {{parse, parse_write, parse_fp, prepared_fp}}, memory: {{before_load: before, after_workload: memory()}}, wasm_linear_memory: hasLinearMemoryMetric ? {{status: "measured", after_init_bytes: linearMemoryAfterInit, after_workload_bytes: mod.wasm_linear_memory_bytes()}} : {{status: "unavailable", reason: "the selected chematic artifact does not export wasm_linear_memory_bytes"}}, accepted_rows: config.smiles.length}});
    return;
  }}
  const initialized = performance.now();
  const rdkit = await withTimeout(window.initRDKitModule({{locateFile: (name) => name.endsWith(".wasm") ? "/rdkit/RDKit_minimal.wasm" : name}}), "RDKit initialization");
  const init_ms = performance.now() - initialized;
  const download_to_ready_ms = performance.now();
  const getMol = (value) => {{ const mol = rdkit.get_mol(value); if (!mol) throw new Error("RDKit rejected input"); return mol; }};
  const parse = time((value) => getMol(value).delete());
  const parse_write = time((value) => {{ const mol = getMol(value); mol.get_smiles(); mol.delete(); }});
  const parse_fp = timeFingerprint((value) => {{ const mol = getMol(value); try {{ return rdkitBitsToPacked(mol.get_morgan_fp(MORGAN_OPTIONS)); }} finally {{ mol.delete(); }} }}, fingerprintRows);
  const preparedMols = fingerprintRows.map((value) => getMol(value));
  let prepared_fp;
  try {{
    prepared_fp = timeFingerprint((mol) => rdkitBitsToPacked(mol.get_morgan_fp(MORGAN_OPTIONS)), preparedMols);
  }} finally {{
    for (const mol of preparedMols) mol.delete();
  }}
  document.querySelector("#result").textContent = JSON.stringify({{arm: config.arm, rdkit_version: rdkit.version(), init_ms, download_to_ready_ms, operations: {{parse, parse_write, parse_fp, prepared_fp}}, memory: {{before_load: before, after_workload: memory()}}, wasm_linear_memory: {{status: "unavailable", reason: "RDKit.js MinimalLib does not expose its WebAssembly.Memory object"}}, accepted_rows: config.smiles.length}});
}};
run().catch((error) => {{ document.querySelector("#result").textContent = JSON.stringify({{error: String(error), stack: error.stack}}); }});
</script>"""


class Handler(BaseHTTPRequestHandler):
    files: dict[str, Path] = {}
    smiles: list[str] = []
    warmup: int = 20
    fingerprint_excluded_smiles: frozenset[str] = frozenset()

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        if parsed.path == "/":
            arm = parse_qs(parsed.query).get("arm", [""])[0]
            if arm not in {"chematic", "rdkit"}:
                self.send_error(400, "arm must be chematic or rdkit")
                return
            body = page_html(
                arm, self.smiles, self.warmup, self.fingerprint_excluded_smiles
            ).encode()
            content_type = "text/html"
        elif parsed.path in self.files:
            body = self.files[parsed.path].read_bytes()
            content_type = (
                "application/wasm"
                if parsed.path.endswith(".wasm")
                else "text/javascript"
            )
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rdkit-package", type=Path, required=True)
    parser.add_argument(
        "--engine", choices=("chromium", "firefox", "webkit"), default="chromium"
    )
    parser.add_argument(
        "--browser",
        type=Path,
        help="explicit Chromium executable; required only with --engine chromium",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--schematic-dir", type=Path, default=ROOT / "demo" / "pkg")
    parser.add_argument(
        "--schematic-tarball",
        type=Path,
        help="optional npm tarball whose extracted package is --schematic-dir; records SHA-256 without inferring registry provenance",
    )
    parser.add_argument(
        "--schematic-package-kind",
        choices=("published", "source_candidate"),
        default="published",
        help="classify the chematic artifact without inferring publication from package metadata",
    )
    parser.add_argument("--schematic-source-revision")
    parser.add_argument("--schematic-source-diff-sha256")
    parser.add_argument(
        "--corpus", type=Path, default=ROOT / "scripts" / "descriptor_census_corpus.smi"
    )
    parser.add_argument("--rows", type=int, default=1000)
    parser.add_argument("--warmup", type=int, default=20)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument(
        "--measure-process-rss",
        action="store_true",
        help="launch Chromium through CDP and sample its process-tree RSS; unsupported for Firefox/WebKit",
    )
    args = parser.parse_args()
    if args.rows <= 0 or args.warmup < 0 or args.repetitions <= 0:
        raise SystemExit("rows/repetitions must be positive and warmup non-negative")
    if args.engine == "chromium" and args.browser is None:
        parser.error("--browser is required with --engine chromium")
    if args.measure_process_rss and args.engine != "chromium":
        parser.error("--measure-process-rss is supported only with --engine chromium")

    corpus_bytes = args.corpus.read_bytes()
    smiles = [
        line.strip() for line in corpus_bytes.decode().splitlines() if line.strip()
    ][: args.rows]
    if len(smiles) != args.rows:
        raise SystemExit(f"corpus has {len(smiles)} usable rows, expected {args.rows}")
    rdkit_dist = args.rdkit_package / "dist"
    Handler.files = {
        "/schematic/chematic_wasm.js": args.schematic_dir / "chematic_wasm.js",
        "/schematic/chematic_wasm_bg.wasm": args.schematic_dir
        / "chematic_wasm_bg.wasm",
        "/rdkit/RDKit_minimal.js": rdkit_dist / "RDKit_minimal.js",
        "/rdkit/RDKit_minimal.wasm": rdkit_dist / "RDKit_minimal.wasm",
    }
    missing = [str(path) for path in Handler.files.values() if not path.is_file()]
    if missing:
        raise SystemExit(f"missing benchmark artifacts: {missing}")
    Handler.smiles, Handler.warmup = smiles, args.warmup
    Handler.fingerprint_excluded_smiles = (
        FINGERPRINT_TYPED_UNSUPPORTED_SMILES.intersection(smiles)
    )
    package_json = args.rdkit_package / "package.json"
    package = json.loads(package_json.read_text(encoding="utf-8"))
    schematic_package_json = args.schematic_dir / "package.json"
    schematic_package = (
        json.loads(schematic_package_json.read_text(encoding="utf-8"))
        if schematic_package_json.is_file()
        else None
    )
    if args.schematic_tarball is not None and not args.schematic_tarball.is_file():
        raise SystemExit(f"schematic tarball does not exist: {args.schematic_tarball}")
    if args.schematic_package_kind == "source_candidate" and (
        not args.schematic_source_revision or not args.schematic_source_diff_sha256
    ):
        parser.error(
            "source candidates require --schematic-source-revision and "
            "--schematic-source-diff-sha256"
        )
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as exc:  # pragma: no cover - host-dependent tool
        server.shutdown()
        raise SystemExit(f"Python Playwright is required: {exc}") from exc
    results: dict[str, list[dict[str, object]]] = {"chematic": [], "rdkit": []}
    port = server.server_address[1]
    try:
        with sync_playwright() as playwright:
            for repetition in range(args.repetitions):
                for arm in (
                    ("chematic", "rdkit")
                    if repetition % 2 == 0
                    else ("rdkit", "chematic")
                ):
                    launched_process = None
                    temp_profile = None
                    monitor = None
                    if args.measure_process_rss:
                        browser, launched_process, temp_profile = (
                            launch_chromium_with_rss(playwright, args.browser)
                        )
                        monitor = ProcessTreeRssMonitor(launched_process.pid)
                        monitor.start()
                    else:
                        launch_options = {"headless": True}
                        if args.browser is not None:
                            launch_options["executable_path"] = str(args.browser)
                        browser = getattr(playwright, args.engine).launch(
                            **launch_options
                        )
                    try:
                        page = browser.new_page()
                        page.goto(
                            f"http://127.0.0.1:{port}/?arm={arm}",
                            wait_until="domcontentloaded",
                        )
                        page.wait_for_function(
                            'document.querySelector("#result").textContent !== "running"',
                            timeout=60_000,
                        )
                        result = json.loads(
                            page.locator("#result").text_content() or "{}"
                        )
                        if result.get("error") or result.get("arm") != arm:
                            raise RuntimeError(f"{arm} browser run failed: {result}")
                        result["process_tree_rss"] = (
                            monitor.finish()
                            if monitor is not None
                            else {
                                "status": "not_measured",
                                "reason": "--measure-process-rss not requested",
                            }
                        )
                        result["repetition"] = repetition
                        results[arm].append(result)
                    finally:
                        browser.close()
                        if launched_process is not None:
                            if launched_process.poll() is None:
                                launched_process.terminate()
                                try:
                                    launched_process.wait(timeout=5)
                                except subprocess.TimeoutExpired:
                                    launched_process.kill()
                                    launched_process.wait(timeout=5)
                            assert temp_profile is not None
                            temp_profile.cleanup()
    finally:
        server.shutdown()
        thread.join()

    # The timed path deliberately materializes the same packed representation in
    # both arms. Reject a run if its aggregate output disagrees before emitting a
    # timing artifact: a fast no-op, a changed bit order, or a comparator option
    # drift must not become a performance result. The separate parity gate still
    # provides the row-level diagnostic for any failure here.
    by_repetition = {
        arm: {int(run["repetition"]): run for run in arm_runs}
        for arm, arm_runs in results.items()
    }
    for repetition in range(args.repetitions):
        left = by_repetition["chematic"].get(repetition)
        right = by_repetition["rdkit"].get(repetition)
        if left is None or right is None:
            raise RuntimeError(f"missing arm result for repetition {repetition}")
        for operation in ("parse_fp", "prepared_fp"):
            left_output = left["operations"][operation]["output"]
            right_output = right["operations"][operation]["output"]
            if left_output != right_output:
                raise RuntimeError(
                    "fingerprint output contract failed for repetition "
                    f"{repetition}, operation {operation}: "
                    f"chematic={left_output}, rdkit={right_output}; "
                    "run the row-level browser parity gate for diagnostics"
                )

    def aggregate(arm: str) -> dict[str, object]:
        runs = results[arm]
        return {
            "runs": len(runs),
            "download_to_ready_ms": summary(
                [float(run["download_to_ready_ms"]) for run in runs]
            ),
            "init_ms": summary([float(run["init_ms"]) for run in runs]),
            "parse_mean_ms": summary(
                [float(run["operations"]["parse"]["mean_ms"]) for run in runs]
            ),
            "parse_write_mean_ms": summary(
                [float(run["operations"]["parse_write"]["mean_ms"]) for run in runs]
            ),
            "parse_fp_mean_ms": summary(
                [float(run["operations"]["parse_fp"]["mean_ms"]) for run in runs]
            ),
            "prepared_fp_mean_ms": summary(
                [float(run["operations"]["prepared_fp"]["mean_ms"]) for run in runs]
            ),
        }

    document = {
        "schema_version": 2,
        "gate": "wasm-vs-official-rdkit-browser-isolated",
        "configuration": {
            "engine": args.engine,
            "browser": str(args.browser) if args.browser else None,
            "rows": args.rows,
            "warmup_rows": args.warmup,
            "repetitions": args.repetitions,
            "timing": "fresh browser process per arm/repetition; no-store local routes",
            "download_to_ready_contract": "navigation start through JavaScript module/script fetch and WebAssembly initialization; excludes the subsequent parse/write/fingerprint workload",
            "process_rss": (
                "summed fresh Chromium process-tree RSS sampled every 50 ms"
                if args.measure_process_rss
                else "not_measured"
            ),
            "operation_contract": "parse; parse+canonical-SMILES write; parse+RDKit-compatible preparation+radius-2/2048-bit fingerprint; prepared-molecule radius-2/2048-bit fingerprint. Both arms perform sanitization-compatible preparation outside prepared_fp timing.",
            "fingerprint_output_contract": "both arms produce and consume a 256-byte LSB-first packed radius-2/2048-bit fingerprint; RDKit's bit string is packed inside the timed operation; a per-run FNV-1a checksum and total set-bit count are retained",
            "fingerprint_typed_unsupported_rows": len(
                Handler.fingerprint_excluded_smiles
            ),
            "fingerprint_typed_unsupported_smiles": sorted(
                Handler.fingerprint_excluded_smiles
            ),
        },
        "corpus": {
            "path": str(args.corpus),
            "sha256": hashlib.sha256(corpus_bytes).hexdigest(),
            "rows": args.rows,
        },
        "artifacts": {
            "schematic_wasm": artifact(
                Handler.files["/schematic/chematic_wasm_bg.wasm"]
            ),
            "schematic_js": artifact(Handler.files["/schematic/chematic_wasm.js"]),
            "schematic_package": (
                {
                    "name": schematic_package.get("name"),
                    "version": schematic_package.get("version"),
                    "kind": args.schematic_package_kind,
                    "source_revision": args.schematic_source_revision,
                    "source_diff_sha256": args.schematic_source_diff_sha256,
                    **npm_package_provenance(
                        args.schematic_dir,
                        schematic_package_json,
                        str(schematic_package.get("name", "")),
                    ),
                    "tarball": (
                        artifact(args.schematic_tarball)
                        if args.schematic_tarball is not None
                        else {"status": "not_supplied"}
                    ),
                }
                if schematic_package is not None
                else {
                    "status": "not_packaged",
                    "detail": "--schematic-dir has no package.json; source artifact only",
                }
            ),
            "rdkit_wasm": artifact(Handler.files["/rdkit/RDKit_minimal.wasm"]),
            "rdkit_js": artifact(Handler.files["/rdkit/RDKit_minimal.js"]),
            "rdkit_package": {
                "name": package.get("name"),
                "version": package.get("version"),
                **npm_package_provenance(
                    args.rdkit_package, package_json, str(package.get("name", ""))
                ),
            },
        },
        "raw_runs": results,
        "aggregate": {"chematic": aggregate("chematic"), "rdkit": aggregate("rdkit")},
        "boundary": (
            "Published package artifact. "
            if args.schematic_package_kind == "published"
            else "Locally built source candidate; not a registry release. "
        )
        + "Both arms receive the same valid SMILES rows for parse and write. The fingerprint operation excludes only its separately recorded typed RDKit coordination-sanitization refusal from both arms, so its timing denominator is explicit rather than silently filtered. Both timed fingerprint paths materialize and consume the same packed-byte representation, but this artifact alone does not establish bit-identical chemistry; use the separate fingerprint compatibility gate for that claim. JS heap is browser-exposed only. chematic linear memory is measured from the active Wasm memory page count; RDKit.js does not expose that object. When requested on Chromium, RSS is sampled as summed process-tree RSS and may double-count shared pages; otherwise it is not measured.",
    }
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(document["aggregate"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
