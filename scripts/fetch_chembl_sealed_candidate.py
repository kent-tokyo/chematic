#!/usr/bin/env python3
"""Acquire a locally held, provenance-recorded ChEMBL candidate cohort.

This does *not* declare a cohort sealed. It records an immutable fetch packet
outside the repository so ``prepare_sealed_accuracy_cohort.py`` can audit it
against exposed development data and requires a later candidate freeze plus
unused-data attestation before using the word "sealed".
"""

from __future__ import annotations

import argparse
import hashlib
import json
import ssl
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
LICENSE_URL = "https://creativecommons.org/licenses/by-sa/3.0/"


def verified_tls_context() -> ssl.SSLContext:
    """Use certifi when a Python framework lacks the macOS certificate store."""
    try:
        import certifi
    except ImportError:
        return ssl.create_default_context()
    return ssl.create_default_context(cafile=certifi.where())


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--offset", type=int, default=100_000)
    parser.add_argument("--source-rows", type=int, default=15_000, help="request this many API records before filtering")
    parser.add_argument("--page-size", type=int, default=100)
    parser.add_argument("--delay-seconds", type=float, default=0.25)
    parser.add_argument("--max-attempts", type=int, default=5)
    parser.add_argument("--retry-delay-seconds", type=float, default=2.0)
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    args = parser.parse_args()
    if (
        args.offset < 0
        or args.source_rows < 100
        or args.page_size <= 0
        or args.max_attempts <= 0
        or args.retry_delay_seconds < 0
        or args.timeout_seconds <= 0
    ):
        parser.error(
            "offset >= 0, source-rows >= 100, positive page-size/max-attempts, "
            "positive timeout-seconds, and non-negative retry-delay-seconds are required"
        )
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    response_cache = output / "response-cache"
    response_cache.mkdir(parents=True, exist_ok=True)
    tls_context = verified_tls_context()
    rows: list[str] = []
    responses: list[dict[str, object]] = []
    for page_offset in range(args.offset, args.offset + args.source_rows, args.page_size):
        params = urlencode({"limit": args.page_size, "offset": page_offset})
        url = f"{API}?{params}"
        request = Request(
            url,
            headers={"Accept": "application/json", "User-Agent": "CheMatic-validation/1.0"},
        )
        cache_path = response_cache / f"limit-{args.page_size}-offset-{page_offset}.json"
        if cache_path.is_file():
            body = cache_path.read_bytes()
        else:
            for attempt in range(1, args.max_attempts + 1):
                try:
                    with urlopen(  # nosec B310: fixed HTTPS endpoint
                        request, timeout=args.timeout_seconds, context=tls_context
                    ) as response:
                        body = response.read()
                    temporary = cache_path.with_suffix(".tmp")
                    temporary.write_bytes(body)
                    temporary.replace(cache_path)
                    break
                except (HTTPError, URLError, TimeoutError) as error:
                    print(
                        f"retry {attempt}/{args.max_attempts} for {url}: {error}",
                        file=sys.stderr,
                    )
                    if attempt == args.max_attempts:
                        raise
                    time.sleep(args.retry_delay_seconds * attempt)
        data = json.loads(body)
        molecules = data.get("molecules")
        if not isinstance(molecules, list):
            raise RuntimeError(f"unexpected ChEMBL response at offset {page_offset}")
        responses.append({"url": url, "sha256": digest(body), "records": len(molecules)})
        for molecule in molecules:
            structures = molecule.get("molecule_structures") or {}
            smiles = structures.get("canonical_smiles")
            if isinstance(smiles, str) and smiles and "." not in smiles:
                rows.append(smiles)
        if len(molecules) < args.page_size:
            break
        time.sleep(args.delay_seconds)
    source_path = output / "source.smi"
    source_path.write_text("".join(f"{smiles}\n" for smiles in rows), encoding="utf-8")
    metadata = {
        "source_url": API,
        "source_release": "ChEMBL REST API snapshot; page URLs and response SHA-256 values retained below",
        "license": "ChEMBL data: CC BY-SA 3.0",
        "license_url": LICENSE_URL,
        "acquired_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_sha256": digest(source_path.read_bytes()),
        "request_offset": args.offset,
        "requested_source_rows": args.source_rows,
        "accepted_single_fragment_rows": len(rows),
        "responses": responses,
        "redistribution": "local evaluation input only; raw source is not committed to this repository",
        "status": "acquired_not_sealed",
    }
    (output / "source-metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": metadata["status"], "accepted_rows": len(rows), "source_sha256": metadata["source_sha256"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
