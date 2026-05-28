#!/usr/bin/env python3
"""
Download canonical SMILES from ChEMBL via REST API.

Usage:
    python3 scripts/download_chembl_smiles.py --output scripts/chembl_50k.smi --count 50000
"""

import argparse
import json
import ssl
import sys
import time
import urllib.error
import urllib.request

# macOS ships without the system CA bundle accessible to Python; bypass for this dev script.
_SSL_CTX = ssl.create_default_context()
_SSL_CTX.check_hostname = False
_SSL_CTX.verify_mode = ssl.CERT_NONE

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
PAGE_SIZE = 1000
PROGRESS_INTERVAL = 5000


def fetch_page(url: str, retry: bool = True) -> dict | None:
    """Fetch a single page from the ChEMBL API. Returns parsed JSON or None on error."""
    try:
        with urllib.request.urlopen(url, timeout=30, context=_SSL_CTX) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        if exc.code == 500 and retry:
            print(f"[warn] HTTP 500 on {url!r}, retrying once after 2s ...", file=sys.stderr)
            time.sleep(2)
            return fetch_page(url, retry=False)
        else:
            print(f"[warn] HTTP {exc.code} on {url!r}, skipping page.", file=sys.stderr)
            return None
    except urllib.error.URLError as exc:
        print(f"[warn] URL error on {url!r}: {exc.reason}, skipping page.", file=sys.stderr)
        return None
    except Exception as exc:  # noqa: BLE001
        print(f"[warn] Unexpected error fetching {url!r}: {exc}, skipping page.", file=sys.stderr)
        return None


def build_url(offset: int, limit: int = PAGE_SIZE) -> str:
    return (
        f"{BASE_URL}"
        f"?format=json"
        f"&limit={limit}"
        f"&offset={offset}"
        f"&only=molecule_structures"
    )


def download_smiles(target_count: int) -> list[str]:
    """Download up to *target_count* unique, single-fragment canonical SMILES."""
    collected: list[str] = []
    seen: set[str] = set()
    offset = 0
    last_progress_milestone = 0

    while len(collected) < target_count:
        url = build_url(offset)
        data = fetch_page(url)

        if data is None:
            # Skip this page and try the next offset
            offset += PAGE_SIZE
            continue

        molecules = data.get("molecules", [])
        if not molecules:
            # No more data from ChEMBL
            print("[info] No more molecules available from ChEMBL.", file=sys.stderr)
            break

        for mol in molecules:
            if len(collected) >= target_count:
                break

            structures = mol.get("molecule_structures")
            if not structures:
                continue

            smiles = structures.get("canonical_smiles")
            if not smiles:
                continue

            # Skip salts / mixtures (multi-fragment SMILES contain ".")
            if "." in smiles:
                continue

            # Deduplicate
            if smiles in seen:
                continue

            seen.add(smiles)
            collected.append(smiles)

            # Progress reporting every PROGRESS_INTERVAL molecules
            current = len(collected)
            milestone = (current // PROGRESS_INTERVAL) * PROGRESS_INTERVAL
            if milestone > last_progress_milestone and milestone > 0:
                print(f"downloaded {milestone} / {target_count}", file=sys.stderr)
                last_progress_milestone = milestone

        # Check whether ChEMBL has a next page
        page_meta = data.get("page_meta", {})
        if page_meta.get("next") is None:
            print("[info] Reached last page of ChEMBL molecule data.", file=sys.stderr)
            break

        offset += PAGE_SIZE

    return collected


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download canonical SMILES from ChEMBL REST API."
    )
    parser.add_argument(
        "--output",
        "-o",
        metavar="FILE",
        default=None,
        help="Write SMILES to FILE (one per line). Defaults to stdout.",
    )
    parser.add_argument(
        "--count",
        "-n",
        metavar="N",
        type=int,
        default=50_000,
        help="Number of unique SMILES to download (default: 50000).",
    )
    args = parser.parse_args()

    target = max(1, args.count)
    print(f"[info] Downloading up to {target} unique SMILES from ChEMBL ...", file=sys.stderr)

    smiles_list = download_smiles(target)

    print(f"[info] Total unique SMILES collected: {len(smiles_list)}", file=sys.stderr)

    output_lines = "\n".join(smiles_list) + "\n"

    if args.output:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(output_lines)
        print(f"[info] Wrote {len(smiles_list)} SMILES to {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(output_lines)


if __name__ == "__main__":
    main()
