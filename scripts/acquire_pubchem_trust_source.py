#!/usr/bin/env python3
"""Acquire the fixed PubChem source universe used by trust holdout v1."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

CID_FIRST = 4097
CID_LAST = 20000
CHUNK = 100
URL_TEMPLATE = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/"
    "property/IsomericSMILES/JSON"
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cid-first", type=int, default=CID_FIRST)
    parser.add_argument("--cid-last", type=int, default=CID_LAST)
    args = parser.parse_args()
    requests = []
    records: dict[int, str] = {}
    if args.cid_first < 1 or args.cid_last < args.cid_first:
        raise SystemExit("invalid CID range")
    cids = list(range(args.cid_first, args.cid_last + 1))
    for offset in range(0, len(cids), CHUNK):
        chunk = cids[offset:offset + CHUNK]
        url = URL_TEMPLATE.format(cids=",".join(map(str, chunk)))
        raw = subprocess.check_output(
            ["curl", "--fail", "--silent", "--show-error", "--location", url],
            timeout=60,
        )
        payload = json.loads(raw)
        found = payload.get("PropertyTable", {}).get("Properties", [])
        for row in found:
            smiles = row.get("SMILES") or row.get("IsomericSMILES")
            if smiles:
                records[int(row["CID"])] = smiles
        requests.append({"cid_first": chunk[0], "cid_last": chunk[-1], "sha256": sha(raw)})
        time.sleep(0.2)
    result = {
        "schema": "chematic.pubchem-holdout-source.v1",
        "acquired_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "url_template": URL_TEMPLATE,
        "cid_first": args.cid_first,
        "cid_last": args.cid_last,
        "requests": requests,
        "records": [{"cid": cid, "smiles": records[cid]} for cid in sorted(records)],
    }
    args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
    print(f"wrote {len(records)} PubChem records to {args.output}")


if __name__ == "__main__":
    main()
