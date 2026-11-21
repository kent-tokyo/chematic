#!/usr/bin/env python3
"""Select a reproducible structural sample from the published COD tar archive.

Selection depends only on member path hashes, never on coordinates or scoring.
Each selected member must contain exactly one SDF record. No archive paths are
extracted to disk. The resulting SDF preserves the supplied coordinate blocks.
"""
import argparse
import hashlib
import heapq
import json
import tarfile
from pathlib import Path


def digest_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--sdf", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--sample-size", type=int, default=12000)
    args = parser.parse_args()
    if args.sample_size < 1:
        parser.error("sample-size must be positive")
    with tarfile.open(args.archive, mode="r:") as archive:
        candidates = (member for member in archive if member.isfile() and member.name.endswith(".sdf"))
        selected = heapq.nsmallest(args.sample_size, candidates, key=lambda m: (hashlib.sha256(m.name.encode()).digest(), m.name))
        if len(selected) != args.sample_size:
            raise ValueError("archive contains fewer SDF members than requested")
        records, rejected = [], []
        with args.sdf.open("wb") as out:
            for member in selected:
                data = archive.extractfile(member).read()
                if data.count(b"$$$$") != 1 or b"\n" not in data:
                    rejected.append(member.name)
                    continue
                # Stable, unique record name without altering the structure block.
                name = Path(member.name).stem
                remainder = data.split(b"\n", 1)[1]
                out.write(name.encode("ascii") + b"\n" + remainder.rstrip() + b"\n")
                records.append(member.name)
    result = {"archive_sha256": digest_file(args.archive), "sample_sdf_sha256": digest_file(args.sdf), "rule": "lowest SHA256 of exact SDF member path, ties by path; preserve single-record structure blocks", "selected_members": len(selected), "written_records": len(records), "rejected_members": rejected, "member_paths": records}
    args.manifest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"selected {len(selected)}, written {len(records)}, rejected {len(rejected)}")


if __name__ == "__main__":
    main()
