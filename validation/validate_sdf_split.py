#!/usr/bin/env python3
"""Original mixed CTAB records, exact raw concatenation and bounded-memory probe."""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
import tracemalloc

import chematic
from rdkit import Chem
import rdkit


def record(index):
    title = f"record_{index}"
    if index % 2:
        block = f"""{title}
  independent    2D
original V3000
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 0 0 0 0
M  V30 2 H 0.9572 0 0 0
M  V30 3 H -0.2399872 0.927297 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 END BOND
M  V30 END CTAB
M  END
"""
    else:
        block = f"""{title}
  independent    2D
original V2000
  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9572    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2400    0.9273    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
"""
    return (
        block
        + f"> <INDEX>\n{index}\n\n> <opaque>\ninline $$$$ value\n未知\n\n> <repeated>\nfirst\n\n> <repeated>\nsecond\n\n$$$$\n"
    ).encode()


def main():
    cases = parts_checked = native_records = rdkit_records = 0
    obabel = "/opt/homebrew/bin/obabel"
    with tempfile.TemporaryDirectory(prefix="chematic-sdf-split-probe-") as folder:
        root = Path(folder)
        path = root / "source.sdf"
        for count in [1, 2, 3, 9, 10, 17, 101, 1003]:
            for newline in [b"\n", b"\r\n"]:
                data = b"".join(record(i) for i in range(count))
                if count % 2:
                    data = data.removesuffix(b"$$$$\n")
                data = data.replace(b"\n", newline)
                path.write_bytes(data)
                for n in sorted({1, min(2, count), min(7, count), min(400, count)}):
                    plan = chematic.plan_sdf_split(path, n)
                    joined = hashlib.sha256()
                    first = 0
                    assert [p.record_count for p in plan.parts] == [
                        count // n + int(i < count % n) for i in range(n)
                    ]
                    for i, part in enumerate(plan.parts):
                        raw = b"".join(plan.iter_part(i, chunk_bytes=173))
                        joined.update(raw)
                        assert (
                            raw
                            == data[part.byte_offset : part.byte_offset + part.byte_len]
                        )
                        native = list(chematic.iter_sdf_str(raw.decode()))
                        oracle = list(
                            Chem.ForwardSDMolSupplier(io.BytesIO(raw), removeHs=False)
                        )
                        assert len(native) == len(oracle) == part.record_count
                        assert [r.name for r in native] == [
                            f"record_{j}"
                            for j in range(first, first + part.record_count)
                        ]
                        assert [r.get("INDEX") for r in native] == [
                            str(j) for j in range(first, first + part.record_count)
                        ]
                        assert all(
                            m is not None and m.GetNumAtoms() == 3 for m in oracle
                        )
                        assert [m.GetProp("_Name") for m in oracle] == [
                            r.name for r in native
                        ]
                        native_records += len(native)
                        rdkit_records += len(oracle)
                        first += part.record_count
                        parts_checked += 1
                    assert (
                        first == count
                        and joined.hexdigest()
                        == plan.source_sha256
                        == hashlib.sha256(data).hexdigest()
                    )
                    cases += 1
        # The upstream command's current behavior is captured, not assumed.
        path.write_bytes(b"".join(record(i) for i in range(7)))
        p = subprocess.run(
            [
                obabel,
                "-isdf",
                str(path),
                "-O",
                str(root / "baseline.sdf"),
                "--splitinto",
                "3",
            ],
            capture_output=True,
            text=True,
        )
        outputs = sorted(f.name for f in root.glob("baseline*.sdf"))
        baseline = {
            "returncode": p.returncode,
            "filenames": outputs,
            "stderr": p.stderr.strip(),
            "requested_parts": 3,
        }
        # External reads of all 3 candidate files, preserving record order.
        plan = chematic.plan_sdf_split(path, 3)
        external_titles = []
        for i, part in enumerate(plan.parts):
            target = root / part.filename
            target.write_bytes(b"".join(plan.iter_part(i)))
            p = subprocess.run(
                [obabel, "-isdf", str(target), "-osmi"],
                capture_output=True,
                text=True,
                check=True,
            )
            external_titles.extend(
                line.split()[-1] for line in p.stdout.splitlines() if line.strip()
            )
        assert external_titles == [f"record_{j}" for j in range(7)]
        # Generate the large source outside the measured region. No native graph
        # allocation participates in this Python allocation measurement.
        large = root / "large.sdf"
        big_record = (
            record(0).removesuffix(b"$$$$\n")
            + b"> <payload>\n"
            + b"x" * 8192
            + b"\n\n$$$$\n"
        )
        with large.open("wb") as file:
            for _ in range(4096):
                file.write(big_record)
        tracemalloc.start()
        bigplan = chematic.plan_sdf_split(large, 3)
        copied = 0
        digest = hashlib.sha256()
        maximum_chunk = 0
        for i in range(3):
            for chunk in bigplan.iter_part(i, chunk_bytes=65536):
                copied += len(chunk)
                digest.update(chunk)
                maximum_chunk = max(maximum_chunk, len(chunk))
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        assert (
            copied == large.stat().st_size
            and digest.hexdigest() == bigplan.source_sha256
        )
        assert maximum_chunk <= 65536 and peak < 8 * 1024 * 1024
        memory = {
            "source_bytes": copied,
            "records": 4096,
            "parts": 3,
            "tracemalloc_peak_bytes": peak,
            "max_yielded_chunk_bytes": maximum_chunk,
            "source_sha256": bigplan.source_sha256,
            "scope": "Python traced allocations, not process RSS",
        }
    print(
        json.dumps(
            {
                "cases": cases,
                "parts_checked": parts_checked,
                "native_record_reads": native_records,
                "rdkit_record_reads": rdkit_records,
                "openbabel_part_reads": 3,
                "openbabel_record_reads": 7,
                "openbabel_baseline": baseline,
                "memory_probe": memory,
                "rdkit_version": rdkit.__version__,
                "openbabel_version": subprocess.run(
                    [obabel, "-V"], capture_output=True, text=True, check=True
                ).stdout.strip(),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
