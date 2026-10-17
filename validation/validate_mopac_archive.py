"""Deterministic archive corpus and Open Babel diagnostics; requires built chematic.

Run with candidate package on PYTHONPATH. Optional upstream zip is downloaded
separately from the issue; it is hashed and inspected, never copied into this repo.
"""

import argparse
import hashlib
import json
import math
import subprocess
import zipfile
from pathlib import Path

import chematic


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream-zip", type=Path)
    args = parser.parse_args()
    corpus_hash = hashlib.sha256()
    n = 0
    cross_reads = 0
    for internal in (False, True):
        for i in range(256):
            charge = i % 2 == 0
            marker = "FINAL GEOMETRY OBTAINED" + (" CHARGE" if charge else "")
            rows = []
            expected = []
            for j, element in enumerate(("O", "C", "H", "H")):
                values = (
                    [1.1 + j * 0.01, 120.0, float(i % 180)]
                    if internal
                    else [j * 1.1 + i * 0.01, math.sin(i + j), math.cos(i - j)]
                )
                refs = (
                    ([0, 0, 0], [1, 0, 0], [2, 1, 0], [2, 1, 3])[j]
                    if internal
                    else None
                )
                row = element + " " + " ".join(f"{v:.12g} +1" for v in values)
                if refs is not None:
                    row += " " + " ".join(map(str, refs))
                if charge:
                    row += " -0.125"
                rows.append(row)
                expected.append([float(f"{v:.12g}") for v in values])
            text = (
                f"Original archive {i}\n {marker}\nPM7\nSynthetic\n\n"
                + "\n".join(rows)
                + "\n\nopaque ancillary section\n"
            )
            if i % 3 == 0:
                text = text.replace("\n", "\r\n")
            corpus_hash.update(text.encode())
            a = chematic.parse_mopac_archive(text)
            assert a["source_text"] == text
            assert len(a["atoms"]) == 4
            assert [row["values"] for row in a["atoms"]] == expected
            assert [row["element"] for row in a["atoms"]] == ["O", "C", "H", "H"]
            assert all(
                row["optimization_flags"] == ["1", "1", "1"] for row in a["atoms"]
            )
            assert [row["references"] for row in a["atoms"]] == (
                [[0, 0, 0], [1, 0, 0], [2, 1, 0], [2, 1, 3]] if internal else [None] * 4
            )
            assert a["atoms"][3]["references"] == ([2, 1, 3] if internal else None)
            assert a["atoms"][3]["partial_charge"] == (-0.125 if charge else None)
            assert "opaque ancillary section" in a["trailing_text"]
            n += 1
            if not internal and i % 16 == 1:
                # Independent coordinate check on the extracted Cartesian deck.
                deck = "PM7\nSynthetic\n\n" + "\n".join(rows) + "\n"
                proc = subprocess.run(
                    ["obabel", "-imopcrt", "-oxyz"],
                    input=deck,
                    text=True,
                    capture_output=True,
                    check=True,
                )
                xyz = proc.stdout.splitlines()
                assert int(xyz[0]) == 4, proc.stderr
                for line, values in zip(xyz[2:], expected):
                    got = list(map(float, line.split()[1:4]))
                    assert max(abs(x - y) for x, y in zip(got, values)) < 1e-4
                cross_reads += 1
    report = {
        "scope": "archive field preservation; no coordinate conversion or MOPAC execution",
        "corpus_records": n,
        "corpus_sha256": corpus_hash.hexdigest(),
        "cartesian_openbabel_cross_reads": cross_reads,
        "openbabel_version": subprocess.check_output(
            ["obabel", "-V"], text=True
        ).strip(),
        "native_extension_sha256": hashlib.sha256(
            Path(chematic.chematic.__file__).read_bytes()
        ).hexdigest(),
    }
    if args.upstream_zip:
        with zipfile.ZipFile(args.upstream_zip) as z:
            raw = z.read("Example_data_set.arc")
        text = raw.decode()
        a = chematic.parse_mopac_archive(text)
        assert len(a["atoms"]) == 4
        assert a["atoms"][3]["references"] == [2, 1, 3]
        assert a["atoms"][3]["values"] == [1.09115836, 121.276097, 180.0]
        assert "3  1    4" in a["trailing_text"]
        proc = subprocess.run(
            ["obabel", "-imopout", "-oxyz"], input=text, text=True, capture_output=True
        )
        report["upstream_attachment"] = {
            "url": "https://github.com/openbabel/openbabel/files/4824357/MOPAC-example.zip",
            "zip_sha256": hashlib.sha256(args.upstream_zip.read_bytes()).hexdigest(),
            "archive_sha256": hashlib.sha256(raw).hexdigest(),
            "atom_count": len(a["atoms"]),
            "openbabel_mopout_returncode": proc.returncode,
            "openbabel_mopout_stdout": proc.stdout,
            "openbabel_mopout_stderr": proc.stderr,
        }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    run()
