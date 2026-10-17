#!/usr/bin/env python3
"""Validate original PDB model corpora and optional issue #2476 source data.

Run with the candidate package on PYTHONPATH. JSON is printed; temporary Open
Babel outputs stay in a disposable directory. No upstream fixture is vendored.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import chematic


def atom(serial, name, element, position):
    x, y, z = position
    return f"ATOM  {serial:5d} {name:4s} HOH A  -2I   {x:8.3f}{y:8.3f}{z:8.3f}  0.50 23.45          {element:>2s}  "


def dataset(count, offset, crlf=False):
    rows = [
        atom(7, " O1 ", "O", (offset, 0, -1.25)),
        atom(10, " H2 ", "H", (offset + 0.95, 0, -1.25)),
        atom(12, " H3 ", "H", (offset - 0.25, 0.91, -1.25)),
    ]
    block = "\n".join(rows) + "\n"
    text = "".join(f"MODEL     {1+2*i:4d}\n{block}ENDMDL\n" for i in range(count))
    text += "CONECT    7   10   12\nEND\n"
    return text.replace("\n", "\r\n") if crlf else text


def check_originals(obabel):
    cases = 0
    source_lines = 0
    crossreads = 0
    with tempfile.TemporaryDirectory(prefix="chematic-pdb-split-probe-") as folder:
        for count in [1, 2, 9, 10, 11, 19, 99, 100]:
            for offset in [-100.125, -1.25, 0, 17.375, 999.125]:
                for crlf in [False, True]:
                    source = dataset(count, offset, crlf)
                    r = chematic.split_pdb_models(source, prefix="trajectory")
                    assert r.source_text == source and not r.losses
                    assert len(r.files) == count
                    assert [f.filename for f in r.files] == sorted(
                        f.filename for f in r.files
                    )
                    width = len(str(count))
                    assert r.files[-1].filename == f"trajectory{count:0{width}d}.pdb"
                    for i, f in enumerate(r.files):
                        assert f.model_num == 1 + 2 * i and f.atom_count == 3
                        expected = [
                            line
                            for line in source.splitlines()
                            if line.startswith("ATOM  ")
                        ][:3]
                        assert [
                            line
                            for line in f.pdb.splitlines()
                            if line.startswith("ATOM  ")
                        ] == expected
                        assert "CONECT    7   10   12\n" in f.pdb
                        assert chematic.select_pdb_atoms(f.pdb)["pdb"] == f.pdb
                        source_lines += 3
                    if obabel and not crlf and offset == 17.375:
                        path = Path(folder) / "candidate.pdb"
                        path.write_text(r.files[-1].pdb)
                        p = subprocess.run(
                            [obabel, "-ipdb", str(path), "-oxyz"],
                            capture_output=True,
                            text=True,
                            check=True,
                        )
                        lines = p.stdout.splitlines()
                        assert lines[0].strip() == "3"
                        assert [line.split()[0] for line in lines[2:5]] == [
                            "O",
                            "H",
                            "H",
                        ]
                        expected_positions = [
                            (offset, 0, -1.25),
                            (offset + 0.95, 0, -1.25),
                            (offset - 0.25, 0.91, -1.25),
                        ]
                        for line, target in zip(lines[2:5], expected_positions):
                            assert (
                                max(
                                    abs(float(a) - b)
                                    for a, b in zip(line.split()[1:4], target)
                                )
                                <= 1e-5
                            )
                        crossreads += 1
                    cases += 1
    return {
        "cases": cases,
        "retained_atom_line_checks": source_lines,
        "openbabel_xyz_crossreads": crossreads,
    }


def check_real(path, obabel):
    source = path.read_text()
    try:
        chematic.split_pdb_models(source)
    except ValueError as error:
        default_refusal = str(error)
    else:
        raise AssertionError("metadata loss must require consent")
    r = chematic.split_pdb_models(source, allow_metadata_loss=True)
    original = [
        line
        for line in source.splitlines()
        if line.startswith(("ATOM  ", "HETATM", "ANISOU", "TER   "))
    ]
    retained = [
        line
        for f in r.files
        for line in f.pdb.splitlines()
        if line.startswith(("ATOM  ", "HETATM", "ANISOU", "TER   "))
    ]
    assert original == retained
    assert [f.model_num for f in r.files] == list(range(1, len(r.files) + 1))
    assert len(r.files) == 30 and all(f.atom_count == 1187 for f in r.files)
    baseline = None
    if obabel:
        with tempfile.TemporaryDirectory(
            prefix="chematic-pdb-split-openbabel-"
        ) as folder:
            subprocess.run(
                [
                    obabel,
                    "-ipdb",
                    str(path),
                    "-O",
                    str(Path(folder) / "split.pdb"),
                    "-m",
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            names = sorted(p.name for p in Path(folder).glob("split*.pdb"))
            assert (
                len(names) == 30
                and "split1.pdb" in names
                and "split01.pdb" not in names
            )
            baseline = {
                "file_count": len(names),
                "lexical_first_12": names[:12],
                "padding_gap_reproduced": True,
            }
    return {
        "source_url": "https://files.rcsb.org/download/1NER.pdb",
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "source_bytes": path.stat().st_size,
        "model_count": len(r.files),
        "atom_counts": [f.atom_count for f in r.files],
        "first_filename": r.files[0].filename,
        "last_filename": r.files[-1].filename,
        "loss_record_count": len(r.losses),
        "atom_and_auxiliary_lines_exact": original == retained,
        "default_refusal": default_refusal,
        "openbabel_baseline": baseline,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path)
    parser.add_argument("--obabel", default=shutil.which("obabel"))
    args = parser.parse_args()
    version = (
        subprocess.run(
            [args.obabel, "-V"], capture_output=True, text=True, check=True
        ).stdout.strip()
        if args.obabel
        else None
    )
    report = {"synthetic": check_originals(args.obabel), "openbabel_version": version}
    if args.pdb:
        report["issue_fixture"] = check_real(args.pdb, args.obabel)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
