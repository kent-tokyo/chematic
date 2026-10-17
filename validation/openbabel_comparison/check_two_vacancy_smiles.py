#!/usr/bin/env python3
"""Compare two-vacancy SMILES recognition with Open Babel 3.2.1.

RDKit 2025.09.3 supplies the isomeric identity and rereads native output.
Each input is passed to a fresh Open Babel process because its SMILES reader
stops on the first unsupported row of a multiline input stream.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(smiles) if smiles else None
    return Chem.MolToSmiles(molecule) if molecule is not None else None


def graph_identity(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(smiles) if smiles else None
    if molecule is None:
        return None
    Chem.RemoveStereochemistry(molecule)
    return Chem.MolToSmiles(molecule)


def openbabel_smiles(binary: Path, source: str) -> tuple[str, int]:
    result = subprocess.run(
        [str(binary), f"-:{source}", "-osmi"], capture_output=True,
        text=True, timeout=10, check=False,
    )
    output = result.stdout.split()[0] if result.stdout.strip() else ""
    return output, result.returncode


def classify(expected: str, output: str) -> str:
    if not output:
        return "rejected"
    actual = identity(output)
    if actual is None:
        return "invalid_output"
    if actual == expected:
        return "strict"
    if graph_identity(output) == graph_identity(expected):
        return "stereo_lost_or_changed"
    return "other_identity"


def chematic_outputs(binary: Path, corpus: Path, count: int) -> dict[int, str]:
    result = subprocess.run(
        [str(binary), str(corpus)], capture_output=True, text=True,
        timeout=120, check=True,
    )
    outputs = {}
    for line in result.stdout.splitlines():
        index_text, smiles = line.split("\t", 1)
        index = int(index_text)
        if index in outputs:
            raise ValueError(f"duplicate CheMatic row {index}")
        outputs[index] = smiles
    if set(outputs) != set(range(count)):
        raise ValueError("CheMatic row indices do not match the input corpus")
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--corpus", type=Path, action="append", required=True)
    parser.add_argument("--chematic", type=Path, required=True)
    parser.add_argument("--openbabel", type=Path,
                        default=Path(shutil.which("obabel") or "obabel"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")
    if Chem.rdBase.rdkitVersion != "2025.09.3":
        parser.error("requires RDKit 2025.09.3")
    openbabel = args.openbabel.resolve(strict=True)
    chematic = args.chematic.resolve(strict=True)
    version = subprocess.run(
        [str(openbabel), "-V"], capture_output=True, text=True, check=True,
    ).stdout.strip()
    if not version.startswith("Open Babel 3.2.1 "):
        parser.error("requires Open Babel 3.2.1")

    control = "F[C@](Cl)(Br)I"
    control_output, _ = openbabel_smiles(openbabel, control)
    if identity(control_output) != identity(control):
        raise ValueError("Open Babel tetrahedral-stereo control failed")

    corpora = []
    for path in args.corpus:
        corpus = path.resolve(strict=True)
        sources = [line.split()[0] for line in corpus.read_text().splitlines()
                   if line.strip()]
        native = chematic_outputs(chematic, corpus, len(sources))
        counts = {engine: Counter() for engine in ("chematic", "openbabel")}
        rows = []
        for index, source in enumerate(sources):
            expected = identity(source)
            if expected is None:
                raise ValueError(f"RDKit rejected input row {index}: {source}")
            babel_output, returncode = openbabel_smiles(openbabel, source)
            outputs = {"chematic": native[index], "openbabel": babel_output}
            classes = {engine: classify(expected, output)
                       for engine, output in outputs.items()}
            for engine, label in classes.items():
                counts[engine][label] += 1
            rows.append({"index": index, "source": source, "expected": expected,
                         "outputs": outputs, "classes": classes,
                         "openbabel_returncode": returncode})
        corpora.append({"path": str(corpus), "sha256": digest(corpus),
                       "count": len(sources), "counts": counts, "rows": rows})

    report = {"scope": "source-built two-vacancy SMILES development comparison",
              "oracle": "RDKit canonical isomeric SMILES read-back",
              "rdkit_version": Chem.rdBase.rdkitVersion,
              "openbabel_version": version, "openbabel_sha256": digest(openbabel),
              "chematic_sha256": digest(chematic),
              "script_sha256": digest(Path(__file__)),
              "tetrahedral_control": {"source": control, "openbabel": control_output},
              "corpora": corpora}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(report, sort_keys=True, separators=(",", ":")) + "\n").encode()
    with args.output.open("wb") as stream:
        with gzip.GzipFile(filename="", fileobj=stream, mode="wb", mtime=0) as zipped:
            zipped.write(payload)
    print(json.dumps({item["path"]: item["counts"] for item in corpora},
                     sort_keys=True))


if __name__ == "__main__":
    main()
