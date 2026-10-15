#!/usr/bin/env python3
"""Generate all distinct-ligand TB4 and OH5 permutations with one vacant site.

Run the pinned RDKit/COSMolKit read-back comparison on the output with
``check_native_smiles_roundtrip.py``. The comparison script gates each source
on RDKit acceptance and reports both engine scores and artifact hashes.
"""

from __future__ import annotations

import argparse
from itertools import permutations
from pathlib import Path


def rows() -> list[str]:
    result = []
    for f, cl, br, i in permutations(("F", "Cl", "Br", "I")):
        for tag in range(1, 21):
            result.append(f"{f}[P@TB{tag}]({cl})({br}){i}")
    for f, cl, br, i, n in permutations(("F", "Cl", "Br", "I", "N")):
        for tag in range(1, 31):
            result.append(f"{f}[S@OH{tag}]({cl})({br})({i}){n}")
    assert len(result) == 4080
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.write_text("\n".join(rows()) + "\n")


if __name__ == "__main__":
    main()
