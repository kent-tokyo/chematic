#!/usr/bin/env python3
"""Generate distinct-ligand TB3 and OH4 SMILES with two vacant sites."""

from __future__ import annotations

import argparse
from itertools import permutations
from pathlib import Path


def rows() -> list[str]:
    result = []
    for a, b, c in permutations(("F", "Cl", "Br")):
        for tag in range(1, 21):
            result.extend((
                f"[P@TB{tag}]({a})({b}){c}",
                f"{a}[P@TB{tag}]({b}){c}",
            ))
    for a, b, c, d in permutations(("F", "Cl", "Br", "I")):
        for tag in range(1, 31):
            result.extend((
                f"[S@OH{tag}]({a})({b})({c}){d}",
                f"{a}[S@OH{tag}]({b})({c}){d}",
            ))
    assert len(result) == 1680
    return result


def square_planar_rows() -> list[str]:
    result = []
    for a, b in permutations(("F", "Cl")):
        for tag in range(1, 4):
            result.extend((f"[Pt@SP{tag}]({a}){b}", f"{a}[Pt@SP{tag}]{b}"))
    assert len(result) == 12
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--square-planar-output", type=Path)
    args = parser.parse_args()
    args.output.write_text("\n".join(rows()) + "\n")
    if args.square_planar_output:
        args.square_planar_output.write_text("\n".join(square_planar_rows()) + "\n")


if __name__ == "__main__":
    main()
