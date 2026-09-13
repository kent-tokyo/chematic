#!/usr/bin/env python3
"""Build the fixed 5-format parser-security corpus with explicit provenance.

Most rows are deliberately small boundary mutations.  They are labelled as
derived rather than misrepresented as byte-identical upstream crash files; the
public issue/reproducer URL and the mutation family remain attached to every
row.  The corpus includes a valid control for each format so rejection-only
implementations cannot pass.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "validation/parser_security_corpus_v1.json"
OPENBABEL_SMILES = "https://raw.githubusercontent.com/openbabel/openbabel/master/test/files/fuzz_regress/cve-2025-10996.smi"
OPENBABEL_MOL = "https://raw.githubusercontent.com/openbabel/openbabel/master/test/files/fuzz_regress/graphsym-nonconvergent.mol"
OPENBABEL_SDF = "https://raw.githubusercontent.com/openbabel/openbabel/master/test/files/fuzz_regress/inchi-maxval-overflow.sdf"
RDKIT_SMARTS = "https://github.com/rdkit/rdkit/issues/8471"
INDIGO_V3000 = "https://github.com/epam/Indigo/issues/3286"


def sha256(payload: str) -> str:
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def record(case_id: str, fmt: str, payload: str, expected: str, source_url: str, cause: str, mutation: str) -> dict[str, str]:
    return {
        "id": case_id,
        "format": fmt,
        "payload": payload,
        "expected_status": expected,
        "source_url": source_url,
        "license": "CheMatic-authored derived boundary fixture; upstream URL is attribution, not redistributed input bytes",
        "affected_version": "upstream issue/reproducer version; chematic current candidate",
        "cause": cause,
        "origin": "derived_from_public_reproducer",
        "mutation": mutation,
        "sha256": sha256(payload),
    }


def v2000(atom: str = "C", bond: str = "", counts: str = "  1  0  0  0  0  0            999 V2000", end: str = "M  END") -> str:
    return f"case\n  chematic\n\n{counts}\n{atom}\n{bond}{end}\n"


def v3000(atom: str = "M  V30 1 C 0 0 0 0", tail: str = "") -> str:
    return "case\n  chematic\n\n  1  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\n" + atom + "\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n" + tail


def main() -> int:
    cases: list[dict[str, str]] = []
    smiles_bad = ["C&", "C@@@", "C\x00", "C\x01", "C\x07", "C\x7f", "C this is junk", "C(", "C1CC", "[C", "C)", "C??", "C%", "C/", "C\\", "C=", "C#", "C[", "C$"]
    cases.append(record("smiles-valid", "smiles", "CCO", "accepted", OPENBABEL_SMILES, "Open Babel CVE-2025-10996 SMILES parser overflow", "valid control"))
    cases.extend(record(f"smiles-{i:02}", "smiles", value, "rejected", OPENBABEL_SMILES, "Open Babel CVE-2025-10996 SMILES parser overflow", "trailing/control/truncation mutation") for i, value in enumerate(smiles_bad, 1))

    smarts_bad = ["CC this is junk", "[#6] this is junk", "[", "[#", "C(", "C)", "C&", "C\x00", "C\x01", "C$", "C%", "C/", "C\\", "[$", "[!", "[R", "[r", "[#6", "C..C"]
    cases.append(record("smarts-valid", "smarts", "[#6]", "accepted", RDKIT_SMARTS, "RDKit #8471 accepted valid prefix followed by junk", "valid control"))
    cases.extend(record(f"smarts-{i:02}", "smarts", value, "rejected", RDKIT_SMARTS, "RDKit #8471 accepted valid prefix followed by junk", "prefix/junk/truncation mutation") for i, value in enumerate(smarts_bad, 1))

    mol_bad = [
        v2000(counts="  X  0  0  0  0  0            999 V2000"), v2000(atom="not an atom"), v2000(end=""),
        v2000(counts=" -1  0  0  0  0  0            999 V2000"), v2000(atom="    0.0 nope 0.0 C"),
        v2000(atom="    0.0 0.0 0.0 Xx"), v2000(counts="  2  0  0  0  0  0            999 V2000"),
        v2000(bond="  1  2  1\n"), "case\n  chematic\n", "not a mol block\n", v2000(end="M  BAD"),
        v2000(atom=""), v2000(counts="  0  1  0  0  0  0            999 V2000"),
        v2000(bond="  x  y  z\n"), v2000(counts="  1  0"), v2000(atom="\x00"),
        v2000(atom="    nan 0.0 0.0 C"), v2000(atom="    inf 0.0 0.0 C"), v2000(end="M  END\n$$$$\ntrailing"),
    ]
    cases.append(record("mol-v2000-valid", "mol_v2000", v2000(atom="    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"), "accepted", OPENBABEL_MOL, "Open Babel graph symmetry fuzz regression", "valid control"))
    cases.extend(record(f"mol-v2000-{i:02}", "mol_v2000", value, "rejected", OPENBABEL_MOL, "Open Babel graph symmetry fuzz regression", "V2000 structural boundary mutation") for i, value in enumerate(mol_bad, 1))

    v3000_bad = ["not V3000\n", "M  V30 COUNTS nope\nM  END\n", v3000(atom="M  V30 1 Xx 0 0 0 0"), v3000(atom="M  V30 1 C nope 0 0 0"), v3000(tail="M V30 4 2 3 2 RXCTR=10\n"), "M  V30 BEGIN CTAB\n", "M  V30 END CTAB\nM END\n", v3000(atom=""), v3000(atom="M V30 1 C 0 0"), v3000(tail="M V30 BEGIN BOND\nM V30 1 1 1 2\n"), v3000(tail="M V30 BEGIN SGROUP\nM V30 1 DAT 0 ATOMS=(1 2)\n"), v3000().replace("END ATOM", "END ATOM\nM V30 END ATOM"), v3000().replace("COUNTS 1 0", "COUNTS 2 1"), "M V30 BEGIN CTAB\nM V30 COUNTS 0 0 0 0 0\n", v3000(atom="M V30 1 C nan 0 0 0"), v3000(atom="M V30 1 C inf 0 0 0"), v3000(tail="\x00"), v3000().replace("M  END", "M BAD"), v3000().replace("M  V30", "M V31", 1)]
    cases.append(record("mol-v3000-valid", "mol_v3000", v3000(), "accepted", INDIGO_V3000, "Indigo #3286 malformed V3000 reaction-center output", "valid control"))
    cases.extend(record(f"mol-v3000-{i:02}", "mol_v3000", value, "rejected", INDIGO_V3000, "Indigo #3286 malformed V3000 reaction-center output", "V3000 section/attribute boundary mutation") for i, value in enumerate(v3000_bad, 1))

    sdf_bad = [value + "$$$$\n" for value in mol_bad]
    cases.append(record("sdf-valid", "sdf", v2000(atom="    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0") + "$$$$\n", "accepted", OPENBABEL_SDF, "Open Babel inchi max-valence SDF fuzz regression", "valid control"))
    cases.extend(record(f"sdf-{i:02}", "sdf", value, "rejected", OPENBABEL_SDF, "Open Babel inchi max-valence SDF fuzz regression", "SDF record boundary mutation") for i, value in enumerate(sdf_bad, 1))
    assert len(cases) == 100
    OUTPUT.write_text(json.dumps({"schema_version": 1, "purpose": "fixed sourced parser-security corpus", "cases": cases}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT.relative_to(ROOT)), "cases": len(cases)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
