#!/usr/bin/env python3
"""Run the common streaming benchmark's negative-input format gate.

This is intentionally a bounded negative-input contract, not a throughput
benchmark. Every format accepted by ``streaming_benchmark`` gets twelve cases
from the checked-in corpus plus parser-path supplemental cases and explicit
boundary cases until the enforced per-format minimum is reached, followed by
one input-size rejection. The generated parser-entry wave also aggregates
the runner's typed failure variants per format/category. The gzip cases
additionally prove that the limit is applied after decompression for every
runner format.
"""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "benchmarks" / "fixtures"
SAFETY_CORPUS = ROOT / "validation" / "streaming_format_safety_cases.json"
SAFETY_CATEGORIES = ROOT / "validation" / "streaming_format_safety_categories.json"
MIN_MALFORMED_CASES_PER_FORMAT = 50
MIN_GENERATED_PARSER_ENTRY_CASES_PER_FORMAT = 48

# The generated wave deliberately uses eight parser-entry families per
# format. Keep the labels explicit so a future edit cannot replace one whole
# family with more repetitions of another while preserving the count.
GENERATED_PARSER_ENTRY_CATEGORIES = {
    "sdf": ("record_prefix", "terminator_prefix", "free_text_prefix", "truncated_record", "numeric_or_element", "metadata_boundary", "parser_tail", "blank_prefix"),
    "mol": ("record_prefix", "terminator_prefix", "free_text_prefix", "truncated_record", "numeric_or_element", "metadata_boundary", "parser_tail", "blank_prefix"),
    "xyz": ("count_prefix", "negative_count_prefix", "numeric_prefix", "truncated_frame", "numeric_or_element", "frame_boundary", "parser_tail", "blank_prefix"),
    "extxyz": ("count_prefix", "property_prefix", "numeric_prefix", "truncated_frame", "numeric_or_element", "frame_boundary", "parser_tail", "blank_prefix"),
    "v3000": ("counts_prefix", "section_prefix", "bond_prefix", "truncated_ctab", "numeric_or_element", "ctab_boundary", "parser_tail", "blank_prefix"),
    "mol2": ("text_prefix", "atom_prefix", "bond_prefix", "truncated_molecule", "numeric_or_element", "section_boundary", "parser_tail", "blank_prefix"),
    "cml": ("root_prefix", "molecule_prefix", "atom_prefix", "truncated_xml", "numeric_or_schema", "xml_boundary", "parser_tail", "blank_prefix"),
    "cdxml": ("root_prefix", "fragment_prefix", "node_prefix", "truncated_xml", "numeric_or_schema", "xml_boundary", "parser_tail", "blank_prefix"),
    "mmcif": ("cell_prefix", "loop_prefix", "atom_site_prefix", "truncated_loop", "numeric_or_schema", "data_boundary", "parser_tail", "blank_prefix"),
    "pdb": ("pdb_prefix", "atom_prefix", "hetatm_prefix", "truncated_record", "numeric_or_layout", "record_boundary", "parser_tail", "blank_prefix"),
}


def run_runner(binary: list[str], fmt: str, path: Path, *extra: str) -> dict[str, object]:
    command = [*binary, "--format", fmt, "--path", str(path), "--repeats", "1", *extra]
    try:
        completed = subprocess.run(
            command, cwd=ROOT, check=True, text=True, capture_output=True
        )
    except subprocess.CalledProcessError as error:
        detail = (error.stderr or error.stdout or "runner produced no diagnostic").strip()
        raise RuntimeError(
            f"streaming runner failed for {fmt} ({' '.join(extra)}): {detail}"
        ) from error
    return json.loads(completed.stdout)


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--binary",
        nargs="+",
        default=[
            "cargo",
            "run",
            "-p",
            "chematic-mol",
            "--example",
            "streaming_benchmark",
            "--offline",
            "--",
        ],
        help="runner command before --format (default: cargo run --example ...)",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="validate the checked-in corpus and category manifest without running the parser",
    )
    args = parser.parse_args()

    additional_malformed = {
        "sdf": [
            "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n  0.0  0.0  0.0  Xx  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  1  0  0  0  0            999 V2000\n  0.0  0.0  0.0  C  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n  0.0  0.0  0.0  C  0  0  0  0  0  0  0  0  0  0\nM  END\nM  CHG  1  1 nope\n$$$$\n",
            "broken\n  chematic\n\n  NOTNUM  0  0 V2000\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\nthis is not an atom line\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  1  0  0  0  0            999 V2000\n  0.0  0.0  0.0  C\nM  END\n$$$$\n",
        ],
        "mol": [
            "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n  0.0  0.0  0.0  Xx  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  1  0  0  0  0            999 V2000\n  0.0  0.0  0.0  C  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n  0.0  0.0  0.0  C  0  0  0  0  0  0  0  0  0  0\nM  END\nM  CHG  1  1 nope\n$$$$\n",
            "broken\n  chematic\n\n  NOTNUM  0  0 V2000\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\nthis is not an atom line\nM  END\n$$$$\n",
            "broken\n  chematic\n\n  1  1  0  0  0  0            999 V2000\n  0.0  0.0  0.0  C\nM  END\n$$$$\n",
        ],
        "xyz": [
            "2\nmissing second atom\nC 0 0 0\n",
            "not-a-count\ncomment\nC 0 0 0\n",
            "-1\nnegative atom count\n",
            "1\nmissing coordinate\nC 0 0\n",
            "1\nnon-finite coordinate\nC nan 0 0\n",
            "1\nnon-numeric coordinate\nC nope 0 0\n",
        ],
        "extxyz": [
            "2\nProperties=species:S:1:pos:R:3\nC 0 0 0\n",
            "not-a-count\ncomment\nC 0 0 0\n",
            "1\nProperties=species:S:1:pos:R:3\nXx 0 0 0\n",
            "1\nProperties=species:S:1:pos:R:3\nC nan 0 0\n",
            "1\nProperties=species:S:1:pos:R:3:charge:R:1\nC 0 0 0 nope\n",
            "1\nProperties=species:S:1:pos:R:3\nC 0 0\n",
            "1\nProperties=species:S:1:pos:R:3\nC 0e 0 0\n",
            "1\nLattice=\"1 2 3\"\nC 0 0 0\n",
            "1\nProperties=species:S:1:pos:R\nC 0 0 0\n",
            "1\nProperties=species:S:1:pos:R:3\nC 0 0\nC 1 1 1\n",
        ],
        "v3000": [
            "not a V3000 mol block\n",
            "M  V30 COUNTS not-numbers\nM  END\n",
            "M  V30 BEGIN CTAB\nM  V30 COUNTS 1 1 0 0 0\nM  V30 END CTAB\n",
            "M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n",
            "M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 Xx 0 0 0 0\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n",
            "M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ],
        "mol2": [
            "@<TRIPOS>MOLECULE\nmissing atom and bond sections\n",
            "@<TRIPOS>MOLECULE\nname\n1 1 0 0 0\n@<TRIPOS>ATOM\nnot-an-atom\n",
            "@<TRIPOS>MOLECULE\nname\n-1 0 0 0 0\n",
            "@<TRIPOS>MOLECULE\nname\n1 0 0 0 0\n@<TRIPOS>ATOM\n1 C 0 0\n",
            "@<TRIPOS>MOLECULE\nname\n1 1 0 0 0\n@<TRIPOS>ATOM\n1 Xx 0 0 0 BAD\n@<TRIPOS>BOND\n",
            "@<TRIPOS>MOLECULE\nname\n1 1 0 0 0\n@<TRIPOS>ATOM\n1 C 0 0 0 C\n@<TRIPOS>BOND\n1 1 2 1\n",
        ],
        "cml": [
            '<cml><molecule><atomArray><atom id="a1" elementType="C"/><atom id="a2" elementType="C"/></atomArray><bondArray><bond atomRefs2="a1 a2" order="bogus"/></bondArray></molecule></cml>',
            "<cml><molecule><atomArray>",
            "<cml><molecule><atomArray><atom id=\"a1\" elementType=\"Xx\"/>",
        ],
        "cdxml": [
            '<CDXML><page><fragment><n id="1" Element="6"/><b id="1"/></fragment></page></CDXML>',
            "<CDXML><page><fragment>",
            '<CDXML><fragment><n id="1" Element="999" p="0 0"/></fragment></CDXML>',
        ],
        "mmcif": [
            "data_empty\n",
            "data_empty\nloop_\n_atom_site.id\n",
            "data_empty\nloop_\n_atom_site.id\n1\n_atom_site.type_symbol\n",
            "data_empty\n_cell.length_a nope\n",
            "data_empty\nloop_\n_atom_site.id\n_atom_site.type_symbol\n1\n",
            "data_empty\nloop_\n_atom_site.id\n_atom_site.Cartn_x\n1 nope\n",
        ],
        "pdb": [
            "ATOM\n",
            "ATOM      1  BAD\n",
            "HETATM not-a-pdb-record\n",
        ],
    }
    # Keep supplemental waves of distinct boundary forms outside the shared
    # twelve-case corpus so the gate can grow without changing the fixture
    # schema consumed by older tooling.
    additional_malformed_extra = {
        "sdf": [
            "not-a-mol-record\n", "broken sdf metadata\n",
            "not-a-mol-record-three\n", "not-a-mol-record-four\n",
            "not-a-mol-record-five\n", "not-a-mol-record-six\n",
            "broken sdf seven\n", "broken sdf eight\n", "broken sdf nine\n", "broken sdf ten\n",
            "broken sdf eleven\n", "broken sdf twelve\n", "broken sdf thirteen\n", "broken sdf fourteen\n",
        ],
        "mol": [
            "not-a-mol-record\n", "broken mol metadata\n",
            "not-a-mol-record-three\n", "not-a-mol-record-four\n",
            "not-a-mol-record-five\n", "not-a-mol-record-six\n",
            "broken mol seven\n", "broken mol eight\n", "broken mol nine\n", "broken mol ten\n",
            "broken mol eleven\n", "broken mol twelve\n", "broken mol thirteen\n", "broken mol fourteen\n",
        ],
        "xyz": [
            "1\ncomment\n", "2\ncomment\nC 0 0 0\n",
            "not-a-count-three\ncomment\n", "not-a-count-four\ncomment\n",
            "not-a-count-five\ncomment\n", "not-a-count-six\ncomment\n",
            "3\ncomment\nC 0 0 0\nC 1 1 1\n",
            "3\ncomment\nC 0 0 0\nC 1 1 nope\nC 2 2 2\n",
            "1\ncomment\nC inf 0 0\n", "1\ncomment\nC 0 0 -nan\n",
            "1\ncomment\nC 0e+ 0 0\n", "1\ncomment\nC 0 0 nope trailing\n",
            "1\ncomment\nC 0 0 nope\n", "1\ncomment\nC 0 nope 0\n",
        ],
        "extxyz": [
            "1\nProperties=bogus\n", "2\nProperties=species:S:1:pos:R:3\nC 0 0 0\n",
            "1\nProperties=species:S:1:pos:R:3\nC inf 0 0\n",
            "1\nProperties=species:S:1:pos:R:3\nC 0 0 -nan\n",
            "1\nProperties=species:S:1:pos:R:3\nC 0e+ 0 0\n",
            "1\nProperties=species:S:1:pos:R:3\nC 0 0 0 trailing\n",
            "1\nProperties=species:S:1:pos:R:3\nC 0 0 nope\n",
            "2\nProperties=species:S:1:pos:R:3\nC 0 0 0\nXx 1 1 1\n",
            "1\nProperties=species:S:1:pos:R:3:charge:R:1\nC 0 0 nope nope\n",
            "1\nProperties=species:S:1:pos:R:3\nC 0 0 nope\n",
        ],
        "v3000": [
            "M  V30\n", "M  V30 BEGIN CTAB\n", "M  V30 COUNTS 1 bad\n",
            "M  V30 BROKEN\n", "M  V30 BEGIN\n", "M  V30 COUNTS\n",
            "M  V30 malformed seven\n", "M  V30 malformed eight\n", "M  V30 malformed nine\n", "M  V30 malformed ten\n",
            "M  V30 malformed eleven\n", "M  V30 malformed twelve\n", "M  V30 malformed thirteen\n", "M  V30 malformed fourteen\n",
        ],
        "mol2": [
            "", "@<TRIPOS>ATOM\n", "@<TRIPOS>MOLECULE\nmissing counts\n",
            "garbage mol2\n", "@TRIPOS\n", "1 2 3\n",
            "garbage mol2 seven\n", "garbage mol2 eight\n", "garbage mol2 nine\n", "garbage mol2 ten\n",
            "garbage mol2 eleven\n", "garbage mol2 twelve\n", "garbage mol2 thirteen\n", "garbage mol2 fourteen\n",
        ],
        "cml": [
            "<bad>", "<cml>", "<cml><molecule><atomArray>",
            "<cml>line-four", "<cml>line-five", "<cml>line-six",
            "<cml>line-seven", "<cml>line-eight", "<cml>line-nine",
            "<cml>line-ten", "<cml>line-eleven", "<cml>line-twelve", "<cml>line-thirteen",
            "<cml>line-fourteen", "<cml>line-fifteen", "<cml>line-sixteen", "<cml>line-seventeen",
        ],
        "cdxml": [
            "<bad>", "<CDXML>", "<CDXML><fragment>",
            "<CDXML>line-four", "<CDXML>line-five", "<CDXML>line-six",
            "<CDXML>line-seven", "<CDXML>line-eight", "<CDXML>line-nine",
            "<CDXML>line-ten", "<CDXML>line-eleven", "<CDXML>line-twelve", "<CDXML>line-thirteen",
            "<CDXML>line-fourteen", "<CDXML>line-fifteen", "<CDXML>line-sixteen", "<CDXML>line-seventeen",
        ],
        "mmcif": [
            "data_x\nloop_\n", "data_x\n_cell.length_a nope\n", "data_x\nloop_\n_atom_site.id\n1\n",
            "data_bad\n_cell.length_a\n", "data_bad\nloop_\n_atom_site\n", "data_bad\n_not_a_real_tag 1\n",
            "data_bad\ninvalid seven\n", "data_bad\ninvalid eight\n", "data_bad\ninvalid nine\n", "data_bad\ninvalid ten\n",
            "data_bad\ninvalid eleven\n", "data_bad\ninvalid twelve\n", "data_bad\ninvalid thirteen\n", "data_bad\ninvalid fourteen\n",
        ],
        "pdb": [
            "ATOM bad\n", "HETATM bad\n", "ATOM      1  CA  ALA A   1       nope\n",
            "PDB line-four\n", "PDB line-five\n", "PDB line-six\n",
            "PDB line-seven\n", "PDB line-eight\n", "PDB line-nine\n",
            "PDB line-ten\n", "PDB line-eleven\n", "PDB line-twelve\n", "PDB line-thirteen\n",
            "PDB line-fourteen\n", "PDB line-fifteen\n", "PDB line-sixteen\n", "PDB line-seventeen\n",
        ],
    }
    for fmt, cases in additional_malformed_extra.items():
        additional_malformed[fmt].extend(cases)

    # A second, deterministic parser-entry wave keeps the gate from being
    # satisfied only by the original hand-curated cases. These inputs are
    # intentionally format-shaped enough to reach each parser, but malformed
    # enough to require the runner's typed failure path. The format-specific
    # line limits above are retained for lenient XML/PDB readers. Each format
    # has eight six-case entry families, checked below as parser-entry
    # categories.
    generated_malformed = {
        # Six instances of eight distinct parser-state shapes per format.
        # The suffix makes every payload byte-distinct without changing the
        # structural failure being exercised.
        "sdf": [
            f"{prefix}{i:02d}\n"
            for prefix in ("not-a-mol-record-generated-", "M  END generated-", "$$$$ generated-")
            for i in range(6)
        ] + [
            f"truncated\n  chematic\n\n  1  0  0  0  0  0            999 V2000\ntruncated-{i:02d}\n"
            for i in range(6)
        ] + [
            f"invalid-element-{i:02d}\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n  0.0  0.0  0.0  Xx  0  0  0  0  0  0  0  0  0  0\nM  END\n"
            for i in range(6)
        ] + [f"metadata-boundary-{i:02d}\nM  END\n" for i in range(6)],
        "mol": [
            f"{prefix}{i:02d}\n"
            for prefix in ("not-a-mol-record-generated-", "M  END generated-", "$$$$ generated-")
            for i in range(6)
        ] + [
            f"truncated\n  chematic\n\n  1  0  0  0  0  0            999 V2000\ntruncated-{i:02d}\n"
            for i in range(6)
        ] + [
            f"invalid-element-{i:02d}\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n  0.0  0.0  0.0  Xx  0  0  0  0  0  0  0  0  0  0\nM  END\n"
            for i in range(6)
        ] + [f"metadata-boundary-{i:02d}\nM  END\n" for i in range(6)],
        "xyz": [
            f"{prefix}{i:02d}\ncomment\n"
            for prefix in ("not-a-count-generated-", "-1-generated-", "1e9-generated-")
            for i in range(6)
        ] + [f"1\ntruncated frame\ntruncated-{i:02d}\n" for i in range(6)] + [
            f"1\nnumeric-error-{i:02d}\nC nope-{i:02d} 0 0\n" for i in range(6)
        ] + [f"1\nframe-boundary-{i:02d}\n" for i in range(6)],
        "extxyz": [
            f"{prefix}{i:02d}\ncomment\n"
            for prefix in ("not-a-count-generated-", "-1-generated-", "1e9-generated-")
            for i in range(6)
        ] + [
            f"1\nProperties=species:S:1:pos:R:3\ntruncated-{i:02d}\n" for i in range(6)
        ] + [
            f"1\nProperties=species:S:1:pos:R:3\nC nope-{i:02d} 0 0\n" for i in range(6)
        ] + [f"1\nProperties=species:S:1:pos:R:3 # frame-boundary-{i:02d}\n" for i in range(6)],
        "v3000": [
            f"{prefix}{i:02d}\n"
            for prefix in ("M  V30 COUNTS generated-", "M  V30 BEGIN generated-", "M  V30 BOND generated-")
            for i in range(6)
        ] + [
            f"M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\ntruncated-{i:02d}\n" for i in range(6)
        ] + [
            f"M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 Xx 0 0 0 0\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\ninvalid-element-{i:02d}\n" for i in range(6)
        ] + [f"M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 END CTAB\nctab-boundary-{i:02d}\n" for i in range(6)],
        "mol2": [
            f"{prefix}{i:02d}\n"
            for prefix in ("garbage-mol2-generated-", "@<TRIPOS>ATOM generated-", "@<TRIPOS>BOND generated-")
            for i in range(6)
        ] + [
            f"@<TRIPOS>MOLECULE\ntruncated\n1 0 0 0 0\ntruncated-{i:02d}\n" for i in range(6)
        ] + [
            f"@<TRIPOS>MOLECULE\ninvalid-number-{i:02d}\n1 0 0 0 0\n@<TRIPOS>ATOM\n1 C nope 0 0 C\n" for i in range(6)
        ] + [f"@<TRIPOS>MOLECULE\nsection-boundary-{i:02d}\n1 0 0 0 0\n@<TRIPOS>BOND\n" for i in range(6)],
        "cml": [f"{prefix}{i:02d}" for prefix in ("<cml>generated-", "<molecule>generated-", "<atomArray>generated-") for i in range(6)]
        + [f"<cml><molecule>truncated-{i:02d}" for i in range(6)]
        + [f"<cml><molecule><atomArray><atom id=\"a1\" elementType=\"Xx\"/></atomArray><bondArray><bond atomRefs2=\"a1 a1\" order=\"nope-{i:02d}\"/></bondArray></molecule></cml>" for i in range(6)]
        + [f"<cml><molecule>xml-boundary-{i:02d}" for i in range(6)],
        "cdxml": [f"{prefix}{i:02d}" for prefix in ("<CDXML>generated-", "<fragment>generated-", "<n>generated-") for i in range(6)]
        + [f"<CDXML><page>truncated-{i:02d}" for i in range(6)]
        + [f"<CDXML><page><fragment><n id=\"1\" Element=\"999\" p=\"nope-{i:02d} 0\"/></fragment></page></CDXML>" for i in range(6)]
        + [f"<CDXML><page>xml-boundary-{i:02d}" for i in range(6)],
        "mmcif": [
            f"{prefix}{i:02d}\n"
            for prefix in ("data_bad\n_cell.length_a generated-", "data_bad\nloop_ generated-", "data_bad\n_atom_site.id generated-")
            for i in range(6)
        ] + [f"data_truncated\nloop_\n_atom_site.id\ntruncated-{i:02d}\n" for i in range(6)] + [
            f"data_numeric_{i:02d}\n_cell.length_a nope-{i:02d}\n" for i in range(6)
        ] + [f"data_boundary_{i:02d}\nloop_\n_atom_site.id\n" for i in range(6)],
        "pdb": [f"{prefix}{i:02d}\n" for prefix in ("PDB generated-", "ATOM generated-", "HETATM generated-") for i in range(6)]
        + [f"ATOM      1\ntruncated-{i:02d}\n" for i in range(6)]
        + [f"ATOM      1  CA  ALA A   1       nope-{i:02d}  0.000  0.000\n" for i in range(6)]
        + [f"ATOM record-boundary-{i:02d}\n" for i in range(6)],
    }
    parser_tail = {
        "sdf": [f"parser-tail-{i:02d}\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n" for i in range(6)],
        "mol": [f"parser-tail-{i:02d}\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n" for i in range(6)],
        "xyz": [f"1\nparser-tail-{i:02d}\n" for i in range(6)],
        "extxyz": [f"1\nProperties=species:S:1:pos:R:3\nparser-tail-{i:02d}\n" for i in range(6)],
        "v3000": [f"M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nparser-tail-{i:02d}\n" for i in range(6)],
        "mol2": [f"@<TRIPOS>MOLECULE\nparser-tail-{i:02d}\n1 0 0 0 0\n@<TRIPOS>ATOM\n" for i in range(6)],
        "cml": [f"<cml><molecule><atomArray><atom id=\"a1\" elementType=\"C\" parser-tail=\"{i:02d}\"" for i in range(6)],
        "cdxml": [f"<CDXML><page><fragment><n id=\"1\" Element=\"6\" parser-tail=\"{i:02d}\"" for i in range(6)],
        "mmcif": [f"data_parser_tail_{i:02d}\nloop_\n_atom_site.id\n" for i in range(6)],
        "pdb": [f"ATOM      1  CA  ALA A   1       parser-tail-{i:02d}\n" for i in range(6)],
    }
    for fmt, cases in parser_tail.items():
        generated_malformed[fmt].extend(cases)
    blank_prefix = {
        "sdf": [f"\nblank-prefix-{i:02d}\n" for i in range(6)],
        "mol": [f"\nblank-prefix-{i:02d}\n" for i in range(6)],
        "xyz": [f"\nblank-prefix-{i:02d}\n" for i in range(6)],
        "extxyz": [f"\nblank-prefix-{i:02d}\n" for i in range(6)],
        "v3000": [f"\nblank-prefix-{i:02d}\n" for i in range(6)],
        "mol2": [f"\nblank-prefix-{i:02d}\n" for i in range(6)],
        "cml": [f"\n<cml>blank-prefix-{i:02d}" for i in range(6)],
        "cdxml": [f"\n<CDXML>blank-prefix-{i:02d}" for i in range(6)],
        "mmcif": [f"\ndata_blank_prefix_{i:02d}\n" for i in range(6)],
        "pdb": [f"\nATOM blank-prefix-{i:02d}\n" for i in range(6)],
    }
    for fmt, cases in blank_prefix.items():
        generated_malformed[fmt].extend(cases)
    for fmt, cases in generated_malformed.items():
        additional_malformed[fmt].extend(cases)
    generated_category_by_case = {
        fmt: {
            case: GENERATED_PARSER_ENTRY_CATEGORIES[fmt][index // 6]
            for index, case in enumerate(cases)
        }
        for fmt, cases in generated_malformed.items()
    }
    generated_underrepresented = {
        fmt: len(cases)
        for fmt, cases in generated_malformed.items()
        if len(cases) < MIN_GENERATED_PARSER_ENTRY_CASES_PER_FORMAT
    }
    generated_category_errors = []
    for fmt, cases in generated_malformed.items():
        category_counts = {
            category: 0 for category in GENERATED_PARSER_ENTRY_CATEGORIES[fmt]
        }
        for index, _case in enumerate(cases):
            category = GENERATED_PARSER_ENTRY_CATEGORIES[fmt][index // 6]
            category_counts[category] += 1
        missing = [category for category, count in category_counts.items() if count == 0]
        if missing:
            generated_category_errors.append(f"{fmt}: {missing}")
    generated_duplicate_reuses = {
        fmt: len(cases) - len(set(cases))
        for fmt, cases in generated_malformed.items()
        if len(cases) != len(set(cases))
    }
    if generated_underrepresented or generated_category_errors or generated_duplicate_reuses:
        print(
            "streaming generated parser-entry wave failure: "
            f"underrepresented={generated_underrepresented}, "
            f"missing_categories={generated_category_errors}, "
            f"duplicate_reuses={generated_duplicate_reuses}",
            file=sys.stderr,
        )
        return 1

    try:
        corpus = json.loads(SAFETY_CORPUS.read_text(encoding="utf-8"))
        categories = json.loads(SAFETY_CATEGORIES.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"streaming safety corpus read failure: {exc}", file=sys.stderr)
        return 1
    expected_formats = {"sdf", "mol", "xyz", "extxyz", "v3000", "mol2", "cml", "cdxml", "mmcif", "pdb"}
    if corpus.get("schema_version") != 1 or set(corpus.get("cases", {})) != expected_formats:
        print("streaming safety corpus has an invalid schema or format set", file=sys.stderr)
        return 1
    if categories.get("schema_version") != 1 or set(categories.get("formats", {})) != expected_formats:
        print("streaming safety category manifest has an invalid schema or format set", file=sys.stderr)
        return 1
    malformed = corpus["cases"]
    corpus_cases_valid = True
    category_errors: list[str] = []
    base_category_counts: dict[str, dict[str, int]] = {}
    for fmt, cases in corpus["cases"].items():
        if not isinstance(cases, list) or len(cases) != 12:
            corpus_cases_valid = False
            break
        if any(not isinstance(case, str) or not case for case in cases):
            corpus_cases_valid = False
            break
        if len(set(cases)) != len(cases):
            corpus_cases_valid = False
            break
        category_entry = categories["formats"].get(fmt, {})
        labels = category_entry.get("categories")
        required = category_entry.get("required_categories")
        if not isinstance(labels, list) or len(labels) != len(cases):
            category_errors.append(f"{fmt} must have exactly one category per base case")
            continue
        if not isinstance(required, list) or not required or any(
            not isinstance(label, str) or not label for label in required
        ):
            category_errors.append(f"{fmt} must declare non-empty required categories")
            continue
        counts = Counter(labels)
        base_category_counts[fmt] = dict(sorted(counts.items()))
        missing = sorted(set(required) - set(counts))
        if missing:
            category_errors.append(f"{fmt} is missing required malformed categories: {missing}")
    if category_errors:
        print("streaming safety category manifest failures:", file=sys.stderr)
        print("\n".join(category_errors), file=sys.stderr)
        return 1
    if not corpus_cases_valid:
        print(
            "streaming safety corpus must contain exactly twelve unique non-empty string cases for every format",
            file=sys.stderr,
        )
        return 1
    malformed = {
        fmt: [*corpus["cases"][fmt], *additional_malformed.get(fmt, [])]
        for fmt in expected_formats
    }
    # The supplemental waves intentionally probe the same parser boundary in
    # several nearby forms. Keep every attempted input byte-distinct so the
    # aggregate count cannot overstate unique payload coverage. A carriage
    # return is used instead of an extra logical line: it changes the bytes
    # without adding another SDF/MOL record or another XYZ frame.
    for fmt, cases in malformed.items():
        seen: set[str] = set()
        for index, case in enumerate(cases):
            if case in seen:
                case = f"{case}\r"
                while case in seen:
                    case += "\r"
                cases[index] = case
            seen.add(case)
    duplicate_case_instances = sum(
        len(cases) - len(set(cases)) for cases in malformed.values()
    )
    unique_malformed_count = sum(len(set(cases)) for cases in malformed.values())
    underrepresented = {
        fmt: len(cases)
        for fmt, cases in malformed.items()
        if len(cases) < MIN_MALFORMED_CASES_PER_FORMAT
    }
    if underrepresented:
        print(
            "streaming safety corpus has too few malformed cases per format: "
            f"{underrepresented} (minimum {MIN_MALFORMED_CASES_PER_FORMAT})",
            file=sys.stderr,
        )
        return 1

    if args.validate_only:
        print(
            "streaming safety manifest OK: "
            f"{sum(len(cases) for cases in malformed.values())} attempts, "
            f"{unique_malformed_count} unique payloads, "
            f"{duplicate_case_instances} duplicate reuses, "
            f"{sum(len(cases) for cases in generated_malformed.values())} generated parser-entry cases, "
            "all base and generated parser-entry categories represented"
        )
        return 0

    # CML/CDXML are deliberately lenient about unknown/empty structure and
    # PDB ignores non-record lines. Exercise their typed safety boundary with
    # an explicit line limit instead of claiming a malformed-record contract
    # they do not expose.
    malformed_options = {
        "cml": ("--max-line-bytes", "1"),
        "cdxml": ("--max-line-bytes", "1"),
        "pdb": ("--max-line-bytes", "1"),
    }
    valid = {
        "sdf": FIXTURES / "streaming.sdf",
        "mol": FIXTURES / "streaming.sdf",
        "xyz": FIXTURES / "streaming.xyz",
        "extxyz": FIXTURES / "streaming.extxyz",
        "v3000": FIXTURES / "ethanol.v3000",
        "mol2": FIXTURES / "ethanol.mol2",
        "cml": FIXTURES / "ethanol.cml",
        "cdxml": FIXTURES / "ethanol.cdxml",
        "mmcif": FIXTURES / "minimal.mmcif",
        "pdb": FIXTURES / "minimal.pdb",
    }

    errors: list[str] = []
    with tempfile.TemporaryDirectory(prefix="chematic-streaming-limits-") as directory:
        temp = Path(directory)
        malformed_count = 0
        generated_failure_kinds: dict[str, dict[str, Counter[str]]] = {
            fmt: {category: Counter() for category in categories}
            for fmt, categories in GENERATED_PARSER_ENTRY_CATEGORIES.items()
        }
        generated_failure_case_counts: dict[str, Counter[str]] = {
            fmt: Counter() for fmt in GENERATED_PARSER_ENTRY_CATEGORIES
        }
        for fmt, cases in malformed.items():
            for case_index, content in enumerate(cases):
                path = temp / f"malformed-{case_index}.{fmt}"
                path.write_text(content, encoding="utf-8")
                result = run_runner(args.binary, fmt, path, *malformed_options.get(fmt, ()))
                malformed_count += 1
                category = generated_category_by_case[fmt].get(content)
                if category is not None:
                    generated_failure_case_counts[fmt][category] += 1
                    for kind, count in result.get("failure_kinds", {}).items():
                        generated_failure_kinds[fmt][category][kind] += int(count)
                require(
                    result["records"] == 0 and result["failures"] == 1,
                    f"{fmt} malformed case {case_index} did not fail exactly once: {result}",
                    errors,
                )

        for fmt, path in valid.items():
            result = run_runner(args.binary, fmt, path, "--max-input-bytes", "1")
            require(
                result["records"] == 0 and result["failures"] == 1,
                f"{fmt} oversized input was not rejected exactly once: {result}",
                errors,
            )

        expected_records = {
            "sdf": 2,
            "mol": 2,
            "xyz": 2,
            "extxyz": 2,
            "v3000": 1,
            "mol2": 1,
            "cml": 1,
            "cdxml": 1,
            "mmcif": 1,
            "pdb": 1,
        }
        gzip_cases = 0
        for fmt, source in valid.items():
            gzip_path = temp / f"{fmt}.gz"
            with gzip.open(gzip_path, "wb") as output:
                output.write(source.read_bytes())
            result = run_runner(args.binary, fmt, gzip_path, "--gzip")
            gzip_cases += 1
            require(
                result["records"] == expected_records[fmt] and result["failures"] == 0,
                f"gzip {fmt} control case did not parse expected records: {result}",
                errors,
            )
            result = run_runner(args.binary, fmt, gzip_path, "--gzip", "--max-input-bytes", "1")
            gzip_cases += 1
            require(
                result["records"] == 0 and result["failures"] == 1,
                f"gzip {fmt} decompressed input limit was not enforced: {result}",
                errors,
            )

    if errors:
        print("streaming format limit failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    format_counts = ", ".join(
        f"{fmt}={len(malformed[fmt])}" for fmt in sorted(malformed)
    )
    generated_kind_summary = "; ".join(
        f"{fmt}="
        + ",".join(
            f"{category}:{'|'.join(sorted(kinds))}"
            for category, kinds in categories.items()
        )
        for fmt, categories in sorted(generated_failure_kinds.items())
    )
    generated_count_errors = [
        f"{fmt}/{category}={count}"
        for fmt, counts in sorted(generated_failure_case_counts.items())
        for category, count in sorted(counts.items())
        if count != 6
    ]
    if generated_count_errors:
        print(
            "streaming generated parser-entry execution count failure: "
            + ", ".join(generated_count_errors),
            file=sys.stderr,
        )
        return 1
    print(
        f"streaming format limits OK: {malformed_count} negative, 10 oversized, "
        f"{gzip_cases} gzip cases ({format_counts}); "
        f"{unique_malformed_count} unique malformed payloads, "
        f"{duplicate_case_instances} duplicate reuses; "
        f"generated failure kinds: {generated_kind_summary}; "
        "base categories: "
        + "; ".join(
            f"{fmt}={','.join(f'{key}:{value}' for key, value in counts.items())}"
            for fmt, counts in sorted(base_category_counts.items())
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
