"""Reproduce closed pybel multi-molecule CML framing on the installed binding."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from xml.etree import ElementTree as ET
from openbabel import pybel, openbabel


def main():
    with TemporaryDirectory(prefix="chematic-cml-pybel-") as directory:
        path = Path(directory) / "output.cml"
        writer = pybel.Outputfile("cml", str(path), overwrite=True)
        for smi in ["CCO", "[NH4+]", "[13CH4]"]:
            writer.write(pybel.readstring("smi", smi))
        writer.close()
        source = path.read_text()
        error = None
        try:
            ET.fromstring(source)
        except ET.ParseError as exc:
            error = str(exc)
        out = Path("validation/results/cml-collection-2026-10-17")
        out.mkdir(parents=True, exist_ok=True)
        (out / "pybel-closed-output.cml").write_text(source)
        print(
            json.dumps(
                {
                    "openbabel_version": openbabel.OBReleaseVersion(),
                    "method": "pybel.Outputfile(cml).write three molecules, followed by close()",
                    "xml_valid": error is None,
                    "xml_error": error,
                    "source_file": "pybel-closed-output.cml",
                    "scope": "original three-molecule synthetic control; no force-field or coordinate generation",
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
