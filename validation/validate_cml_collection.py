"""CML framing/opaque payload checks with lxml, Indigo and Open Babel."""

import hashlib
import json
from pathlib import Path
import subprocess

from indigo import Indigo
from lxml import etree
import chematic

NS = "http://www.xml-cml.org/schema"


def tree(node, is_root=True):
    # Expanded names and literal attribute/text values, independently parsed.
    attrs = {k: v for k, v in node.attrib.items() if not (is_root and k == "id")}
    return (
        node.tag,
        tuple(sorted(attrs.items())),
        node.text,
        tuple((tree(c, False), c.tail) for c in node),
    )


def main():
    indigo = Indigo()
    smiles = [
        "CCO",
        "[NH4+]",
        "[13CH4]",
        "[Na+].[Cl-]",
        "C1CCCCC1",
        "c1ccncc1",
        "CC(=O)[O-]",
        "[2H]O[2H]",
        "C[C@H](O)F",
        "F/C=C/F",
        "F/C=C\\F",
        "N#CC(=O)O",
    ]
    sources = []
    for lane in ["native", "indigo-stereo"]:
        for smi in smiles:
            if lane == "native":
                block = chematic.from_smiles(smi).to_cml()
            else:
                # Indigo emits an unnamespaced cml envelope. Author a standalone
                # namespaced source block from its generated molecular payload.
                original = indigo.loadMolecule(smi).cml()
                payload = original.split("<molecule>", 1)[1].split("</molecule>", 1)[0]
                block = f'<molecule xmlns="{NS}">' + payload + "</molecule>"
            sources.append(block)
    cases = records_checked = 0
    checksum = hashlib.sha256()
    max_atoms = 0
    for n in [1, 2, 7, 24]:
        for sample in range(8):
            blocks = [sources[(sample * 3 + i) % len(sources)] for i in range(n)]
            checksum.update("".join(blocks).encode())
            output = chematic.write_cml_collection(blocks)
            root = etree.fromstring(output.encode())
            assert root.tag == f"{{{NS}}}cml" and len(root) == n
            records = chematic.read_cml_collection(output)
            models = list(indigo.iterateCML(indigo.loadString(output)))
            assert len(records) == len(models) == n
            assert len({r.molecule_id for r in records}) == n
            for i, (block, record, model) in enumerate(zip(blocks, records, models)):
                expected = indigo.loadMolecule(block)
                assert model.canonicalSmiles() == expected.canonicalSmiles()
                assert model.countAtoms() == expected.countAtoms()
                assert model.countBonds() == expected.countBonds()
                assert tree(etree.fromstring(record.molecule_xml.encode())) == tree(
                    etree.fromstring(block.encode())
                )
                assert record.source_text == output and record.index == i
                max_atoms = max(max_atoms, model.countAtoms())
                records_checked += 1
            # Reader output remains safely re-assemblable, with supplied IDs.
            again = chematic.write_cml_collection([r.molecule_xml for r in records])
            assert [
                m.canonicalSmiles() for m in indigo.iterateCML(indigo.loadString(again))
            ] == [m.canonicalSmiles() for m in models]
            cases += 1
    observations = []
    for n in [1, 2, 3, 7, 12]:
        blocks = [sources[i % 12] for i in range(n)]
        output = chematic.write_cml_collection(blocks)
        r = subprocess.run(
            ["/opt/homebrew/bin/obabel", "-icml", "-ocan"],
            input=output,
            text=True,
            capture_output=True,
        )
        rows = [line.split("\t")[0] for line in r.stdout.splitlines() if line.strip()]
        observations.append(
            {
                "molecules": n,
                "returned_rows": len(rows),
                "exit_code": r.returncode,
                "xml_diagnostic": "parser error" in r.stderr.lower(),
                "stderr": r.stderr,
            }
        )
    generated = subprocess.run(
        ["/opt/homebrew/bin/obabel", "-ismi", "-ocml"],
        input="CCO\n[NH4+]\n[13CH4]\n",
        text=True,
        capture_output=True,
    )
    control = subprocess.run(
        ["/opt/homebrew/bin/obabel", "-icml", "-ocan"],
        input=generated.stdout,
        text=True,
        capture_output=True,
    )
    Path(
        "validation/results/cml-collection-2026-10-17/obabel-cli-output.cml"
    ).write_text(generated.stdout)
    own_records = chematic.read_cml_collection(generated.stdout)
    assert len(own_records) == 3
    cli_control = {
        "writer_command": "obabel -ismi -ocml",
        "reader_command": "obabel -icml -ocan",
        "writer_exit_code": generated.returncode,
        "reader_exit_code": control.returncode,
        "input_molecules": 3,
        "candidate_returned_records": len(own_records),
        "xml_molecule_elements": len(etree.fromstring(generated.stdout.encode())),
        "returned_rows": len(
            [line for line in control.stdout.splitlines() if line.strip()]
        ),
        "stderr": control.stderr,
    }
    print(
        json.dumps(
            {
                "indigo_version": indigo.version(),
                "lxml_version": etree.LXML_VERSION,
                "collection_cases": cases,
                "molecular_records_checked": records_checked,
                "input_sha256": checksum.hexdigest(),
                "xml_payload": "expanded element/attribute/text trees equal; all internal source spelling retained separately",
                "graph_reference": "Indigo input-CML vs collection canonical SMILES/atom/bond counts, including the Indigo-stereo source lane; native writer stereo omissions are not repaired",
                "reassembly_cases": cases,
                "openbabel_version": subprocess.check_output(
                    ["/opt/homebrew/bin/obabel", "-V"], text=True
                ).strip(),
                "openbabel_observations": observations,
                "openbabel_cli_self_read_control": cli_control,
                "scope": "CML framing and opaque molecular subtree retention, not complete schema/chemistry validation or full native CML projection parity",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
