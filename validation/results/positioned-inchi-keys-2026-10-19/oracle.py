"""Check every successful key directly against RDKit's input-string key API."""

from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess

import chematic as ch
import chematic.chematic as native
import rdkit
from rdkit import Chem
from rdkit.Chem import inchi


def main():
    paths = [
        Path("scripts/descriptor_census_corpus.smi"),
        Path("scripts/chembl_accuracy_corpus_4999.smi"),
    ]
    smiles = [
        row.split()[0] for path in paths for row in path.read_text().splitlines()[:32]
    ]
    smiles += [
        "O",
        "CCO",
        "C",
        "CC(=O)O",
        "[Na+].[Cl-]",
        "[13CH3]CO",
        "F/C=C/F",
        "N[C@@H](C)O",
    ]
    identifiers = []
    for value in smiles:
        molecule = Chem.MolFromSmiles(value)
        assert molecule is not None
        identifier = inchi.MolToInchi(molecule)
        assert identifier.startswith("InChI=1S/")
        identifiers.append(identifier)
    inputs = []
    for identifier in identifiers:
        inputs.extend([identifier, "not-an-inchi", identifier])
    inputs += ["InChI=1/CH4/h1H4", "InChI=1S/CH4", "InChI=1S/999999999C"]
    records = ch.inchi_keys_checked(iter(inputs))
    assert len(records) == len(inputs)
    statuses = Counter()
    checked = 0
    declines = []
    for index, (source, record) in enumerate(zip(inputs, records)):
        assert record.index == index and record.source_text == source
        statuses[record.status] += 1
        if record.status == "ok":
            # Reference keys hash the ORIGINAL input identifier, without graph readback.
            expected = inchi.InchiToInchiKey(source)
            assert expected and record.key == expected, (
                index,
                source,
                record.key,
                expected,
            )
            checked += 1
        else:
            assert record.key is None and record.code and record.message
            if source in identifiers:
                declines.append(
                    {
                        "index": index,
                        "inchi": source,
                        "status": record.status,
                        "code": record.code,
                    }
                )
    for index in range(0, 3 * len(identifiers), 3):
        assert records[index].key == records[index + 2].key
        assert records[index + 1].status == "invalid"
    water, ethanol = "InChI=1S/H2O/h1H2", "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    diagnostic_inputs = [water, "InChI=1S/not-a-formula", ethanol]
    probe = subprocess.run(
        ["obabel", "-iinchi", "-oinchikey", "-e"],
        input="\n".join(diagnostic_inputs) + "\n",
        capture_output=True,
        text=True,
        check=False,
        timeout=15,
    )
    summary = {
        "rdkit": rdkit.__version__,
        "valid_identifier_rows": 2 * len(identifiers),
        "input_rows": len(inputs),
        "output_rows": len(records),
        "statuses": dict(statuses),
        "successful_keys_compared": checked,
        "key_mismatches": 0,
        "valid_declines": declines,
        "corpus_sha256": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths
        },
        "native_sha256": hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        "openbabel_diagnostic": {
            "version": subprocess.run(
                ["obabel", "-V"], capture_output=True, text=True
            ).stdout.strip(),
            "input_rows": 3,
            "output_lines": probe.stdout.splitlines(),
            "returncode": probe.returncode,
            "stderr": probe.stderr,
        },
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
