#!/usr/bin/env python3
"""Compare canonical SMILES strings without conflating representation and identity.

The engines may intentionally choose different canonical spellings.  The report
therefore records exact-string equality separately from semantic equality, using
RDKit to re-canonicalize each emitted spelling.  Open Babel is optional and
Indigo is recorded as unavailable unless it is available in the environment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone

import chematic
from rdkit import Chem, RDLogger
from rdkit import rdBase

RDLogger.DisableLog("rdApp.*")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def obabel_canonical_many(smiles_list: list[str], executable: str) -> list[str | None]:
    """Run one Open Babel process so corpus-scale probes do not fork per row."""
    try:
        result = subprocess.run(
            [executable, "-ismi", "-osmi", "--canonical"],
            input=("\n".join(smiles_list) + "\n").encode(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=30,
        )
    except subprocess.TimeoutExpired:
        return [None] * len(smiles_list)
    if result.returncode != 0:
        return [None] * len(smiles_list)
    lines = [line.strip() for line in result.stdout.decode(errors="replace").splitlines() if line.strip()]
    if len(lines) != len(smiles_list):
        return [None] * len(smiles_list)
    return [line.split("\t", 1)[0].strip() or None for line in lines]


def rdkit_canonical(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol is not None else None


def obabel_version(executable: str) -> str | None:
    result = subprocess.run([executable, "-V"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    line = result.stdout.decode(errors="replace").splitlines()
    return line[0].strip() if line else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="scripts/descriptor_census_corpus.smi")
    parser.add_argument("--limit", type=int, default=200)
    parser.add_argument("--output", default="validation/results/canonical-cross-engine-v1.0.13.json")
    args = parser.parse_args()
    input_path = os.path.join(ROOT, args.input)
    output_path = os.path.join(ROOT, args.output)
    raw = open(input_path, encoding="utf-8").read().encode()
    smiles_list = [line.strip() for line in raw.decode().splitlines() if line.strip()][:args.limit]
    obabel = shutil.which("obabel")
    obabel_outputs = obabel_canonical_many(smiles_list, obabel) if obabel else [None] * len(smiles_list)
    rows = []
    for index, smiles in enumerate(smiles_list):
        rd_native = rdkit_canonical(smiles)
        if rd_native is None:
            rows.append({"input_index": index, "status": "rdkit_parse_failed"})
            continue
        try:
            schematic = chematic.from_smiles(smiles).smiles
        except Exception as error:
            rows.append({"input_index": index, "status": "chematic_parse_failed", "error": str(error)})
            continue
        engines = {"rdkit": rd_native, "chematic": schematic}
        if obabel:
            engines["openbabel"] = obabel_outputs[index]
        semantic = {}
        inchis = {}
        native_mol = Chem.MolFromSmiles(rd_native)
        for name, emitted in engines.items():
            semantic[name] = rdkit_canonical(emitted) if emitted else None
            if emitted and native_mol is not None:
                emitted_mol = Chem.MolFromSmiles(emitted)
                if emitted_mol is not None:
                    try:
                        inchis[name] = Chem.MolToInchi(emitted_mol) == Chem.MolToInchi(native_mol)
                    except Exception:
                        inchis[name] = None
        rows.append({
            "input_index": index,
            "status": "ok",
            "engines": engines,
            "rdkit_recanonicalized": semantic,
            "exact_string_match": {
                name: emitted == rd_native for name, emitted in engines.items() if emitted is not None
            },
            "semantic_match_to_rdkit": {
                name: canonical == rd_native for name, canonical in semantic.items() if canonical is not None
            },
            "inchi_match_to_rdkit": inchis,
        })

    ok = [row for row in rows if row["status"] == "ok"]
    summary = {
        "schema_version": 1,
        "target_version": "1.0.13",
        "scope": "bounded canonical spelling and semantic comparison",
        "corpus": {"path": args.input, "sha256": hashlib.sha256(raw).hexdigest(), "rows_requested": len(smiles_list)},
        "engines": {
            "chematic": {"status": "measured", "version": "1.0.13"},
            "rdkit": {"status": "measured", "version": rdBase.rdkitVersion},
            "openbabel": {"status": "measured" if obabel else "not_installed", "executable": obabel,
                          "version": obabel_version(obabel) if obabel else None},
            "indigo": {"status": "not_installed"},
        },
        "counts": {
            "rows_ok": len(ok),
            "parse_or_engine_failures": len(rows) - len(ok),
            "schematic_exact_string_matches": sum(row["exact_string_match"].get("chematic", False) for row in ok),
            "schematic_semantic_matches": sum(row["semantic_match_to_rdkit"].get("chematic", False) for row in ok),
            "openbabel_exact_string_matches": sum(row["exact_string_match"].get("openbabel", False) for row in ok),
            "openbabel_semantic_matches": sum(row["semantic_match_to_rdkit"].get("openbabel", False) for row in ok),
            "openbabel_inchi_matches": sum(row["inchi_match_to_rdkit"].get("openbabel", False) for row in ok),
        },
        "rows": rows,
        "not_claimed": ["Indigo parity", "canonical string parity across algorithms", "full corpus coverage"],
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    print(json.dumps(summary["counts"], indent=2))
    print(f"wrote {os.path.relpath(output_path, ROOT)}")


if __name__ == "__main__":
    main()
