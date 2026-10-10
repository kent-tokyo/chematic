#!/usr/bin/env python3
"""Common-schema and persistent benchmark adapter for chematic, RDKit, Indigo.

The adapter deliberately exposes only operations with an equivalent public API
in all three engines.  Values are normalized before comparison; the benchmark
server returns an output digest so a timing cannot silently omit work.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

SMARTS = ["[#6]", "[#7]", "[#8]", "c1ccccc1", "[C,c](=O)[O,N]", "[N;H1,H2,H3;+0]"]
OPERATIONS = [
    "canonical_stable",
    "formula",
    "molecular_weight_centi_da",
    "tpsa_milli",
    "hba",
    "hbd",
    "rings",
    "mol_roundtrip",
    *[f"smarts:{query}" for query in SMARTS],
]


def _commit() -> str | None:
    try:
        root = Path(__file__).resolve().parents[2]
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _formula(value: str, charge: int = 0) -> str:
    compact = "".join(value.split())
    if compact.startswith("[") and "]" in compact:
        close = compact.index("]")
        compact = compact[1:close] + compact[close + 1:]
    if charge and not compact.endswith(("+", "-")):
        magnitude = "" if abs(charge) == 1 else str(abs(charge))
        compact += magnitude + ("+" if charge > 0 else "-")
    return compact


def _digest(values: list[Any]) -> str:
    payload = json.dumps(values, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode()).hexdigest()


@dataclass
class Engine:
    name: str
    version: str
    parse: Callable[[str], Any]
    operations: dict[str, Callable[[Any], Any]]
    source_commit: str | None = None


def rdkit_engine() -> Engine:
    from rdkit import Chem, rdBase
    from rdkit.Chem import Descriptors
    from rdkit.Chem import rdMolDescriptors as D

    queries = {query: Chem.MolFromSmarts(query) for query in SMARTS}

    def canonical(mol):
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

    def mol_roundtrip(mol):
        reread = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))
        return reread is not None and canonical(reread) == canonical(mol)

    operations: dict[str, Callable[[Any], Any]] = {
        "canonical_stable": lambda mol: canonical(Chem.MolFromSmiles(canonical(mol))) == canonical(mol),
        "formula": lambda mol: D.CalcMolFormula(mol),
        "molecular_weight_centi_da": lambda mol: round(Descriptors.MolWt(mol) * 100),
        "tpsa_milli": lambda mol: round(D.CalcTPSA(mol) * 1000),
        "hba": D.CalcNumHBA,
        "hbd": D.CalcNumHBD,
        "rings": D.CalcNumRings,
        "mol_roundtrip": mol_roundtrip,
    }
    operations.update({
        f"smarts:{query}": (lambda compiled: lambda mol: mol.HasSubstructMatch(compiled))(compiled)
        for query, compiled in queries.items()
    })
    return Engine("rdkit", rdBase.rdkitVersion, Chem.MolFromSmiles, operations)


def chematic_engine() -> Engine:
    import chematic

    queries = {query: chematic.compile_smarts(query) for query in SMARTS}

    def canonical(mol):
        return mol.rdkit_smiles

    def mol_roundtrip(mol):
        reread = chematic.from_mol_block(mol.rdkit_mol_block_2d())
        return canonical(reread) == canonical(mol)

    operations: dict[str, Callable[[Any], Any]] = {
        "canonical_stable": lambda mol: canonical(chematic.from_smiles(canonical(mol))) == canonical(mol),
        "formula": lambda mol: mol.formula,
        "molecular_weight_centi_da": lambda mol: round(mol.rdkit_mw * 100),
        "tpsa_milli": lambda mol: round(mol.rdkit_tpsa * 1000),
        "hba": lambda mol: mol.rdkit_hba,
        "hbd": lambda mol: mol.hbd,
        "rings": lambda mol: mol.num_rings,
        "mol_roundtrip": mol_roundtrip,
    }
    operations.update({
        f"smarts:{query}": (lambda compiled: lambda mol: compiled.matches(mol))(compiled)
        for query, compiled in queries.items()
    })
    return Engine(
        "chematic",
        getattr(chematic, "__version__", "unknown"),
        chematic.from_smiles,
        operations,
        _commit(),
    )


def indigo_engine() -> Engine:
    from indigo import Indigo

    indigo = Indigo()
    queries = {query: indigo.loadSmarts(query) for query in SMARTS}

    def canonical(mol):
        return mol.canonicalSmiles()

    def canonical_stable(mol):
        text = canonical(mol)
        return canonical(indigo.loadMolecule(text)) == text

    def mol_roundtrip(mol):
        reread = indigo.loadMolecule(mol.molfile())
        return canonical(reread) == canonical(mol)

    operations: dict[str, Callable[[Any], Any]] = {
        "canonical_stable": canonical_stable,
        "formula": lambda mol: _formula(
            mol.grossFormula(), sum(atom.charge() for atom in mol.iterateAtoms())
        ),
        "molecular_weight_centi_da": lambda mol: round(mol.molecularWeight() * 100),
        "tpsa_milli": lambda mol: round(mol.tpsa() * 1000),
        "hba": lambda mol: mol.numHydrogenBondAcceptors(),
        "hbd": lambda mol: mol.numHydrogenBondDonors(),
        "rings": lambda mol: mol.countSSSR(),
        "mol_roundtrip": mol_roundtrip,
    }

    def matches(query):
        def apply(mol):
            matcher = indigo.substructureMatcher(mol)
            return matcher.match(query) is not None
        return apply

    operations.update({f"smarts:{query}": matches(compiled) for query, compiled in queries.items()})
    return Engine("indigo", indigo.version(), indigo.loadMolecule, operations)


def load_engine(name: str) -> Engine:
    return {"rdkit": rdkit_engine, "chematic": chematic_engine, "indigo": indigo_engine}[name]()


def operation_result(engine: Engine, molecule: Any, operation: str) -> dict[str, Any]:
    try:
        return {"status": "ok", "value": engine.operations[operation](molecule)}
    except Exception as exc:  # competitor exception classes are not shared
        return {"status": "error", "error": f"{type(exc).__name__}: {exc}"}


def fresh_operation_result(engine: Engine, smiles: str, operation: str) -> dict[str, Any]:
    """Evaluate one operation without inheriting another operation's mutations."""
    try:
        molecule = engine.parse(smiles)
        if molecule is None:
            raise ValueError("parser returned no molecule")
    except Exception as exc:
        return {"status": "error", "error": f"{type(exc).__name__}: {exc}"}
    return operation_result(engine, molecule, operation)


def emit_results(engine: Engine, corpus: Path) -> None:
    corpus_hash = hashlib.sha256(corpus.read_bytes()).hexdigest()
    for line in corpus.read_text().splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        record = {
            "schema_version": 1,
            "engine": engine.name,
            "engine_version": engine.version,
            "source_commit": engine.source_commit,
            "corpus_sha256": corpus_hash,
            "id": row["id"],
            "smiles": row["smiles"],
            "status": "ok",
            "operations": {},
        }
        try:
            molecule = engine.parse(row["smiles"])
            if molecule is None:
                raise ValueError("parser returned no molecule")
        except Exception as exc:
            record["status"] = "parse_error"
            record["operations"]["parse"] = {
                "status": "parse_error",
                "error": f"{type(exc).__name__}: {exc}",
            }
        else:
            record["operations"] = {
                operation: fresh_operation_result(engine, row["smiles"], operation)
                for operation in OPERATIONS
            }
        print(json.dumps(record, sort_keys=True, separators=(",", ":")))


def benchmark(engine: Engine, request: dict[str, Any]) -> dict[str, Any]:
    smiles = request["smiles"]
    operation = request["operation"]
    lane = request["lane"]
    iterations = int(request.get("iterations", 1))
    if operation != "parse" and operation not in engine.operations:
        raise ValueError(f"unsupported operation: {operation}")
    prepared = None
    if lane == "prepared":
        prepared = [engine.parse(text) for text in smiles]
    elif lane != "pipeline":
        raise ValueError(f"unsupported lane: {lane}")

    values: list[Any] = []
    errors = 0
    start = time.perf_counter_ns()
    for _ in range(iterations):
        for index, text in enumerate(smiles):
            try:
                molecule = engine.parse(text) if prepared is None else prepared[index]
                value = True if operation == "parse" else engine.operations[operation](molecule)
                values.append(value)
            except Exception as exc:
                errors += 1
                values.append(["error", type(exc).__name__])
    elapsed = time.perf_counter_ns() - start
    return {
        "engine": engine.name,
        "engine_version": engine.version,
        "operation": operation,
        "lane": lane,
        "rows": len(smiles) * iterations,
        "elapsed_ns": elapsed,
        "errors": errors,
        "digest": _digest(values),
    }


def serve(engine: Engine) -> None:
    for line in sys.stdin:
        try:
            request = json.loads(line)
            command = request.get("command")
            if command == "metadata":
                response = {
                    "engine": engine.name,
                    "engine_version": engine.version,
                    "operations": ["parse", *OPERATIONS],
                }
            elif command == "benchmark":
                response = benchmark(engine, request)
            else:
                raise ValueError(f"unknown command: {command}")
        except Exception as exc:
            response = {"error": f"{type(exc).__name__}: {exc}"}
        print(json.dumps(response, sort_keys=True, separators=(",", ":")), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", choices=["rdkit", "chematic", "indigo"], required=True)
    parser.add_argument("--corpus", type=Path)
    parser.add_argument("--server", action="store_true")
    args = parser.parse_args()
    if args.server == (args.corpus is not None):
        parser.error("choose exactly one of --server or --corpus")
    engine = load_engine(args.engine)
    serve(engine) if args.server else emit_results(engine, args.corpus)


if __name__ == "__main__":
    main()
