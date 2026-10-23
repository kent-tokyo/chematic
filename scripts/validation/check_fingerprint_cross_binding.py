#!/usr/bin/env python3
"""Compare one fingerprint operation across Rust, Python, and Node/WASM.

Historical reports remain frozen. Choose an explicit output path for each new
run so a current source build cannot silently overwrite older evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CORPUS = ROOT / "scripts/fixtures/chembl_accuracy_corpus_4999.smi"


@dataclass(frozen=True)
class Profile:
    rust_example: str
    node_mode: str
    python_method: str
    field: str
    contract: str
    boundary: tuple[str, ...]
    detail_index: int | None = None
    exact_rows: bool = False


PROFILES = {
    "ecfp4": Profile(
        "ecfp4_binding_dump", "ecfp4", "ecfp4", "ecfp4_hex",
        "same chematic source SMILES, native ECFP4 2048-bit LSB-first bytes",
        ("This is native chematic binding parity, not an RDKit fingerprint-accuracy claim.",
         "RDKit-compatible ECFP4 and sparse/explanation parity remain separate gates.",
         "The existing held-out parity ledger remains unchanged."), exact_rows=True,
    ),
    "rdkit-ecfp4": Profile(
        "rdkit_ecfp4_binding_dump", "rdkit-ecfp4", "rdkit_ecfp4", "ecfp4_hex",
        "same chematic source SMILES, RDKit-compatible ECFP4 2048-bit LSB-first bytes; failures compared by typed status",
        ("This is cross-binding parity for chematic's RDKit-compatible implementation, not an independent RDKit oracle run.",
         "Sparse/count/bitInfo parity remains a separate contract."),
    ),
    "rdkit-ecfp4-bitinfo": Profile(
        "rdkit_ecfp4_bitinfo_binding_dump", "rdkit-ecfp4-bitinfo", "rdkit_ecfp4_detail",
        "folded_bit_info", "same chematic source SMILES, RDKit-compatible ECFP4 folded bit -> atom/radius provenance",
        ("This checks folded bitInfo binding parity for chematic's RDKit-compatible implementation, not an independent RDKit oracle run.",
         "Raw identifier bitInfo remains a separate provenance gate."), detail_index=3,
    ),
    "rdkit-ecfp4-raw-bitinfo": Profile(
        "rdkit_ecfp4_raw_bitinfo_binding_dump", "rdkit-ecfp4-raw-bitinfo", "rdkit_ecfp4_detail",
        "raw_bit_info", "same chematic source SMILES, RDKit-compatible ECFP4 raw identifier -> atom/radius provenance",
        ("This checks raw bitInfo binding parity for chematic's RDKit-compatible implementation, not an independent RDKit oracle run.",), detail_index=2,
    ),
    "rdkit-ecfp4-sparse": Profile(
        "rdkit_ecfp4_sparse_binding_dump", "rdkit-ecfp4-sparse", "rdkit_ecfp4_detail",
        "sparse_counts", "same chematic source SMILES, RDKit-compatible ECFP4 raw sparse identifier/count pairs",
        ("This checks sparse-count binding parity for chematic's RDKit-compatible implementation, not an independent RDKit oracle run.",
         "Raw/folded bitInfo provenance remains a separate gate."), detail_index=1,
    ),
    "rdkit-path": Profile(
        "rdkit_path_binding_dump", "rdkit-path", "path_fp", "path_hex",
        "same chematic source SMILES, RDKit-compatible path 2048-bit LSB-first bytes",
        ("This is cross-binding parity for chematic's RDKit-compatible path implementation, not an independent RDKit oracle run.",
         "Native topo_path remains a separate operation."),
    ),
    "rdkit-torsion": Profile(
        "rdkit_torsion_binding_dump", "rdkit-torsion", "rdkit_torsion_fp", "torsion_hex",
        "same chematic source SMILES, RDKit-compatible hashed torsion 2048-bit LSB-first bytes",
        ("This is cross-binding parity for chematic's RDKit-compatible implementation, not an independent RDKit oracle run.",),
    ),
    "topo-path": Profile(
        "topo_path_binding_dump", "topo-path", "topo_path_fp", "topo_path_hex",
        "same chematic source SMILES, native topo_path 2048-bit LSB-first bytes",
        ("This is native chematic binding parity, not an RDKit topological fingerprint accuracy claim.",
         "RDKit-compatible path parity is a separate operation."), exact_rows=True,
    ),
    "torsion": Profile(
        "torsion_binding_dump", "torsion", "torsion_fp", "torsion_hex",
        "same chematic source SMILES, native torsion 2048-bit LSB-first bytes",
        ("This is native chematic binding parity, not an RDKit torsion accuracy claim.",
         "RDKit-compatible torsion remains a separate operation."), exact_rows=True,
    ),
}


def workspace_version() -> str:
    import re

    cargo = (ROOT / "Cargo.toml").read_text(encoding="utf-8")
    match = re.search(r'(?ms)^\[workspace\.package\]\s*$.*?^version\s*=\s*"([^"]+)"\s*$', cargo)
    if not match:
        raise ValueError("workspace version is missing from Cargo.toml")
    return match.group(1)


def runtime_provenance(version: str) -> dict:
    """Refuse stale bindings before labeling a report with the workspace version."""
    chematic = importlib.import_module("chematic")
    python_version = getattr(chematic, "__version__", None)
    if python_version != version:
        raise ValueError(f"Python chematic version {python_version!r} does not match workspace {version}")
    package = ROOT / "crates/chematic-wasm/pkg-node/package.json"
    node_version = json.loads(package.read_text(encoding="utf-8")).get("version")
    if node_version != version:
        raise ValueError(f"Node/WASM package version {node_version!r} does not match workspace {version}")
    javascript = package.parent / "chematic_wasm.js"
    wasm = package.parent / "chematic_wasm_bg.wasm"
    module = Path(chematic.__file__).resolve()
    extension = next(module.parent.glob("chematic*.so"), None)
    return {
        "python_version": python_version,
        "python_module": str(module),
        "python_extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest() if extension else None,
        "node_package_version": node_version,
        "node_package_sha256": hashlib.sha256(package.read_bytes()).hexdigest(),
        "node_javascript_sha256": hashlib.sha256(javascript.read_bytes()).hexdigest(),
        "node_wasm_sha256": hashlib.sha256(wasm.read_bytes()).hexdigest(),
    }


def payload(molecule: object, profile: Profile) -> object:
    value = getattr(molecule, profile.python_method)()
    if profile.detail_index is None:
        return value.hex()
    entries = value[profile.detail_index]
    if profile.field == "sparse_counts":
        return [[int(identifier), int(count)] for identifier, count in sorted(entries.items())]
    return [
        [int(identifier), [[int(atom), int(radius)] for atom, radius in sorted(locations)]]
        for identifier, locations in sorted(entries.items())
    ]


def python_records(corpus: str, profile: Profile) -> list[dict]:
    import chematic

    rows = []
    for index, raw in enumerate(corpus.splitlines()):
        smiles = raw.strip()
        if not smiles:
            continue
        try:
            molecule = chematic.from_smiles(smiles)
            rows.append({"index": index, "smiles": smiles, "status": "ok", profile.field: payload(molecule, profile)})
        except Exception as error:
            rows.append({"index": index, "smiles": smiles, "status": "error", "error": str(error)})
    return rows


def external_records(command: list[str], corpus: str) -> list[dict]:
    output = subprocess.run(command, cwd=ROOT, input=corpus, text=True, capture_output=True, check=True).stdout
    return [json.loads(line) for line in output.splitlines() if line.strip()]


def validate_rows(rows: list[dict], inputs: list[tuple[int, str]], field: str, lane: str) -> None:
    if len(rows) != len(inputs):
        raise ValueError(f"{lane}: {len(rows)} rows, expected {len(inputs)}")
    for position, (row, (index, smiles)) in enumerate(zip(rows, inputs, strict=True)):
        if not isinstance(row, dict) or row.get("index") != index or row.get("smiles") != smiles:
            raise ValueError(f"{lane}: input identity mismatch at position {position}")
        if row.get("status") not in {"ok", "error"}:
            raise ValueError(f"{lane}: missing or invalid status at position {position}")
        if row["status"] == "ok" and field not in row:
            raise ValueError(f"{lane}: missing {field} at position {position}")


def compare(left: list[dict], right: list[dict], profile: Profile) -> dict:
    examples = []
    matches = 0
    for position, (first, second) in enumerate(zip(left, right, strict=True)):
        same = first == second if profile.exact_rows else (
            first["index"] == second["index"]
            and first["smiles"] == second["smiles"]
            and first["status"] == second["status"]
            and (first["status"] != "ok" or first[profile.field] == second[profile.field])
        )
        if same:
            matches += 1
        elif len(examples) < 10:
            examples.append({"index": position, "left": first, "right": second})
    return {"rows": len(left), "matches": matches, "mismatches": len(left) - matches, "examples": examples}


def build_report(profile: Profile, corpus_path: Path, corpus: str, version: str) -> dict:
    inputs = [(index, raw.strip()) for index, raw in enumerate(corpus.splitlines()) if raw.strip()]
    bindings = {
        "rust": external_records(["cargo", "run", "-p", "chematic-chem", "--release", "--offline", "--example", profile.rust_example], corpus),
        "python": python_records(corpus, profile),
        "node_wasm": external_records(["node", "scripts/diagnostics/binding_dump.mjs", profile.node_mode], corpus),
    }
    for lane, rows in bindings.items():
        validate_rows(rows, inputs, profile.field, lane)
    pairwise = {
        "rust_vs_python": compare(bindings["rust"], bindings["python"], profile),
        "rust_vs_node_wasm": compare(bindings["rust"], bindings["node_wasm"], profile),
        "python_vs_node_wasm": compare(bindings["python"], bindings["node_wasm"], profile),
    }
    try:
        corpus_name = str(corpus_path.relative_to(ROOT))
    except ValueError:
        corpus_name = str(corpus_path)
    return {
        "schema_version": 1,
        "target_version": version,
        "contract": profile.contract,
        "corpus": {"path": corpus_name, "rows": len(inputs), "sha256": hashlib.sha256(corpus_path.read_bytes()).hexdigest()},
        "binding_status_counts": {lane: {status: sum(row["status"] == status for row in rows) for status in ("ok", "error")} for lane, rows in bindings.items()},
        "pairwise": pairwise,
        "gate_passed": all(result["mismatches"] == 0 for result in pairwise.values()),
        "boundary": list(profile.boundary),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=tuple(PROFILES), required=True)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replace", action="store_true", help="allow replacing an existing output")
    args = parser.parse_args()
    source = args.corpus.resolve()
    destination = args.output.resolve()
    if destination.exists() and not args.replace:
        parser.error(f"output already exists: {destination}; choose a new path or pass --replace")
    version = workspace_version()
    provenance = runtime_provenance(version)
    report = build_report(PROFILES[args.profile], source, source.read_text(encoding="utf-8"), version)
    report["runtime_provenance"] = provenance
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(destination), "rows": report["corpus"]["rows"], "gate_passed": report["gate_passed"]}))
    return 0 if report["gate_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
