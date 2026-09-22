#!/usr/bin/env python3
"""Run the published chematic wheel through the fixed RDKit 3D corpus.

This runner is intentionally Python-only: it must be launched with an isolated
environment containing the registry wheel under test.  It never imports the
checkout's Rust extension.  Output is JSONL compatible with the existing
pipeline-v2/RDKit common geometry scorer.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RANDOM_SEED = 20260801
MAX_ATTEMPTS = 8
FORCE_FIELD_MAX_ITERATIONS = 200
TOTAL_TIMEOUT_MS = 20_000
BEST_OF_N = 10

ARM_FORCE_FIELDS = {
    "chematic_pipeline_v2_no_ff": "none",
    "chematic_pipeline_v2_uff_only": "uff_only",
    "chematic_pipeline_v2_mmff94_strict": "mmff94_bond_angle_strict",
}
STEREO_SAFE_ARM_FORCE_FIELDS = {
    "chematic_pipeline_v2_uff_only_stereo_safe": "uff_only",
    "chematic_pipeline_v2_mmff94_strict_stereo_safe": "mmff94_bond_angle_strict",
}
BEST_OF_N_ARM = "chematic_pipeline_v2_uff_best_of_10"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_rows(tiers: str) -> list[dict[str, object]]:
    selected = set(tiers)
    rows: list[dict[str, object]] = []
    for tier, relative in (
        ("A", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json"),
        ("B", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"),
    ):
        if tier not in selected:
            continue
        manifest = json.loads((ROOT / relative).read_text(encoding="utf-8"))
        for molecule in manifest["molecules"]:
            rows.append({"tier": tier, **molecule})
    return rows


def pipeline_config(
    chematic: object, force_field: str, *, stereo_safe: bool = False
) -> object:
    # Match the established source-vs-RDKit benchmark contract: do not let the
    # engine under test judge its own stereo output, and do not turn a missing
    # ring-torsion application into a coverage failure.  The identical external
    # scorer evaluates geometry/stereo for both engines after this run.
    common = dict(
        force_field=force_field,
        ring_torsion_policy="diagnostic_only",
        fail_on_unevaluable_stereo=False,
        embed_seed=RANDOM_SEED,
        max_attempts=MAX_ATTEMPTS,
        use_exp_torsions=True,
        use_small_ring_torsions=True,
        use_macrocycle_torsions=True,
        use_macrocycle_14_bounds=True,
        include_legacy_torsion_heuristic=False,
        force_field_max_iterations=FORCE_FIELD_MAX_ITERATIONS,
        gate_mmff94_torsion_oop=False,
        gate_mmff94_stretch_bend=False,
        total_timeout_ms=TOTAL_TIMEOUT_MS,
    )
    if stereo_safe:
        return chematic.PipelineV2Config.stereo_safe(**common)  # type: ignore[attr-defined]
    return chematic.PipelineV2Config.safe(  # type: ignore[attr-defined]
        **common,
        stereo_policy="ignore",
        enforce_chirality=False,
        expand_implicit_h_through_pipeline=False,
    )


def run_pipeline(mol: object, config: object) -> dict[str, object]:
    started = time.perf_counter_ns()
    try:
        result = mol.embed_pipeline_v2(config)  # type: ignore[attr-defined]
    except Exception as error:  # retain every typed or unexpected failure
        return {
            "status": "typed_failure",
            "failure_cause": f"{type(error).__name__}: {error}",
            "elapsed_ms": (time.perf_counter_ns() - started) / 1_000_000,
        }
    elapsed_ms = (time.perf_counter_ns() - started) / 1_000_000
    validation = result.get("final_validation")
    force_field = result.get("force_field")
    return {
        "status": "success",
        "elapsed_ms": elapsed_ms,
        "coords": result["coords"],
        "force_field": force_field,
        "final_validation": validation,
        "pipeline_stage_ms": result.get("elapsed_ms_by_stage"),
    }


def run_best_of_n(mol: object, chematic: object, config: object) -> dict[str, object]:
    ensemble_config = chematic.EnsembleV2Config(
        per_conformer=config,
        count=BEST_OF_N,
        base_seed=RANDOM_SEED,
        rmsd_threshold=0.5,
        use_symmetric_rmsd_pruning=True,
        ensemble_timeout_ms=TOTAL_TIMEOUT_MS * BEST_OF_N,
    )
    started = time.perf_counter_ns()
    try:
        result = mol.conformer_ensemble_v2(ensemble_config)
    except Exception as error:
        return {
            "status": "typed_failure",
            "failure_cause": f"{type(error).__name__}: {error}",
            "elapsed_ms": (time.perf_counter_ns() - started) / 1_000_000,
        }
    elapsed_ms = (time.perf_counter_ns() - started) / 1_000_000
    conformers = result.get("conformers", [])
    if not conformers:
        return {
            "status": "typed_failure",
            "failure_cause": "no_conformer_kept",
            "elapsed_ms": elapsed_ms,
            "termination": result.get("termination"),
            "attempts": result.get("attempts"),
        }
    provenance = result.get("conformer_provenance", [])
    return {
        "status": "success",
        "elapsed_ms": elapsed_ms,
        "coords": conformers[0],
        "force_field": "uff",
        "conformers_kept": len(conformers),
        "conformer_provenance": provenance,
        "termination": result.get("termination"),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata-output", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--tiers", choices=("A", "B", "AB"), default="AB")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--count", type=int)
    parser.add_argument(
        "--arms",
        nargs="+",
        choices=(*ARM_FORCE_FIELDS, *STEREO_SAFE_ARM_FORCE_FIELDS, BEST_OF_N_ARM),
        default=[*ARM_FORCE_FIELDS, BEST_OF_N_ARM],
    )
    parser.add_argument("--wall-budget-seconds", type=float, default=3300.0)
    parser.add_argument(
        "--package-kind",
        choices=("published", "source_candidate"),
        default="published",
    )
    parser.add_argument("--source-revision")
    parser.add_argument("--source-diff-sha256")
    args = parser.parse_args()
    if args.start < 0 or (args.count is not None and args.count <= 0):
        parser.error("--start must be non-negative and --count must be positive")
    if args.wall_budget_seconds <= 0:
        parser.error("--wall-budget-seconds must be positive")
    if not args.wheel.is_file():
        parser.error(f"wheel not found: {args.wheel}")
    if args.package_kind == "source_candidate" and (
        not args.source_revision or not args.source_diff_sha256
    ):
        parser.error(
            "source candidates require --source-revision and --source-diff-sha256"
        )
    return args


def main() -> int:
    args = parse_args()
    import chematic

    engine = (
        "chematic_public_python"
        if args.package_kind == "published"
        else "chematic_source_candidate_python"
    )

    all_rows = load_rows(args.tiers)
    selected = all_rows[args.start :]
    if args.count is not None:
        selected = selected[: args.count]
    configs = {
        arm: pipeline_config(chematic, force_field)
        for arm, force_field in ARM_FORCE_FIELDS.items()
        if arm in args.arms
    }
    configs.update(
        {
            arm: pipeline_config(chematic, force_field, stereo_safe=True)
            for arm, force_field in STEREO_SAFE_ARM_FORCE_FIELDS.items()
            if arm in args.arms
        }
    )
    best_config = (
        pipeline_config(chematic, "uff_only") if BEST_OF_N_ARM in args.arms else None
    )
    session_started = time.monotonic()
    emitted = 0
    termination = "completed"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as stream:
        for row_index, row in enumerate(selected, start=args.start):
            if time.monotonic() - session_started >= args.wall_budget_seconds:
                termination = "wall_budget_exhausted"
                break
            smiles = str(row["smiles"])
            try:
                mol = chematic.from_smiles(smiles)
            except Exception as error:
                result_rows = [
                    {
                        "arm": arm,
                        "status": "parse_failure",
                        "failure_cause": f"{type(error).__name__}: {error}",
                        "elapsed_ms": 0.0,
                    }
                    for arm in args.arms
                ]
            else:
                result_rows = [
                    {
                        "arm": arm,
                        **(
                            run_best_of_n(mol, chematic, best_config)
                            if arm == BEST_OF_N_ARM
                            else run_pipeline(mol, configs[arm])
                        ),
                    }
                    for arm in args.arms
                ]
            for result in result_rows:
                output_row = {
                    "tier": row["tier"],
                    "row_index": row_index,
                    "name": row["name"],
                    "smiles": smiles,
                    "primary_category": row.get("primary_category", "unknown"),
                    "engine": engine,
                    **result,
                }
                stream.write(json.dumps(output_row, separators=(",", ":")) + "\n")
                stream.flush()
                emitted += 1

    metadata = {
        "schema_version": 1,
        "benchmark": "public-package-3d-vs-rdkit"
        if args.package_kind == "published"
        else "source-candidate-3d-vs-rdkit",
        "engine": engine,
        "package": {
            "name": "chematic",
            "distribution_version": importlib.metadata.version("chematic"),
            "runtime_version": getattr(chematic, "__version__", None),
            "module_path": chematic.__file__,
            "wheel": {
                "path": str(args.wheel),
                "bytes": args.wheel.stat().st_size,
                "sha256": sha256(args.wheel),
            },
            "kind": args.package_kind,
            "source_revision": args.source_revision,
            "source_diff_sha256": args.source_diff_sha256,
        },
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version,
            "python_executable": sys.executable,
        },
        "configuration": {
            "random_seed": RANDOM_SEED,
            "max_attempts": MAX_ATTEMPTS,
            "force_field_max_iterations": FORCE_FIELD_MAX_ITERATIONS,
            "total_timeout_ms": TOTAL_TIMEOUT_MS,
            "best_of_n": BEST_OF_N,
            "stereo_policy": {
                "comparable_arms": "ignore_external_common_scorer_is_authoritative",
                "stereo_safe_arms": "repair_and_verify_with_chirality_and_explicit_h",
            },
            "ring_torsion_policy": "diagnostic_only",
            "enforce_chirality": False,
            "expand_implicit_h_through_pipeline": False,
            "tiers": args.tiers,
            "start": args.start,
            "requested_count": args.count,
            "arms": args.arms,
            "wall_budget_seconds": args.wall_budget_seconds,
        },
        "result": {
            "termination": termination,
            "selected_molecules": len(selected),
            "emitted_rows": emitted,
            "output": str(args.output),
            "output_sha256": sha256(args.output),
            "elapsed_seconds": time.monotonic() - session_started,
        },
        "collected_at_utc": datetime.now(timezone.utc).isoformat(),
        "boundary": (
            "Published registry wheel only. "
            if args.package_kind == "published"
            else "Locally built source candidate; not a registry release. "
        )
        + "Parsing is outside per-arm timing. Comparable arms use stereo_policy=ignore and treat the identical external scorer as authoritative. Explicit *_stereo_safe diagnostic arms use the package's bundled repair_and_verify + chirality + explicit-H configuration. All arms keep ring_torsion_policy=diagnostic_only and retain coordinates and failures for the external scorer.",
    }
    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_output.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata["result"], indent=2))
    return 0 if termination == "completed" else 3


if __name__ == "__main__":
    raise SystemExit(main())
