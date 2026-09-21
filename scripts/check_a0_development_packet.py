#!/usr/bin/env python3
"""Validate the public A0 development packet without reading sealed raw data."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def relative_file(value: object) -> Path:
    if not isinstance(value, str):
        raise ValueError("evidence path must be a string")
    path = (ROOT / value).resolve()
    if ROOT not in path.parents or not path.is_file():
        raise ValueError(f"evidence path is outside the repository or missing: {value}")
    return path


def validate(packet_path: Path) -> list[str]:
    packet = load(packet_path)
    if not isinstance(packet, dict):
        return ["packet must be an object"]
    errors: list[str] = []
    if packet.get("status") != "development_packet_not_adoption":
        errors.append("packet must declare development-only status")
    if "not an unused-data evaluation" not in str(packet.get("scope", "")):
        errors.append("packet scope must forbid unused-data interpretation")

    candidates = packet.get("sealed_candidates")
    if not isinstance(candidates, list) or len(candidates) != 3:
        errors.append("packet must preserve exactly three sealed candidate summaries")
    else:
        for entry in candidates:
            if not isinstance(entry, dict):
                errors.append("sealed candidate entry must be an object")
                continue
            try:
                path = relative_file(entry.get("path"))
            except ValueError as error:
                errors.append(str(error))
                continue
            if digest(path) != entry.get("sha256"):
                errors.append(f"sealed summary digest changed: {path}")
            summary = load(path)
            if not isinstance(summary, dict) or summary.get("status") != entry.get("expected_status"):
                errors.append(f"sealed summary status changed: {path}")
            if "raw_rows" in summary or "cases" in summary:
                errors.append(f"sealed summary must not embed raw rows: {path}")

    classification = packet.get("nonsealed_classification")
    if not isinstance(classification, dict):
        errors.append("nonsealed classification is missing")
    else:
        try:
            path = relative_file(classification.get("path"))
            report = load(path)
            if digest(path) != classification.get("sha256"):
                errors.append("nonsealed classification digest changed")
            if not isinstance(report, dict) or report.get("status") != "development_classification":
                errors.append("classification must remain development-only")
            if report.get("rows") != classification.get("expected_rows"):
                errors.append("classification row count changed")
            expected = classification.get("contract_differences", {})
            observed = {
                field: record.get("mismatches")
                for field, record in report.get("summary", {}).items()
                if isinstance(record, dict) and record.get("mismatches")
            }
            if observed != expected:
                errors.append("classification mismatch disposition changed")
        except (ValueError, OSError, json.JSONDecodeError) as error:
            errors.append(f"cannot validate nonsealed classification: {error}")

    tpsa_probe = packet.get("tpsa_atom_type_probe")
    if not isinstance(tpsa_probe, dict):
        errors.append("TPSA atom-type probe is missing")
    else:
        try:
            path = relative_file(tpsa_probe.get("path"))
            report = load(path)
            if digest(path) != tpsa_probe.get("sha256"):
                errors.append("TPSA atom-type probe digest changed")
            if not isinstance(report, dict) or report.get("status") != "development_regression":
                errors.append("TPSA atom-type probe must remain development-only")
            if report.get("rows") != tpsa_probe.get("expected_rows"):
                errors.append("TPSA atom-type probe row count changed")
            if report.get("strict_matches") != tpsa_probe.get("expected_strict_matches"):
                errors.append("TPSA atom-type probe is not strict-green")
            if report.get("strict_tolerance") != tpsa_probe.get("strict_tolerance"):
                errors.append("TPSA atom-type probe tolerance changed")
            if report.get("failures") != []:
                errors.append("TPSA atom-type probe contains failures")
        except (ValueError, OSError, json.JSONDecodeError) as error:
            errors.append(f"cannot validate TPSA atom-type probe: {error}")

    tpsa_corpora = packet.get("tpsa_public_corpus_parity")
    if not isinstance(tpsa_corpora, dict):
        errors.append("TPSA public-corpus parity is missing")
    else:
        try:
            path = relative_file(tpsa_corpora.get("path"))
            report = load(path)
            if digest(path) != tpsa_corpora.get("sha256"):
                errors.append("TPSA public-corpus parity digest changed")
            if not isinstance(report, dict) or report.get("status") != "development_regression":
                errors.append("TPSA public-corpus parity must remain development-only")
            if report.get("rows") != tpsa_corpora.get("expected_rows"):
                errors.append("TPSA public-corpus row count changed")
            if report.get("strict_matches") != tpsa_corpora.get("expected_strict_matches"):
                errors.append("TPSA public-corpus parity is not strict-green")
            if report.get("strict_tolerance") != tpsa_corpora.get("strict_tolerance"):
                errors.append("TPSA public-corpus tolerance changed")
            if report.get("mismatches") != 0:
                errors.append("TPSA public-corpus parity contains mismatches")
        except (ValueError, OSError, json.JSONDecodeError) as error:
            errors.append(f"cannot validate TPSA public-corpus parity: {error}")

    sealed_tpsa = packet.get("sealed_tpsa_resolution")
    if not isinstance(sealed_tpsa, dict):
        errors.append("sealed TPSA resolution is missing")
    else:
        try:
            path = relative_file(sealed_tpsa.get("path"))
            report = load(path)
            if digest(path) != sealed_tpsa.get("sha256"):
                errors.append("sealed TPSA resolution digest changed")
            if not isinstance(report, dict) or report.get("status") != "passed":
                errors.append("sealed TPSA resolution is not passed")
            if report.get("candidate", {}).get("tag") != sealed_tpsa.get("candidate_tag"):
                errors.append("sealed TPSA candidate tag changed")
            if report.get("cohort", {}).get("rows") != sealed_tpsa.get("expected_rows"):
                errors.append("sealed TPSA row count changed")
            if report.get("accounting") != {"parsed": sealed_tpsa.get("expected_rows"), "parse_failures": 0}:
                errors.append("sealed TPSA accounting is incomplete")
            tpsa = report.get("tpsa", {})
            if tpsa.get("strict_matches") != sealed_tpsa.get("expected_strict_matches"):
                errors.append("sealed TPSA is not strict-green")
            if tpsa.get("strict_tolerance") != sealed_tpsa.get("strict_tolerance"):
                errors.append("sealed TPSA tolerance changed")
            if tpsa.get("mismatches") != 0:
                errors.append("sealed TPSA contains mismatches")
            if "raw_rows" in report or "cases" in report:
                errors.append("sealed TPSA summary must not embed raw rows")
        except (ValueError, OSError, json.JSONDecodeError) as error:
            errors.append(f"cannot validate sealed TPSA resolution: {error}")

    binding = packet.get("binding_impact")
    if not isinstance(binding, dict):
        errors.append("binding impact is missing")
    else:
        try:
            path = relative_file(binding.get("path"))
            report = load(path)
            if digest(path) != binding.get("sha256"):
                errors.append("binding impact digest changed")
            if not isinstance(report, dict) or report.get("corpus", {}).get("rows") != binding.get("expected_rows"):
                errors.append("binding impact row count changed")
            for name in binding.get("bindings", []):
                counts = report.get("binding_status_counts", {}).get(name, {})
                if counts != {"ok": binding.get("expected_rows"), "error": 0}:
                    errors.append(f"binding impact is incomplete for {name}")
        except (ValueError, OSError, json.JSONDecodeError) as error:
            errors.append(f"cannot validate binding impact: {error}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path, nargs="?", default=ROOT / "validation/a0-development-packet.json")
    args = parser.parse_args()
    errors = validate(args.packet)
    if errors:
        print("A0 development packet: BLOCKED")
        print("\n".join(f"- {error}" for error in errors))
        return 1
    print("A0 development packet: PASS (three rejected sealed summaries, sealed TPSA 8,000/8,000, 52-row descriptor classification, 63-row TPSA atom-type probe, 10,000-row TPSA public-corpus parity, 7,737-row cross-binding impact)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
