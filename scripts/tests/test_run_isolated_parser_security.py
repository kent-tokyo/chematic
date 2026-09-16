import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).parents[2]
MODULE_PATH = ROOT / "scripts/run_isolated_parser_security.py"
SPEC = importlib.util.spec_from_file_location("isolated_parser_security", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def case(case_id: str, fmt: str, expected_status: str) -> dict[str, str]:
    valid = {
        "smiles": "C",
        "smarts": "[#6]",
        "mol_v2000": "methane\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n",
        "mol_v3000": "M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n",
        "sdf": "methane\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n",
    }
    payload = valid[fmt] if expected_status == "accepted" else f"broken-{fmt}-{case_id}"
    return {
        "id": case_id,
        "format": fmt,
        "payload": payload,
        "source_url": "https://example.invalid/public-reproducer",
        "license": "test-only",
        "affected_version": "test",
        "cause": "test fixture",
        "sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "expected_status": expected_status,
    }


def corpus() -> dict[str, list[dict[str, str]]]:
    cases = []
    for fmt in MODULE.FORMATS:
        cases.append(case(f"{fmt}-valid", fmt, "accepted"))
        cases.extend(case(f"{fmt}-invalid-{index}", fmt, "rejected") for index in range(19))
    return {"cases": cases}


def test_load_cases_requires_20_sourced_cases_per_format_and_valid_control(tmp_path):
    path = tmp_path / "corpus.json"
    path.write_text(json.dumps(corpus()), encoding="utf-8")
    loaded = MODULE.load_cases(path)
    assert len(loaded) == 100


def test_load_cases_rejects_digest_drift(tmp_path):
    document = corpus()
    document["cases"][0]["sha256"] = "0" * 64
    path = tmp_path / "corpus.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    try:
        MODULE.load_cases(path)
    except ValueError as error:
        assert "digest mismatch" in str(error)
    else:
        raise AssertionError("digest drift must fail closed")


def test_linux_network_isolation_uses_a_privileged_fresh_namespace(tmp_path):
    command = MODULE.isolated_command(["runner"], True, "smiles", tmp_path / "input")
    assert command[:5] == ["sudo", "--non-interactive", "unshare", "--net", "--"]
    assert command[5:] == ["runner", "--format", "smiles", "--input", str(tmp_path / "input")]


def test_diagnostic_command_does_not_claim_network_isolation(tmp_path):
    command = MODULE.isolated_command(["runner"], False, "smiles", tmp_path / "input")
    assert command == ["runner", "--format", "smiles", "--input", str(tmp_path / "input")]


def test_peak_rss_gate_rejects_missing_boolean_and_over_budget_values():
    assert MODULE.peak_rss_within_limit(256 * 1024, 256)
    assert MODULE.peak_rss_within_limit(0, 256)
    assert not MODULE.peak_rss_within_limit(None, 256)
    assert not MODULE.peak_rss_within_limit(True, 256)
    assert not MODULE.peak_rss_within_limit(256 * 1024 + 1, 256)
