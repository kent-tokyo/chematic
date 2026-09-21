import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).parents[2]
VALIDATOR = ROOT / "scripts/check_a0_development_packet.py"
PACKET = ROOT / "validation/a0-development-packet.json"


def run(packet: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(VALIDATOR), str(packet)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )


def test_a0_development_packet_accepts_current_evidence():
    result = run(PACKET)
    assert result.returncode == 0, result.stdout


def test_a0_development_packet_rejects_relabelled_sealed_summary(tmp_path):
    packet = json.loads(PACKET.read_text(encoding="utf-8"))
    packet["sealed_candidates"][0]["expected_status"] = "accepted"
    path = tmp_path / "packet.json"
    path.write_text(json.dumps(packet), encoding="utf-8")
    result = run(path)
    assert result.returncode != 0
    assert "sealed summary status changed" in result.stdout


def test_a0_development_packet_rejects_incomplete_tpsa_probe(tmp_path):
    packet = json.loads(PACKET.read_text(encoding="utf-8"))
    packet["tpsa_atom_type_probe"]["expected_strict_matches"] -= 1
    path = tmp_path / "packet.json"
    path.write_text(json.dumps(packet), encoding="utf-8")
    result = run(path)
    assert result.returncode != 0
    assert "TPSA atom-type probe is not strict-green" in result.stdout


def test_a0_development_packet_rejects_incomplete_tpsa_public_corpus(tmp_path):
    packet = json.loads(PACKET.read_text(encoding="utf-8"))
    packet["tpsa_public_corpus_parity"]["expected_strict_matches"] -= 1
    path = tmp_path / "packet.json"
    path.write_text(json.dumps(packet), encoding="utf-8")
    result = run(path)
    assert result.returncode != 0
    assert "TPSA public-corpus parity is not strict-green" in result.stdout
