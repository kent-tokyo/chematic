import importlib.util
import json
import sys
from pathlib import Path
from urllib.error import HTTPError


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "fetch_chembl_sealed_candidate.py"
SPEC = importlib.util.spec_from_file_location("fetch_chembl_sealed_candidate", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FakeResponse:
    def __init__(self, body: bytes):
        self.body = body

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False

    def read(self) -> bytes:
        return self.body


def test_fetch_retries_then_reuses_cached_response(tmp_path, monkeypatch):
    body = json.dumps(
        {
            "molecules": [
                {"molecule_structures": {"canonical_smiles": "CCO"}},
            ]
        }
    ).encode()
    calls = 0

    def flaky_urlopen(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise HTTPError("https://example.invalid", 500, "test", {}, None)
        return FakeResponse(body)

    output = tmp_path / "candidate"
    argv = [
        str(SCRIPT),
        "--output-dir",
        str(output),
        "--offset",
        "0",
        "--source-rows",
        "100",
        "--page-size",
        "100",
        "--delay-seconds",
        "0",
        "--max-attempts",
        "2",
        "--retry-delay-seconds",
        "0",
    ]
    monkeypatch.setattr(MODULE, "urlopen", flaky_urlopen)
    monkeypatch.setattr(MODULE, "verified_tls_context", lambda: None)
    monkeypatch.setattr(sys, "argv", argv)
    assert MODULE.main() == 0
    assert calls == 2
    assert (output / "source.smi").read_text() == "CCO\n"

    monkeypatch.setattr(
        MODULE,
        "urlopen",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("cache was ignored")),
    )
    assert MODULE.main() == 0
    assert calls == 2

