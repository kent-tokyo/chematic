from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/changed_public_crates.py"


def load_module():
    spec = importlib.util.spec_from_file_location("changed_public_crates", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_public_crates_only_returns_library_targets() -> None:
    module = load_module()
    crates = dict(module.public_crates())

    assert "chematic-cli" not in crates
    assert "chematic-core" in crates
    for path in crates.values():
        manifest = ROOT / path / "Cargo.toml"
        source = ROOT / path / "src/lib.rs"
        data = module.tomllib.loads(manifest.read_text())
        assert "lib" in data or source.is_file()
