from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/changed_public_crates.py"


def load_module():
    spec = importlib.util.spec_from_file_location("changed_public_crates", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PublicCratesTests(unittest.TestCase):
    def test_public_crates_only_returns_library_targets(self) -> None:
        module = load_module()
        crates = dict(module.public_crates())

        self.assertNotIn("chematic-cli", crates)
        self.assertIn("chematic-core", crates)
        for path in crates.values():
            manifest = ROOT / path / "Cargo.toml"
            source = ROOT / path / "src/lib.rs"
            data = module.tomllib.loads(manifest.read_text())
            self.assertTrue("lib" in data or source.is_file())
