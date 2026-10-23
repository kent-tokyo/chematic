"""Unit tests for release-to-release binding compatibility helpers."""

import unittest

from scripts.checks.check_binding_api_compat import (
    dts_surface,
    json_contract_surface,
    public_python_api,
)


class BindingApiCompatibilityTests(unittest.TestCase):
    def test_python_surface_keeps_least_restrictive_overload(self) -> None:
        surface = public_python_api(
            """
from typing import overload
@overload
def parse(text: str) -> int: ...
@overload
def parse(text: str, strict: bool = ...) -> int: ...
"""
        )
        self.assertEqual(surface["parse"], 1)

    def test_wasm_surface_includes_exported_class_members(self) -> None:
        surface = dts_surface(
            """
export function parse(text: string): Mol;
export class Mol {
  readonly atoms: number;
  smiles(): string;
}
"""
        )
        self.assertEqual(surface, {"parse", "Mol", "Mol.atoms", "Mol.smiles"})

    def test_json_contract_surface_tracks_contract_and_bindings(self) -> None:
        surface = json_contract_surface(
            '{"operations":[{"id":"parse","contract":"fixtures",'
            '"bindings":["wasm","rust"]}]}'
        )
        self.assertEqual(surface, {("parse", "fixtures", ("rust", "wasm"))})


if __name__ == "__main__":
    unittest.main()
