"""The reaction wheel gate must authenticate the loaded native extension."""

from types import SimpleNamespace
from zipfile import ZipFile

import pytest

from scripts.reaction_python_checked_provenance_gate import source_artifact


def test_source_artifact_checks_installed_extension_bytes(tmp_path):
    package = tmp_path / "site-packages" / "chematic"
    package.mkdir(parents=True)
    init = package / "__init__.py"
    init.write_text("", encoding="utf-8")
    extension = package / "chematic.cpython-313-test.so"
    extension.write_bytes(b"compiled native extension")
    wheel = tmp_path / "chematic-1.0.31-py3-none-any.whl"
    with ZipFile(wheel, "w") as archive:
        archive.write(extension, f"chematic/{extension.name}")
    module = SimpleNamespace(__file__=str(init))

    artifact = source_artifact(module, wheel)
    assert artifact["kind"] == "release_profile_source_wheel"
    assert artifact["published"] is False
    assert artifact["wheel"] == wheel.name
    assert len(artifact["wheel_sha256"]) == 64

    extension.write_bytes(b"different installed extension")
    with pytest.raises(ValueError, match="does not match"):
        source_artifact(module, wheel)


def test_candidates_must_cover_the_pinned_fixtures(tmp_path):
    import json

    from scripts.reaction_python_checked_provenance_gate import candidates_rows

    cases = [{"id": "a"}, {"id": "b"}]
    hashes = {"base_sha256": "x", "strata_sha256": "y"}
    row = {"status": "products", "sets": [[{"smiles": "CO", "atom_sources": [[0, 0], [0, 1]],
                                            "template_map_numbers": [None, 1]}]]}

    def write(**changes):
        doc = {"schema": "python-checked-reaction-candidates/v1", "fixtures": hashes,
               "artifact": {"kind": "published_wheel", "published": True},
               "rows": {"a": row, "b": {"status": "no_match", "sets": []}}, **changes}
        path = tmp_path / "candidates.json"
        path.write_text(json.dumps(doc), encoding="utf-8")
        return path

    artifact, rows = candidates_rows(write(), cases, hashes)
    assert artifact["published"] is True and len(artifact["candidates_sha256"]) == 64
    assert [r["status"] for r in rows] == ["products", "no_match"]
    with pytest.raises(ValueError, match="different fixtures"):
        candidates_rows(write(fixtures={"base_sha256": "z", "strata_sha256": "y"}), cases, hashes)
    with pytest.raises(ValueError, match="do not cover"):
        candidates_rows(write(rows={"a": row}), cases, hashes)
    with pytest.raises(ValueError, match="RDKit-readable"):
        bad = {"status": "products", "sets": [[{"smiles": "C(", "atom_sources": [[0, 0]],
                                                "template_map_numbers": [None]}]]}
        candidates_rows(write(rows={"a": bad, "b": row}), cases, hashes)
