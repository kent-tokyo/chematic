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
