"""Shared fixtures for chematic pytest suite."""
import pytest

pytest_plugins: list[str] = []


@pytest.fixture(scope="session")
def aspirin():
    import chematic
    return chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")


@pytest.fixture(scope="session")
def caffeine():
    import chematic
    return chematic.from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")


@pytest.fixture(scope="session")
def ethanol():
    import chematic
    return chematic.from_smiles("CCO")


@pytest.fixture(scope="session")
def benzene():
    import chematic
    return chematic.from_smiles("c1ccccc1")


@pytest.fixture(scope="session")
def biphenyl():
    import chematic
    return chematic.from_smiles("c1ccc(-c2ccccc2)cc1")
