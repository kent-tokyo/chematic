"""
chematic.rdkit_compat — RDKit API compatibility layer for chematic.

Provides a drop-in subset of the RDKit Python API so that common RDKit
workflows can run with minimal changes in environments where RDKit is
unavailable (WASM, serverless, pure-Rust services, no-conda Python envs).

Usage::

    from chematic import rdkit_compat as Chem
    from chematic.rdkit_compat import Descriptors, rdMolDescriptors, DataStructs

    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    print(Descriptors.MolWt(mol))            # 180.16
    print(rdMolDescriptors.CalcTPSA(mol))    # 63.6
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
    print(DataStructs.TanimotoSimilarity(fp, fp))  # 1.0

Known differences vs RDKit:

- ``GetNumAtoms()`` returns heavy-atom count regardless of ``onlyHeavy``.
  Use ``mol.add_hydrogens()`` before calling if you need H-inclusive counts.
- ``GetMorganFingerprintAsBitVect`` uses chematic's ECFP algorithm; Morgan
  bit patterns differ from RDKit's due to different hash implementations.
- ``SanitizeMol`` / ``Kekulize`` are no-ops — chematic sanitizes on parse.
- ``SDMolSupplier.__len__`` performs an O(n) scan; avoid in hot loops.
- Atom/bond object access (``GetAtomWithIdx``, ``GetBonds``, ``GetPropNames``)
  is not supported in this compatibility layer.
- ``MolFromSmarts`` returns a lightweight query object; SMARTS pattern
  matching is delegated to ``chematic.smarts_match``.
"""

from __future__ import annotations

__all__ = [
    "Mol",
    "MolFromSmiles", "MolToSmiles", "MolFromMolBlock", "MolFromMolFile",
    "MolToMolBlock", "SanitizeMol", "Kekulize", "AddHs", "RemoveHs",
    "MolFromSmarts",
    "SDMolSupplier", "SDWriter",
    "Descriptors", "rdMolDescriptors", "DataStructs",
]


import chematic as _ch


# ---------------------------------------------------------------------------
# Mol wrapper
# ---------------------------------------------------------------------------

class Mol:
    """Thin wrapper around ``chematic.Mol`` with RDKit-compatible methods."""

    __slots__ = ("_mol",)

    def __init__(self, mol: _ch.Mol) -> None:
        self._mol = mol

    # -- atom / bond counts --------------------------------------------------

    def GetNumAtoms(self, onlyHeavy: bool = True) -> int:
        # ponytail: always heavy-atom count; AddHs not tracked in this wrapper
        return self._mol.heavy_atoms

    def GetNumHeavyAtoms(self) -> int:
        return self._mol.heavy_atoms

    # -- SMILES / InChI ------------------------------------------------------

    def GetSmiles(self) -> str:
        return self._mol.smiles

    # -- Substructure --------------------------------------------------------

    def HasSubstructMatch(self, query) -> bool:
        """SMARTS string or another Mol (matched by SMILES as query)."""
        if isinstance(query, _SmartsQuery):
            return _ch.smarts_match(query._pattern, self._mol)
        if isinstance(query, Mol):
            return _ch.smarts_match(query._mol.smiles, self._mol)
        if isinstance(query, str):
            return _ch.smarts_match(query, self._mol)
        return False

    def GetSubstructMatches(self, query) -> list:
        if isinstance(query, _SmartsQuery):
            return _ch.smarts_find(query._pattern, self._mol)
        if isinstance(query, Mol):
            return _ch.smarts_find(query._mol.smiles, self._mol)
        if isinstance(query, str):
            return _ch.smarts_find(query, self._mol)
        return []

    # -- MOL block -----------------------------------------------------------

    def ToMolBlock(self) -> str:
        return self._mol.to_mol_block()

    # -- repr ----------------------------------------------------------------

    def __repr__(self) -> str:
        return f"Mol({self._mol.smiles!r})"


class _SmartsQuery:
    """Internal: stores a raw SMARTS string for substructure matching."""
    __slots__ = ("_pattern",)

    def __init__(self, pattern: str) -> None:
        self._pattern = pattern


# ---------------------------------------------------------------------------
# Module-level functions (Chem.*)
# ---------------------------------------------------------------------------

def MolFromSmiles(smi: str, sanitize: bool = True) -> Mol | None:
    """Parse SMILES, return ``Mol`` or ``None`` on failure."""
    try:
        return Mol(_ch.from_smiles(smi))
    except Exception:
        return None


def MolToSmiles(
    mol: Mol,
    isomericSmiles: bool = True,
    kekuleSmiles: bool = False,
    canonical: bool = True,
) -> str:
    """Return canonical SMILES for *mol*."""
    return mol._mol.smiles


def MolFromMolBlock(
    block: str,
    sanitize: bool = True,
    removeHs: bool = True,
) -> Mol | None:
    """Parse a MOL/SDF block string."""
    try:
        return Mol(_ch.from_mol_block(block))
    except Exception:
        return None


def MolFromMolFile(
    filename: str,
    sanitize: bool = True,
    removeHs: bool = True,
) -> Mol | None:
    """Read the first molecule from a MOL file."""
    try:
        text = open(filename).read()
        return MolFromMolBlock(text, sanitize=sanitize, removeHs=removeHs)
    except Exception:
        return None


def MolToMolBlock(
    mol: Mol,
    includeStereo: bool = True,
    confId: int = -1,
) -> str:
    """Return a MOL block string for *mol*."""
    return mol._mol.to_mol_block()


def MolFromSmarts(pattern: str) -> _SmartsQuery:
    """Parse a SMARTS pattern for use with ``HasSubstructMatch``."""
    return _SmartsQuery(pattern)


def SanitizeMol(mol: Mol) -> None:
    """No-op: chematic sanitizes on parse."""


def Kekulize(mol: Mol, clearAromaticFlags: bool = False) -> None:
    """No-op: chematic kekulizes on demand internally."""


def AddHs(mol: Mol, addCoords: bool = False) -> Mol:
    """Return a new ``Mol`` with explicit hydrogens added."""
    return Mol(mol._mol.add_hydrogens())


def RemoveHs(mol: Mol) -> Mol:
    """Return a new ``Mol`` with explicit hydrogens removed."""
    return Mol(mol._mol.remove_hydrogens())


# ---------------------------------------------------------------------------
# SDMolSupplier — streaming SDF reader
# ---------------------------------------------------------------------------

class SDMolSupplier:
    """Iterate over molecules in an SDF file, yielding ``Mol | None``.

    .. note::
       ``__len__`` performs an O(n) file scan. Avoid in hot loops.
    """

    def __init__(
        self,
        filename: str,
        sanitize: bool = True,
        removeHs: bool = True,
    ) -> None:
        self._filename = filename

    def __iter__(self):
        for rec in _ch.iter_sdf(self._filename):
            if rec.mol is not None:
                yield Mol(rec.mol)
            else:
                yield None

    def __len__(self) -> int:
        # ponytail: O(n) scan; document limitation
        return sum(1 for _ in _ch.iter_sdf(self._filename))


# ---------------------------------------------------------------------------
# SDWriter — streaming SDF writer
# ---------------------------------------------------------------------------

class SDWriter:
    """Write molecules to an SDF file, one record at a time.

    Supports context-manager protocol::

        with SDWriter("out.sdf") as w:
            w.write(mol)
    """

    def __init__(self, filename: str) -> None:
        self._f = open(filename, "w")

    def write(self, mol: Mol, confId: int = -1) -> None:
        """Append *mol* as an SDF record (MOL block + ``$$$$``)."""
        block = mol._mol.to_mol_block()
        self._f.write(block)
        if not block.endswith("\n"):
            self._f.write("\n")
        self._f.write("$$$$\n")

    def flush(self) -> None:
        self._f.flush()

    def close(self) -> None:
        self._f.close()

    def __enter__(self) -> "SDWriter":
        return self

    def __exit__(self, *args) -> None:
        self.close()


# ---------------------------------------------------------------------------
# Descriptors namespace
# ---------------------------------------------------------------------------

class Descriptors:
    """Mirrors ``rdkit.Chem.Descriptors``."""

    @staticmethod
    def MolWt(mol: Mol) -> float:
        return mol._mol.mw

    @staticmethod
    def ExactMolWt(mol: Mol) -> float:
        return mol._mol.exact_mass

    @staticmethod
    def MolLogP(mol: Mol) -> float:
        return mol._mol.logp

    @staticmethod
    def TPSA(mol: Mol) -> float:
        return mol._mol.tpsa

    @staticmethod
    def NumHDonors(mol: Mol) -> int:
        return mol._mol.hbd

    @staticmethod
    def NumHAcceptors(mol: Mol) -> int:
        return mol._mol.hba

    @staticmethod
    def NumRotatableBonds(mol: Mol) -> int:
        return mol._mol.rotatable_bonds

    @staticmethod
    def NumHeavyAtoms(mol: Mol) -> int:
        return mol._mol.heavy_atoms

    @staticmethod
    def FractionCSP3(mol: Mol) -> float:
        return mol._mol.fsp3

    @staticmethod
    def MolMR(mol: Mol) -> float:
        return mol._mol.molar_refractivity


# ---------------------------------------------------------------------------
# rdMolDescriptors namespace
# ---------------------------------------------------------------------------

class rdMolDescriptors:
    """Mirrors ``rdkit.Chem.rdMolDescriptors``."""

    @staticmethod
    def CalcTPSA(mol: Mol, includeSandP: bool = True) -> float:
        return mol._mol.tpsa

    @staticmethod
    def CalcNumHBA(mol: Mol) -> int:
        return mol._mol.hba

    @staticmethod
    def CalcNumHBD(mol: Mol) -> int:
        return mol._mol.hbd

    @staticmethod
    def CalcNumAtomStereoCenters(mol: Mol) -> int:
        return mol._mol.num_stereocenters

    @staticmethod
    def CalcNumRotatableBonds(mol: Mol) -> int:
        return mol._mol.rotatable_bonds

    @staticmethod
    def CalcNumHeavyAtoms(mol: Mol) -> int:
        return mol._mol.heavy_atoms

    @staticmethod
    def CalcExactMolWt(mol: Mol) -> float:
        return mol._mol.exact_mass

    @staticmethod
    def CalcFractionCSP3(mol: Mol) -> float:
        return mol._mol.fsp3

    @staticmethod
    def CalcNumAromaticRings(mol: Mol) -> int:
        return mol._mol.aromatic_ring_count

    @staticmethod
    def CalcNumRings(mol: Mol) -> int:
        return mol._mol.aromatic_ring_count + mol._mol.num_aliphatic_rings

    @staticmethod
    def CalcNumSpiroAtoms(mol: Mol) -> int:
        return mol._mol.num_spiro_atoms

    @staticmethod
    def CalcNumBridgeheadAtoms(mol: Mol) -> int:
        return mol._mol.num_bridgehead_atoms

    @staticmethod
    def CalcNumAmideBonds(mol: Mol) -> int:
        return mol._mol.num_amide_bonds

    @staticmethod
    def CalcNumHeteroatoms(mol: Mol) -> int:
        return mol._mol.num_heteroatoms

    @staticmethod
    def GetMorganFingerprintAsBitVect(
        mol: Mol,
        radius: int,
        nBits: int = 2048,
        useChirality: bool = False,
    ) -> bytes:
        """Return ECFP fingerprint as bytes.

        .. note::
           Bit patterns differ from RDKit's Morgan fingerprints due to
           different hash implementations. Use for similarity ranking, not
           for cross-library bit-level comparison.
        """
        # ponytail: radius→ECFP variant; nBits other than 2048 not supported
        if radius <= 2:
            return mol._mol.ecfp4()
        else:
            return mol._mol.ecfp6()


# ---------------------------------------------------------------------------
# DataStructs namespace
# ---------------------------------------------------------------------------

class DataStructs:
    """Mirrors ``rdkit.DataStructs``."""

    @staticmethod
    def TanimotoSimilarity(fp1: bytes, fp2: bytes) -> float:
        return _ch.tanimoto(fp1, fp2)

    @staticmethod
    def DiceSimilarity(fp1: bytes, fp2: bytes) -> float:
        return _ch.dice_similarity(fp1, fp2)
