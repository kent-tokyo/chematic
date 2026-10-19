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

Common Morgan imports from both ``rdkit.Chem.AllChem`` and
``rdkit.Chem.rdFingerprintGenerator`` are also supported. Unsupported options
raise explicitly rather than silently changing the fingerprint algorithm.

Known differences vs RDKit:

- ``GetMorganFingerprintAsBitVect`` uses chematic's ECFP algorithm; Morgan
  bit patterns differ from RDKit's due to different hash implementations.
- ``SanitizeMol`` / ``Kekulize`` are no-ops — chematic sanitizes on parse.
- ``SDMolSupplier.__len__`` performs an O(n) scan; avoid in hot loops.
- ``Mol``'s atom/bond objects (``GetAtomWithIdx``, ``GetBonds``) are
  read-only; use ``RWMol`` for structure editing (``AddAtom``/``AddBond``/
  ``RemoveAtom``/``RemoveBond``/``GetMol``) — a smaller subset than RDKit's
  ``RWMol`` (no mid-edit atom/bond iteration; call ``GetMol()`` first).
- ``MolFromSmarts`` returns a lightweight query object; SMARTS pattern
  matching is delegated to ``chematic.smarts_match``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
    from ._numpy_fingerprint import _BitVectorSource

__all__ = [
    "Mol", "Atom", "Bond", "BondType", "RingInfo", "RWMol",
    "ExplicitBitVect", "LayeredFingerprint", "LayeredAtomCounts", "LayeredBitMask", "rdMIF",
    "PatternFingerprint", "PatternAtomCounts", "PatternBitMask",
    "MolFromSmiles", "MolToSmiles", "MolFromMolBlock", "MolFromMolFile",
    "MolToMolBlock", "SanitizeMol", "Kekulize", "AddHs", "RemoveHs",
    "MolFromSmarts", "RenumberAtoms", "GetAdjacencyMatrix", "GetDistanceMatrix",
    "CanonSmiles", "AllChem", "rdFingerprintGenerator",
    "MolFromMrvBlock", "MolFromMrvFile", "MolToMrvBlock", "MolToMrvFile",
    "SDMolSupplier", "SDWriter", "SmilesMolSupplier", "SmilesWriter",
    "TDTMolSupplier", "TDTWriter",
    "Descriptors", "rdMolDescriptors", "DataStructs", "pyAvalonTools", "rdAlignment", "rdInfoTheory", "rdMetricMatrixCalc", "MetricBitVector", "Scoring", "RankedScoreRow",
]


from ._point_alignment import rdAlignment
from ._info_theory import rdInfoTheory
from ._graph_matrices import (
    AdjacencyMatrix as _AdjacencyMatrix,
    DistanceMatrix as _DistanceMatrix,
    adjacency_matrix as _adjacency_matrix,
    distance_matrix as _distance_matrix,
)
from ._metric_matrix import MetricBitVector as MetricBitVector, rdMetricMatrixCalc as rdMetricMatrixCalc
from ._scoring_enrichment import RankedScoreRow as RankedScoreRow, Scoring as Scoring

import chematic as _ch

from typing import Optional, Sequence, SupportsIndex, cast
from ._layered_counts import (
    LayeredAtomCounts as LayeredAtomCounts,
    LayeredBitMask as LayeredBitMask,
    LayeredResult as _LayeredResult,
    prepare as _prepare_layered,
)
from ._coulomb_field import Coulomb as _Coulomb
from ._coulomb_dielectric import CoulombDielectric as _CoulombDielectric
from ._pattern_inputs import (
    PatternAtomCounts as PatternAtomCounts,
    PatternBitMask as PatternBitMask,
    prepare as _prepare_pattern,
)


_MISSING = object()


# ---------------------------------------------------------------------------
# ExplicitBitVect — RDKit-compatible bit vector
# ---------------------------------------------------------------------------

class ExplicitBitVect:
    """RDKit-compatible bit vector backed by a bytearray (LSB-first)."""

    __slots__ = ("_nbits", "_bits")

    def __init__(self, nBits: int) -> None:
        self._nbits = nBits
        self._bits = bytearray((nBits + 7) // 8)

    def GetNumBits(self) -> int:
        return self._nbits

    def __len__(self) -> int:
        return self._nbits

    def GetNumOnBits(self) -> int:
        return sum(byte.bit_count() for byte in self._bits)

    def GetBit(self, i: int) -> bool:
        if i < 0 or i >= self._nbits:
            raise IndexError(i)
        return bool(self._bits[i >> 3] & (1 << (i & 7)))

    def SetBit(self, i: int) -> None:
        if i < 0 or i >= self._nbits:
            raise IndexError(i)
        self._bits[i >> 3] |= 1 << (i & 7)

    def GetOnBits(self) -> list:
        return [i for i in range(self._nbits) if self.GetBit(i)]

    def ToBitString(self) -> str:
        return "".join("1" if self.GetBit(i) else "0" for i in range(self._nbits))

    def _to_bytes(self) -> bytes:
        return bytes(self._bits)

    @classmethod
    def _from_bytes(cls, data: bytes, nBits: int) -> "ExplicitBitVect":
        bv = cls(nBits)
        n = min(len(data), len(bv._bits))
        bv._bits[:n] = data[:n]
        return bv

    def __repr__(self) -> str:
        return f"ExplicitBitVect({self._nbits})"


def LayeredFingerprint(
    mol: Mol,
    layerFlags: SupportsIndex = 0xFFFFFFFF,
    minPath: SupportsIndex = 1,
    maxPath: SupportsIndex = 7,
    fpSize: SupportsIndex = 2048,
    atomCounts: Optional[LayeredAtomCounts] = None,
    setOnlyBits: Optional[LayeredBitMask] = None,
    branchedPaths: bool = True,
    fromAtoms: Optional[Sequence[SupportsIndex]] = None,
    *,
    maxVisits: SupportsIndex = 1_000_000,
    maxPaths: SupportsIndex = 100_000,
) -> ExplicitBitVect:
    """Compute Layered bits and update checked caller-owned u32 counters.

    Args:
        mol: A compatibility Mol containing the indexed native graph.
        atomCounts: Optional plain list or writable, non-overlapping uint32
            ndarray covering the atoms. Seeds accumulate with u32 wrapping.
        setOnlyBits: Optional binary bit-vector mask with fpSize bits.
        maxVisits: Native enumeration visit budget; exhaustion returns no result.
        maxPaths: Native selected-path budget; exhaustion returns no result.

    Returns:
        An ExplicitBitVect. All calculation and conversion finishes before the
        caller's counters are updated; trailing counter entries stay unchanged.

    Raises:
        TypeError: Wrong input types, including bool in integer fields.
        ValueError: Invalid bounds, unsupported graphs or exhausted budgets.
    """
    if not isinstance(mol, Mol):
        raise TypeError("mol must be a chematic.rdkit_compat.Mol")
    options, counters = _prepare_layered(
        mol.GetNumAtoms(),
        layerFlags,
        minPath,
        maxPath,
        fpSize,
        atomCounts,
        setOnlyBits,
        branchedPaths,
        fromAtoms,
        maxVisits,
        maxPaths,
    )
    result = cast(_LayeredResult, mol._mol.rdkit_layered_config(**options))
    fingerprint = ExplicitBitVect._from_bytes(result["bits"], options["nbits"])
    counters.commit(result["atom_counts"])
    return fingerprint


def PatternFingerprint(
    mol: Mol,
    fpSize: SupportsIndex = 2048,
    atomCounts: Optional[PatternAtomCounts] = None,
    setOnlyBits: Optional[PatternBitMask] = None,
    tautomerFingerprints: bool = False,
) -> ExplicitBitVect:
    """Return the native Pattern fingerprint with checked compatibility inputs.

    Args:
        mol: A concrete compatibility Mol, at most 1024 atoms/4096 bonds.
        fpSize: Power of two from 1 through 2048. Smaller widths fold the native
            2048-bit hashes exactly. Other widths are rejected.
        atomCounts: Optional plain u32-valued list or writable non-overlapping
            1D uint32 ndarray covering all atoms (an empty destination is also
            accepted). Validated but unchanged, matching pinned RDKit behavior.
        setOnlyBits: Optional structural binary vector of fpSize bits. Validated
            but does not filter output, matching pinned RDKit behavior.
        tautomerFingerprints: Only False is supported by the native algorithm.

    Returns:
        An independent ExplicitBitVect. The existing concrete-molecule native
        Pattern algorithm and its chemistry compatibility limits are retained.

    Raises:
        TypeError: Invalid molecule, option, vector or counter type.
        ValueError: Unsupported width/mode, bad counts/mask or graph size limit.

    All options are validated before native calculation. This does not create
    atom contribution counts; use LayeredFingerprint for its documented counter
    behavior. Query molecules and MolBundle overloads are outside this profile.
    """
    if not isinstance(mol, Mol):
        raise TypeError("mol must be a concrete chematic.rdkit_compat.Mol")
    width = _prepare_pattern(
        mol.GetNumAtoms(), mol.GetNumBonds(), fpSize, atomCounts,
        setOnlyBits, tautomerFingerprints,
    )
    raw = mol._mol.rdkit_pattern_fp()
    return _fold_bits(raw, width)


def _fp_bytes(fp) -> bytes:
    return fp._to_bytes() if isinstance(fp, ExplicitBitVect) else fp


def _fold_bits(raw: bytes, nbits: int) -> "ExplicitBitVect":
    """Fold a 2048-bit fingerprint to *nbits* via modulo (RDKit-style nBits)."""
    src = ExplicitBitVect._from_bytes(raw, 2048)
    bv = ExplicitBitVect(nbits)
    for i in src.GetOnBits():
        bv.SetBit(i % nbits)
    return bv


# ---------------------------------------------------------------------------
# BondType constants
# ---------------------------------------------------------------------------

class BondType:
    """RDKit-compatible bond type constants."""
    SINGLE   = "SINGLE"
    DOUBLE   = "DOUBLE"
    TRIPLE   = "TRIPLE"
    AROMATIC = "AROMATIC"
    OTHER    = "OTHER"


# ---------------------------------------------------------------------------
# Atom wrapper
# ---------------------------------------------------------------------------

class Atom:
    """Read-only wrapper around a chematic atom (RDKit-compatible)."""

    __slots__ = ("_mol", "_idx")

    def __init__(self, mol: "Mol", idx: int) -> None:
        self._mol = mol
        self._idx = idx

    def _d(self):
        # (symbol, atomic_num, formal_charge, is_aromatic, implicit_h, degree, is_in_ring)
        return self._mol._get_atom_cache()[self._idx]

    def GetIdx(self) -> int:           return self._idx
    def GetSymbol(self) -> str:        return self._d()[0]
    def GetAtomicNum(self) -> int:     return self._d()[1]
    def GetFormalCharge(self) -> int:  return self._d()[2]
    def GetIsAromatic(self) -> bool:   return self._d()[3]
    def GetTotalNumHs(self) -> int:    return self._d()[4]
    def GetDegree(self) -> int:        return self._d()[5]
    def GetTotalDegree(self) -> int:   return self._d()[5] + self._d()[4]
    def IsInRing(self) -> bool:        return self._d()[6]

    def IsInRingSize(self, n: int) -> bool:
        ri = self._mol._get_ring_info()
        return any(len(r) == n and self._idx in r for r in ri.AtomRings())

    def __repr__(self) -> str:
        return f"Atom({self._idx}, {self._d()[0]})"


# ---------------------------------------------------------------------------
# Bond wrapper
# ---------------------------------------------------------------------------

_BOND_TYPE_TO_DOUBLE = {"SINGLE": 1.0, "DOUBLE": 2.0, "TRIPLE": 3.0, "AROMATIC": 1.5}


class Bond:
    """Read-only wrapper around a chematic bond (RDKit-compatible)."""

    __slots__ = ("_mol", "_idx")

    def __init__(self, mol: "Mol", idx: int) -> None:
        self._mol = mol
        self._idx = idx

    def _d(self):
        # (atom1_idx, atom2_idx, bond_type_str, is_aromatic)
        return self._mol._get_bond_cache()[self._idx]

    def GetIdx(self) -> int:           return self._idx
    def GetBeginAtomIdx(self) -> int:  return self._d()[0]
    def GetEndAtomIdx(self) -> int:    return self._d()[1]
    def GetBeginAtom(self) -> Atom:    return Atom(self._mol, self._d()[0])
    def GetEndAtom(self) -> Atom:      return Atom(self._mol, self._d()[1])
    def GetBondType(self) -> str:      return self._d()[2]
    def GetIsAromatic(self) -> bool:   return self._d()[3]

    def GetBondTypeAsDouble(self) -> float:
        return _BOND_TYPE_TO_DOUBLE.get(self._d()[2], 0.0)

    def GetOtherAtomIdx(self, idx: int) -> int:
        a1, a2 = self._d()[0], self._d()[1]
        if idx == a1: return a2
        if idx == a2: return a1
        raise ValueError(f"Atom {idx} is not an endpoint of bond {self._idx}")

    def IsInRing(self) -> bool:
        return self._mol._get_ring_info().NumBondRings(self._idx) > 0

    def __repr__(self) -> str:
        d = self._d()
        return f"Bond({d[0]}-{d[1]}, {d[2]})"


# ---------------------------------------------------------------------------
# RingInfo wrapper
# ---------------------------------------------------------------------------

class RingInfo:
    """RDKit-compatible ring information wrapper.

    .. note::
       Uses SSSR. Fused-ring SSSR decomposition may differ from RDKit's
       in degenerate cases (e.g. indolizine).
    """

    __slots__ = ("_atom_rings", "_bond_rings", "_atom_rc", "_bond_rc")

    def __init__(self, mol: "Mol") -> None:
        atom_rings_raw = mol._mol.sssr_atom_rings
        self._atom_rings = tuple(tuple(r) for r in atom_rings_raw)

        bond_cache = mol._get_bond_cache()
        bond_lookup = {
            (min(b[0], b[1]), max(b[0], b[1])): i
            for i, b in enumerate(bond_cache)
        }

        bond_rings = []
        for ring in atom_rings_raw:
            n = len(ring)
            br = []
            for i in range(n):
                key = (min(ring[i], ring[(i + 1) % n]), max(ring[i], ring[(i + 1) % n]))
                if key in bond_lookup:
                    br.append(bond_lookup[key])
            bond_rings.append(tuple(br))
        self._bond_rings = tuple(bond_rings)

        n_atoms = len(mol._get_atom_cache())
        ac = [0] * n_atoms
        for ring in self._atom_rings:
            for a in ring:
                ac[a] += 1
        self._atom_rc = tuple(ac)

        n_bonds = len(bond_cache)
        bc = [0] * n_bonds
        for ring in self._bond_rings:
            for b in ring:
                bc[b] += 1
        self._bond_rc = tuple(bc)

    def NumRings(self) -> int:
        return len(self._atom_rings)

    def AtomRings(self):
        return self._atom_rings

    def BondRings(self):
        return self._bond_rings

    def NumAtomRings(self, i: int) -> int:
        if i < 0 or i >= len(self._atom_rc):
            raise IndexError(i)
        return self._atom_rc[i]

    def NumBondRings(self, i: int) -> int:
        if i < 0 or i >= len(self._bond_rc):
            raise IndexError(i)
        return self._bond_rc[i]


# ---------------------------------------------------------------------------
# Mol wrapper
# ---------------------------------------------------------------------------

class Mol:
    """Thin wrapper around ``chematic.Mol`` with RDKit-compatible methods."""

    __slots__ = ("_mol", "_atom_cache", "_bond_cache", "_ring_info_cache")

    def __init__(self, mol: _ch.Mol) -> None:
        self._mol = mol
        self._atom_cache = None
        self._bond_cache = None
        self._ring_info_cache = None

    # -- lazy caches for atom/bond tables ------------------------------------

    def _get_atom_cache(self):
        if self._atom_cache is None:
            self._atom_cache = self._mol.atom_table
        return self._atom_cache

    def _get_bond_cache(self):
        if self._bond_cache is None:
            self._bond_cache = self._mol.bond_table
        return self._bond_cache

    def _get_ring_info(self) -> "RingInfo":
        if self._ring_info_cache is None:
            self._ring_info_cache = RingInfo(self)
        return self._ring_info_cache

    def GetRingInfo(self) -> "RingInfo":
        return self._get_ring_info()

    # -- atom / bond counts --------------------------------------------------

    def GetNumAtoms(self, onlyHeavy: bool = False) -> int:
        return self._mol.heavy_atoms if onlyHeavy else self._mol.atom_count

    def GetNumHeavyAtoms(self) -> int:
        return self._mol.heavy_atoms

    def GetNumBonds(self) -> int:
        return len(self._get_bond_cache())

    # -- atom / bond iteration -----------------------------------------------

    def GetAtoms(self):
        """Iterate over all heavy atoms as ``Atom`` objects."""
        return (Atom(self, i) for i in range(self.GetNumAtoms()))

    def GetBonds(self):
        """Iterate over all bonds as ``Bond`` objects."""
        return (Bond(self, i) for i in range(self.GetNumBonds()))

    def GetAtomWithIdx(self, i: int) -> "Atom":
        if i < 0 or i >= self.GetNumAtoms():
            raise IndexError(i)
        return Atom(self, i)

    def GetBondWithIdx(self, i: int) -> "Bond":
        if i < 0 or i >= self.GetNumBonds():
            raise IndexError(i)
        return Bond(self, i)

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

    def _smarts_find(self, query) -> list:
        if isinstance(query, _SmartsQuery):
            return _ch.smarts_find(query._pattern, self._mol)
        if isinstance(query, Mol):
            return _ch.smarts_find(query._mol.smiles, self._mol)
        if isinstance(query, str):
            return _ch.smarts_find(query, self._mol)
        return []

    def GetSubstructMatch(self, query) -> tuple:
        """First substructure match as a tuple of atom indices (``()`` if none)."""
        matches = self._smarts_find(query)
        return tuple(matches[0]) if matches else ()

    def GetSubstructMatches(self, query, uniquify: bool = True) -> tuple:
        """All substructure matches as a tuple of atom-index tuples."""
        matches = self._smarts_find(query)
        if uniquify:
            seen = set()
            unique = []
            for m in matches:
                key = frozenset(m)
                if key not in seen:
                    seen.add(key)
                    unique.append(tuple(m))
            return tuple(unique)
        return tuple(tuple(m) for m in matches)

    # -- MOL block -----------------------------------------------------------

    def ToMolBlock(self) -> str:
        return self._mol.to_mol_block()

    # -- SD properties (RDKit-compatible) ------------------------------------

    def GetProp(self, key: str, autoConvert: bool = False, default=_MISSING):
        if default is _MISSING:
            return self._mol.GetProp(key, autoConvert)
        return self._mol.GetProp(key, autoConvert, default)

    def SetProp(self, key: str, val: str) -> None:
        self._mol.SetProp(key, val)

    def HasProp(self, key: str) -> bool:
        return self._mol.HasProp(key)

    def GetPropsAsDict(self) -> dict:
        return self._mol.GetPropsAsDict()

    def GetPropNames(self) -> list:
        return self._mol.GetPropNames()

    def ClearProp(self, key: str) -> None:
        self._mol.ClearProp(key)

    def GetName(self) -> str:
        """Return the molecule title, or ``""`` when ``_Name`` is unset."""
        return self._mol.GetName()

    def SetName(self, name: str) -> None:
        """Set the molecule title through RDKit's conventional ``_Name`` property."""
        self._mol.SetName(name)

    def SetIntProp(self, key: str, val: int) -> None:
        self._mol.SetIntProp(key, val)

    def GetIntProp(self, key: str, default=_MISSING) -> int:
        if default is _MISSING:
            return self._mol.GetIntProp(key)
        return self._mol.GetIntProp(key, default)

    def SetUnsignedProp(self, key: str, val: int) -> None:
        self._mol.SetUnsignedProp(key, val)

    def GetUnsignedProp(self, key: str, default=_MISSING) -> int:
        if default is _MISSING:
            return self._mol.GetUnsignedProp(key)
        return self._mol.GetUnsignedProp(key, default)

    def SetDoubleProp(self, key: str, val: float) -> None:
        self._mol.SetDoubleProp(key, val)

    def GetDoubleProp(self, key: str, default=_MISSING) -> float:
        if default is _MISSING:
            return self._mol.GetDoubleProp(key)
        return self._mol.GetDoubleProp(key, default)

    def SetBoolProp(self, key: str, val: bool) -> None:
        self._mol.SetBoolProp(key, val)

    def GetBoolProp(self, key: str, default=_MISSING) -> bool:
        if default is _MISSING:
            return self._mol.GetBoolProp(key)
        return self._mol.GetBoolProp(key, default)

    # -- repr ----------------------------------------------------------------

    def __repr__(self) -> str:
        return f"Mol({self._mol.smiles!r})"


# ---------------------------------------------------------------------------
# RWMol — editable molecule
# ---------------------------------------------------------------------------


class RWMol:
    """RDKit-compatible editable molecule (``AddAtom``/``AddBond``/``RemoveAtom``/``RemoveBond``).

    .. note::
       Supports only the most common ``RWMol`` operations. Unlike RDKit,
       atom/bond iteration and queries (``GetAtoms``, ``GetAtomWithIdx``, …)
       are not available mid-edit — call :meth:`GetMol` first to get a
       read-only :class:`Mol` for those.
    """

    __slots__ = ("_rw",)

    def __init__(self, mol: "Mol | None" = None) -> None:
        self._rw = _ch.RWMol(mol._mol if mol is not None else None)

    def Clear(self) -> None:
        """chematic extension: reset graph, title and properties; preserve snapshots."""
        self._rw.Clear()

    def AddAtom(self, atom) -> int:
        """Add an atom given an atomic number (``int``), element symbol
        (``str``), or any object with ``GetAtomicNum()`` (e.g. an existing
        ``Atom`` from another ``Mol``). Returns the new atom's index."""
        if isinstance(atom, bool):
            raise TypeError(f"invalid atom spec: {atom!r}")
        if isinstance(atom, int):
            atomic_num = atom
        elif isinstance(atom, str):
            atomic_num = _ch.element_atomic_number(atom)
        else:
            atomic_num = atom.GetAtomicNum()
        return self._rw.AddAtom(atomic_num)

    def AddBond(self, beginAtomIdx: int, endAtomIdx: int, order: str = BondType.SINGLE) -> int:
        """Add a bond; returns the molecule's bond count after adding
        (RDKit's convention — not the new bond's own index)."""
        return self._rw.AddBond(beginAtomIdx, endAtomIdx, order)

    def RemoveAtom(self, idx: int) -> None:
        self._rw.RemoveAtom(idx)

    def RemoveBond(self, beginAtomIdx: int, endAtomIdx: int) -> None:
        self._rw.RemoveBond(beginAtomIdx, endAtomIdx)

    def GetNumAtoms(self) -> int:
        return self._rw.GetNumAtoms()

    def GetNumBonds(self) -> int:
        return self._rw.GetNumBonds()

    def GetName(self) -> str:
        return self._rw.GetName()

    def SetName(self, name: str) -> None:
        self._rw.SetName(name)

    def GetMol(self) -> "Mol":
        """Snapshot the current state into an independent, read-only ``Mol``."""
        return Mol(self._rw.GetMol())

    def __repr__(self) -> str:
        return f"RWMol({self.GetNumAtoms()} atoms, {self.GetNumBonds()} bonds)"


class _SmartsQuery:
    """Internal: stores a raw SMARTS string for substructure matching."""
    __slots__ = ("_pattern",)

    def __init__(self, pattern: str) -> None:
        self._pattern = pattern


# ---------------------------------------------------------------------------
# Module-level functions (Chem.*)
# ---------------------------------------------------------------------------

def MolFromSmiles(smi: str, sanitize: bool = True) -> Mol | None:
    """Parse SMILES, return ``Mol`` or ``None`` on failure.

    With ``sanitize=True`` (RDKit default), aromaticity is perceived so that
    Kekulé input (e.g. ``C1=CC2=CC=CC=C2C=C1``) matches aromatic SMARTS like
    ``c`` and round-trips to aromatic canonical SMILES, as in RDKit.
    """
    try:
        mol = _ch.from_smiles(smi)
        if sanitize:
            mol = mol.apply_aromaticity("rdkit_like")
        return Mol(mol)
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


def CanonSmiles(smiles: str, **kwargs) -> str:
    """RDKit-compatible convenience wrapper for canonical SMILES."""
    return MolToSmiles(MolFromSmiles(smiles), **kwargs)


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


def MolFromMrvBlock(
    block: str,
    sanitize: bool = True,
    removeHs: bool = True,
) -> Mol | None:
    """Parse a ChemAxon Marvin (MRV) block string.

    S-groups, polymers, reactions, multicenter bonds, query atoms/bonds,
    R-groups, enhanced stereo groups, and embedded/compressed data are
    deliberately unsupported and yield ``None``, same as any other parse
    failure -- see the ``chematic_mol::mrv`` Rust module docs for the
    full scope boundary.
    """
    try:
        return Mol(_ch.from_mrv_block(block))
    except Exception:
        return None


def MolFromMrvFile(
    filename: str,
    sanitize: bool = True,
    removeHs: bool = True,
) -> Mol | None:
    """Read the first molecule from an MRV file."""
    try:
        text = open(filename).read()
        return MolFromMrvBlock(text, sanitize=sanitize, removeHs=removeHs)
    except Exception:
        return None


def MolToMrvBlock(
    mol: Mol,
    includeStereo: bool = True,
    confId: int = -1,
    kekulize: bool = True,
    prettyPrint: bool = False,
) -> str:
    """Return an MRV block string for *mol*.

    ``confId`` must be ``-1`` -- chematic's ``Mol`` has no multi-conformer
    slot to select from (raises ``NotImplementedError`` otherwise, never
    silently ignored). ``prettyPrint`` is not implemented (output is
    always single-line; raises ``NotImplementedError`` if ``True``).
    """
    if confId != -1:
        raise NotImplementedError("MolToMrvBlock: confId is not supported")
    if prettyPrint:
        raise NotImplementedError("MolToMrvBlock: prettyPrint is not supported")
    return mol._mol.to_mrv_block(kekulize=kekulize, include_stereo=includeStereo)


def MolToMrvFile(
    mol: Mol,
    filename: str,
    includeStereo: bool = True,
    confId: int = -1,
    kekulize: bool = True,
    prettyPrint: bool = False,
) -> None:
    """Write *mol* to an MRV file."""
    block = MolToMrvBlock(
        mol,
        includeStereo=includeStereo,
        confId=confId,
        kekulize=kekulize,
        prettyPrint=prettyPrint,
    )
    with open(filename, "w") as f:
        f.write(block)


def MolFromSmarts(pattern: str) -> _SmartsQuery:
    """Parse a SMARTS pattern for use with ``HasSubstructMatch``."""
    return _SmartsQuery(pattern)


def SanitizeMol(mol: Mol) -> None:
    """No-op: chematic sanitizes on parse."""


def Kekulize(mol: Mol, clearAromaticFlags: bool = False) -> None:
    """No-op: chematic kekulizes on demand internally."""


def RenumberAtoms(mol: Mol, newOrder) -> Mol:
    """Copy with newOrder[new_index] = old_index, preserving stereo and properties.

    Bond indices stay fixed. External coordinate/atom arrays must be reordered
    separately; this wrapper does not store conformers.
    """
    return Mol(mol._mol.renumber_atoms(newOrder))


def GetAdjacencyMatrix(
    mol: Mol,
    useBO: bool = False,
    emptyVal: SupportsIndex = 0,
    force: bool = False,
    prefix: str = "",
) -> _AdjacencyMatrix:
    """Return a checked int32 binary or float64 bond-weighted NumPy matrix.

    Args:
        mol: A compatibility Mol. Native preparation must retain its indices.
        useBO: Use endpoint bond contributions (aromatic 1.5, dative asymmetric).
        emptyVal: Only integer zero is supported; nonzero byte fills are refused.
        force: Both bool values calculate a fresh independent result.
        prefix: Only the empty string is supported; there is no matrix cache.

    Raises:
        TypeError: Wrong Mol, flag, integer or identifier types.
        ValueError: Native unsupported graph, changed indices, nonzero fill,
            or input beyond 1024 atoms / native adjacency bond limits.
        NotImplementedError: Nonempty cache prefix.
    """
    if not isinstance(mol, Mol):
        raise TypeError("mol must be a chematic.rdkit_compat.Mol")
    return _adjacency_matrix(mol._mol, useBO, emptyVal, force, prefix)


def GetDistanceMatrix(
    mol: Mol,
    useBO: bool = False,
    useAtomWts: bool = False,
    force: bool = False,
    prefix: str = "",
) -> _DistanceMatrix:
    """Return a checked float64 NumPy topological distance matrix.

    Args:
        mol: A compatibility Mol. Native preparation must retain its indices.
        useBO: Reciprocal bond-order edge lengths, aromatic edges 2/3.
        useAtomWts: Set diagonal values to 6 / atomic number.
        force: Both bool values calculate a fresh independent result.
        prefix: Only the empty string is supported; there is no matrix cache.

    Disconnected distances use the native RDKit profile's 1e8 sentinel.
    Empty molecules return shape (0, 0). Native model errors or more than
    1024 input atoms raise ValueError; unsupported cache prefixes raise
    NotImplementedError, and wrong argument types raise TypeError.
    """
    if not isinstance(mol, Mol):
        raise TypeError("mol must be a chematic.rdkit_compat.Mol")
    return _distance_matrix(mol._mol, useBO, useAtomWts, force, prefix)


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

    SD properties are attached to each molecule and accessible via
    ``GetProp`` / ``GetPropsAsDict``.

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
        self._sanitize = sanitize
        self._removeHs = removeHs

    def __iter__(self):
        # Use chematic.SDMolSupplier which propagates SD props onto mol.props
        sup = _ch.SDMolSupplier(
            self._filename,
            sanitize=self._sanitize,
            removeHs=self._removeHs,
            strictParsing=False,
        )
        for mol in sup:
            yield Mol(mol) if mol is not None else None

    def __len__(self) -> int:
        # ponytail: O(n) scan; document limitation
        return sum(1 for _ in _ch.iter_sdf(self._filename))

    def __getitem__(self, i: int):
        """Random access by record index (RDKit-compatible).

        .. note:: O(i) — iterates from the start. Avoid in tight loops.
        """
        if i < 0:
            raise IndexError(i)
        for j, mol in enumerate(self):
            if j == i:
                return mol
        raise IndexError(i)

    def __enter__(self) -> "SDMolSupplier":
        return self

    def __exit__(self, *args) -> None:
        pass


# ---------------------------------------------------------------------------
# ForwardSDMolSupplier — binary stream SDF reader
# ---------------------------------------------------------------------------

class ForwardSDMolSupplier:
    """Forward-only SDF iterator over a binary ``read(size)`` source.

    The source is retained and never closed. Missing/non-callable ``read``
    raises ValueError at construction; read exceptions propagate unchanged.
    ``sanitize`` and ``removeHs`` are compatibility no-ops. Text-mode streams
    are rejected on first read. No random access or restart is provided.
    """

    def __init__(self, fileobj: object, sanitize: bool = True,
                 removeHs: bool = True, strictParsing: bool = True) -> None:
        self._supplier = _ch.ForwardSDMolSupplier(
            fileobj, sanitize=sanitize, removeHs=removeHs,
            strictParsing=strictParsing,
        )

    def __iter__(self):
        return self

    def __next__(self):
        mol = next(self._supplier)
        return Mol(mol) if mol is not None else None


# ---------------------------------------------------------------------------
# SDWriter — streaming SDF writer
# ---------------------------------------------------------------------------

class SDWriter:
    """Write molecules to an SDF file, one record at a time.

    SD properties set on the molecule are written as SD data fields.
    Use ``SetProps`` to restrict which fields are written.

    Supports context-manager protocol::

        with SDWriter("out.sdf") as w:
            w.write(mol)
    """

    def __init__(self, filename: str) -> None:
        self._writer = _ch.SDWriter(filename)

    def write(self, mol: Mol, confId: int = -1) -> None:
        """Append *mol* as an SDF record including its SD properties."""
        self._writer.write(mol._mol)

    def SetProps(self, props) -> None:
        """Restrict which SD properties are written (list of key names)."""
        self._writer.SetProps(list(props))

    def SetKekulize(self, val: bool) -> None:
        """No-op: chematic does not need explicit kekulization."""
        self._writer.SetKekulize(val)

    def SetForceV3000(self, val: bool) -> None:
        """No-op: V3000 output is not yet supported."""
        self._writer.SetForceV3000(val)

    def flush(self) -> None:
        self._writer.flush()

    def close(self) -> None:
        self._writer.close()

    def __enter__(self) -> "SDWriter":
        return self

    def __exit__(self, *args) -> None:
        self.close()


# ---------------------------------------------------------------------------
# SmilesMolSupplier — read molecules from a SMILES file
# ---------------------------------------------------------------------------

class SmilesMolSupplier:
    """Read molecules from a SMILES table file (RDKit-compatible).

    Thin wrapper over :class:`chematic.SmilesMolSupplier`, the Rust-backed
    streaming reader in ``chematic-mol::smiles_table`` (a source-cited port
    of RDKit's own ``SmilesMolSupplier`` semantics) -- this class does not
    implement any parsing itself, unlike an earlier version of this wrapper
    which loaded the whole file into memory in pure Python.

    Yields ``Mol | None`` (``None`` on a malformed row), like RDKit.

    .. note::
       ``__len__``/``__getitem__`` perform an O(n) re-scan from the start
       each time (the underlying reader is forward-only streaming, matching
       ``chematic-mol``'s own documented divergence from RDKit's seek-based
       ``SmilesMolSupplier``). Avoid in hot loops.
    """

    def __init__(
        self,
        filename: str,
        delimiter: str = " ",
        smilesColumn: int = 0,
        nameColumn: int = 1,
        titleLine: bool = True,
        sanitize: bool = True,
    ) -> None:
        self._filename = filename
        self._delimiter = delimiter
        self._smiles_col = smilesColumn
        self._name_col = nameColumn
        self._title_line = titleLine
        self._sanitize = sanitize

    def __iter__(self):
        sup = _ch.SmilesMolSupplier(
            self._filename,
            delimiter=self._delimiter,
            smilesColumn=self._smiles_col,
            nameColumn=self._name_col,
            titleLine=self._title_line,
            sanitize=self._sanitize,
        )
        for mol in sup:
            yield Mol(mol) if mol is not None else None

    def __len__(self) -> int:
        return sum(1 for _ in self)

    def __getitem__(self, i: int):
        """Random access by record index (RDKit-compatible).

        .. note:: O(i) -- iterates from the start. Avoid in tight loops.
        """
        if i < 0:
            raise IndexError(i)
        for j, mol in enumerate(self):
            if j == i:
                return mol
        raise IndexError(i)


# ---------------------------------------------------------------------------
# SmilesWriter — write molecules as a SMILES/CSV file
# ---------------------------------------------------------------------------

class SmilesWriter:
    """Write molecules as a delimited SMILES table file (RDKit-compatible).

    Thin wrapper over :class:`chematic.SmilesWriter`. Columns: SMILES + name
    + any properties set via ``SetProps``. Supports the context-manager
    protocol.

    ``isomericSmiles=False`` and ``kekuleSmiles=True`` raise
    ``NotImplementedError`` -- chematic has no non-isomeric or Kekule-form
    SMILES writer at present.
    """

    def __init__(
        self,
        filename: str,
        delimiter: str = " ",
        nameHeader: str = "Name",
        includeHeader: bool = True,
        isomericSmiles: bool = True,
        kekuleSmiles: bool = False,
    ) -> None:
        self._writer = _ch.SmilesWriter(
            filename,
            delimiter=delimiter,
            nameHeader=nameHeader,
            includeHeader=includeHeader,
            isomericSmiles=isomericSmiles,
            kekuleSmiles=kekuleSmiles,
        )

    def write(self, mol: Mol) -> None:
        self._writer.write(mol._mol)

    def SetProps(self, props) -> None:
        """Restrict (and order) which properties are written as extra columns."""
        self._writer.SetProps(list(props))

    def flush(self) -> None:
        self._writer.flush()

    def close(self) -> None:
        self._writer.close()

    def __enter__(self) -> "SmilesWriter":
        return self

    def __exit__(self, *args) -> None:
        self.close()


# ---------------------------------------------------------------------------
# TDTMolSupplier / TDTWriter — Daylight Tagged Data file I/O
# ---------------------------------------------------------------------------

class TDTMolSupplier:
    """Read molecules from a Daylight TDT file (RDKit-compatible).

    Thin wrapper over :class:`chematic.TDTMolSupplier`. ``nameRecord``
    defaults to ``"NAME"`` here, not RDKit's own ``""`` default -- see
    ``chematic.TDTMolSupplier``'s doc comment for why. ``confId2D``/
    ``confId3D`` >= 0 raise ``NotImplementedError`` (chematic's Python
    ``Mol`` wrapper has no coordinate-carrying slot yet).
    """

    def __init__(
        self,
        filename: str,
        nameRecord: str = "NAME",
        confId2D: int = -1,
        confId3D: int = -1,
        sanitize: bool = True,
    ) -> None:
        self._filename = filename
        self._name_record = nameRecord
        self._conf_id_2d = confId2D
        self._conf_id_3d = confId3D
        self._sanitize = sanitize

    def __iter__(self):
        sup = _ch.TDTMolSupplier(
            self._filename,
            nameRecord=self._name_record,
            confId2D=self._conf_id_2d,
            confId3D=self._conf_id_3d,
            sanitize=self._sanitize,
        )
        for mol in sup:
            yield Mol(mol) if mol is not None else None

    def __len__(self) -> int:
        return sum(1 for _ in self)

    def __getitem__(self, i: int):
        """Random access by record index (RDKit-compatible).

        .. note:: O(i) -- iterates from the start. Avoid in tight loops.
        """
        if i < 0:
            raise IndexError(i)
        for j, mol in enumerate(self):
            if j == i:
                return mol
        raise IndexError(i)


class TDTWriter:
    """Write molecules to a Daylight TDT file (RDKit-compatible).

    Thin wrapper over :class:`chematic.TDTWriter`. Only name + properties
    round-trip through this binding at present -- see
    ``chematic.TDTWriter``'s doc comment.
    """

    def __init__(self, filename: str) -> None:
        self._writer = _ch.TDTWriter(filename)

    def write(self, mol: Mol) -> None:
        self._writer.write(mol._mol)

    def SetProps(self, props) -> None:
        """Restrict (and order) which properties are written."""
        self._writer.SetProps(list(props))

    def SetWriteNames(self, val: bool) -> None:
        self._writer.SetWriteNames(val)

    def SetWrite2D(self, val: bool) -> None:
        """No visible effect at present -- see the class doc comment."""
        self._writer.SetWrite2D(val)

    def SetNumDigits(self, n: int) -> None:
        self._writer.SetNumDigits(n)

    def NumMols(self) -> int:
        return self._writer.NumMols()

    def flush(self) -> None:
        self._writer.flush()

    def close(self) -> None:
        self._writer.close()

    def __enter__(self) -> "TDTWriter":
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
        return mol._mol.rdkit_hba

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
        useFeatures: bool = False,
        useBondTypes: bool = True,
        bitInfo=None,
        **kwargs,
    ) -> "ExplicitBitVect":
        """Return ECFP (or FCFP, with ``useFeatures=True``) fingerprint as ExplicitBitVect.

        .. note::
           Bit patterns differ from RDKit's Morgan fingerprints due to
           different hash implementations. Use for similarity ranking, not
           for cross-library bit-level comparison. ``nBits`` other than 2048
           is produced by modulo folding the internal 2048-bit fingerprint
           and is likewise not RDKit bit-exact. When ``bitInfo`` is a dict it
           is filled with ``{bit: ((atom_idx, radius), ...)}`` — shape- and
           origin-consistent with chematic's own bits, but not RDKit-identical.
           ``useFeatures=True`` uses chematic's pharmacophore feature-class
           (FCFP) invariants instead of plain atomic properties; it does not
           support ``useChirality=True`` (chirality is not tracked in the
           feature-class atom invariant).
        """
        if useFeatures and useChirality:
            raise NotImplementedError(
                "useFeatures=True with useChirality=True is not supported"
            )
        if not useBondTypes:
            raise NotImplementedError("useBondTypes=False is not supported")
        if kwargs:
            raise TypeError(f"Unsupported keyword arguments: {sorted(kwargs)}")
        if nBits <= 0:
            raise ValueError(f"nBits must be positive, got {nBits}")
        if bitInfo is not None:
            if useChirality:
                raise NotImplementedError(
                    "bitInfo with useChirality is not supported"
                )
            if useFeatures:
                raw, info = mol._mol.fcfp_bitinfo(radius)  # 2048-bit fp + dict
            else:
                raw, info = mol._mol.ecfp_bitinfo(radius)  # 2048-bit fp + dict
            bitInfo.clear()
            if nBits == 2048:
                for bit, pairs in info.items():
                    bitInfo[bit] = tuple(tuple(p) for p in pairs)
                return ExplicitBitVect._from_bytes(bytes(raw), 2048)
            # Fold both fp and bitInfo keys with the same modulo.
            for bit, pairs in info.items():
                fb = bit % nBits
                bitInfo[fb] = bitInfo.get(fb, ()) + tuple(tuple(p) for p in pairs)
            return _fold_bits(bytes(raw), nBits)
        if useChirality:
            if radius > 2:
                raise NotImplementedError(
                    "useChirality=True is only supported for radius≤2 (ECFP4)"
                )
            raw = mol._mol.ecfp4_chiral()
        elif useFeatures:
            raw = mol._mol.fcfp4() if radius <= 2 else mol._mol.fcfp6()
        elif radius <= 2:
            raw = mol._mol.ecfp4()
        else:
            raw = mol._mol.ecfp6()
        raw = bytes(raw)
        if nBits == 2048:
            return ExplicitBitVect._from_bytes(raw, 2048)
        # ponytail: modulo fold of the internal 2048-bit fp; not RDKit bit-exact
        return _fold_bits(raw, nBits)


class AllChem:
    """Common ``rdkit.Chem.AllChem`` entry points backed by chematic.

    The namespace intentionally exposes only operations with a faithful
    chematic implementation.  3D embedding and force-field routines are not
    silently approximated here; callers should use chematic's explicit 3D API.
    """

    @staticmethod
    def GetMorganFingerprintAsBitVect(mol: Mol, radius: int, nBits: int = 2048,
                                      useChirality: bool = False,
                                      useFeatures: bool = False, **kwargs):
        return rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, radius, nBits=nBits, useChirality=useChirality,
            useFeatures=useFeatures, **kwargs
        )


class _MorganGenerator:
    """Small stable-API analogue of RDKit's MorganGenerator."""

    __slots__ = ("_radius", "_nBits", "_includeChirality", "_useFeatures")

    def __init__(self, radius: int, nBits: int, includeChirality: bool,
                 useFeatures: bool):
        if radius < 0:
            raise ValueError("radius must be non-negative")
        if nBits <= 0:
            raise ValueError("fpSize must be positive")
        self._radius = radius
        self._nBits = nBits
        self._includeChirality = includeChirality
        self._useFeatures = useFeatures

    def GetFingerprint(self, mol: Mol, additionalOutput=None) -> ExplicitBitVect:
        if additionalOutput is not None:
            raise NotImplementedError("AdditionalOutput is not supported")
        return rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, self._radius, nBits=self._nBits,
            useChirality=self._includeChirality, useFeatures=self._useFeatures,
        )


class rdFingerprintGenerator:
    """Common ``rdkit.Chem.rdFingerprintGenerator`` Morgan entry point."""

    @staticmethod
    def GetMorganGenerator(radius: int = 3, countSimulation: bool = False,
                           includeChirality: bool = False,
                           useBondTypes: bool = True,
                           onlyNonzeroInvariants: bool = False,
                           includeRedundantEnvironments: bool = False,
                           fpSize: int = 2048, **kwargs):
        if countSimulation:
            raise NotImplementedError("countSimulation is not supported")
        if not useBondTypes:
            raise NotImplementedError("useBondTypes=False is not supported")
        if onlyNonzeroInvariants or includeRedundantEnvironments:
            raise NotImplementedError(
                "onlyNonzeroInvariants/includeRedundantEnvironments are not supported"
            )
        if kwargs:
            raise TypeError(f"Unsupported keyword arguments: {sorted(kwargs)}")
        return _MorganGenerator(radius, fpSize, includeChirality, False)


# ---------------------------------------------------------------------------
# pyAvalonTools namespace (mirrors rdkit.Avalon.pyAvalonTools import path)
# ---------------------------------------------------------------------------

class pyAvalonTools:
    """Mirrors ``rdkit.Avalon.pyAvalonTools``."""

    @staticmethod
    def GetAvalonFP(
        mol: Mol,
        nBits: int = 512,
        isQuery: bool = False,
        resetVect: bool = False,
        bitFlags: int = 15761407,
    ) -> "ExplicitBitVect":
        """Return RDKit's Avalon fingerprint as ExplicitBitVect.

        Bit-identical to ``rdkit.Avalon.pyAvalonTools.GetAvalonFP`` for the
        default ``bitFlags`` with ``isQuery=False`` (see
        :meth:`chematic.Mol.rdkit_avalon_fp`). As in RDKit, when ``nBits`` is
        not a multiple of 8 only the first ``nBits // 8`` bytes carry bits.
        """
        if isQuery or bitFlags != 15761407:
            raise NotImplementedError(
                "only the default bitFlags with isQuery=False are supported"
            )
        if nBits <= 0:
            raise ValueError(f"nBits must be positive, got {nBits}")
        raw = bytes(mol._mol.rdkit_avalon_fp(nBits))
        raw += bytes((nBits + 7) // 8 - len(raw))
        return ExplicitBitVect._from_bytes(raw, nBits)


# ---------------------------------------------------------------------------
# DataStructs namespace
# ---------------------------------------------------------------------------

class DataStructs:
    """Mirrors selected ``rdkit.DataStructs`` operations."""

    @staticmethod
    def FoldToTargetDensity(fp, density=0.3, minLength=64):
        """OR-fold with the corrected RDKit #9715 stopping conditions.

        Supports this adapter's ExplicitBitVect, up to 1048576 bits. The
        result is a new independent vector, including when no fold is needed.
        Finite density in [0, 1] and positive minLength are required.
        """
        if not isinstance(fp, ExplicitBitVect):
            raise TypeError("FoldToTargetDensity requires ExplicitBitVect")
        if fp.GetNumBits() > 1048576:
            raise ValueError("fold density width exceeds limit 1048576")
        result = _ch.fold_to_target_density(
            fp._to_bytes(), density, minLength, bit_width=fp.GetNumBits()
        )
        return ExplicitBitVect._from_bytes(result["fingerprint"], result["bit_width"])

    @staticmethod
    def TanimotoSimilarity(fp1, fp2) -> float:
        return _ch.tanimoto(_fp_bytes(fp1), _fp_bytes(fp2))

    @staticmethod
    def DiceSimilarity(fp1, fp2) -> float:
        return _ch.dice_similarity(_fp_bytes(fp1), _fp_bytes(fp2))

    @staticmethod
    def BulkTanimotoSimilarity(fp, fps) -> list:
        """Return Tanimoto similarity of *fp* against each fingerprint in *fps*."""
        b = _fp_bytes(fp)
        return [_ch.tanimoto(b, _fp_bytes(f)) for f in fps]

    @staticmethod
    def ConvertToNumpyArray(
        fp: _BitVectorSource, dest: np.ndarray | None = None,
    ) -> np.ndarray:
        """Replace every destination element with the corresponding binary bit.

        Sources provide GetNumBits()/GetOnBits(), including explicit/sparse
        RDKit vectors. A supplied dest must be a plain writable 1D NumPy array
        of exactly that width, with boolean/integer/float dtype and disjoint
        elements. Strided/reversed and non-native byte-order arrays are allowed.
        No resizing occurs. Validation/source errors leave dest unchanged.

        Width and raw on-bit IDs are each limited to 1048576; duplicate IDs are
        idempotent. Invalid shapes/ranges/limits raise ValueError; wrong types
        (including bool widths/IDs) raise TypeError. Source call/iterator errors
        propagate. At most one excess ID is consumed to detect the input cap.
        Arbitrary source callbacks have no timeout.

        Returns the supplied dest, or a new independent int8 array if omitted.
        Optional dest and the returned array are chematic extensions to RDKit's
        in-place-only API. Count fingerprints and whole-RDKit parity are outside
        this binary projection contract.
        """
        from ._numpy_fingerprint import _convert

        return _convert(fp, dest)
        return _convert(fp, dest)


class rdMIF:
    """Checked explicit-charge vacuum and local-solvation interaction fields."""

    Coulomb = _Coulomb
    CoulombDielectric = _CoulombDielectric
