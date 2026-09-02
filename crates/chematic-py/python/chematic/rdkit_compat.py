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

- ``GetNumAtoms()`` returns heavy-atom count regardless of ``onlyHeavy``.
  Use ``mol.add_hydrogens()`` before calling if you need H-inclusive counts.
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

__all__ = [
    "Mol", "Atom", "Bond", "BondType", "RingInfo", "RWMol",
    "ExplicitBitVect",
    "MolFromSmiles", "MolToSmiles", "MolFromMolBlock", "MolFromMolFile",
    "MolToMolBlock", "SanitizeMol", "Kekulize", "AddHs", "RemoveHs",
    "MolFromSmarts",
    "CanonSmiles", "AllChem", "rdFingerprintGenerator",
    "MolFromMrvBlock", "MolFromMrvFile", "MolToMrvBlock", "MolToMrvFile",
    "SDMolSupplier", "SDWriter", "SmilesMolSupplier", "SmilesWriter",
    "TDTMolSupplier", "TDTWriter",
    "Descriptors", "rdMolDescriptors", "DataStructs", "pyAvalonTools",
]


import chematic as _ch


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

    def GetNumAtoms(self, onlyHeavy: bool = True) -> int:
        # ponytail: always heavy-atom count; AddHs not tracked in this wrapper
        return self._mol.heavy_atoms

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

    def GetProp(self, key: str) -> str:
        return self._mol.GetProp(key)

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

    def SetIntProp(self, key: str, val: int) -> None:
        self._mol.SetIntProp(key, val)

    def SetDoubleProp(self, key: str, val: float) -> None:
        self._mol.SetDoubleProp(key, val)

    def SetBoolProp(self, key: str, val: bool) -> None:
        self._mol.SetBoolProp(key, val)

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
    def GetAvalonFP(mol: Mol, nBits: int = 512, **kwargs) -> "ExplicitBitVect":
        """Return an Avalon-style structural fingerprint as ExplicitBitVect.

        .. note::
           Bit patterns differ from RDKit's Avalon fingerprint (chematic uses
           a broad atom/bond/ring/path feature mix hashed with FNV-1a; RDKit
           uses its own C++ Avalon toolkit implementation). Use for
           similarity ranking, not cross-library bit-level comparison.
        """
        if kwargs:
            raise TypeError(f"Unsupported keyword arguments: {sorted(kwargs)}")
        if nBits <= 0:
            raise ValueError(f"nBits must be positive, got {nBits}")
        raw = bytes(mol._mol.avalon_fp())
        if nBits == 2048:
            return ExplicitBitVect._from_bytes(raw, 2048)
        # ponytail: modulo fold of the internal 2048-bit fp; not RDKit bit-exact
        return _fold_bits(raw, nBits)


# ---------------------------------------------------------------------------
# DataStructs namespace
# ---------------------------------------------------------------------------

class DataStructs:
    """Mirrors ``rdkit.DataStructs``."""

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
    def ConvertToNumpyArray(fp, dest=None):
        """Fill *dest* (or a new int8 array) with the bits of *fp*.

        Unlike RDKit's in-place-only signature, *dest* is optional; when
        omitted a new ``numpy.int8`` array of length ``fp.GetNumBits()`` is
        returned.
        """
        import numpy as np

        n = fp.GetNumBits()
        if dest is None:
            dest = np.zeros(n, dtype=np.int8)
        for i in fp.GetOnBits():
            dest[i] = 1
        return dest
