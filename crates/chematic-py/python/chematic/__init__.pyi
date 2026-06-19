"""
chematic — Pure-Rust cheminformatics for Python.

Quick start::

    import chematic

    mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
    print(mol.mw, mol.logp, mol.tpsa)
    print(mol.admet())
    fp = mol.ecfp4()
    sim = chematic.tanimoto(fp, chematic.from_smiles("c1ccccc1").ecfp4())
"""

from __future__ import annotations

from typing import Iterator, Optional

import numpy as np
from numpy import ndarray

__version__: str

# ---------------------------------------------------------------------------
# Mol
# ---------------------------------------------------------------------------

class Mol:
    """A parsed molecule.

    Create with :func:`from_smiles`, :func:`from_mol_block`, or :func:`from_inchi`.
    """

    # -- Identity ------------------------------------------------------------

    @property
    def smiles(self) -> str:
        """Canonical SMILES string."""
        ...

    @property
    def formula(self) -> str:
        """Molecular formula in Hill notation (C first, H second, then alphabetical)."""
        ...

    @property
    def heavy_atoms(self) -> int:
        """Number of heavy atoms (does not count implicit H)."""
        ...

    @property
    def inchi(self) -> str:
        """Non-standard InChI string (pure-Rust approximation). Use ``standard_inchi`` for IUPAC-compliant output."""
        ...

    @property
    def inchikey(self) -> str:
        """Non-standard InChIKey (pure-Rust approximation). Use ``standard_inchikey`` for IUPAC-compliant output."""
        ...

    @property
    def standard_inchi(self) -> str:
        """Standard IUPAC InChI string via the vendored InChI C library (v1.07.5).
        Requires the ``native-inchi`` Cargo feature. Raises ``RuntimeError`` on failure."""
        ...

    @property
    def standard_inchikey(self) -> str:
        """Standard IUPAC InChIKey (27 characters) via the vendored InChI C library.
        Requires the ``native-inchi`` Cargo feature. Raises ``RuntimeError`` on failure."""
        ...

    @property
    def iupac_name(self) -> str:
        """IUPAC systematic name. Returns an empty string for unsupported structures."""
        ...

    # -- Core physicochemical descriptors ------------------------------------

    @property
    def mw(self) -> float:
        """Average molecular weight (Da)."""
        ...

    @property
    def exact_mass(self) -> float:
        """Monoisotopic (exact) mass."""
        ...

    @property
    def logp(self) -> float:
        """Crippen–Wildman LogP."""
        ...

    @property
    def tpsa(self) -> float:
        """Topological polar surface area (Å²)."""
        ...

    @property
    def qed(self) -> float:
        """Quantitative Estimate of Drug-likeness [0, 1]."""
        ...

    @property
    def hbd(self) -> int:
        """Hydrogen bond donors."""
        ...

    @property
    def hba(self) -> int:
        """Hydrogen bond acceptors."""
        ...

    @property
    def rotatable_bonds(self) -> int:
        """Number of rotatable bonds."""
        ...

    @property
    def fsp3(self) -> float:
        """Fraction of sp3 carbons (Fsp3)."""
        ...

    @property
    def sa_score(self) -> float:
        """Synthetic Accessibility Score [1–10]; lower = easier to synthesize."""
        ...

    @property
    def molar_refractivity(self) -> float:
        """Wildman–Crippen molar refractivity."""
        ...

    @property
    def formal_charge(self) -> int:
        """Sum of formal charges."""
        ...

    @property
    def esol(self) -> float:
        """ESOL estimated aqueous solubility (log mol/L). Negative = less soluble."""
        ...

    # -- Ring / stereo counts ------------------------------------------------

    @property
    def ring_count(self) -> int:
        """Total number of rings (SSSR count)."""
        ...

    @property
    def aromatic_ring_count(self) -> int:
        """Number of aromatic rings."""
        ...

    @property
    def num_stereocenters(self) -> int:
        """Number of assigned stereocenters (R/S)."""
        ...

    # -- Drug-likeness rules -------------------------------------------------

    @property
    def lipinski_passes(self) -> bool:
        """True if Lipinski's Rule of Five passes (MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10)."""
        ...

    @property
    def veber_passes(self) -> bool:
        """True if Veber's oral bioavailability criteria pass (TPSA ≤ 140, RotBonds ≤ 10)."""
        ...

    @property
    def pains_passes(self) -> bool:
        """True if no PAINS structural alerts are present."""
        ...

    @property
    def ghose_passes(self) -> bool:
        """True if Ghose drug-likeness criteria pass (MW 160–480, LogP −0.4–5.6, MR 40–130, atoms 20–70)."""
        ...

    @property
    def egan_passes(self) -> bool:
        """True if Egan absorption criteria pass (TPSA ≤ 131.6, LogP ≤ 5.88)."""
        ...

    @property
    def reos_passes(self) -> bool:
        """True if REOS drug-likeness filter passes."""
        ...

    @property
    def brenk_passes(self) -> bool:
        """True if no Brenk structural alerts are present."""
        ...

    # -- pKa and ADMET -------------------------------------------------------

    def pka(self) -> dict[str, Optional[float]]:
        """pKa prediction.

        Returns a dict with keys ``most_acidic`` and ``most_basic``
        (float or ``None`` when no such site is found).

        Example::

            p = mol.pka()
            print(p["most_acidic"])  # e.g. 3.49
        """
        ...

    def admet(self) -> dict[str, float | bool | str]:
        """ADMET profile.

        Returns a dict with keys:
        ``bbb`` (bool), ``bbb_score`` (float),
        ``caco2`` (float, logPCaco2),
        ``herg_risk`` (float, 0–1),
        ``cyp3a4_risk`` (float, 0–1),
        ``ames_risk`` (float, 0–1),
        ``ppb`` (float, plasma protein binding %),
        ``clearance`` (str: ``"Low"`` / ``"Medium"`` / ``"High"``),
        ``gi_absorbed`` (bool), ``bbb_penetrant`` (bool).
        """
        ...

    def boiled_egg(self) -> dict[str, float | bool]:
        """Predict GI absorption and BBB penetration (Daina & Zoete 2016).

        Uses LogP and TPSA thresholds:

        - **GI absorbed** (egg-white): LogP ≤ 5.88 and TPSA ≤ 131.6 Å²
        - **BBB penetrant** (egg-yolk): LogP ∈ [−0.3, 6.1] and TPSA ≤ 71.1 Å²

        Returns a dict with keys:
        ``gi_absorbed`` (bool), ``bbb_penetrant`` (bool),
        ``logp`` (float), ``tpsa`` (float).

        Example::

            egg = mol.boiled_egg()
            print(egg["gi_absorbed"])    # True
            print(egg["bbb_penetrant"])  # False
        """
        ...

    # -- All descriptors in one call ----------------------------------------

    def descriptors(self) -> dict[str, float | int | bool | None]:
        """Return all scalar descriptors as a dict (70+ keys).

        Useful for building a Pandas DataFrame row::

            import pandas as pd
            df = pd.DataFrame([chematic.from_smiles(s).descriptors() for s in smiles_list])
        """
        ...

    # -- Fingerprints --------------------------------------------------------

    def ecfp4(self) -> bytes:
        """ECFP4 fingerprint as bytes (256 bytes = 2048 bits, LSB-first)."""
        ...

    def ecfp6(self) -> bytes:
        """ECFP6 fingerprint as bytes (256 bytes = 2048 bits, LSB-first)."""
        ...

    def ecfp4_chiral(self) -> bytes:
        """ECFP4 with chirality fingerprint as bytes (256 bytes = 2048 bits)."""
        ...

    def fcfp4(self) -> bytes:
        """FCFP4 (functional-class ECFP4) fingerprint as bytes (256 bytes = 2048 bits)."""
        ...

    def atom_pair_fp(self) -> bytes:
        """Atom-pair fingerprint as bytes (256 bytes = 2048 bits, LSB-first)."""
        ...

    def torsion_fp(self) -> bytes:
        """Topological torsion fingerprint as bytes (256 bytes = 2048 bits, LSB-first)."""
        ...

    def maccs(self) -> bytes:
        """MACCS 166-bit keys as bytes (21 bytes, LSB-first)."""
        ...

    def ecfp4_numpy(self) -> ndarray:
        """ECFP4 fingerprint as a numpy array of shape ``(2048,)``, dtype ``uint8`` (0/1)."""
        ...

    def maccs_numpy(self) -> ndarray:
        """MACCS 166-bit fingerprint as a numpy array of shape ``(166,)``, dtype ``uint8`` (0/1)."""
        ...

    # -- Visualization -------------------------------------------------------

    def svg(self) -> str:
        """2D SVG depiction as a string.

        In Jupyter notebooks::

            from IPython.display import SVG
            SVG(mol.svg())
        """
        ...

    def svg_highlighted(
        self,
        atom_indices: list[int],
        color: str = "#FFFF00",
    ) -> str:
        """2D SVG depiction with highlighted atoms.

        Args:
            atom_indices: Zero-based atom indices to highlight.
            color: CSS color string (default ``"#FFFF00"`` yellow).

        Example::

            svg = mol.svg_highlighted([0, 1, 2])
            svg = mol.svg_highlighted([0], color="#FF0000")
        """
        ...

    # -- Transformations -----------------------------------------------------

    def standardize(self) -> Mol:
        """Return the standardized molecule (largest fragment, charges neutralized, tautomer canonicalized)."""
        ...

    def scaffold(self) -> Mol:
        """Return the Bemis–Murcko scaffold as a new Mol."""
        ...

    def generic_scaffold(self) -> Mol:
        """Return the generic Murcko scaffold (all atoms replaced with carbons, all bonds single)."""
        ...

    def canonical_tautomer(self) -> Mol:
        """Return the canonical tautomer as a new Mol."""
        ...

    def enumerate_tautomers(self) -> list[Mol]:
        """Return all tautomers as a list of Mol objects."""
        ...

    def enumerate_stereoisomers(self) -> list[Mol]:
        """Return all stereoisomers as a list of Mol objects.

        Returns an empty list when the molecule has more than 6 unspecified
        stereocenters (combinatorial explosion guard). Use ``mol.num_stereocenters``
        to distinguish this from having no centers.
        """
        ...

    def add_hydrogens(self) -> Mol:
        """Return a copy with all implicit hydrogens made explicit."""
        ...

    def remove_hydrogens(self) -> Mol:
        """Return a copy with all explicit hydrogen atoms removed."""
        ...

    def remove_stereo(self) -> Mol:
        """Return a copy with all stereochemistry assignments removed."""
        ...

    def remove_isotopes(self) -> Mol:
        """Return a copy with all isotope labels removed."""
        ...

    def largest_fragment(self) -> Mol:
        """Return the largest covalently connected fragment."""
        ...

    def neutralize(self) -> Mol:
        """Return a charge-neutralized copy."""
        ...

    def brics_fragments(self) -> list[Mol]:
        """Fragment the molecule using BRICS rules.

        When no BRICS-breakable bonds are found, returns a list containing the
        original molecule (not an empty list).
        """
        ...

    # -- Dunder methods ------------------------------------------------------

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

# ---------------------------------------------------------------------------
# Module-level functions
# ---------------------------------------------------------------------------

def from_smiles(smiles: str) -> Mol:
    """Parse a SMILES string and return a :class:`Mol`.

    Raises:
        ValueError: on invalid SMILES.

    Example::

        mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    """
    ...

def from_mol_block(block: str) -> Mol:
    """Parse a MOL/SDF block and return a :class:`Mol`.

    Raises:
        ValueError: on parse failure.
    """
    ...

def from_inchi(inchi: str) -> Mol:
    """Parse an InChI string and return a :class:`Mol`.

    Raises:
        ValueError: on parse failure.

    Example::

        mol = chematic.from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
    """
    ...

def is_valid_smiles(smiles: str) -> bool:
    """Return True if the SMILES can be parsed without error."""
    ...

def tanimoto(a: bytes, b: bytes) -> float:
    """Tanimoto similarity between two fingerprint byte arrays.

    Works with any equal-length ``bytes`` objects (ECFP4, ECFP6, MACCS, …).
    Returns a value in [0.0, 1.0].

    Raises:
        ValueError: if ``a`` and ``b`` have different lengths.

    Example::

        sim = chematic.tanimoto(mol1.ecfp4(), mol2.ecfp4())
    """
    ...

def smarts_match(smarts: str, mol: Mol) -> bool:
    """Test whether a SMARTS pattern matches a molecule.

    Raises:
        ValueError: on invalid SMARTS.

    Example::

        if chematic.smarts_match("[OH]", mol):
            print("has hydroxyl")
    """
    ...

def smarts_find(smarts: str, mol: Mol) -> list[list[int]]:
    """Return all substructure matches of a SMARTS pattern in a molecule.

    Each inner list contains atom indices (in query-atom order).
    Returns an empty list when there are no matches.

    Raises:
        ValueError: on invalid SMARTS.

    Example::

        matches = chematic.smarts_find("[OH]", mol)
        # → [[3], [7], ...]   (one list per match; each element is a mol atom index)
    """
    ...

def depict_grid(mols: list[Mol], cols: int) -> str:
    """Render a list of molecules as a grid SVG.

    Example::

        svg = chematic.depict_grid([mol1, mol2, mol3], cols=3)
    """
    ...

def run_smirks(smirks: str, reactants: list[Mol]) -> list[list[Mol]]:
    """Apply a SMIRKS reaction template to a list of reactant molecules.

    Returns a list of product sets; each product set is a list of :class:`Mol`.

    Raises:
        ValueError: on SMIRKS parse failure or reactant count mismatch.

    Example::

        products = chematic.run_smirks("[OH:1]>>[O-:1]", [mol])
        # → [[product_mol], ...]
    """
    ...

def find_mcs(mols: list[Mol]) -> Optional[Mol]:
    """Find the Maximum Common Substructure (MCS) of a list of molecules.

    Returns the MCS as a :class:`Mol`, or ``None`` when there is no common substructure.

    Example::

        mcs = chematic.find_mcs([mol1, mol2])
        if mcs:
            print(mcs.smiles)
    """
    ...

# ---------------------------------------------------------------------------
# SimilarityIndex (MHFP LSH)
# ---------------------------------------------------------------------------

class SimilarityIndex:
    """Locality-sensitive hashing index for fast approximate similarity search.

    Uses MinHash fingerprints (MHFP). Much faster than brute-force Tanimoto
    for large libraries (10k+ molecules).

    Example::

        idx = chematic.SimilarityIndex.from_smiles(smiles_list)
        hits = idx.search("c1ccccc1", threshold=0.7, k=10)
        # hits → [(index, score), ...]
    """

    def __new__(cls, num_hashes: int = 128) -> SimilarityIndex:
        """Create an empty index.

        Args:
            num_hashes: Number of hash functions (default 128). Higher = more accurate.
        """
        ...

    @staticmethod
    def from_smiles(smiles: list[str]) -> SimilarityIndex:
        """Build an index from a list of SMILES strings.

        Invalid SMILES are silently skipped.
        """
        ...

    def add(self, smiles: str) -> int:
        """Add a molecule and return its index in the library."""
        ...

    def search(
        self,
        query: str,
        threshold: float = 0.7,
        k: Optional[int] = None,
    ) -> list[tuple[int, float]]:
        """Search for similar molecules.

        Args:
            query: SMILES string of the query molecule.
            threshold: Minimum similarity (default 0.7).
            k: Maximum number of results (``None`` = all above threshold).

        Returns:
            List of ``(index, similarity)`` tuples, sorted by descending similarity.
        """
        ...

    def get_smiles(self, index: int) -> str:
        """Return the SMILES stored at the given index."""
        ...

    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...

# ---------------------------------------------------------------------------
# SDF I/O
# ---------------------------------------------------------------------------

class SdfRecord:
    """One record from an SDF file."""

    mol: Mol
    name: str

    @property
    def smiles(self) -> str:
        """Canonical SMILES of the molecule."""
        ...

    def properties(self) -> dict[str, str]:
        """Return all SD data fields as a dict."""
        ...

    def get(self, key: str) -> Optional[str]:
        """Get one SD property by name, or ``None`` if absent."""
        ...

    def __repr__(self) -> str: ...

class SdfIter:
    """Lazy iterator over SDF records."""

    def __iter__(self) -> Iterator[SdfRecord]: ...
    def __next__(self) -> SdfRecord: ...
    def __len__(self) -> int: ...

def iter_sdf(path: str) -> SdfIter:
    """Lazily iterate over records in an SDF file.

    Example::

        for record in chematic.iter_sdf("library.sdf"):
            print(record.name, record.mol.mw)
    """
    ...

def iter_sdf_str(content: str) -> SdfIter:
    """Lazily iterate over records in an SDF string (already loaded into memory)."""
    ...

# ---------------------------------------------------------------------------
# bulk submodule
# ---------------------------------------------------------------------------

class bulk:
    """Parallel batch operations for processing large molecule libraries.

    All functions use Rayon for multi-core parallelism.
    """

    @staticmethod
    def parse(smiles: list[str]) -> list[Optional[Mol]]:
        """Parse a list of SMILES in parallel. Returns ``None`` for invalid entries."""
        ...

    @staticmethod
    def ecfp4(smiles: list[str]) -> ndarray:
        """Compute ECFP4 fingerprints in parallel.

        Returns:
            numpy array of shape ``(N, 2048)``, dtype ``uint8``.
            Invalid SMILES produce all-zero rows.
        """
        ...

    @staticmethod
    def maccs(smiles: list[str]) -> ndarray:
        """Compute MACCS 166-bit fingerprints in parallel.

        Returns:
            numpy array of shape ``(N, 166)``, dtype ``uint8``.
        """
        ...

    @staticmethod
    def descriptors(smiles: list[str]) -> list[dict[str, float | int | bool | None]]:
        """Compute 55+ descriptors for each SMILES in parallel.

        Returns:
            List of descriptor dicts (one per molecule). Failed parses produce empty dicts.
        """
        ...

    @staticmethod
    def tanimoto(smiles_a: list[str], smiles_b: list[str]) -> ndarray:
        """Compute pairwise ECFP4 Tanimoto similarity matrix.

        Returns:
            numpy array of shape ``(M, N)``, dtype ``float32``.
        """
        ...

    @staticmethod
    def tanimoto_search(query: str, smiles: list[str]) -> ndarray:
        """Compute ECFP4 Tanimoto similarity of one query against a library.

        Returns:
            numpy array of shape ``(N,)``, dtype ``float32``.
        """
        ...
