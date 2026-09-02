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

from typing import Any, Iterable, Iterator, Optional

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

    def iupac_name_stereo(self) -> str:
        """IUPAC name with CIP stereochemistry prefix (e.g. ``"(R)-butan-2-ol"``).

        Returns an empty string for structures outside the IUPAC naming scope.
        """
        ...

    def random_smiles(self, seed: int) -> str:
        """Generate a random (non-canonical) SMILES string.

        ``seed`` controls the atom traversal order deterministically.
        Useful for ML data augmentation (SMILES enumeration).

        Example::

            smi = mol.random_smiles(42)   # e.g. "OCC" for ethanol
        """
        ...

    def random_smiles_n(self, n: int, seed: int = 0) -> list[str]:
        """Generate up to ``n`` unique random SMILES strings.

        Example::

            variants = mol.random_smiles_n(n=5, seed=42)
        """
        ...

    def are_atoms_equivalent(self, a: int, b: int) -> bool:
        """Return ``True`` if atoms ``a`` and ``b`` are symmetry-equivalent (same Morgan rank).

        Example::

            assert mol.are_atoms_equivalent(0, 1)  # in benzene all C are equivalent
        """
        ...

    @property
    def recap_breakable_bond_count(self) -> int:
        """Number of RECAP-breakable bonds (C–N, C–O, C–S single bonds)."""
        ...

    def to_pdbqt(
        self,
        coords: list[tuple[float, float, float]],
        charges: list[float],
        name: str = "LIG",
    ) -> str:
        """Write molecule to AutoDock PDBQT format (rigid body, no torsion tree).

        Args:
            coords: list of (x, y, z) tuples in Å, one per heavy atom.
            charges: partial charges (float) per heavy atom. Use
                     ``chematic_ff.gasteiger_charges()`` or MMFF94 BCI charges
                     for best docking accuracy.
            name: ligand name in the REMARK header (e.g. ``"LIG"``).

        Returns:
            str: PDBQT-format string.
        """
        ...

    def minimize_uff(
        self,
        coords: list[list[float]],
        max_iter: int = 500,
    ) -> dict[str, object]:
        """Minimise geometry using the Universal Force Field (UFF, Rappé 1992).

        UFF covers all elements including metals, unlike MMFF94 which is limited
        to common organic/heteroatoms.

        Args:
            coords: list of ``[x, y, z]`` lists (Å) — initial 3D coordinates.
            max_iter: maximum steepest-descent iterations.

        Returns:
            dict with keys ``coords`` (list[list[float]]), ``energy`` (float,
            kcal/mol), ``iterations`` (int), ``converged`` (bool), ``sound``
            (bool — all-finite coordinates and no bond stretched past a sane
            covalent-bond length; independent of ``converged``, check this
            before trusting a result).
        """
        ...

    def to_mol_block(self) -> str:
        """Serialize to MDL MOL V2000 format (without 3D coordinates).

        Equivalent to RDKit's ``Chem.MolToMolBlock(mol)``.
        Use :meth:`to_mol2` for Tripos format or :meth:`to_pdb` for PDB with 3D.

        Example::

            block = mol.to_mol_block()
            m2 = chematic.from_mol_block(block)  # round-trip
        """
        ...

    def to_mol_block_2d(
        self,
        coords: list[list[float]],
        name: Optional[str] = None,
    ) -> str:
        """Serialize to MDL MOL V2000 format preserving 2D layout coordinates.

        Each element of ``coords`` is an ``[x, y]`` pair in Å.
        Designed for round-tripping with :func:`from_mol_block_with_coords`.

        Example::

            mol, name, coords_2d = chematic.from_mol_block_with_coords(block)
            new_block = mol.to_mol_block_2d(coords_2d, name=name)
        """
        ...

    def to_mol_v3000(
        self,
        coords: list[list[float]],
        name: Optional[str] = None,
    ) -> str:
        """Serialize to MDL MOL V3000 format with 2D layout coordinates.

        V3000 supports >999 atoms and extended atom/bond features.
        Accepts the same ``[[x, y], ...]`` coordinate format as :meth:`to_mol_block_2d`.
        Pass an empty list for zero coordinates.

        Equivalent to RDKit ``Chem.MolToV3KMolBlock(mol)``.

        Example::

            block = mol.to_mol_v3000(coords_2d, name="my_mol")
        """
        ...

    def to_cml(self, coords: Optional[list[list[float]]] = None) -> str:
        """Serialize this molecule to Chemical Markup Language (CML) XML.

        ``coords``: optional list of ``[x, y]`` pairs (Å) — one per heavy atom.
        If omitted or ``None``, no coordinate attributes are written.

        Equivalent to RDKit ``Chem.MolToCMLBlock(mol)``.

        Example::

            cml = mol.to_cml()
            cml_with_layout = mol.to_cml(coords_2d)
        """
        ...

    def to_mol2(self) -> str:
        """Serialize to Tripos MOL2 format string (SYBYL MOL2).

        The output contains the mandatory @<TRIPOS>MOLECULE, @<TRIPOS>ATOM,
        and @<TRIPOS>BOND sections.  Atom coordinates are all zero when no
        3D geometry is available; use :func:`chematic.generate_3d` first if
        you need 3D coordinates.

        Returns:
            str: A complete MOL2-format string.

        Example::

            mol = chematic.from_smiles("CCO")
            with open("ethanol.mol2", "w") as f:
                f.write(mol.to_mol2())

            # With 3D coordinates
            mol3d = chematic.generate_3d("CCO")   # returns Mol with coords
            with open("ethanol_3d.mol2", "w") as f:
                f.write(mol3d.to_mol2())
        """
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
    def ring_system_count(self) -> int:
        """Number of distinct connected ring systems.

        Two SSSR rings form the same system when they share at least one atom
        (fused, bridged, or spiro).  Differs from :attr:`ring_count` which
        counts SSSR rings individually.

        Example::

            naphthalene = chematic.from_smiles("c1ccc2ccccc2c1")
            naphthalene.ring_system_count  # 1 (two fused rings = one system)
            biphenyl = chematic.from_smiles("c1ccc(-c2ccccc2)cc1")
            biphenyl.ring_system_count     # 2 (two independent benzene rings)
        """
        ...

    @property
    def hba_count_lipinski(self) -> int:
        """Lipinski (1997) HBA count — total number of N and O heavy atoms.

        The original Rule-of-Five HBA definition: count all N and O atoms
        regardless of hybridisation or substitution.  For the chemically more
        accurate Ertl (2000) definition use :attr:`hba`.

        Example::

            caffeine = chematic.from_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
            caffeine.hba_count_lipinski  # 5 (2 O + 3 N)
        """
        ...

    @property
    def fraction_rotatable_bonds(self) -> float:
        """Fraction of heavy atoms that are rotatable bonds (0.0–1.0).

        ``fraction_rotatable_bonds = rotatable_bond_count / heavy_atom_count``.
        Returns 0.0 for rigid or acyclic molecules with no rotatable bonds.

        Example::

            benzene = chematic.from_smiles("c1ccccc1")
            benzene.fraction_rotatable_bonds  # 0.0
        """
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

    def mcf_passes(self) -> bool:
        """Med-Chem Friendly (MCF) composite filter.

        Returns ``True`` when all of the following hold:
        no PAINS alerts, no Brenk alerts, Lipinski Ro5, and Veber oral bioavailability.
        Mirrors the "MCF" concept in the medchemfilters Python library.

        Example::

            caffeine = chematic.from_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
            caffeine.mcf_passes  # True
        """
        ...

    def ro3_passes(self) -> bool:
        """True if Rule of Three criteria pass (Congreve 2003).

        MW ≤ 300, LogP ≤ 3, HBD ≤ 3, HBA ≤ 3, RotBonds ≤ 3.
        Used for fragment-based drug discovery (FBDD) library screening.
        """
        ...

    def lead_like_passes(self) -> bool:
        """True if lead-like criteria pass (Oprea 2001).

        MW ≤ 450, LogP −3.5–4.5, RotBonds ≤ 10, RingCount 1–4.
        Lead-like compounds have lower MW/LogP than drugs, leaving room
        for optimisation-related property increases.
        """
        ...

    def pfizer_3_75_passes(self) -> bool:
        """True if compound is NOT in the Pfizer 3/75 high-metabolic-liability zone.

        The danger zone is ``LogP > 3 AND TPSA < 75``; compounds there have
        higher CYP3A4 metabolic clearance risk (Leeson & Springthorpe 2007).
        Returns ``True`` = safe (not in danger zone).

        Example::

            ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
            ibuprofen.pfizer_3_75_passes  # False (LogP≈3.8, TPSA≈37)
        """
        ...

    def cns_mpo_score(self) -> float:
        """CNS Multi-Parameter Optimisation (MPO) score (Wager 2010), range 0–6.

        Combines desirability functions for cLogP, cLogD (pH 7.4), MW, TPSA,
        HBD, and pKa (most basic site). Higher scores indicate better CNS
        drug-like properties. Scores ≥ 4 are generally considered CNS-appropriate.

        Example::

            mol = chematic.from_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C")  # caffeine
            mol.cns_mpo_score  # ≥ 3.0 (small, low HBD)
        """
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

    def layered_fp_layers(self) -> list[bytes]:
        """Layered fingerprint decomposed into 7 individual layers.

        Each layer is a 2048-bit (256-byte) fingerprint encoding progressively
        more structural detail:

        - Layer 0: raw atom types (element, H count, charge)
        - Layer 1: + bond orders
        - Layer 2: + aromaticity
        - Layer 3: + ring membership
        - Layer 4: + is-ring-bond
        - Layer 5: + stereochemistry
        - Layer 6: all features combined

        Returns a list of 7 ``bytes`` objects compatible with :func:`tanimoto`,
        :func:`dice_similarity`, etc.
        Equivalent to RDKit's ``Chem.LayeredFingerprint(mol, layerFlags=0x7F)``.

        Example::

            layers = mol.layered_fp_layers()
            sim = chematic.tanimoto(layers[3], other.layered_fp_layers()[3])
        """
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

        In Jupyter notebooks, use :meth:`_repr_svg_` (automatic) or::

            from IPython.display import SVG
            SVG(mol.svg())
        """
        ...

    def depict_data(self) -> dict:
        """Structured 2D depiction data (atoms + bonds with layout coordinates).

        Use this instead of :meth:`svg` when you want to drive your own
        renderer (e.g. matplotlib, a custom canvas) rather than parse SVG.

        Returns:
            dict with keys:
            ``atoms`` (list of dicts: ``idx``, ``element`` (symbol string),
            ``x``, ``y``, ``label`` (``None`` when suppressed), ``color``
            (CSS hex string), ``charge``) and ``bonds`` (list of dicts:
            ``idx``, ``atom1``, ``atom2``, ``kind`` — one of ``"Single"``,
            ``"Double"``, ``"Triple"``, ``"Aromatic"``, ``"Up"``, ``"Down"``).

        Example::

            mol = chematic.from_smiles("CCO")
            data = mol.depict_data()
            for atom in data["atoms"]:
                print(atom["element"], atom["x"], atom["y"])
        """
        ...

    def _repr_svg_(self) -> str:
        """Jupyter auto-display hook — renders the molecule automatically in a cell.

        Just write ``mol`` in a Jupyter cell and the 2D structure appears.
        No ``IPython.display.SVG(...)`` wrapper needed.
        """
        ...

    def has_substructure(self, smarts: str) -> bool:
        """Return True if this molecule matches the SMARTS pattern.

        Raises ``ValueError`` for invalid SMARTS.

        Example::

            mol = chematic.from_smiles("CC(=O)O")
            mol.has_substructure("[OH]")   # True
        """
        ...

    def find_matches(self, smarts: str) -> list[list[int]]:
        """Return atom-index lists for all SMARTS matches in this molecule.

        Each inner list contains sorted atom indices for one match.
        Returns an empty list when there are no matches.
        Raises ``ValueError`` for invalid SMARTS.

        Example::

            mol = chematic.from_smiles("OCC(=O)O")
            mol.find_matches("[OH]")          # [[0], [4]]
            mol.find_matches("[CX3](=O)[OH]") # [[1, 2, 4]]
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

    def connected_components(self) -> list[Mol]:
        """Split this molecule into its connected components (fragments).

        Returns one :class:`Mol` per connected component.
        A fully connected molecule returns a single-element list.

        Equivalent to RDKit ``Chem.GetMolFrags(mol, asMols=True)``.

        Example::

            mol = chematic.from_smiles("CC.[NH3]")
            parts = mol.connected_components()   # [Mol("CC"), Mol("N")]
            assert len(parts) == 2
        """
        ...

    def is_same_as(self, other: Mol) -> bool:
        """Return True if this molecule and other represent the same chemical structure.

        Uses canonical SMILES comparison (reliable after fix for issue #14).
        Equivalent to :func:`chematic.are_identical`.

        Example::

            m1 = chematic.from_smiles("CC(=O)O")
            m2 = chematic.from_smiles("OC(C)=O")
            assert m1.is_same_as(m2)    # True — same acetic acid
            assert not m1.is_same_as(chematic.from_smiles("CCO"))
        """
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

    # -- Extended fingerprints -----------------------------------------------

    def map4(self) -> list[int]:
        """MAP4 fingerprint (Minervini 2020) — 1024 u32 hash values.

        Use :func:`tanimoto_map4` for similarity (not :func:`tanimoto`).

        Example::

            sim = chematic.tanimoto_map4(mol1.map4(), mol2.map4())
        """
        ...

    def map4_numpy(self) -> ndarray:
        """MAP4 fingerprint as a numpy array of shape ``(1024,)`` dtype ``uint32``."""
        ...

    def erg(self) -> bytes:
        """Extended Reduced Graph (ERG) fingerprint as bytes (256 bytes = 2048 bits)."""
        ...

    # -- Extended descriptors ------------------------------------------------

    def logd(self, ph: float = 7.4) -> float:
        """LogD at a given pH — accounts for ionization of acids/bases.

        Default pH is 7.4 (physiological). More relevant than LogP for ADMET.
        """
        ...

    def logd_profile(self) -> list[tuple[float, float]]:
        """LogD profile — list of ``(pH, LogD)`` pairs from pH 0 to 14 (28 steps)."""
        ...

    def mqn(self) -> list[int]:
        """Molecular Quantum Numbers — 42-element topological descriptor vector.

        Returns a list of 42 integers. Reference: Ertl et al., *J. Chem. Inf. Model.* 2009.
        """
        ...

    def logp_per_atom(self) -> list[float]:
        """Per-atom Crippen LogP contributions — one float per heavy atom, in atom order."""
        ...

    def tpsa_per_atom(self) -> list[float]:
        """Per-atom TPSA contributions (Ertl 2000) — one float per heavy atom, in atom order.

        Only N, O, S, P atoms have non-zero contributions. ``sum(mol.tpsa_per_atom()) == mol.tpsa``.

        Example::

            aspirin = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
            ta = aspirin.tpsa_per_atom()
            assert len(ta) == aspirin.heavy_atoms
        """
        ...

    def logp_map_svg(self) -> str:
        """2D SVG with atoms coloured by their Crippen LogP contribution.

        Positive (lipophilic) atoms → blue; negative (hydrophilic) atoms → red;
        zero-contribution atoms → white (no tint).

        Example::

            svg = mol.logp_map_svg()
            with open("logp_map.svg", "w") as f:
                f.write(svg)
        """
        ...

    def tpsa_map_svg(self) -> str:
        """2D SVG with atoms coloured by their TPSA contribution.

        TPSA-contributing atoms (N, O, S, P) → blue; carbon and halogens → white.

        Example::

            svg = mol.tpsa_map_svg()
        """
        ...

    def similarity_map_svg(self, weights: list[float]) -> str:
        """2D SVG with atoms coloured by custom per-atom weights.

        ``weights``: one float per heavy atom. Positive → blue, negative → red, zero → white.
        Weights are normalised to the maximum absolute value before colour mapping.

        Example::

            weights = mol.logp_per_atom()
            svg = mol.similarity_map_svg(weights)
        """
        ...

    def isotope_distribution(self) -> list[tuple[float, float]]:
        """Isotopic distribution — list of ``(mass, relative_intensity)`` pairs.

        The highest-intensity peak is normalised to 1.0.
        """
        ...

    # -- Chemical analysis ---------------------------------------------------

    def functional_groups(self) -> list[dict]:
        """Identify functional groups (Ertl 2017 algorithm).

        Returns a list of dicts, each with:
          ``atom_indices`` (list of int), ``atom_types`` (str, e.g. ``"N,O"``).
        """
        ...

    def scaffold_network(self) -> list[str]:
        """Schuffenhauer scaffold parents — list of SMILES (outermost scaffold first)."""
        ...

    # -- 3D generation -------------------------------------------------------

    def generate_3d(self) -> list[list[float]]:
        """Generate 3D coordinates (distance geometry + DREIDING minimization).

        Returns a list of ``[x, y, z]`` lists (Å), one per heavy atom.
        Use the returned coords with :meth:`whim`, :meth:`getaway`,
        :meth:`mmff94_energy_breakdown`, :meth:`to_pdb`, etc.

        Example::

            coords = mol.generate_3d()
            pdb = mol.to_pdb(coords)
        """
        ...

    def conformer_ensemble(
        self,
        n: int,
        rmsd_threshold: float = 0.5,
    ) -> list[list[list[float]]]:
        """Generate multiple conformers with RMSD-based pruning.

        Returns a list of coordinate arrays — each is a ``[[x,y,z], ...]`` list.

        Args:
            n: Number of conformers to attempt.
            rmsd_threshold: Minimum RMSD (Å) between retained conformers.
        """
        ...

    # -- 3D descriptors ------------------------------------------------------

    def whim(self, coords: list[list[float]]) -> list[float]:
        """WHIM 3D descriptors (Todeschini & Gramatica 1997).

        ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
        Returns a flat list of floats (shape/symmetry descriptors).
        """
        ...

    def getaway(self, coords: list[list[float]]) -> list[float]:
        """GETAWAY 3D descriptors (Consonni et al. 2002).

        ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
        """
        ...

    def autocorr_3d(self, coords: list[list[float]]) -> list[float]:
        """3D autocorrelation descriptors.

        ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
        """
        ...

    def spectrophores(
        self,
        coords: list[list[float]],
        normalize: str = "none",
    ) -> list[float]:
        """Spectrophores 3D fingerprint — 48-element vector.

        Encodes the 3D electrostatic, lipophilic, aromatic, and H-bond
        character of the molecule's surface.  Requires 3D coordinates.

        Args:
            coords: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
            normalize: ``"none"`` (default), ``"zscore"``, or ``"l2"``.

        Returns:
            List of 48 floats (4 properties × 12 probe positions).

        Example::

            coords = mol.generate_3d()
            fp = mol.spectrophores(coords)
            sim = chematic.tanimoto_spectrophores(fp1, fp2)

        Reference: Silicos-it Spectrophores (patent expired 2024).
        """
        ...

    # -- 3D file I/O ---------------------------------------------------------

    def to_pdb(self, coords: list[list[float]]) -> str:
        """Write this molecule's 3D structure to PDB format.

        ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
        """
        ...

    def to_xyz(self, coords: list[list[float]], comment: str = "") -> str:
        """Write this molecule's 3D structure to XYZ format.

        ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
        """
        ...

    def sasa_per_atom_3d(self, coords: list[list[float]]) -> list[float]:
        """Per-atom Solvent-Accessible Surface Area (Å²) from explicit 3D coordinates.

        Uses the Shrake-Rupley algorithm (probe 1.4 Å, 100 sphere points/atom).
        Unlike :meth:`sasa_per_atom` (which generates DG coords internally), this method
        accepts coordinates from :meth:`generate_3d`, :meth:`minimize_mmff94`, etc.

        Example::

            coords = mol.generate_3d()
            per_atom = mol.sasa_per_atom_3d(coords)
            total = sum(per_atom)
        """
        ...

    # -- Force field analysis ------------------------------------------------

    def mmff94_total_energy(self, coords: list[list[float]]) -> float:
        """Total MMFF94 force field energy in kcal/mol for the given 3D coordinates.

        Returns ``0.0`` if MMFF94 typing fails. Complements :meth:`mmff94_energy_breakdown`.

        Example::

            coords = mol.generate_3d()
            e = mol.mmff94_total_energy(coords)
        """
        ...

    def mmff94_atom_types(self) -> list[str]:
        """Per-atom MMFF94 force field type names.

        Returns one string per heavy atom (e.g. ``"C_sp3"``, ``"C=O"``, ``"N_amide"``).
        Returns an empty list if MMFF94 typing is not supported for this molecule.
        Equivalent to RDKit ``AllChem.MMFFGetMoleculeProperties(mol).GetMMFFAtomType(i)``.
        """
        ...

    def mmff94_charges_3d(self, coords: list[list[float]]) -> list[float]:
        """MMFF94 partial charges incorporating 3D polarization effects.

        Requires 3D coordinates (from :meth:`generate_3d`).
        More accurate than :meth:`mmff94_charges` (2D topology only) for polar molecules.
        Returns one float per heavy atom.

        Raises:
            ValueError: if MMFF94 typing fails or coords length mismatches atom count.

        Example::

            coords = mol.generate_3d()
            charges = mol.mmff94_charges_3d(coords)
        """
        ...

    def mmff94_energy_breakdown(
        self, coords: list[list[float]]
    ) -> dict[str, float]:
        """MMFF94 energy breakdown for given 3D coordinates.

        Returns a dict with keys: ``bond``, ``angle``, ``stretch_bend``,
        ``torsion``, ``oop``, ``vdw``, ``electrostatic``, ``total`` (kcal/mol).

        Raises:
            ValueError: for atoms not parameterised by MMFF94.
        """
        ...

    def embed_pipeline_v2(self, config: "PipelineV2Config") -> dict[str, object]:
        """Run the opt-in v2 embedding pipeline.

        Torsion-knowledge-aware distance geometry + stereo verification/repair
        + policy-gated force field, applied directly to this ``Mol``'s own
        atom order (never canonicalizes/reparses first -- coordinates and
        every atom/bond index in the result correspond 1:1 to this ``Mol``'s
        existing atom/bond tables).

        Returns a dict with keys: ``coords``, ``embed_stats``,
        ``bound_adjustment_report``, ``torsion_knowledge_report``,
        ``ring_torsion_evidence``, ``torsion_optimization_report``,
        ``stereo_before``, ``stereo_repair``, ``stereo_after_repair``,
        ``force_field``, ``final_stereo``, ``final_validation``,
        ``elapsed_ms_by_stage``.

        Raises:
            PipelineV2Error: on pipeline failure. ``error.diagnostics`` carries
                the same per-stage partial evidence as a Rust
                ``PipelineV2Failure`` -- ``diagnostics["last_known_coords"]``
                is diagnostic only, never a usable result.
        """
        ...

    def mmff94_torsion_scan(
        self,
        coords: list[list[float]],
        atom_i: int,
        atom_j: int,
        atom_k: int,
        atom_l: int,
        steps: int = 36,
    ) -> list[tuple[float, float]]:
        """Scan a torsion dihedral over 360°.

        Returns a list of ``(angle_deg, energy_kcal)`` pairs.
        The molecule is not modified.
        """
        ...

    # -- Sprint 2: additional descriptors & fingerprints --------------------

    @property
    def xlogp3(self) -> float:
        """XLogP3 — alternative logP (Wang et al. 2000)."""
        ...

    def xlogp3_per_atom(self) -> list[float]:
        """Per-atom XLogP3 contributions — one float per heavy atom."""
        ...

    def autocorr_2d(self) -> list[float]:
        """2D autocorrelation descriptors (mordred/Dragon compatible)."""
        ...

    @property
    def hall_kier_alpha(self) -> float:
        """Hall-Kier alpha — correction term for kappa shape indices."""
        ...

    @property
    def zagreb_m2(self) -> int:
        """Second Zagreb index M2 — Σ(deg(a) × deg(b)) over all heavy-atom bonds.

        Complements :attr:`zagreb_m1` (Σ deg(v)²). Both measure molecular branching;
        M2 is edge-based while M1 is vertex-based.

        Example::

            ethane = chematic.from_smiles("CC")
            ethane.zagreb_m2  # 1 (one bond, degree-1 × degree-1)
            benzene = chematic.from_smiles("c1ccccc1")
            benzene.zagreb_m2  # 24 (6 bonds × (2×2))
        """
        ...

    def peoe_vsa(self) -> list[float]:
        """PEOE VSA — van der Waals surface area bins by partial charge."""
        ...

    def slogp_vsa(self) -> list[float]:
        """SLogP VSA — van der Waals surface area bins by SLogP contribution."""
        ...

    def smr_vsa(self) -> list[float]:
        """SMR VSA — van der Waals surface area bins by molar refractivity."""
        ...

    def estate_vsa(self) -> list[float]:
        """EState VSA — van der Waals surface area bins by E-state index."""
        ...

    def usrcat(self) -> list[float]:
        """USRCAT — 42-element topological shape/pharmacophore descriptor."""
        ...

    def to_sdf_with_charges(
        self,
        coords: list[list[float]],
        charges: list[float],
    ) -> str:
        """Write molecule as 3D SDF with ``> <PARTIAL_CHARGES>`` property."""
        ...

    def pharmacophore_fp(self) -> bytes:
        """Pharmacophore 2D fingerprint (256 bytes = 2048 bits). Compatible with :func:`tanimoto`."""
        ...

    def mhfp(self) -> list[int]:
        """MHFP — 128 u64 hash values. Use :func:`tanimoto_mhfp` for similarity."""
        ...

    def mhfp_config(
        self,
        radius: int = 2,
        num_hashes: int = 128,
        seed: int = 0,
    ) -> list[int]:
        """MinHash fingerprint with custom parameters.

        Returns a list of ``num_hashes`` unsigned 64-bit integers.
        Use :func:`tanimoto_mhfp` for similarity comparison.

        Args:
            radius: circular subgraph radius (default 2, like ECFP4).
            num_hashes: fingerprint length in hash slots (default 128).
            seed: hash seed for reproducibility (default 0).

        Example::

            fp = mol.mhfp_config(radius=3, num_hashes=256)
            sim = chematic.tanimoto_mhfp(fp, other.mhfp_config(radius=3, num_hashes=256))
        """
        ...

    # -- Sprint 3: charges, structure analysis, 3D shape, conformer tools ----

    def gasteiger_charges(self) -> list[float]:
        """Gasteiger–Marsili partial charges — one float per heavy atom.

        Use with :meth:`to_pdbqt` to complete the docking pipeline::

            coords = mol.generate_3d()
            charges = mol.gasteiger_charges()
            pdbqt = mol.to_pdbqt(coords, charges, "LIG")
        """
        ...

    def cip_stereo(self) -> list[dict]:
        """CIP stereochemistry assignments — list of ``{"atom_idx": int, "descriptor": str}`` dicts.

        ``descriptor`` is ``"R"``, ``"S"``, ``"E"``, or ``"Z"``.
        """
        ...

    def pains_alerts(self) -> list[str]:
        """Names of PAINS structural alerts matched by this molecule."""
        ...

    def brenk_alerts(self) -> list[str]:
        """Names of Brenk instability/toxicity alerts matched by this molecule."""
        ...

    @property
    def bbb_score(self) -> float:
        """Blood-brain barrier penetration score (0 = low, 1 = high penetration).

        See :meth:`admet` for the full profile.
        """
        ...

    @property
    def caco2(self) -> float:
        """Predicted Caco-2 intestinal permeability (nm/s). Higher = better oral absorption."""
        ...

    @property
    def herg_risk(self) -> float:
        """hERG cardiac toxicity risk score (0 = low, 1 = high risk)."""
        ...

    @property
    def cyp3a4_risk(self) -> float:
        """CYP3A4 inhibition risk score (0 = low, 1 = high)."""
        ...

    def ames_alerts(self) -> list[str]:
        """Names of Ames mutagenicity SMARTS alerts matched by this molecule.

        An empty list means the molecule passes all Ames filters
        (equivalent to :attr:`ames_passes` returning ``True``).
        Complements :meth:`pains_alerts` and :meth:`brenk_alerts`.

        Example::

            alerts = mol.ames_alerts()   # ["aromatic_amine", ...]
        """
        ...

    def named_functional_groups(self) -> list[str]:
        """Named functional groups detected — list of group names (e.g. ``"carboxyl"``)."""
        ...

    def pmi(self, coords: list[list[float]]) -> list[float]:
        """Principal Moments of Inertia — ``[PMI1, PMI2, PMI3]`` (ascending order)."""
        ...

    def npr(self, coords: list[list[float]]) -> list[float]:
        """Normalised Principal Moments — ``[NPR1, NPR2]``. Values in [0, 1]."""
        ...

    def asphericity(self, coords: list[list[float]]) -> float:
        """Asphericity — deviation from a perfect sphere. Range [0, 1]."""
        ...

    def eccentricity(self, coords: list[list[float]]) -> float:
        """Eccentricity — elongation measure. Range [0, 1]."""
        ...

    def radius_of_gyration(self, coords: list[list[float]]) -> float:
        """Radius of gyration (Å)."""
        ...

    def plane_of_best_fit(self, coords: list[list[float]]) -> float:
        """Plane of Best Fit (PBF) — deviation from least-squares plane (Å)."""
        ...

    def generate_3d_etkdg(self) -> list[list[float]]:
        """Generate 3D coordinates using the ETKDG algorithm (higher quality than :meth:`generate_3d`)."""
        ...

    def get_dihedral(
        self,
        coords: list[list[float]],
        i: int,
        j: int,
        k: int,
        l: int,
    ) -> Optional[float]:
        """Measure dihedral angle i–j–k–l in degrees. Returns ``None`` if collinear."""
        ...

    def set_dihedral(
        self,
        coords: list[list[float]],
        i: int,
        j: int,
        k: int,
        l: int,
        angle_deg: float,
    ) -> list[list[float]]:
        """Set dihedral i–j–k–l to ``angle_deg`` degrees. Returns new coords."""
        ...

    def get_bond_length(self, coords: list[list[float]], i: int, j: int) -> float:
        """Measure bond length between atoms i and j (Å)."""
        ...

    def get_bond_angle(self, coords: list[list[float]], i: int, j: int, k: int) -> float:
        """Measure bond angle i–j–k in degrees (j is the central atom)."""
        ...

    # -- Sprint 9 additions --------------------------------------------------

    def estate_indices(self) -> list[float]:
        """Per-atom E-state electrotopological indices (Kier & Hall).

        Returns one float per heavy atom. Equivalent to RDKit's
        ``EState.EState.EStateIndices(mol)``.

        Example::

            idx = mol.estate_indices()
        """
        ...

    def minimize_mmff94(self, coords: list[list[float]]) -> list[list[float]]:
        """Minimize 3D coordinates with the MMFF94 force field.

        Complements :meth:`minimize_uff`. Returns minimized coords as ``[[x,y,z], ...]``.

        Example::

            coords = mol.generate_3d()
            minimized = mol.minimize_mmff94(coords)
        """
        ...

    def minimize_dreiding(self, coords: list[list[float]]) -> list[list[float]]:
        """Minimize 3D coordinates with the DREIDING force field.

        Returns minimized coordinates as ``[[x, y, z], ...]``.
        Complements :meth:`minimize_mmff94` and :meth:`minimize_uff`.

        Example::

            coords = mol.generate_3d()
            minimized = mol.minimize_dreiding(coords)
        """
        ...

    @property
    def num_unspecified_stereocenters(self) -> int:
        """Number of stereocenters with unspecified (unknown) configuration.

        Equivalent to RDKit's ``rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)``.
        """
        ...

    def whim_getaway(self, coords: list[list[float]]) -> list[float]:
        """Combined WHIM + GETAWAY 3D descriptor vector.

        Equivalent to concatenating :meth:`whim` and :meth:`getaway`.
        Useful for single-call 3D featurisation pipelines (mordred compatible).

        Example::

            coords = mol.generate_3d()
            vec = mol.whim_getaway(coords)
        """
        ...

    # -- Sprint 17: ames_alerts, clearance_score, mr_per_atom, mmff94_charges_3d -

    @property
    def clearance_score(self) -> float:
        """Predicted hepatic clearance score (raw float 0.0–1.0).

        Lower = slower clearance; higher = faster.
        Complements :attr:`clearance_class` which returns ``"Low"``/``"Medium"``/``"High"``.
        Useful for ML regression targets.
        """
        ...

    def mr_per_atom(self) -> list[float]:
        """Per-atom molar refractivity contributions.

        Returns one float per heavy atom. ``sum(mol.mr_per_atom()) ≈ mol.molar_refractivity``.

        Example::

            mr = mol.mr_per_atom()
        """
        ...

    def topological_distance_matrix(self) -> list[list[int]]:
        """Topological (graph) distance matrix for all heavy atoms.

        Entry ``[i][j]`` is the shortest path length in bonds between heavy atom
        ``i`` and heavy atom ``j``. Diagonal entries are 0.  Disconnected atoms
        get ``2147483647`` (``u32::MAX``). Row/column order follows atom-insertion
        order (same as :meth:`hybridization_per_atom` and other per-atom vectors).

        Useful for scaffold graph analysis (ScaffoldGraph-style topology),
        molecular topology descriptors, and custom fingerprinting.

        Example::

            propane = chematic.from_smiles("CCC")
            dm = propane.topological_distance_matrix()
            # [[0, 1, 2], [1, 0, 1], [2, 1, 0]]
            assert dm[0][2] == 2   # C1 to C3 = 2 bonds
        """
        ...

    def hybridization_per_atom(self) -> list[int]:
        """Per-atom hybridization state as integers.

        Returns one value per heavy atom (index = atom order):
        ``1`` = sp, ``2`` = sp2, ``3`` = sp3, ``0`` = other/wildcard.

        Rules:
        - Aromatic atom → 2 (sp2)
        - Has triple bond → 1 (sp)
        - Has double bond → 2 (sp2)
        - Otherwise → 3 (sp3)

        Useful for scaffold modification (PromptSMILES), fragment building
        (BuildAMol), and custom QSAR feature generation.

        Example::

            mol = chematic.from_smiles("CC=O")   # ethanol-like, acetaldehyde
            mol.hybridization_per_atom()  # [3, 2, 2] (CH3=sp3, C=sp2, O=sp2)
        """
        ...

    def formal_charge_per_atom(self) -> list[int]:
        """Per-atom formal charge — one ``int`` per heavy atom.

        All values are 0 for neutral molecules.
        Charged atoms (e.g. ``[NH4+]``, ``[O-]``) have non-zero entries.

        Example::

            mol = chematic.from_smiles("[NH4+]")
            mol.formal_charge_per_atom()  # [1, 0, 0, 0, 0]  (N+, 4 H)
        """
        ...

    def implicit_hcount_per_atom(self) -> list[int]:
        """Per-atom implicit hydrogen count — one ``int`` per heavy atom.

        Counts the number of implicit (non-explicit) H atoms attached to each
        heavy atom. Consistent with ``sum(mol.implicit_hcount_per_atom()) ≈``
        total implicit H count.

        Useful for building 3D structures, atom-level featurization (BuildAMol),
        and ML model inputs.

        Example::

            mol = chematic.from_smiles("CC")   # ethane
            mol.implicit_hcount_per_atom()     # [3, 3]
        """
        ...

    # -- Sprint 10: element/bond counts, ring topology, ERG vec, canonical ----

    @property
    def num_fluorines(self) -> int:
        """Number of fluorine atoms."""
        ...

    @property
    def num_chlorines(self) -> int:
        """Number of chlorine atoms."""
        ...

    @property
    def num_bromines(self) -> int:
        """Number of bromine atoms."""
        ...

    @property
    def num_iodines(self) -> int:
        """Number of iodine atoms."""
        ...

    @property
    def num_phosphorus(self) -> int:
        """Number of phosphorus atoms."""
        ...

    @property
    def num_heteroatoms(self) -> int:
        """Number of heteroatoms (non-C, non-H heavy atoms). Equivalent to RDKit ``CalcNumHeteroatoms``."""
        ...

    @property
    def num_carbons(self) -> int:
        """Number of carbon atoms."""
        ...

    @property
    def num_nitrogens(self) -> int:
        """Number of nitrogen atoms."""
        ...

    @property
    def num_oxygens(self) -> int:
        """Number of oxygen atoms."""
        ...

    @property
    def num_sulfurs(self) -> int:
        """Number of sulfur atoms."""
        ...

    @property
    def num_hydrogens(self) -> int:
        """Total hydrogen count (implicit + explicit)."""
        ...

    @property
    def num_amide_bonds(self) -> int:
        """Number of amide bonds (–C(=O)–N–). Equivalent to RDKit ``CalcNumAmideBonds``."""
        ...

    @property
    def num_ester_bonds(self) -> int:
        """Number of ester bonds (–C(=O)–O–)."""
        ...

    @property
    def num_spiro_atoms(self) -> int:
        """Number of spiro atoms. Equivalent to RDKit ``CalcNumSpiroAtoms``."""
        ...

    @property
    def num_bridgehead_atoms(self) -> int:
        """Number of bridgehead atoms. Equivalent to RDKit ``CalcNumBridgeheadAtoms``."""
        ...

    @property
    def num_aromatic_heterocycles(self) -> int:
        """Number of aromatic heterocyclic rings."""
        ...

    @property
    def num_aliphatic_heterocycles(self) -> int:
        """Number of aliphatic (non-aromatic) heterocyclic rings."""
        ...

    @property
    def num_saturated_heterocycles(self) -> int:
        """Number of saturated heterocyclic rings."""
        ...

    @property
    def num_aliphatic_rings(self) -> int:
        """Number of aliphatic (non-aromatic) carbocyclic rings."""
        ...

    def morgan_fp_counts(self, radius: int = 2) -> dict[int, int]:
        """Count-based Morgan (ECFP) fingerprint — maps substructure hash → occurrence count.

        Unlike :meth:`ecfp4` (bit vector), preserves multiplicity.
        Equivalent to RDKit ``GetMorganFingerprint(mol, radius)``.

        Example::

            counts = mol.morgan_fp_counts(radius=2)
        """
        ...

    def pharmacophore_feature_counts(self) -> list[int]:
        """Count of each pharmacophore feature type: ``[donor, acceptor, aromatic, hydrophobic, positive, negative]``."""
        ...

    def mmff94_charges_typed(self) -> list[float]:
        """MMFF94 partial charges using the atom-type BCI model.

        A more accurate alternative to :meth:`mmff94_charges` (element-pair BCI).
        Returns one float per heavy atom.
        """
        ...

    def erg_vec(self) -> list[float]:
        """ERG (Extended Reduced Graph) continuous feature vector (length 315).

        Use with :func:`cosine_erg_vec` or :func:`tanimoto_erg_vec`.

        Example::

            v1 = mol1.erg_vec()
            v2 = mol2.erg_vec()
            sim = chematic.cosine_erg_vec(v1, v2)
        """
        ...

    def morgan_ranks(self) -> list[int]:
        """Morgan canonical rank per heavy atom. Equivalent to RDKit ``CanonicalRankAtoms``."""
        ...

    def canonical_atom_order(self) -> list[int]:
        """Canonical atom permutation — maps original index to canonical position."""
        ...

    def equivalent_atom_classes(self) -> list[int]:
        """Equivalent atom class IDs — atoms with same ID are symmetry-equivalent."""
        ...

    # -- Sprint 12: saturated rings, zwitterion, remove_salts, invert_stereocenter

    @property
    def num_saturated_rings(self) -> int:
        """Number of saturated (sp³) rings. Equivalent to RDKit ``CalcNumSaturatedRings``."""
        ...

    def has_zwitterion(self) -> bool:
        """Return ``True`` if the molecule has simultaneous positive and negative charges."""
        ...

    def normalize_zwitterion(self) -> "Mol":
        """Normalize zwitterion to neutral form via proton transfer."""
        ...

    def remove_salts(self) -> "Mol":
        """Remove salt fragments using the built-in salt catalog."""
        ...

    def invert_stereocenter(self, atom_idx: int) -> "Mol":
        """Invert stereochemistry at atom ``atom_idx`` (flip wedge/dash bonds → enantiomer)."""
        ...

    # -- Sprint 11: topo descriptors, ring perception, stereo validation, pharma

    @property
    def kappa1(self) -> float:
        """Hall-Kier κ₁ shape index — molecular size vs linear chain."""
        ...

    @property
    def kappa2(self) -> float:
        """Hall-Kier κ₂ shape index — branching degree."""
        ...

    @property
    def kappa3(self) -> float:
        """Hall-Kier κ₃ shape index — centrality of branching."""
        ...

    @property
    def wiener_index(self) -> float:
        """Wiener index — sum of topological distances between all heavy atom pairs."""
        ...

    @property
    def bertz_ct(self) -> float:
        """Bertz complexity index — information-theoretic graph complexity."""
        ...

    @property
    def chi0(self) -> float:
        """Zero-order path connectivity index χ⁰ (Kier & Hall)."""
        ...

    @property
    def chi1(self) -> float:
        """First-order path connectivity index χ¹."""
        ...

    @property
    def chi2(self) -> float:
        """Second-order path connectivity index χ²."""
        ...

    @property
    def chi3(self) -> float:
        """Third-order path connectivity index χ³."""
        ...

    @property
    def chi4(self) -> float:
        """Fourth-order path connectivity index χ⁴."""
        ...

    @property
    def chi0v(self) -> float:
        """Zero-order valence connectivity index χ⁰ᵥ."""
        ...

    @property
    def chi1v(self) -> float:
        """First-order valence connectivity index χ¹ᵥ."""
        ...

    @property
    def chi2v(self) -> float:
        """Second-order valence connectivity index χ²ᵥ."""
        ...

    @property
    def chi3v(self) -> float:
        """Third-order valence connectivity index χ³ᵥ."""
        ...

    @property
    def chi4v(self) -> float:
        """Fourth-order valence connectivity index χ⁴ᵥ."""
        ...

    def ring_membership(self) -> list[list[int]]:
        """SSSR ring membership per atom.

        Returns a list of N lists (one per heavy atom). Each inner list contains
        the 0-based SSSR ring indices to which that atom belongs.

        Example::

            membership = mol.ring_membership()
            ring_idxs = membership[atom_i]
        """
        ...

    def ring_sizes_for_atom(self, atom_idx: int) -> list[int]:
        """Ring sizes of all SSSR rings containing ``atom_idx``. Empty for acyclic atoms."""
        ...

    def is_fused_ring_system(self) -> bool:
        """Return ``True`` if the molecule has a fused ring system (rings sharing an edge)."""
        ...

    def ring_families(self) -> list[dict]:
        """Classify the ring systems (families) in this molecule by topology.

        Returns a list of dicts, one per connected ring system, each with:

        - ``kind``: ``"simple"`` | ``"fused"`` | ``"spiro"`` | ``"bridged"``
        - ``atom_indices``: list of heavy-atom indices in the ring system
        - ``ring_count``: number of SSSR rings in this family

        Example::

            for fam in mol.ring_families():
                print(fam['kind'], fam['ring_count'])

            # naphthalene → [{'kind': 'fused', 'ring_count': 2, 'atom_indices': [...]}]
        """
        ...

    def validate_stereo(self) -> list[dict]:
        """Validate stereochemistry; return list of errors (empty = consistent).

        Each error dict has keys ``atom_idx`` (int) and ``kind`` (str):
        ``"ImpossibleCenter"``, ``"ConflictingWedges"``, or ``"RedundantStereo"``.
        """
        ...

    def stereo_completeness(self) -> dict[str, int]:
        """Summarise stereocenters: ``{specified, unspecified, total_centers}``."""
        ...

    def stereo_from_coords(self, coords: list[list[float]]) -> list[dict]:
        """Perceive stereochemistry (R/S and E/Z) from 3D coordinates.

        Returns a list of dicts, one per assigned stereocentre:

        - ``atom_idx``: heavy-atom index of the chiral centre or E/Z bond atom
        - ``code``: ``"R"``, ``"S"``, ``"E"``, or ``"Z"``

        An empty list means no stereocentres could be assigned.

        Note:
            Only assigns R/S for atoms with four heavy-atom neighbours (no
            implicit H). Chiral centres with an implicit H (e.g. amino acids)
            are not assigned by this function.

        Equivalent to RDKit ``Chem.AssignStereochemistryFrom3D(mol)``.

        Example::

            mol = chematic.from_smiles("[C@@](N)(C)(F)Cl")
            coords = mol.generate_3d()
            for a in mol.stereo_from_coords(coords):
                print(a['atom_idx'], a['code'])   # 0 R or 0 S
        """
        ...

    def stereo_from_2d_coords(self, coords: list[list[float]]) -> list[dict]:
        """Perceive stereochemistry (R/S and E/Z) from 2D layout coordinates.

        Returns a list of dicts, one per assigned stereocentre:

        - ``atom_idx``: heavy-atom index of the chiral centre or E/Z bond atom
        - ``code``: ``"R"``, ``"S"``, ``"E"``, or ``"Z"``

        Coordinates are typically obtained from :func:`from_mol_block_with_coords`.
        For 3D-coordinate-based assignment use :meth:`stereo_from_coords`.

        Example::

            mol, name, coords_2d = chematic.from_mol_block_with_coords(block)
            for a in mol.stereo_from_2d_coords(coords_2d):
                print(a['atom_idx'], a['code'])
        """
        ...

    def pharmacophore_features(self) -> list[dict]:
        """Detect pharmacophore features (Donor, Acceptor, Aromatic, Hydrophobic, Positive, Negative).

        Each feature dict has keys:
        - ``"type"`` (str): feature type
        - ``"atom_idx"`` (int): primary atom index
        - ``"neighbor_indices"`` (list[int]): secondary atoms

        Example::

            feats = mol.pharmacophore_features()
            donors = [f for f in feats if f["type"] == "Donor"]
        """
        ...

    # -- Dunder methods ------------------------------------------------------

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

# ---------------------------------------------------------------------------
# Module-level functions
# ---------------------------------------------------------------------------

def tanimoto_erg(mol1: Mol, mol2: Mol) -> float:
    """Tanimoto similarity between two molecules using ERG fingerprints.

    Convenience alternative to ``tanimoto_erg_vec(m1.erg_vec(), m2.erg_vec())``.

    Example::

        sim = chematic.tanimoto_erg(mol1, mol2)
    """
    ...

def tanimoto_matrix(fps_a: list[bytes], fps_b: list[bytes]) -> list[list[float]]:
    """Compute an M×N Tanimoto similarity matrix.

    Returns M rows, each of length N: ``result[i][j] = Tanimoto(fps_a[i], fps_b[j])``.
    All fingerprints must have the same byte length.

    Example::

        matrix = chematic.tanimoto_matrix(
            [m.ecfp4() for m in queries],
            [m.ecfp4() for m in library],
        )
    """
    ...

def nearest_neighbors_from_fp(
    query_fp: bytes,
    db_fps: list[bytes],
    k: int = 10,
) -> list[tuple[int, float]]:
    """Find top-K nearest neighbors from precomputed fingerprint byte arrays.

    More efficient than :func:`top_k_similar_fp` when ``db_fps`` is reused
    across multiple queries (fingerprints computed only once).

    Returns list of ``(index, tanimoto_score)`` tuples, descending by score.

    Example::

        db_fps = [mol.ecfp4() for mol in library]   # compute once
        hits = chematic.nearest_neighbors_from_fp(query.ecfp4(), db_fps, k=10)
    """
    ...

def tanimoto_slice(query: bytes, db: list[bytes]) -> list[float]:
    """Compute Tanimoto similarity of one fingerprint against a list of fingerprints.

    All byte arrays must be the same length (e.g., all from :meth:`Mol.ecfp4`).
    More efficient than repeated :func:`tanimoto` calls for virtual screening.

    Example::

        db_fps = [mol.ecfp4() for mol in library]
        scores = chematic.tanimoto_slice(query.ecfp4(), db_fps)
    """
    ...

def tanimoto_reaction_fp(rxn1: str, rxn2: str) -> float:
    """Tanimoto similarity between two reactions using reaction fingerprints.

    Parses both reaction SMILES and computes fingerprint similarity.

    Raises:
        ValueError: on invalid reaction SMILES.

    Example::

        sim = chematic.tanimoto_reaction_fp("CC>>CO", "CC>>CN")
    """
    ...

def named_pattern(name: str) -> Optional[str]:
    """Look up a built-in named SMARTS pattern by name.

    Returns the SMARTS string if ``name`` is known, ``None`` otherwise.

    Available names include: ``"donor"``, ``"donor_strict"``,
    ``"acceptor"``, ``"acceptor_strict"``, ``"aromatic"``,
    ``"aromatic_ring"``, ``"hydrophobic"``, ``"positive"``, ``"negative"``.

    Example::

        if smarts := chematic.named_pattern("donor"):
            hits = chematic.smarts_find(smarts, mol)
    """
    ...

def parse_smi_file(content: str) -> list[tuple[Mol, str]]:
    """Parse a ``.smi`` file into ``(Mol, name)`` pairs.

    Each line is ``SMILES[<tab>name]``. Lines with invalid SMILES, blank lines,
    and ``#``-comment lines are silently skipped.
    Equivalent to RDKit's ``Chem.SmilesMolSupplier``.

    Example::

        records = chematic.parse_smi_file(open("library.smi").read())
        for mol, name in records:
            print(name, mol.mw)
    """
    ...

def write_smi_file(records: list[tuple[Mol, str]]) -> str:
    """Write ``(Mol, name)`` pairs to ``.smi`` format.

    Format: ``SMILES<TAB>name<NEWLINE>`` per record (name omitted if empty).
    Equivalent to RDKit's ``Chem.SmilesWriter``.

    Example::

        text = chematic.write_smi_file([(mol1, "cpd1"), (mol2, "cpd2")])
    """
    ...

def atom_color(atomic_num: int) -> str:
    """CSS color string for an element by atomic number.

    Returns CPK/standard coloring (e.g. ``"#FF0000"`` for oxygen).

    Example::

        print(chematic.atom_color(8))   # "#FF0000"
    """
    ...

def atom_color_rgb(atomic_num: int) -> tuple[int, int, int]:
    """RGB color triple for an element by atomic number.

    Returns the same color as :func:`atom_color` as a ``(R, G, B)`` tuple (0–255).

    Example::

        r, g, b = chematic.atom_color_rgb(8)   # (255, 0, 0)
    """
    ...

def from_smiles(smiles: str) -> Mol:
    """Parse a SMILES string and return a :class:`Mol`.

    Raises:
        ValueError: on invalid SMILES.

    Example::

        mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    """
    ...

def from_cxsmiles(s: str) -> tuple[Mol, dict]:
    """Parse a CXSMILES string and return the molecule with CX metadata.

    Returns a 2-tuple ``(mol, cx)`` where ``cx`` is a dict with:

    - ``atom_labels``: list of atom label strings (or ``None`` per atom)
    - ``atom_props``: list of ``{"atom_idx", "key", "value"}`` dicts
    - ``atom_radicals``: list of radical class integers (or ``None`` per atom)

    CXSMILES without a CX extension block behaves like :func:`from_smiles`
    (all CX fields contain ``None`` entries).

    Raises:
        ValueError: on parse failure.

    Example::

        mol, cx = chematic.from_cxsmiles("CC |$R1;R2$|")
        print(cx['atom_labels'])   # ['R1', 'R2']
    """
    ...

def from_condensed(formula: str) -> Optional[Mol]:
    """Parse a condensed molecular formula (e.g., ``"CH3OH"``, ``"C6H12O6"``) into a :class:`Mol`.

    Returns ``None`` if the formula is unknown or cannot be resolved to a unique structure.
    Equivalent to chempy's condensed formula support.

    Example::

        mol = chematic.from_condensed("CH3OH")   # methanol
        if mol:
            print(mol.smiles)  # CO
    """
    ...

def from_cml(cml_str: str) -> Mol:
    """Parse a Chemical Markup Language (CML) string and return a :class:`Mol`.

    Raises:
        ValueError: on parse failure.
    """
    ...

def from_cdxml(cdxml_str: str) -> Mol:
    """Parse a ChemDraw XML (CDXML) string and return a :class:`Mol`.

    Raises:
        ValueError: on parse failure.
    """
    ...

def from_mol_v3000(block: str) -> Mol:
    """Parse an MDL MOL V3000 block and return a :class:`Mol`.

    Raises:
        ValueError: on parse failure.
    """
    ...

def from_mol_v3000_with_coords(block: str) -> tuple[Mol, str, list[list[float]]]:
    """Parse a MDL MOL V3000 block and return the molecule with its 2D layout coordinates.

    Returns a 3-tuple ``(mol, name, coords_2d)`` identical to
    :func:`from_mol_block_with_coords` but for V3000 input.

    Raises:
        ValueError: on parse failure.

    Example::

        mol, name, coords_2d = chematic.from_mol_v3000_with_coords(block)
        new_block = mol.to_mol_v3000(coords_2d, name=name)
    """
    ...

def from_mol_v3000_with_diagnostics(
    block: str,
) -> tuple[Mol, str, list[list[float]], list[dict[str, object]]]:
    """Parse a MDL MOL V3000 block, returning stereo-perception diagnostics
    alongside the molecule.

    Same shape as :func:`from_mol_block_with_diagnostics` but for V3000 input.

    Raises:
        ValueError: on parse failure.
    """
    ...

def from_mol_block(block: str) -> Mol:
    """Parse a MOL/SDF block and return a :class:`Mol`.

    Raises:
        ValueError: on parse failure.
    """
    ...

def from_mol_block_with_coords(block: str) -> tuple[Mol, str, list[list[float]]]:
    """Parse a MDL MOL V2000 block and return the molecule with its 2D layout coordinates.

    Returns a 3-tuple ``(mol, name, coords_2d)`` where:

    - ``mol``: :class:`Mol` object
    - ``name``: molecule name from the MOL header (may be empty)
    - ``coords_2d``: list of ``[x, y]`` pairs (one per heavy atom, Å)

    Use :func:`from_mol_block` if you only need the molecule graph.
    Use this function to preserve 2D layout for round-tripping via :meth:`Mol.to_mol_block_2d`.

    Raises:
        ValueError: on parse failure.

    Example::

        mol, name, coords_2d = chematic.from_mol_block_with_coords(block)
        new_block = mol.to_mol_block_2d(coords_2d, name=name)
    """
    ...

def from_mol_block_with_diagnostics(
    block: str,
) -> tuple[Mol, str, list[list[float]], list[dict[str, object]]]:
    """Parse a MDL MOL V2000 block, returning stereo-perception diagnostics
    alongside the molecule.

    Returns a 4-tuple ``(mol, name, coords_2d, stereo_diagnostics)``.
    ``stereo_diagnostics`` is a list of ``{"atom_idx": int, "reason": str}``
    dicts, one per wedge/hash center that could not be resolved — ``reason``
    is one of ``"contradictory_wedges"``, ``"missing_coordinate"``,
    ``"degenerate_geometry"``, or ``"unsupported_coordination"``. Empty
    unless a wedge/hash bond was present at some center and got rejected; an
    atom with no wedge/hash bond at all never produces an entry.

    Local tetrahedral parity (``Atom.chirality``) is always perceived
    automatically — this function differs from
    :func:`from_mol_block_with_coords` only in also surfacing *why* any
    center was rejected.

    Raises:
        ValueError: on parse failure.
    """
    ...

def parse_sdf_with_coords(text: str) -> list[tuple[Mol, str, list[list[float]]]]:
    """Parse a multi-record SDF string and return all molecules with their 2D layout coordinates.

    Returns a list of 3-tuples ``(mol, name, coords_2d)`` — one per SDF record.
    Invalid records are silently skipped (same behaviour as :func:`iter_sdf`).

    This is the batch equivalent of :func:`from_mol_block_with_coords`.

    Example::

        with open("library.sdf") as f:
            records = chematic.parse_sdf_with_coords(f.read())
        for mol, name, coords_2d in records:
            new_block = mol.to_mol_block_2d(coords_2d, name=name)
    """
    ...

def from_mol2(mol2_str: str) -> Mol:
    """Parse a Tripos MOL2 string and return a :class:`Mol`.

    Reads the mandatory @<TRIPOS>MOLECULE, @<TRIPOS>ATOM, and @<TRIPOS>BOND
    sections.  Atom coordinates are parsed but not stored on the :class:`Mol`
    object (use :func:`generate_3d` to generate or reload coordinates).

    Args:
        mol2_str: Complete MOL2-format string (not a file path).

    Returns:
        Mol: Parsed molecule.

    Raises:
        ValueError: If the MOL2 string is malformed or missing required sections.

    Example::

        with open("ligand.mol2") as f:
            mol = chematic.from_mol2(f.read())
        print(mol.mw)
        print(mol.formula)

        # Round-trip
        mol = chematic.from_smiles("c1ccccc1")
        mol2_str = mol.to_mol2()
        mol2 = chematic.from_mol2(mol2_str)
        assert mol2.atom_count() == mol.atom_count()
    """
    ...

def from_pdbqt(pdbqt_str: str) -> Mol:
    """Parse an AutoDock PDBQT string and return a :class:`Mol`.

    Only the molecular graph is extracted; 3D coordinates and partial charges
    are discarded.

    Raises:
        ValueError: on parse failure.

    Example::

        with open("ligand.pdbqt") as f:
            mol = chematic.from_pdbqt(f.read())
        print(mol.mw)
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

def is_valid_smarts(smarts: str) -> bool:
    """Return ``True`` if ``smarts`` is a valid SMARTS pattern.

    Mirrors :func:`is_valid_smiles` for SMARTS validation.
    Useful for validating user-supplied patterns before calling
    :func:`smarts_match` or :func:`smarts_find`.

    Example::

        chematic.is_valid_smarts("c1ccccc1")  # True
        chematic.is_valid_smarts("[invalid")  # False
        chematic.is_valid_smarts("[#6]-[#7]") # True
    """
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

def from_pdb(pdb_str: str) -> tuple[Mol, list[list[float]]]:
    """Parse a PDB string and return ``(Mol, coords)``.

    ``coords`` is a list of ``[x, y, z]`` lists (Å), one per heavy atom.
    Bond information is inferred from inter-atom distances.

    Raises:
        ValueError: when no ATOM/HETATM records are found.

    Example::

        mol, coords = chematic.from_pdb(open("ligand.pdb").read())
        print(mol.formula, len(coords))
    """
    ...

def from_xyz(xyz_str: str) -> tuple[Mol, list[list[float]]]:
    """Parse an XYZ string and return ``(Mol, coords)``.

    ``coords`` is a list of ``[x, y, z]`` lists (Å), one per heavy atom.

    Raises:
        ValueError: on parse failure.

    Example::

        mol, coords = chematic.from_xyz(open("molecule.xyz").read())
    """
    ...

def from_extxyz(text: str) -> dict:
    """Parse an Extended XYZ (extxyz) frame and return a dict describing it.

    A plain XYZ file (free-form comment, no ``Lattice=``/``Properties=``)
    parses too, with ``lattice=None`` and empty ``properties``/``info``.

    Returns:
        dict: ``{"mol": Mol, "coords": list[list[float]],
        "lattice": list[float] | None, "properties": dict, "info": dict}``.

    Raises:
        ValueError: on malformed input.

    Example::

        result = chematic.from_extxyz(open("frame.xyz").read())
        forces = result["properties"].get("forces")
    """
    ...

def from_extxyz_all(text: str) -> list[dict]:
    """Parse every frame of a multi-frame extxyz trajectory.

    See :func:`from_extxyz` for the shape of each returned dict.

    Raises:
        ValueError: on the first parse failure.
    """
    ...

def to_extxyz(
    mol: Mol,
    coords: list[list[float]],
    lattice: Optional[list[float]] = None,
    properties: Optional[dict] = None,
    info: Optional[dict] = None,
) -> str:
    """Write a molecule + coordinates as an Extended XYZ (extxyz) frame.

    ``properties`` is ``dict[str, list[list[float]]]`` (real-valued per-atom
    columns only, e.g. ``{"forces": [[fx, fy, fz], ...]}``). ``info`` is
    ``dict[str, str]`` of extra frame metadata (e.g. ``{"energy": "-76.4"}``).

    Raises:
        ValueError: if `coords`' length doesn't match `mol`'s atom count.
    """
    ...

def parse_mmcif(
    text: str,
    max_input_bytes: Optional[int] = None,
    max_atoms: Optional[int] = None,
    max_line_len: Optional[int] = None,
) -> dict:
    """Parse an mmCIF (macromolecular CIF, ``_atom_site.*`` loop) string.

    Returns:
        dict: ``{"mol": Mol, "coords": list[list[float]], "atoms":
        list[dict], "cell": dict | None, "space_group": str | None,
        "unhandled_columns": list[str]}``.

    Raises:
        ValueError: on parse failure.
    """
    ...

def write_mmcif(
    atoms: list[dict],
    cell: Optional[dict] = None,
    space_group: Optional[str] = None,
    data_block_name: str = "chematic",
) -> str:
    """Write an mmCIF file from atom records (same dict shape as
    :func:`parse_mmcif`'s ``"atoms"``).
    """
    ...

def parse_pqr(
    text: str,
    max_input_bytes: Optional[int] = None,
    max_atoms: Optional[int] = None,
    max_line_len: Optional[int] = None,
) -> dict:
    """Parse a PQR (PDB-like ATOM/HETATM + per-atom charge/radius) string.

    Returns:
        dict: ``{"mol": Mol, "coords": list[list[float]], "atoms": list[dict]}``.

    Raises:
        ValueError: on parse failure.
    """
    ...

def write_pqr(atoms: list[dict]) -> str:
    """Write a PQR file from atom records (same dict shape as
    :func:`parse_pqr`'s ``"atoms"``); ``element`` is inferred via
    :func:`infer_element` when omitted.
    """
    ...

def infer_element(group_pdb: str, res_name: str, atom_name: str) -> Optional[str]:
    """Infer an element symbol from a PQR/PDB atom name."""
    ...

def parse_orca_input(text: str) -> dict:
    """Parse an ORCA input file (``.inp``).

    Returns:
        dict: ``{"comments": list[str], "keywords": list[str], "blocks":
        list[dict], "coords": dict | None}``.

    Raises:
        ValueError: on parse failure.
    """
    ...

def write_orca_input(input: dict) -> str:
    """Write an ORCA input file from a dict, same shape as
    :func:`parse_orca_input`'s return value.
    """
    ...

def parse_orca_output(text: str) -> dict:
    """Parse an ORCA output file (``.out``/``.log``).

    Returns:
        dict: ``{"charge": int | None, "multiplicity": int | None,
        "final_energy_hartree": float | None, "trajectory": list[dict],
        "frequencies_cm1": list[float], "termination": dict,
        "optimization_convergence": str}``.

    Raises:
        ValueError: on a non-finite value or oversized input.
    """
    ...

def parse_qcschema_molecule(text: str) -> dict:
    """Parse a QCSchema ``Molecule`` JSON document.

    Returns every QCSchema ``Molecule`` field plus ``"mol"``/``"coords"``
    convenience keys; :func:`write_qcschema_molecule` strips both before
    re-serializing, so this dict round-trips through it directly.

    Raises:
        ValueError: on malformed JSON or a schema violation.
    """
    ...

def write_qcschema_molecule(molecule: dict) -> str:
    """Serialize a QCSchema ``Molecule`` dict to canonical JSON text."""
    ...

def chematic_to_qc_molecule(
    mol: Mol,
    coords: list[list[float]],
    molecular_charge: float = 0.0,
    molecular_multiplicity: int = 1,
) -> dict:
    """Convert a chematic ``Mol`` + coordinates into a QCSchema ``Molecule`` dict."""
    ...

def qc_molecule_to_chematic(molecule: dict) -> tuple[Mol, list[list[float]], float, int]:
    """Convert a QCSchema ``Molecule`` dict into ``(Mol, coords,
    molecular_charge, molecular_multiplicity)``.
    """
    ...

def parse_atomic_input(text: str) -> dict:
    """Parse a QCSchema ``AtomicInput`` JSON document.

    Returns every ``AtomicInput`` field plus ``"mol"``/``"coords"``
    convenience keys built from ``"molecule"``.

    Raises:
        ValueError: on malformed JSON or a schema violation.
    """
    ...

def write_atomic_input(input: dict) -> str:
    """Serialize an ``AtomicInput`` dict to canonical QCSchema JSON text."""
    ...

def parse_atomic_result(text: str) -> dict:
    """Parse a QCSchema ``AtomicResult`` JSON document.

    Raises:
        ValueError: on malformed JSON or a schema violation.
    """
    ...

def write_atomic_result(result: dict) -> str:
    """Serialize an ``AtomicResult`` dict to canonical QCSchema JSON text."""
    ...

def parse_lammps_data(text: str, atom_style: str) -> dict:
    """Parse a LAMMPS data file (``read_data`` command format).

    ``atom_style`` is one of ``"atomic"``/``"charge"``/``"molecular"``/``"full"``.

    Returns:
        dict: ``{"counts": dict[str, int], "atom_style": str, "box": dict,
        "masses": list[dict], "atoms": list[dict], "velocities": list[dict],
        "bonds": list[dict], "unparsed_sections": list[tuple[str, str]]}``.

    Raises:
        ValueError: on parse failure.
    """
    ...

def write_lammps_data(data: dict) -> str:
    """Write a LAMMPS data file from a dict, same shape as
    :func:`parse_lammps_data`'s return value.
    """
    ...

def box_bounds_to_true(
    bound_lo: tuple[float, float, float],
    bound_hi: tuple[float, float, float],
    tilt: Optional[tuple[float, float, float]] = None,
) -> dict:
    """Convert a dump file's ``ITEM: BOX BOUNDS`` values into the true box."""
    ...

def true_to_box_bounds(
    box_dict: dict,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Inverse of :func:`box_bounds_to_true`."""
    ...

def parse_lammps_dump_frame(text: str) -> LammpsDumpFrame:
    """Parse a single LAMMPS dump frame.

    Raises:
        ValueError: on parse failure.
    """
    ...

def parse_lammps_dump_all(text: str) -> list[LammpsDumpFrame]:
    """Parse every frame of a LAMMPS dump/trajectory file.

    Raises:
        ValueError: on the first parse failure.
    """
    ...

def write_lammps_dump_frame(frame: LammpsDumpFrame) -> str: ...
def write_lammps_trajectory(frames: list[LammpsDumpFrame]) -> str: ...
def tanimoto_map4(a: list[int], b: list[int]) -> float:
    """Estimate MAP4 Tanimoto similarity between two MAP4 fingerprints.

    ``a`` and ``b`` must be lists of 1024 integers as returned by :meth:`Mol.map4`.
    Uses position-wise matching (not bitwise AND/OR), so do NOT use :func:`tanimoto`.

    Returns:
        float in [0, 1].

    Example::

        sim = chematic.tanimoto_map4(mol1.map4(), mol2.map4())
    """
    ...

def butina_cluster(smiles: list[str], cutoff: float = 0.65) -> list[list[int]]:
    """Butina clustering — group molecules by ECFP4 Tanimoto similarity.

    Returns a list of clusters; each cluster is a list of SMILES indices (centroid first).
    Clusters are sorted by size (largest first). Invalid SMILES are silently skipped.

    Args:
        smiles: list of SMILES strings.
        cutoff: Tanimoto similarity threshold (default 0.65).

    Example::

        clusters = chematic.butina_cluster(smiles, 0.65)
        for c in clusters:
            print(f"centroid: {smiles[c[0]]}, size: {len(c)}")
    """
    ...

def maxmin_picks(smiles: list[str], n: int) -> list[int]:
    """MaxMin diversity picking — select ``n`` maximally diverse molecules.

    Returns a list of indices into ``smiles``, in selection order.
    Uses ECFP4 Tanimoto distance. Invalid SMILES are silently skipped.

    Example::

        picks = chematic.maxmin_picks(smiles, 100)
        diverse_set = [smiles[i] for i in picks]
    """
    ...

def from_rxn_file(text: str) -> str:
    """Parse a MDL RXN V2000 file and return the canonical reaction SMILES.

    Raises:
        ValueError: on parse failure.

    Example::

        with open("reaction.rxn") as f:
            rxn_smiles = chematic.from_rxn_file(f.read())
        ae = chematic.atom_economy(rxn_smiles)
    """
    ...

def to_rxn_file(reaction_smiles: str) -> str:
    """Convert a reaction SMILES string to MDL RXN V2000 format.

    Raises:
        ValueError: on invalid reaction SMILES.

    Example::

        block = chematic.to_rxn_file("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
        with open("output.rxn", "w") as f:
            f.write(block)
    """
    ...

def atom_economy(reaction_smiles: str) -> float:
    """Atom economy of a reaction (green chemistry metric).

    atom_economy = MW(desired products) / MW(all reactants) × 100.

    Example::

        ae = chematic.atom_economy("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
        print(f"{ae:.1f}%")
    """
    ...

def balance_check(reaction_smiles: str) -> dict[str, object]:
    """Check whether a reaction SMILES is atom-balanced.

    Returns a dict with keys ``balanced`` (bool) and ``diff`` (list of str).

    Example::

        result = chematic.balance_check("C+O>>CO")
        print(result["balanced"], result["diff"])
    """
    ...

def enumerate_library(
    smirks: str,
    fragment_sets: list[list[str]],
    max_size: int = 1_000_000,
) -> list[str]:
    """Enumerate a combinatorial library from a SMIRKS template and fragment SMILES.

    Args:
        smirks: Reaction SMIRKS template.
        fragment_sets: List of SMILES lists (one per reactant slot).
        max_size: Maximum library size.

    Returns:
        List of product SMILES strings.

    Example::

        products = chematic.enumerate_library(
            "[C:1]C(=O)Cl.[N:2]>>[C:1]C(=O)[N:2]",
            [["c1ccccc1", "CC"], ["N", "CN"]],
        )
    """
    ...

def enumerate_library_2way(
    smirks: str,
    scaffolds: list[str],
    building_blocks: list[str],
    max_size: int = 1_000_000,
) -> list[str]:
    """Enumerate a 2-fragment combinatorial library (scaffold × building block).

    Convenience alternative to ``enumerate_library(smirks, [scaffolds, building_blocks])``.
    The most common combinatorial chemistry pattern.

    Args:
        smirks: Reaction SMIRKS template.
        scaffolds: SMILES for the first reactant slot (scaffolds).
        building_blocks: SMILES for the second reactant slot (building blocks).
        max_size: Maximum library size (default 1,000,000).

    Returns:
        List of product SMILES strings.

    Example::

        products = chematic.enumerate_library_2way(
            "[C:1]C(=O)Cl.[N:2]>>[C:1]C(=O)[N:2]",
            scaffolds=["c1ccccc1C(=O)Cl", "CC(=O)Cl"],
            building_blocks=["N", "CN"],
        )
    """
    ...

def enumerate_library_3way(
    smirks: str,
    scaffolds: list[str],
    r1_set: list[str],
    r2_set: list[str],
    max_size: int = 1_000_000,
) -> list[str]:
    """Enumerate a 3-fragment combinatorial library (scaffold × R1 × R2).

    Convenience alternative to ``enumerate_library(smirks, [scaffolds, r1_set, r2_set])``.
    Covers scaffold-decoration with two variable positions.

    Args:
        smirks: Reaction SMIRKS template.
        scaffolds: SMILES for the scaffold slot.
        r1_set: SMILES for the R1 slot.
        r2_set: SMILES for the R2 slot.
        max_size: Maximum library size (default 1,000,000).

    Returns:
        List of product SMILES strings.

    Example::

        products = chematic.enumerate_library_3way(
            "[C:1]C(=O)Cl.[N:2]>>[C:1]C(=O)[N:2]",
            scaffolds=["CC(=O)Cl"],
            r1_set=["N", "CN"],
            r2_set=["c1ccccc1", "CC"],
        )
    """
    ...

def dice_similarity(a: bytes, b: bytes) -> float:
    """Dice similarity between two fingerprint byte arrays.

    Dice = 2|A∩B| / (|A| + |B|). Compatible with ECFP4, MACCS, ERG, etc.

    Example::

        sim = chematic.dice_similarity(mol1.ecfp4(), mol2.ecfp4())
    """
    ...

def tversky_similarity(a: bytes, b: bytes, alpha: float, beta: float) -> float:
    """Tversky similarity between two fingerprint byte arrays.

    Tversky(α, β) = |A∩B| / (α|A\\B| + β|B\\A| + |A∩B|).

    - α=β=0.5 → Dice;  α=β=1.0 → Tanimoto
    - α=0, β=1 → recall-biased (useful for sub-structure queries)

    Example::

        sim = chematic.tversky_similarity(query.ecfp4(), target.ecfp4(), 0.0, 1.0)
    """
    ...

def tanimoto_mhfp(a: list[int], b: list[int]) -> float:
    """Estimate MHFP Tanimoto similarity between two MHFP fingerprints.

    ``a`` and ``b`` must be lists of 128 u64 values as returned by :meth:`Mol.mhfp`.

    Example::

        sim = chematic.tanimoto_mhfp(mol1.mhfp(), mol2.mhfp())
    """
    ...

def tanimoto_spectrophores(a: list[float], b: list[float]) -> float:
    """Tanimoto-like similarity between two Spectrophores fingerprints.

    Uses the USR formula ``S = 1 / (1 + mean|a − b|)``, returning values in (0, 1].
    Both vectors must have the same length (typically 48).

    Example::

        coords1 = mol1.generate_3d()
        coords2 = mol2.generate_3d()
        fp1 = mol1.spectrophores(coords1)
        fp2 = mol2.spectrophores(coords2)
        sim = chematic.tanimoto_spectrophores(fp1, fp2)
    """
    ...

def cosine_erg_vec(a: list[float], b: list[float]) -> float:
    """Cosine similarity between two ERG continuous feature vectors.

    Both ``a`` and ``b`` must have length 315 (from :meth:`Mol.erg_vec`).
    Returns a value in [0, 1].

    Example::

        sim = chematic.cosine_erg_vec(mol1.erg_vec(), mol2.erg_vec())
    """
    ...

def tanimoto_erg_vec(a: list[float], b: list[float]) -> float:
    """Tanimoto similarity between two ERG continuous feature vectors.

    Both ``a`` and ``b`` must have length 315 (from :meth:`Mol.erg_vec`).
    Returns a value in [0, 1].

    Example::

        sim = chematic.tanimoto_erg_vec(mol1.erg_vec(), mol2.erg_vec())
    """
    ...

def top_k_similar_fp(
    query: str,
    smiles: list[str],
    k: int = 10,
    fp: Optional[str] = None,
) -> list[tuple[int, float]]:
    """Find top-K similar molecules using a selectable fingerprint type.

    Args:
        query: SMILES of the query molecule.
        smiles: Library of SMILES to search.
        k: Number of results to return.
        fp: Fingerprint type — ``"ecfp4"`` (default), ``"ecfp6"``,
            ``"ecfp4_chiral"``, ``"fcfp4"``, ``"maccs"``, ``"topo_path"``.

    Returns:
        List of ``(index, tanimoto_score)`` tuples, descending by score.

    Example::

        results = chematic.top_k_similar_fp("c1ccccc1", smiles, k=5, fp="maccs")
    """
    ...

def align_coords(
    probe: list[list[float]],
    reference: list[list[float]],
) -> tuple[list[list[float]], float]:
    """Align ``probe`` onto ``reference`` using the Kabsch algorithm.

    Both lists must have the same number of ``[x,y,z]`` entries (atom correspondence
    is assumed). Returns ``(aligned_coords, rmsd)``.

    Example::

        aligned, rmsd = chematic.align_coords(mol.generate_3d(), ref_coords)
    """
    ...

def rmsd(coords_a: list[list[float]], coords_b: list[list[float]]) -> float:
    """RMSD between two sets of paired 3D coordinates **without** alignment.

    Example::

        r = chematic.rmsd(mol.generate_3d(), ref_coords)
    """
    ...

def depict_grid(mols: list[Mol], cols: int) -> str:
    """Render a list of molecules as a grid SVG.

    Example::

        svg = chematic.depict_grid([mol1, mol2, mol3], cols=3)
    """
    ...

def query_reaction(reaction_smiles: str, smarts: str) -> bool:
    """Check whether a reaction SMILES matches a reaction SMARTS pattern.

    Returns ``True`` if matched. Raises ``ValueError`` on invalid input.

    Example::

        matched = chematic.query_reaction("CC>>CO", "[C:1]>>[C:1]O")
    """
    ...

def batch_query_reactions(reactions: list[str], smarts: str) -> dict[str, object]:
    """Query a list of reaction SMILES against a single SMARTS pattern.

    Returns a dict with keys:
    - ``total`` (int): reactions processed
    - ``matching`` (int): reactions that matched
    - ``match_pct`` (float): match percentage (0–100)
    - ``matches`` (list[tuple[int, bool]]): per-reaction (original_index, matched) pairs

    Invalid SMILES are silently skipped. Raises ``ValueError`` on invalid SMARTS.

    Example::

        r = chematic.batch_query_reactions(["CC>>CO", "c1ccccc1>>c1ccccc1N"], "[C:1]>>[C:1]O")
        print(r["matching"], "/", r["total"])
    """
    ...

def tanimoto_pharmacophore_3d(a: bytes, b: bytes) -> float:
    """Tanimoto similarity between two 3D pharmacophore fingerprints.

    Both ``a`` and ``b`` must be byte arrays from :meth:`Mol.pharmacophore_fp_3d`.
    Returns a value in [0, 1].

    Example::

        fp1 = mol1.pharmacophore_fp_3d(coords1)
        fp2 = mol2.pharmacophore_fp_3d(coords2)
        sim = chematic.tanimoto_pharmacophore_3d(fp1, fp2)
    """
    ...

def reaction_svg(reaction_smiles: str) -> str:
    """Render a reaction SMILES as an SVG diagram.

    Returns an SVG string showing reactants → products with an arrow.
    Equivalent to RDKit's ``Draw.ReactionToImage(rxn)``.

    Raises:
        ValueError: on invalid reaction SMILES.

    Example::

        svg = chematic.reaction_svg("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
    """
    ...

def rgroup_decompose(
    scaffold_smarts: str, mols: list[Mol]
) -> list[dict[str, object] | None]:
    """R-group decomposition — split molecules into a scaffold core and R-group substituents.

    ``scaffold_smarts``: SMARTS pattern defining the common scaffold.
    ``mols``: list of :class:`Mol` objects to decompose.

    Returns a list parallel to ``mols``.  Each element is either:
    - A dict with keys ``mol_idx`` (int), ``core`` (str, scaffold SMILES with ``[*]``
      at attachment points), and ``R1``, ``R2``, … (str, R-group SMILES).
    - ``None`` if the scaffold did not match that molecule.

    R-group numbering is 1-based and determined by ascending core-atom index.

    Raises:
        ValueError: if ``scaffold_smarts`` is an invalid SMARTS pattern.

    Example::

        mols = [chematic.from_smiles(s) for s in ["CCc1ccccc1", "CCCc1ccccc1", "Nc1ccccc1"]]
        results = chematic.rgroup_decompose("c1ccccc1", mols)
        for r in results:
            if r is not None:
                print(r["R1"])  # e.g. "[*]CC", "[*]CCC", "[*]N"
    """
    ...

def similarity_map_svg(mol: Mol, weights: list[float]) -> str:
    """Render a molecule SVG with atoms coloured by per-atom weights.

    ``mol``: molecule to render.
    ``weights``: list of floats, one per heavy atom.
    Positive → blue tint, negative → red tint, zero → white (no tint).
    Weights are normalised to the maximum absolute value before colour mapping.

    Example::

        weights = mol.logp_per_atom()
        svg = chematic.similarity_map_svg(mol, weights)
    """
    ...

def activity_cliffs(
    mols: list[Mol],
    activities: list[float],
    sim_threshold: float = 0.65,
    cliff_delta: float = 2.0,
) -> list[dict[str, object]]:
    """Detect activity cliffs in a set of molecules with known activity values.

    An activity cliff is a structurally similar pair with a large activity gap —
    a classic signal of SAR sensitivity, as used in MolScore and mol-eval.

    ``mols``: list of :class:`Mol` objects.
    ``activities``: list of floats (one per mol, e.g. pIC50 values).
    ``sim_threshold``: minimum ECFP4 Tanimoto similarity (default 0.65).
    ``cliff_delta``: minimum ``|activity_i − activity_j|`` to be a cliff (default 2.0).

    Returns a list of dicts sorted by ``similarity`` descending:
      ``mol_a_idx`` (int), ``mol_b_idx`` (int), ``similarity`` (float), ``activity_delta`` (float).

    Example::

        mols = [chematic.from_smiles(s) for s in ["c1ccccc1", "Cc1ccccc1"]]
        cliffs = chematic.activity_cliffs(mols, [5.0, 8.5], sim_threshold=0.0, cliff_delta=2.0)
        # [{"mol_a_idx": 0, "mol_b_idx": 1, "similarity": 0.xx, "activity_delta": 3.5}]
    """
    ...

def parse_formula(formula: str) -> dict[str, int]:
    """Parse a Hill-notation molecular formula string into an element count dict.

    Mirrors the API of PyPI libraries **chemparse** and **chemformula**.

    Supported syntax:
      - Simple formulas: ``"H2O"``, ``"C6H12O6"``
      - Parentheses with multipliers: ``"Ca(OH)2"`` → ``{"Ca":1,"O":2,"H":2}``
      - SMILES-style brackets: ``"[NH4]+"`` → ``{"N":1,"H":4}``
      - Trailing charge signs are ignored: ``"NH4+"`` → same as ``"NH4"``

    Raises:
        ValueError: on empty formula or unbalanced parentheses.

    Example::

        chematic.parse_formula("C6H12O6")  # {"C": 6, "H": 12, "O": 6}
        chematic.parse_formula("Ca(OH)2")  # {"Ca": 1, "O": 2, "H": 2}
        chematic.parse_formula("[NH4]+")   # {"N": 1, "H": 4}
    """
    ...

def scaffold_network_counts(smiles: list[str]) -> dict[str, object]:
    """Compute scaffold network statistics across a molecule library.

    Returns a dict with three parallel lists:
    - ``scaffolds``: canonical SMILES of each unique scaffold
    - ``counts``: how many input molecules contain each scaffold
    - ``parents``: index of the parent scaffold, or ``None`` for root

    Invalid SMILES are silently skipped.

    Example::

        result = chematic.scaffold_network_counts(smiles_list)
        for smi, n in zip(result["scaffolds"], result["counts"]):
            print(smi, n)
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

def find_mcs(
    mols: list[Mol],
    match_bonds: bool = True,
    min_atoms: int = 1,
    timeout_ms: Optional[int] = None,
    ring_matches_ring_only: bool = False,
    complete_rings_only: bool = False,
    atom_compare: str = "elements",
    bond_compare: str = "order_or_aromatic",
    match_chiral_tag: bool = False,
    match_charge: bool = False,
    match_isotope: bool = False,
    maximize_bonds: bool = True,
) -> Optional[Mol]:
    """Find the Maximum Common Substructure (MCS) of a list of molecules.

    Returns the MCS as a :class:`Mol`, or ``None`` when there is no common substructure.
    If ``timeout_ms`` is reached before the search finishes, returns the best result
    found so far -- indistinguishable here from an exhaustive result; use
    :func:`find_mcs_checked` when that distinction matters.

    ``atom_compare`` is one of ``"elements"`` (default), ``"any_heavy_atom"``, or ``"any"``.
    ``bond_compare`` is one of ``"order_or_aromatic"`` (default) or ``"any"``.

    Example::

        mcs = chematic.find_mcs([mol1, mol2])
        if mcs:
            print(mcs.smiles)

        # Ring-aware scaffold extraction
        scaffold = chematic.find_mcs(mols, ring_matches_ring_only=True, complete_rings_only=True)
    """
    ...

def find_mcs_checked(
    mols: list[Mol],
    match_bonds: bool = True,
    min_atoms: int = 1,
    timeout_ms: Optional[int] = None,
    ring_matches_ring_only: bool = False,
    complete_rings_only: bool = False,
    atom_compare: str = "elements",
    bond_compare: str = "order_or_aromatic",
    match_chiral_tag: bool = False,
    match_charge: bool = False,
    match_isotope: bool = False,
    maximize_bonds: bool = True,
) -> tuple[Optional[Mol], bool]:
    """Like :func:`find_mcs`, but also reports whether ``timeout_ms`` was reached
    before the search finished exhaustively.

    Returns ``(mcs, was_timed_out)``: ``mcs`` is the MCS as a :class:`Mol` (or ``None``
    if there is no common substructure), and ``was_timed_out`` is ``True`` if the search
    was cut off before proving ``mcs`` optimal.

    Example::

        mcs, timed_out = chematic.find_mcs_checked(mols, timeout_ms=500)
        if timed_out:
            print("warning: MCS may not be optimal")
    """
    ...

# ---------------------------------------------------------------------------
# PipelineV2Config / PipelineV2Error (Mol.embed_pipeline_v2)
# ---------------------------------------------------------------------------

class PipelineV2Config:
    """Configuration for :meth:`Mol.embed_pipeline_v2`.

    Every field is required -- there is no hidden default, matching the Rust
    ``PipelineV2Config``'s own deliberate lack of a ``Default`` impl (force
    field / stereo / ring-torsion policy are judgment calls). Use
    :meth:`safe` for a convenience constructor that still requires those
    three policies explicitly.

    ``stereo_policy``: one of ``"ignore"``, ``"verify_only"``,
    ``"repair_and_verify"``.

    ``ring_torsion_policy``: one of ``"fail_closed"``, ``"diagnostic_only"``.

    ``force_field_policy``: one of ``"mmff94_bond_angle_strict"``,
    ``"mmff94_with_uff_fallback"``, ``"uff_only"``, ``"dreiding"``, ``"none"``.
    """

    def __init__(
        self,
        embed_seed: int,
        max_attempts: int,
        embed_timeout_ms: Optional[int],
        use_exp_torsions: bool,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle_14_bounds: bool,
        include_legacy_torsion_heuristic: bool,
        stereo_policy: str,
        fail_on_unevaluable_stereo: bool,
        force_field_policy: str,
        force_field_max_iterations: int,
        gate_mmff94_torsion_oop: bool,
        gate_mmff94_stretch_bend: bool,
        ring_torsion_policy: str,
        total_timeout_ms: Optional[int],
        enforce_chirality: bool = False,
        expand_implicit_h_through_pipeline: bool = False,
    ) -> None: ...
    @staticmethod
    def safe(
        force_field: str,
        stereo_policy: str,
        ring_torsion_policy: str,
        fail_on_unevaluable_stereo: bool = False,
        embed_seed: int = ...,
        max_attempts: int = 8,
        embed_timeout_ms: Optional[int] = None,
        use_exp_torsions: bool = False,
        use_small_ring_torsions: bool = False,
        use_macrocycle_torsions: bool = False,
        use_macrocycle_14_bounds: bool = False,
        include_legacy_torsion_heuristic: bool = False,
        force_field_max_iterations: int = 200,
        gate_mmff94_torsion_oop: bool = False,
        gate_mmff94_stretch_bend: bool = False,
        total_timeout_ms: Optional[int] = None,
        enforce_chirality: bool = False,
        expand_implicit_h_through_pipeline: bool = False,
    ) -> "PipelineV2Config":
        """Convenience constructor.

        ``force_field``, ``stereo_policy``, and ``ring_torsion_policy`` are
        still required, explicit arguments -- never hidden defaults --
        everything else takes a conservative default (every torsion-knowledge
        flag off, ``fail_on_unevaluable_stereo=False``, no timeouts).
        """
        ...
    @staticmethod
    def stereo_safe(
        force_field: str,
        ring_torsion_policy: str,
        fail_on_unevaluable_stereo: bool = False,
        embed_seed: int = ...,
        max_attempts: int = 8,
        embed_timeout_ms: Optional[int] = None,
        use_exp_torsions: bool = False,
        use_small_ring_torsions: bool = False,
        use_macrocycle_torsions: bool = False,
        use_macrocycle_14_bounds: bool = False,
        include_legacy_torsion_heuristic: bool = False,
        force_field_max_iterations: int = 200,
        gate_mmff94_torsion_oop: bool = False,
        gate_mmff94_stretch_bend: bool = False,
        total_timeout_ms: Optional[int] = None,
    ) -> "PipelineV2Config":
        """Convenience constructor for the "stereo-safe" configuration (issue
        #291/#383): sets ``stereo_policy="repair_and_verify"``,
        ``enforce_chirality=True``, and
        ``expand_implicit_h_through_pipeline=True`` together -- the exact
        combination measured to correctly handle ring-fused declared
        stereocenters (e.g. testosterone, cholesterol) that
        ``enforce_chirality`` alone cannot repair. Prefer this over setting
        those three individually via :meth:`safe`/the constructor: they only
        work correctly as a set, and forgetting one silently falls back to a
        configuration issue #291 measured as unsound for that molecule
        class. ``force_field``/``ring_torsion_policy`` are still required,
        explicit arguments; everything else takes the same conservative
        defaults :meth:`safe` does.
        """
        ...

    embed_seed: int
    max_attempts: int
    embed_timeout_ms: Optional[int]
    use_exp_torsions: bool
    use_small_ring_torsions: bool
    use_macrocycle_torsions: bool
    use_macrocycle_14_bounds: bool
    include_legacy_torsion_heuristic: bool
    stereo_policy: str
    fail_on_unevaluable_stereo: bool
    force_field_policy: str
    force_field_max_iterations: int
    gate_mmff94_torsion_oop: bool
    gate_mmff94_stretch_bend: bool
    ring_torsion_policy: str
    total_timeout_ms: Optional[int]
    enforce_chirality: bool
    """When True, each embedding attempt is checked against declared E/Z and
    tetrahedral stereo; violations are repaired where the raw embedder's own
    repair can reach them. Compatible with every ``stereo_policy`` value,
    including ``"repair_and_verify"`` (validated as of issue #291 Step A --
    composing the two repair mechanisms measurably improves correctness with
    no observed regressions). Default ``False``, matching the Rust
    ``EmbedParameters`` default -- existing callers are unaffected."""

    expand_implicit_h_through_pipeline: bool
    """Issue #291/#383: for declared stereocenters whose only non-ring
    substituent is an implicit H (ring-fused steroid-like centers such as
    testosterone/cholesterol), run the whole pipeline on a temporary
    ``add_hydrogens``-expanded copy of the molecule instead of the original,
    then map the result back onto the original atom count before returning.
    Requires ``enforce_chirality=True`` (raises :class:`PipelineV2Error`
    otherwise). Prefer :meth:`stereo_safe` over setting this flag alone --
    it only works correctly combined with ``stereo_policy="repair_and_verify"``
    and ``enforce_chirality=True``, which ``stereo_safe`` sets together so a
    caller can't set one but forget another. Default ``False``, matching the
    Rust ``PipelineV2Config`` default -- existing callers are unaffected."""

    def __repr__(self) -> str: ...

class PipelineV2Error(ValueError):
    """A failed :meth:`Mol.embed_pipeline_v2` call.

    ``.diagnostics`` carries the same per-stage partial evidence a Rust
    caller sees on ``PipelineV2Failure`` --
    ``diagnostics["last_known_coords"]`` is diagnostic only, never a usable
    result (see ``diagnostics["coords_are_diagnostic_only"]``).
    """

    diagnostics: dict[str, object]

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

    def stereo_diagnostics(self) -> list[dict[str, object]]:
        """Rejected wedge/hash stereocenters for this record.

        A list of ``{"atom_idx": int, "reason": str}`` dicts — see
        :func:`from_mol_block_with_diagnostics` for the reason vocabulary.
        Empty unless a wedge/hash bond was present at some center and got
        rejected.
        """
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
# Crystal (chematic-crystal bindings) — Lattice / PeriodicStructure / Site
# ---------------------------------------------------------------------------

class CifSymmetryStatus:
    """How a CIF's declared symmetry relates to a PeriodicStructure's sites.

    ``is_p1=False`` means the source CIF declared symmetry beyond P1 that
    chematic's CIF adapter did not expand — the structure's ``sites`` are
    only the asymmetric unit, not a full unit cell.
    """

    is_p1: bool
    space_group_name: Optional[str]
    operation_count: int

    def __repr__(self) -> str: ...

class PeriodicNeighbor:
    """One periodic neighbor relationship from ``PeriodicStructure.neighbors()``."""

    center_index: int
    neighbor_index: int
    image: tuple[int, int, int]
    displacement: tuple[float, float, float]
    distance: float

    def __repr__(self) -> str: ...

class Site:
    """A periodic site: one or more ``(element_symbol, occupancy)`` species
    (more than one models disorder), a fractional position, and an optional
    label.
    """

    def __new__(
        cls,
        species: list[tuple[str, float]],
        fractional: tuple[float, float, float],
        label: Optional[str] = None,
    ) -> Site:
        """Construct a validated site.

        Raises:
            ValueError: unknown element symbol, non-finite/negative
                occupancy, occupancy sum over 1.0 (+ tolerance), empty
                species list, or non-finite fractional position.
        """
        ...

    @property
    def species(self) -> list[tuple[str, float]]: ...
    @property
    def fractional(self) -> tuple[float, float, float]: ...
    @property
    def label(self) -> Optional[str]: ...
    def __repr__(self) -> str: ...

class Lattice:
    """A validated 3x3 lattice matrix (rows = lattice vectors a, b, c)."""

    @staticmethod
    def from_matrix(matrix: list[list[float]]) -> Lattice: ...
    @staticmethod
    def from_parameters(
        a: float, b: float, c: float, alpha: float, beta: float, gamma: float
    ) -> Lattice: ...
    @staticmethod
    def cubic(a: float) -> Lattice: ...
    @staticmethod
    def orthorhombic(a: float, b: float, c: float) -> Lattice: ...
    @property
    def matrix(self) -> ndarray: ...
    @property
    def inverse_matrix(self) -> ndarray: ...
    @property
    def reciprocal_matrix(self) -> ndarray: ...
    @property
    def volume(self) -> float: ...
    @property
    def lengths(self) -> tuple[float, float, float]: ...
    @property
    def angles_degrees(self) -> tuple[float, float, float]: ...
    def frac_to_cart(self, point: tuple[float, float, float]) -> tuple[float, float, float]: ...
    def cart_to_frac(self, point: tuple[float, float, float]) -> tuple[float, float, float]: ...
    def __repr__(self) -> str: ...

class PeriodicStructure:
    """A periodic structure: a Lattice plus an ordered list of Sites.

    Immutable by convention — ``wrap_into_cell()``/``make_supercell()``
    return a new ``PeriodicStructure`` rather than mutating in place.

    Example::

        s = chematic.PeriodicStructure.from_cif(cif_text)
        s.lattice.volume
        s.cartesian_positions()
        s.neighbors(cutoff=3.0)
        s.make_supercell((2, 2, 2)).to_cif()
    """

    def __new__(cls, lattice: Lattice, sites: list[Site]) -> PeriodicStructure: ...
    @staticmethod
    def from_cif(text: str) -> PeriodicStructure: ...
    @staticmethod
    def from_poscar(text: str) -> PeriodicStructure: ...
    @property
    def lattice(self) -> Lattice: ...
    @property
    def sites(self) -> list[Site]: ...
    def site_count(self) -> int: ...
    def cartesian_positions(self) -> ndarray: ...
    def fractional_positions(self) -> ndarray: ...
    def neighbors(self, cutoff: float) -> list[PeriodicNeighbor]: ...
    def make_supercell(self, mult: tuple[int, int, int]) -> PeriodicStructure: ...
    def wrap_into_cell(self) -> PeriodicStructure: ...
    @property
    def symmetry_status(self) -> Optional[CifSymmetryStatus]: ...
    @property
    def formula(self) -> str: ...
    def to_cif(self) -> str: ...
    def to_poscar(self) -> str: ...
    def __repr__(self) -> str: ...

# ---------------------------------------------------------------------------
# VolumetricGrid (Gaussian Cube / OpenDX)
# ---------------------------------------------------------------------------

class VolumetricGrid:
    """A scalar field on a regular 3D grid, shared by Gaussian Cube and OpenDX.

    ``values`` is a flat numpy array in row-major, third-axis-fastest order:
    ``index = (i * shape[1] + j) * shape[2] + k``.
    """

    def __new__(
        cls,
        origin: tuple[float, float, float],
        axes: list[list[float]],
        shape: tuple[int, int, int],
        values: list[float],
        atoms: list[tuple[str, float, tuple[float, float, float]]] = [],
        units: str = "angstrom",
    ) -> VolumetricGrid: ...
    @staticmethod
    def from_cube(
        text: str,
        max_input_bytes: Optional[int] = None,
        max_atoms: Optional[int] = None,
        max_grid_points: Optional[int] = None,
    ) -> VolumetricGrid: ...
    @staticmethod
    def from_opendx(
        text: str,
        max_input_bytes: Optional[int] = None,
        max_grid_points: Optional[int] = None,
    ) -> VolumetricGrid: ...
    def to_cube(self) -> str: ...
    def to_opendx(self) -> str:
        """Raises ``ValueError`` for a Bohr-units grid or a grid with atoms."""
        ...
    def to_opendx_lossy(self) -> str: ...
    @property
    def origin(self) -> tuple[float, float, float]: ...
    @property
    def axes(self) -> ndarray: ...
    @property
    def shape(self) -> tuple[int, int, int]: ...
    @property
    def values(self) -> ndarray: ...
    @property
    def values_3d(self) -> ndarray:
        """``values`` reshaped to ``self.shape`` (``(nx, ny, nz)``)."""
        ...
    @property
    def units(self) -> str: ...
    @property
    def atoms(self) -> list[tuple[str, float, tuple[float, float, float]]]: ...
    def point_count(self) -> int: ...
    def checked_index(self, i: int, j: int, k: int) -> Optional[int]: ...
    def get(self, i: int, j: int, k: int) -> Optional[float]: ...
    def to_molecule(self) -> tuple[Mol, list[list[float]]]: ...
    def __repr__(self) -> str: ...

# ---------------------------------------------------------------------------
# LammpsDumpFrame
# ---------------------------------------------------------------------------

class LammpsDumpFrame:
    """One frame of a LAMMPS dump/trajectory file."""

    def __new__(
        cls,
        timestep: int,
        box_bounds: dict,
        column_names: list[str],
        rows: list[list[float]],
        boundary_flags: tuple[str, str, str] = ("pp", "pp", "pp"),
        num_atoms: Optional[int] = None,
    ) -> LammpsDumpFrame: ...
    @property
    def timestep(self) -> int: ...
    @property
    def num_atoms(self) -> int: ...
    @property
    def box_bounds(self) -> dict: ...
    @property
    def boundary_flags(self) -> tuple[str, str, str]: ...
    @property
    def column_names(self) -> list[str]: ...
    @property
    def rows(self) -> ndarray: ...
    def column_index(self, name: str) -> Optional[int]: ...
    def column(self, name: str) -> Optional[list[float]]: ...
    def cartesian_positions(self) -> Optional[ndarray]: ...
    def __repr__(self) -> str: ...

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
    def descriptors_array(
        smiles: list[str],
        columns: list[str],
    ) -> dict[str, ndarray]:
        """Compute descriptors and return selected columns as numpy arrays.

        Faster than ``descriptors()`` + ``pd.DataFrame()`` for column-oriented access
        because it avoids per-molecule Python dict allocation.

        Args:
            smiles: List of SMILES strings. Invalid entries are silently skipped.
            columns: Descriptor column names to return (e.g. ``["mw", "logp", "tpsa"]``).
                Float columns use ``float64``; bool columns use ``bool``; optional float
                columns (``"pka_acid"``, ``"pka_base"``) use ``float64`` with ``NaN`` for None.

        Returns:
            Dict mapping column name to 1-D numpy array.

        Raises:
            ValueError: If any column name is unknown.

        Example::

            result = chematic.bulk.descriptors_array(smiles, ["mw", "logp", "tpsa"])
            df = pd.DataFrame(result)          # fast, no per-molecule dict
            mw = result["mw"]                  # numpy.ndarray, dtype float64
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
    def map4(smiles: list[str]) -> ndarray:
        """Compute MAP4 fingerprints in parallel.

        Returns:
            numpy array of shape ``(N, 1024)``, dtype ``uint32``.
            Invalid SMILES are silently skipped.

        Note:
            Use :func:`chematic.tanimoto_map4` for similarity between MAP4 fps,
            not :func:`chematic.tanimoto` (which is bitwise).
        """
        ...

    @staticmethod
    def tanimoto_search(query: str, smiles: list[str]) -> ndarray:
        """Compute ECFP4 Tanimoto similarity of one query against a library.

        Returns:
            numpy array of shape ``(N,)``, dtype ``float32``.
        """
        ...

# ---------------------------------------------------------------------------
# Convenience functions (pure Python layer)
# ---------------------------------------------------------------------------

def from_smiles_list(
    smiles: Iterable[str],
    /,
    *,
    skip_invalid: bool = True,
) -> list[Mol]:
    """Parse a list of SMILES strings into Mol objects.

    Runs in parallel (Rayon). Invalid SMILES are silently dropped by default.

    Args:
        smiles: Iterable of SMILES strings.
        skip_invalid: If True (default), drop invalid entries.
                      If False, keep ``None`` for invalid entries.

    Returns:
        List of :class:`Mol` objects.

    Example::

        mols = chematic.from_smiles_list(["CCO", "c1ccccc1", "INVALID"])
        # → [<Mol CCO>, <Mol c1ccccc1>]
    """
    ...

def convert_format(
    text: str,
    input_format: str,
    output_format: str,
    /,
    *,
    coords: Optional[list[list[float]]] = None,
    charges: Optional[list[float]] = None,
    name: str = "LIG",
    comment: str = "",
) -> str: ...

def descriptors_df(smiles: Iterable[str]) -> Any:
    """Compute 55+ descriptors for a list of SMILES and return a DataFrame.

    Requires pandas (``pip install pandas``). Runs in parallel via Rayon.

    Args:
        smiles: Iterable of SMILES strings. Invalid entries are skipped.

    Returns:
        ``pd.DataFrame`` with one row per valid molecule and 55+ descriptor
        columns (mw, logp, tpsa, hbd, hba, qed, sa_score, pains_passes, …).

    Example::

        df = chematic.descriptors_df(["CCO", "c1ccccc1", "CC(=O)O"])
        df[["mw", "logp", "tpsa"]].head()
    """
    ...


def screen(
    smiles: Union[str, list[str]],
    profile: str = "druglike",
    filters: Optional[list[str]] = None,
) -> list[dict]:
    """Screen compounds against a preset or custom filter profile.

    Args:
        smiles: One or more SMILES strings.
        profile: Preset profile — "druglike" (default), "fragment", or "leadlike".
            Ignored when *filters* is provided.
        filters: Explicit filter list (overrides *profile*). Supported values:
            "lipinski", "veber", "pains", "brenk", "egan", "ghose", "ro3",
            "lead_like", "reos", "mcf", "ames", "pfizer_3_75", "qed", "sa_score".

    Returns:
        One dict per SMILES with fields ``smiles``, ``valid``, ``mw``, ``logp``,
        ``tpsa``, ``hbd``, ``hba``, ``qed``, ``sa_score``, one ``<filter>_pass``
        bool per requested filter, and ``overall_pass``.

    Example::

        results = chematic.screen(smiles_list, profile="druglike")
        df = pd.DataFrame(results)
        passing = df[df.overall_pass]
    """
    ...
