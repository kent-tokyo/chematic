//! RDKit-compatible 3D scalar descriptors ("Phase G1").
//!
//! This module is structurally separate from [`crate::shape_descriptors`]. That
//! module already ships functions with the *same names* (`pmi`, `npr1`,
//! `radius_of_gyration`, `asphericity`, `eccentricity`, `plane_of_best_fit`, …) —
//! but a source audit of RDKit itself (pinned commit
//! `8afba32ec539dcb2369bc84549d802aca3f7eb39`) showed several of those are
//! **not** RDKit's formula despite the matching name:
//!
//! - `shape_descriptors::eccentricity` computes `sqrt(1 - PMI1/PMI3)`; RDKit's
//!   `eccentricity` computes `sqrt(PMI3² - PMI1²) / PMI3`
//!   (`Code/GraphMol/Descriptors/PMI.cpp:188-202`) — a different tensor-ratio
//!   entirely, not just a different weighting.
//! - `shape_descriptors::asphericity` computes `PMI3 - (PMI1+PMI2)/2` from the
//!   **inertia** tensor; RDKit's `asphericity` computes
//!   `0.5·Σ(tᵢ−tⱼ)² / (Σtᵢ)²` from the **gyration** (mass-normalized
//!   covariance) tensor — a different tensor AND a different formula
//!   (`PMI.cpp:204-222`).
//! - `shape_descriptors::plane_of_best_fit` uses **heavy atoms only** and an
//!   **RMS** deviation; RDKit's `PBF` uses **all atoms including explicit H**
//!   (`ignoreHs` is hardcoded `false`) and a **mean absolute** deviation
//!   (`Code/GraphMol/Descriptors/PBF.cpp:114-152`).
//! - `shape_descriptors` has no `InertialShapeFactor` or `SpherocityIndex` at
//!   all.
//! - `descriptors_3d`'s WHIM/GETAWAY/RDF/AutoCorr3D mass-weight atoms with
//!   `Element::atomic_mass()`, which is documented (see
//!   `chematic-core/src/element.rs`) as the **monoisotopic** mass used for CIP
//!   rule 4 — not RDKit's **average atomic weight** (e.g. Cl: 34.9689 u
//!   monoisotopic vs 35.453 u average). Out of scope for this module (G1 does
//!   not touch WHIM/GETAWAY/RDF/AutoCorr3D), but the same root cause applies
//!   there and is why this module keeps its own weight table below instead of
//!   reusing `atomic_mass()`.
//!
//! Per this crate's accuracy-wave mandate, existing `shape_descriptors`/
//! `descriptors_3d` functions are left **byte-identical** (see
//! `tests/g1_regression_byte_identical.rs`), and RDKit parity is offered
//! **only** through the new `rdkit_*`-prefixed functions in this module —
//! mirroring the same never-silently-interchanged pattern `chematic-fp` uses
//! for `ecfp4()` vs. `rdkit_morgan_ecfp4_experimental()`.
//!
//! # Scope: caller-supplied coordinates only
//!
//! RDKit-compatible 3D descriptors describe supplied coordinates.
//! They do not certify the quality of the supplied conformer.
//!
//! For publication-quality conformers, docking preparation, or complex
//! macrocycles, generate and validate the conformer ensemble with a mature
//! 3D toolkit such as RDKit before passing coordinates to chematic.
//!
//! This module never calls [`crate::dg::generate_coords`] (or any other
//! chematic conformer generator) internally as a fallback — a missing or
//! malformed [`Coords3D`] is always a typed [`RdkitDescriptorError`], never a
//! silently-substituted geometry.
//!
//! # Explicit-H contract
//!
//! RDKit's own C++ implementation hardcodes `ignoreHs = false` for every G1
//! descriptor (see `computeInertiaTerms`/`computeCovarianceTerms` call sites in
//! `PMI.cpp`, and `PBF.cpp`'s unconditional `mol.getNumAtoms()` loop) — RDKit
//! always uses *every atom in the conformer*, heavy and H alike. chematic's
//! `Molecule`/`Coords3D` have no separate "implicit H" storage (see the crate
//! root docs): every atom you want counted — H included — must be an explicit
//! [`chematic_core::Atom`] with a matching [`Coords3D`] entry. To match RDKit
//! numerically, freeze coordinates from an **explicit-H** RDKit
//! conformer (post `Chem.AddHs`) and build the chematic `Molecule` with the
//! same explicit H atoms, in the same order. Passing a heavy-atom-only
//! `Molecule` is legal (no error) but will generally **not** agree with an
//! RDKit reference value computed on the explicit-H conformer, since the atom
//! sets differ.
//!
//! # Conformer index/count contract
//!
//! chematic has no multi-conformer storage: a [`Coords3D`] is exactly one
//! coordinate set, analogous to selecting a single RDKit `confId`. There is no
//! "conformer index" parameter on these functions — a caller with multiple
//! RDKit conformers calls once per conformer, each with its own frozen
//! `Coords3D` (see the multi-conformer `ibuprofen` fixture entries in
//! `validation/rdkit_3d_g1_fixtures.json`, one JSON record per RDKit `confId`).
//! There is consequently no distinct "missing conformer" error variant either:
//! an unpopulated/all-zero `Coords3D` (e.g. `Coords3D::new_zeroed` never
//! actually filled in) is geometrically indistinguishable from a molecule
//! whose atoms are all genuinely coincident at the origin, and both are
//! caught by the same [`RdkitDescriptorError::DegenerateGeometry`] path.
//!
//! # Degenerate-geometry handling (descriptor by descriptor)
//!
//! RDKit itself guards several of these formulas against division by a
//! near-zero quantity by returning a silent placeholder `0.0` (see the
//! `PMI.cpp` line citations in each function's doc comment below). That
//! placeholder is easy to mistake for a real "this molecule is a sphere"
//! answer. This module applies the **same** RDKit-sourced epsilon thresholds
//! to detect the identical guard conditions, but returns
//! `Err(RdkitDescriptorError::DegenerateGeometry)` instead of `Ok(0.0)` —
//! this is an intentional, documented departure from RDKit's silent-fallback
//! behavior, required by this PR's no-silent-fallback policy. It does not
//! change numeric agreement on any well-conditioned molecule (see the fixture
//! parity test): only inputs RDKit's own code already flags as numerically
//! degenerate are affected.
//!
//! Not every descriptor has such a guard:
//! - **PMI1/PMI2/PMI3/radius_of_gyration** are never degenerate-guarded: a
//!   single atom (or a fully coincident point cloud) has a genuinely
//!   well-defined value of `0.0` for all four (a point has zero moment of
//!   inertia / zero radius of gyration — that's correct physics, not an
//!   error). Only `ZeroAtoms`/`AtomCoordCountMismatch`/`NonFiniteCoordinate`/
//!   `UnsupportedElement` apply.
//! - **NPR1/NPR2/InertialShapeFactor/Eccentricity/Asphericity/
//!   SpherocityIndex** divide by a principal moment that can be legitimately
//!   near-zero (linear or perfectly flat geometry) — guarded per RDKit's own
//!   thresholds, see each function's doc comment.
//! - **PBF** additionally guards a case RDKit itself does *not*: if a point
//!   cloud is (near-)collinear, the "smallest-eigenvalue eigenvector" used as
//!   the plane normal is not unique (the eigenspace is 2-D), so RDKit's own
//!   solver would return an arbitrary, numerically unstable normal rather
//!   than erroring. chematic detects this (second-smallest gyration
//!   eigenvalue also near zero) and returns `DegenerateGeometry` — stricter
//!   than upstream RDKit, called out explicitly since it is a deliberate
//!   divergence.
//!
//! # Macrocycle handling
//!
//! [`RdkitDescriptorError`] is not raised for macrocycles — the geometric
//! formulas apply just as well to a big ring as to a small one. Instead, per
//! this PR's requirement to never return a "normal-looking" success value for
//! conformers this module cannot vouch for, every successful result is
//! wrapped in [`DescriptorValue`], which carries a [`MacrocycleStatus`]
//! alongside the number: [`MacrocycleStatus::Macrocyclic`] whenever the input
//! `Molecule` contains an SSSR ring at least [`MACROCYCLE_RING_THRESHOLD`]
//! atoms large, [`MacrocycleStatus::NotMacrocyclic`] otherwise. That
//! threshold (9) is RDKit's own `minMacrocycleRingSize`
//! (`Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp:38`, pinned
//! commit) — the ring size at which RDKit's *own* ETKDG embedder switches to
//! macrocycle-specific torsion sampling, i.e. RDKit's own signal that
//! "ordinary small-molecule conformer generation assumptions may not hold
//! here". The status does not change the returned value or block
//! computation — it flags that this module cannot certify the *conformer's*
//! quality for such rings (see the publication-quality-conformer boundary
//! note above).
//!
//! Ring perception needs bond connectivity, which a chematic `Molecule` is
//! not guaranteed to have (atoms-only construction is legal). Rather than
//! silently treating "no bonds to check" the same as "checked, no
//! macrocycle" — which would be exactly the kind of normal-looking-but-wrong
//! success value this section opened by ruling out —
//! [`MacrocycleStatus::Unscreenable`] is a third, explicit state for a
//! `Molecule` with atoms but zero bonds. Molecules parsed from
//! SMILES/SDF/MOL (the normal construction path) always carry bonds and so
//! never hit this state; it is reachable only via direct, bond-free
//! `MoleculeBuilder` construction from raw coordinates.

use crate::coords::Coords3D;
use crate::shape_descriptors::jacobi3;
use chematic_core::{AtomIdx, Molecule};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// SSSR ring size at/above which a ring is flagged as (potentially)
/// macrocyclic. Sourced from RDKit's own `minMacrocycleRingSize` constant,
/// `Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp:38` (pinned commit
/// `8afba32ec539dcb2369bc84549d802aca3f7eb39`) — the same threshold RDKit's
/// ETKDG embedder itself uses to switch to macrocycle-specific torsion
/// sampling.
pub const MACROCYCLE_RING_THRESHOLD: usize = 9;

/// RDKit-sourced degenerate-geometry threshold for NPR1/NPR2
/// (`PMI.cpp:116,128`).
const NPR_EPS: f64 = 1e-8;

/// RDKit-sourced degenerate-geometry threshold for InertialShapeFactor
/// (`PMI.cpp:181`), Eccentricity (`PMI.cpp:196`), Asphericity (`PMI.cpp:213`),
/// and SpherocityIndex (`PMI.cpp:232`). Also reused (chematic-specific, see
/// module docs) for PBF's collinearity guard.
const SHAPE_EPS: f64 = 1e-4;

/// Average atomic weights (u), one entry per atomic number `1..=118`
/// (index `z - 1`). Sourced from a live oracle query —
/// `rdkit.Chem.GetPeriodicTable().GetAtomicWeight(z)` against
/// `rdkit==2026.03.3` — not from chematic-core's `Element::atomic_mass()`,
/// which is the *monoisotopic* mass (see module docs). Regenerate via the
/// one-liner documented in `scripts/gen_rdkit_3d_g1_fixtures.py`.
#[rustfmt::skip]
static RDKIT_ATOMIC_WEIGHTS: [f64; 118] = [
    1.008, 4.003, 6.941, 9.012, 10.812, 12.011, 14.007, 15.999, 18.998, 20.180,
    22.990, 24.305, 26.982, 28.086, 30.974, 32.067, 35.453, 39.948, 39.098, 40.078,
    44.956, 47.867, 50.944, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.390,
    69.723, 72.610, 74.922, 78.960, 79.904, 83.800, 85.468, 87.620, 88.906, 91.224,
    92.906, 95.940, 98.000, 101.070, 102.906, 106.420, 107.868, 112.412, 114.818, 118.711,
    121.760, 127.600, 126.904, 131.290, 132.905, 137.328, 138.906, 140.116, 140.908, 144.240,
    145.000, 150.360, 151.964, 157.250, 158.925, 162.500, 164.930, 167.260, 168.934, 173.040,
    174.967, 178.490, 180.948, 183.840, 186.207, 190.230, 192.217, 195.078, 196.967, 200.590,
    204.383, 207.200, 208.980, 209.000, 210.000, 222.000, 223.000, 226.000, 227.000, 232.038,
    231.036, 238.029, 237.000, 244.000, 243.000, 247.000, 247.000, 251.000, 252.000, 257.000,
    258.000, 259.000, 262.000, 267.000, 268.000, 269.000, 270.000, 269.000, 278.000, 281.000,
    281.000, 285.000, 284.000, 289.000, 288.000, 293.000, 292.000, 294.000,
];

fn rdkit_atomic_weight(atomic_number: u8) -> Option<f64> {
    if atomic_number == 0 {
        return None;
    }
    RDKIT_ATOMIC_WEIGHTS
        .get(atomic_number as usize - 1)
        .copied()
}

// ---------------------------------------------------------------------------
// Error / result types
// ---------------------------------------------------------------------------

/// Typed failure modes for every `rdkit_*` function in this module. No
/// function in this module ever falls back to a silently-substituted
/// geometry or an unflagged placeholder value — every failure surfaces here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RdkitDescriptorError {
    /// `mol.atom_count() == 0`: nothing to compute a shape from.
    ZeroAtoms,
    /// `coords.atom_count() != mol.atom_count()`: the caller-supplied
    /// coordinate array does not describe this molecule.
    AtomCoordCountMismatch { atoms: usize, coords: usize },
    /// A coordinate component was NaN or +/-Infinity.
    NonFiniteCoordinate { atom_index: usize },
    /// The atom's element has no entry in chematic's RDKit-sourced atomic
    /// weight table (atomic number 0 or out of the supported 1..=118 range).
    /// In practice unreachable for real molecules — every real element is
    /// covered — but checked explicitly rather than assumed unreachable.
    UnsupportedElement {
        atom_index: usize,
        atomic_number: u8,
    },
    /// Fewer atoms than the descriptor requires to be meaningful (currently
    /// only PBF, which RDKit itself defines as needing >= 4 atoms;
    /// `Code/GraphMol/Descriptors/PBF.cpp:117`).
    TooFewAtoms { required: usize, found: usize },
    /// The input geometry is numerically degenerate for this specific
    /// descriptor (near-zero variance along an axis the formula divides by).
    /// See the module docs' "Degenerate-geometry handling" section for which
    /// descriptors this applies to and the exact RDKit-sourced threshold
    /// used for each.
    DegenerateGeometry { reason: &'static str },
}

impl std::fmt::Display for RdkitDescriptorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RdkitDescriptorError::ZeroAtoms => {
                write!(f, "rdkit 3d descriptor: molecule has zero atoms")
            }
            RdkitDescriptorError::AtomCoordCountMismatch { atoms, coords } => write!(
                f,
                "rdkit 3d descriptor: atom count ({atoms}) != coordinate count ({coords})"
            ),
            RdkitDescriptorError::NonFiniteCoordinate { atom_index } => write!(
                f,
                "rdkit 3d descriptor: non-finite coordinate at atom index {atom_index}"
            ),
            RdkitDescriptorError::UnsupportedElement {
                atom_index,
                atomic_number,
            } => write!(
                f,
                "rdkit 3d descriptor: atom {atom_index} has unsupported atomic number {atomic_number} (no RDKit atomic-weight entry)"
            ),
            RdkitDescriptorError::TooFewAtoms { required, found } => write!(
                f,
                "rdkit 3d descriptor: needs >= {required} atoms, found {found}"
            ),
            RdkitDescriptorError::DegenerateGeometry { reason } => {
                write!(f, "rdkit 3d descriptor: degenerate geometry: {reason}")
            }
        }
    }
}

impl std::error::Error for RdkitDescriptorError {}

/// Size (atom count) of the largest SSSR ring found, when it meets or
/// exceeds [`MACROCYCLE_RING_THRESHOLD`]. See the module docs' "Macrocycle
/// handling" section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MacrocycleWarning {
    pub max_ring_size: usize,
}

/// Result of screening a `Molecule` for a macrocyclic ring. Deliberately a
/// 3-state enum rather than `Option<MacrocycleWarning>`: ring perception
/// ([`chematic_perception::find_sssr`]) needs bond connectivity, and a
/// `Molecule` built with atoms but zero bonds (legal in chematic's data
/// model -- e.g. coordinates loaded without bond perception) would otherwise
/// silently report "not a macrocycle" when it in fact was never actually
/// checked. `Option::None` cannot distinguish "screened, no macrocycle
/// found" from "couldn't screen at all" -- exactly the kind of
/// normal-looking-but-wrong success value this PR's failure policy forbids,
/// so this type makes the distinction explicit instead.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MacrocycleStatus {
    /// Bond connectivity was available; no SSSR ring reached the threshold.
    NotMacrocyclic,
    /// Bond connectivity was available; at least one ring reached the
    /// threshold.
    Macrocyclic(MacrocycleWarning),
    /// `mol` has atoms but zero bonds, so ring perception is meaningless --
    /// macrocycle status could not be determined. Callers who build a
    /// `Molecule` purely from raw coordinates (no bond perception step) will
    /// always see this; run bond perception (or parse from SMILES/SDF/MOL,
    /// which carry bonds) first if a real macrocycle answer is needed.
    Unscreenable,
}

/// A successfully computed descriptor value. Deliberately not a bare `f64`:
/// per this PR's requirement not to return a "normal-looking" success value
/// for conformers this module cannot vouch for, every `Ok` result carries an
/// explicit macrocycle screening status alongside the number.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DescriptorValue {
    pub value: f64,
    pub macrocycle_status: MacrocycleStatus,
}

/// Screen `mol` for a macrocyclic ring, per [`MACROCYCLE_RING_THRESHOLD`].
/// Exposed publicly so callers can check this independently of any
/// particular descriptor call. See [`MacrocycleStatus`] for why this is not
/// a plain `Option`.
pub fn detect_macrocycle_status(mol: &Molecule) -> MacrocycleStatus {
    if mol.atom_count() > 1 && mol.bonds().next().is_none() {
        return MacrocycleStatus::Unscreenable;
    }
    let rings = chematic_perception::find_sssr(mol);
    let max_ring_size = rings.rings().iter().map(|r| r.len()).max().unwrap_or(0);
    if max_ring_size >= MACROCYCLE_RING_THRESHOLD {
        MacrocycleStatus::Macrocyclic(MacrocycleWarning { max_ring_size })
    } else {
        MacrocycleStatus::NotMacrocyclic
    }
}

// ---------------------------------------------------------------------------
// Shared geometry plumbing
// ---------------------------------------------------------------------------

/// Validate `mol`/`coords` and collect `(atomic_weight, [x,y,z])` per atom.
///
/// Shared by every function in this module regardless of whether it actually
/// needs mass weighting (SpherocityIndex/PBF don't) — a single validation
/// path keeps the error contract identical across all eleven descriptors,
/// and an "unsupported element" is such a narrow edge case (atomic number 0
/// or > 118; never occurs for a real molecule) that splitting the path in
/// two isn't worth the complexity.
fn collect_atoms(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<Vec<(f64, [f64; 3])>, RdkitDescriptorError> {
    let n = mol.atom_count();
    if n == 0 {
        return Err(RdkitDescriptorError::ZeroAtoms);
    }
    if coords.atom_count() != n {
        return Err(RdkitDescriptorError::AtomCoordCountMismatch {
            atoms: n,
            coords: coords.atom_count(),
        });
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let p = coords.get(idx);
        if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
            return Err(RdkitDescriptorError::NonFiniteCoordinate { atom_index: i });
        }
        let atomic_number = mol.atom(idx).element.atomic_number();
        let weight =
            rdkit_atomic_weight(atomic_number).ok_or(RdkitDescriptorError::UnsupportedElement {
                atom_index: i,
                atomic_number,
            })?;
        out.push((weight, [p.x, p.y, p.z]));
    }
    Ok(out)
}

fn centroid(atoms: &[(f64, [f64; 3])], mass_weighted: bool) -> [f64; 3] {
    let mut c = [0.0f64; 3];
    let mut wsum = 0.0f64;
    for &(m, p) in atoms {
        let w = if mass_weighted { m } else { 1.0 };
        wsum += w;
        c[0] += w * p[0];
        c[1] += w * p[1];
        c[2] += w * p[2];
    }
    [c[0] / wsum, c[1] / wsum, c[2] / wsum]
}

/// Unnormalized mass-weighted inertia tensor about the (possibly
/// mass-weighted) centroid. Mirrors RDKit's `computeInertiaTerms`
/// (`Code/GraphMol/MolTransforms/MolTransforms.cpp:129-155`) term for term.
fn inertia_tensor(atoms: &[(f64, [f64; 3])], mass_weighted: bool) -> [[f64; 3]; 3] {
    let com = centroid(atoms, mass_weighted);
    let (mut xx, mut xy, mut xz, mut yy, mut yz, mut zz) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for &(m, p) in atoms {
        let w = if mass_weighted { m } else { 1.0 };
        let (x, y, z) = (p[0] - com[0], p[1] - com[1], p[2] - com[2]);
        xx += w * (y * y + z * z);
        yy += w * (x * x + z * z);
        zz += w * (x * x + y * y);
        xy -= w * x * y;
        xz -= w * x * z;
        yz -= w * y * z;
    }
    [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]
}

/// Mass-normalized gyration (covariance) tensor about the (possibly
/// mass-weighted) centroid — divided by total weight, unlike the inertia
/// tensor above. Mirrors RDKit's `computeCovarianceTerms` with
/// `normalize=true` (`Code/GraphMol/MolTransforms/MolTransforms.cpp:72-109`).
fn gyration_tensor(atoms: &[(f64, [f64; 3])], mass_weighted: bool) -> [[f64; 3]; 3] {
    let com = centroid(atoms, mass_weighted);
    let (mut xx, mut xy, mut xz, mut yy, mut yz, mut zz) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut wsum = 0.0f64;
    for &(m, p) in atoms {
        let w = if mass_weighted { m } else { 1.0 };
        wsum += w;
        let (x, y, z) = (p[0] - com[0], p[1] - com[1], p[2] - com[2]);
        xx += w * x * x;
        xy += w * x * y;
        xz += w * x * z;
        yy += w * y * y;
        yz += w * y * z;
        zz += w * z * z;
    }
    [
        [xx / wsum, xy / wsum, xz / wsum],
        [xy / wsum, yy / wsum, yz / wsum],
        [xz / wsum, yz / wsum, zz / wsum],
    ]
}

/// Ascending inertia-tensor eigenvalues (PMI1 <= PMI2 <= PMI3), mass-weighted
/// (RDKit's `useAtomicMasses = true` default for every PMI-family
/// descriptor: `Code/GraphMol/Descriptors/PMI.h:23,30,36,41,46`).
fn pmi_moments(atoms: &[(f64, [f64; 3])]) -> [f64; 3] {
    let (evals, _) = jacobi3(inertia_tensor(atoms, true));
    [evals[0].max(0.0), evals[1].max(0.0), evals[2].max(0.0)]
}

/// Ascending gyration-tensor eigenvalues, optionally mass-weighted.
fn gyration_moments_evecs(
    atoms: &[(f64, [f64; 3])],
    mass_weighted: bool,
) -> ([f64; 3], [[f64; 3]; 3]) {
    let (evals, evecs) = jacobi3(gyration_tensor(atoms, mass_weighted));
    (
        [evals[0].max(0.0), evals[1].max(0.0), evals[2].max(0.0)],
        evecs,
    )
}

fn wrap(mol: &Molecule, value: f64) -> DescriptorValue {
    DescriptorValue {
        value,
        macrocycle_status: detect_macrocycle_status(mol),
    }
}

// ---------------------------------------------------------------------------
// Public API — PMI / NPR
// ---------------------------------------------------------------------------

/// Smallest principal moment of inertia (mass-weighted; Da·Å²).
///
/// RDKit: `Code/GraphMol/Descriptors/PMI.cpp:133-141` (`PMI1`), default
/// `useAtomicMasses = true` (`PMI.h:35-37`). Never degenerate-guarded: a
/// single atom or fully coincident point cloud legitimately has `PMI1 = 0`.
pub fn rdkit_pmi1(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<DescriptorValue, RdkitDescriptorError> {
    let atoms = collect_atoms(mol, coords)?;
    Ok(wrap(mol, pmi_moments(&atoms)[0]))
}

/// Middle principal moment of inertia (mass-weighted; Da·Å²).
///
/// RDKit: `Code/GraphMol/Descriptors/PMI.cpp:142-150` (`PMI2`).
pub fn rdkit_pmi2(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<DescriptorValue, RdkitDescriptorError> {
    let atoms = collect_atoms(mol, coords)?;
    Ok(wrap(mol, pmi_moments(&atoms)[1]))
}

/// Largest principal moment of inertia (mass-weighted; Da·Å²).
///
/// RDKit: `Code/GraphMol/Descriptors/PMI.cpp:151-159` (`PMI3`).
pub fn rdkit_pmi3(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<DescriptorValue, RdkitDescriptorError> {
    let atoms = collect_atoms(mol, coords)?;
    Ok(wrap(mol, pmi_moments(&atoms)[2]))
}

/// Normalized principal moments ratio 1 = PMI1 / PMI3.
///
/// RDKit: `Code/GraphMol/Descriptors/PMI.cpp:109-120` (`NPR1`). RDKit guards
/// `PMI3 < 1e-8` by returning `0.0`; chematic returns
/// `Err(DegenerateGeometry)` for that case instead (see module docs).
pub fn rdkit_npr1(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<DescriptorValue, RdkitDescriptorError> {
    let atoms = collect_atoms(mol, coords)?;
    let [p1, _, p3] = pmi_moments(&atoms);
    if p3 < NPR_EPS {
        return Err(RdkitDescriptorError::DegenerateGeometry {
            reason: "PMI3 < 1e-8 (RDKit PMI.cpp:116): NPR1 is 0/0-undefined, not a real ratio",
        });
    }
    Ok(wrap(mol, p1 / p3))
}

/// Normalized principal moments ratio 2 = PMI2 / PMI3.
///
/// RDKit: `Code/GraphMol/Descriptors/PMI.cpp:121-132` (`NPR2`).
pub fn rdkit_npr2(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<DescriptorValue, RdkitDescriptorError> {
    let atoms = collect_atoms(mol, coords)?;
    let [_, p2, p3] = pmi_moments(&atoms);
    if p3 < NPR_EPS {
        return Err(RdkitDescriptorError::DegenerateGeometry {
            reason: "PMI3 < 1e-8 (RDKit PMI.cpp:128): NPR2 is 0/0-undefined, not a real ratio",
        });
    }
    Ok(wrap(mol, p2 / p3))
}

// ---------------------------------------------------------------------------
// Public API — gyration-tensor descriptors
// ---------------------------------------------------------------------------

/// Mass-weighted radius of gyration (Å): `sqrt(t1 + t2 + t3)`, `tᵢ` the
/// gyration-tensor eigenvalues.
///
/// RDKit: `Code/GraphMol/Descriptors/PMI.cpp:161-171`
/// (`radiusOfGyration`). Never degenerate-guarded (a sum-of-squares square
/// root is always well-defined, including `0.0` for a single atom).
pub fn rdkit_radius_of_gyration(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<DescriptorValue, RdkitDescriptorError> {
    let atoms = collect_atoms(mol, coords)?;
    let ([t1, t2, t3], _) = gyration_moments_evecs(&atoms, true);
    Ok(wrap(mol, (t1 + t2 + t3).sqrt()))
}

/// Inertial shape factor = PMI2 / (PMI1 * PMI3) (mass-weighted).
///
/// RDKit: `Code/GraphMol/Descriptors/PMI.cpp:173-187`
/// (`inertialShapeFactor`). RDKit guards `PMI1 < 1e-4 || PMI3 < 1e-4` by
/// returning `0.0`; chematic returns `Err(DegenerateGeometry)` instead.
pub fn rdkit_inertial_shape_factor(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<DescriptorValue, RdkitDescriptorError> {
    let atoms = collect_atoms(mol, coords)?;
    let [p1, p2, p3] = pmi_moments(&atoms);
    if p1 < SHAPE_EPS || p3 < SHAPE_EPS {
        return Err(RdkitDescriptorError::DegenerateGeometry {
            reason: "PMI1 < 1e-4 or PMI3 < 1e-4 (RDKit PMI.cpp:181): planar/no-coordinate degeneracy",
        });
    }
    Ok(wrap(mol, p2 / (p1 * p3)))
}

/// Molecular eccentricity = `sqrt(PMI3² − PMI1²) / PMI3` (mass-weighted).
///
/// RDKit: `Code/GraphMol/Descriptors/PMI.cpp:188-202` (`eccentricity`). RDKit
/// guards `PMI3 < 1e-4 || (PMI3² − PMI1²) < 1e-4` by returning `0.0`; chematic
/// returns `Err(DegenerateGeometry)` instead.
pub fn rdkit_eccentricity(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<DescriptorValue, RdkitDescriptorError> {
    let atoms = collect_atoms(mol, coords)?;
    let [p1, _, p3] = pmi_moments(&atoms);
    if p3 < SHAPE_EPS || (p3 * p3 - p1 * p1) < SHAPE_EPS {
        return Err(RdkitDescriptorError::DegenerateGeometry {
            reason: "PMI3 < 1e-4 or (PMI3^2 - PMI1^2) < 1e-4 (RDKit PMI.cpp:196): near-spherical/no-coordinate degeneracy",
        });
    }
    Ok(wrap(mol, (p3 * p3 - p1 * p1).sqrt() / p3))
}

/// Molecular asphericity = `0.5 * Σ(tᵢ−tⱼ)² / (Σtᵢ)²` from the gyration
/// tensor eigenvalues (mass-weighted).
///
/// RDKit: `Code/GraphMol/Descriptors/PMI.cpp:204-222` (`asphericity`). RDKit
/// guards `t3 < 1e-4` by returning `0.0`; chematic returns
/// `Err(DegenerateGeometry)` instead.
pub fn rdkit_asphericity(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<DescriptorValue, RdkitDescriptorError> {
    let atoms = collect_atoms(mol, coords)?;
    let ([t1, t2, t3], _) = gyration_moments_evecs(&atoms, true);
    if t3 < SHAPE_EPS {
        return Err(RdkitDescriptorError::DegenerateGeometry {
            reason: "gyration t3 < 1e-4 (RDKit PMI.cpp:213): no-coordinate degeneracy",
        });
    }
    let denom = t1 + t2 + t3;
    let value = 0.5 * ((t1 - t2).powi(2) + (t1 - t3).powi(2) + (t2 - t3).powi(2)) / (denom * denom);
    Ok(wrap(mol, value))
}

/// Spherocity index = `3 * t1 / (t1+t2+t3)` from the **unweighted** gyration
/// tensor eigenvalues.
///
/// RDKit: `Code/GraphMol/Descriptors/PMI.cpp:223-238`
/// (`spherocityIndex`) — `useAtomicMasses` is hardcoded `false` (line 225),
/// unlike every other PMI-family descriptor. RDKit guards `t3 < 1e-4` by
/// returning `0.0`; chematic returns `Err(DegenerateGeometry)` instead.
pub fn rdkit_spherocity_index(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<DescriptorValue, RdkitDescriptorError> {
    let atoms = collect_atoms(mol, coords)?;
    let ([t1, t2, t3], _) = gyration_moments_evecs(&atoms, false);
    if t3 < SHAPE_EPS {
        return Err(RdkitDescriptorError::DegenerateGeometry {
            reason: "gyration t3 < 1e-4 (RDKit PMI.cpp:232): no-coordinate degeneracy",
        });
    }
    Ok(wrap(mol, 3.0 * t1 / (t1 + t2 + t3)))
}

/// Plane of Best Fit (PBF): mean **absolute** perpendicular distance of every
/// atom (heavy + H) from the least-squares plane (Å), unweighted.
///
/// RDKit: `Code/GraphMol/Descriptors/PBF.cpp:74-152` (`PBF`/`getBestFitPlane`).
/// RDKit requires `>= 4` atoms (`PBF.cpp:117`, else returns `0.0`); chematic
/// returns `Err(TooFewAtoms)` instead. RDKit does **not** guard against a
/// (near-)collinear point cloud, where the best-fit-plane normal (the
/// smallest-eigenvalue eigenvector of the unweighted gyration tensor) is not
/// unique — chematic detects this (second-smallest gyration eigenvalue also
/// `< 1e-4`) and returns `Err(DegenerateGeometry)`, a deliberate,
/// stricter-than-upstream-RDKit addition (see module docs).
pub fn rdkit_pbf(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<DescriptorValue, RdkitDescriptorError> {
    let atoms = collect_atoms(mol, coords)?;
    let n = atoms.len();
    if n < 4 {
        return Err(RdkitDescriptorError::TooFewAtoms {
            required: 4,
            found: n,
        });
    }

    let ([_g1, g2, _g3], evecs) = gyration_moments_evecs(&atoms, false);
    if g2 < SHAPE_EPS {
        return Err(RdkitDescriptorError::DegenerateGeometry {
            reason: "second-smallest gyration eigenvalue < 1e-4: point cloud is (near-)collinear, plane normal is not unique (chematic-added guard, RDKit itself does not check this)",
        });
    }

    // Smallest-eigenvalue eigenvector (ascending index 0) = plane normal.
    let normal = [evecs[0][0], evecs[1][0], evecs[2][0]];
    let denom = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();

    // Unweighted centroid (matches RDKit's `getBestFitPlane` with `weights = nullptr`).
    let origin = centroid(&atoms, false);
    let d = -(normal[0] * origin[0] + normal[1] * origin[1] + normal[2] * origin[2]);

    let mut sum_abs = 0.0f64;
    for &(_, p) in &atoms {
        let dist = (p[0] * normal[0] + p[1] * normal[1] + p[2] * normal[2] + d).abs() / denom;
        sum_abs += dist;
    }
    Ok(wrap(mol, sum_abs / n as f64))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coords::Point3;
    use chematic_core::{Atom, Element, MoleculeBuilder};

    fn mol_from_symbols(symbols: &[u8]) -> Molecule {
        let mut b = MoleculeBuilder::new();
        for &z in symbols {
            b.add_atom(Atom::new(Element::from_atomic_number(z).unwrap()));
        }
        b.build()
    }

    fn coords_from(points: &[[f64; 3]]) -> Coords3D {
        let mut c = Coords3D::new_zeroed(points.len());
        for (i, p) in points.iter().enumerate() {
            c.set(AtomIdx(i as u32), Point3::new(p[0], p[1], p[2]));
        }
        c
    }

    // --- Typed failure modes ------------------------------------------------

    #[test]
    fn zero_atoms_is_typed_error() {
        let mol = mol_from_symbols(&[]);
        let coords = coords_from(&[]);
        assert_eq!(
            rdkit_pmi1(&mol, &coords),
            Err(RdkitDescriptorError::ZeroAtoms)
        );
        assert_eq!(
            rdkit_pbf(&mol, &coords),
            Err(RdkitDescriptorError::ZeroAtoms)
        );
    }

    #[test]
    fn atom_coord_count_mismatch_is_typed_error() {
        let mol = mol_from_symbols(&[6, 6]);
        let coords = coords_from(&[[0.0, 0.0, 0.0]]); // only 1, mol has 2
        assert_eq!(
            rdkit_pmi1(&mol, &coords),
            Err(RdkitDescriptorError::AtomCoordCountMismatch {
                atoms: 2,
                coords: 1
            })
        );
    }

    #[test]
    fn nan_coordinate_is_typed_error() {
        let mol = mol_from_symbols(&[6, 6]);
        let coords = coords_from(&[[0.0, 0.0, 0.0], [1.5, f64::NAN, 0.0]]);
        assert_eq!(
            rdkit_pmi1(&mol, &coords),
            Err(RdkitDescriptorError::NonFiniteCoordinate { atom_index: 1 })
        );
    }

    #[test]
    fn infinite_coordinate_is_typed_error() {
        let mol = mol_from_symbols(&[6, 6]);
        let coords = coords_from(&[[0.0, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0]]);
        assert_eq!(
            rdkit_pbf(&mol, &coords),
            Err(RdkitDescriptorError::NonFiniteCoordinate { atom_index: 1 })
        );
    }

    #[test]
    fn single_atom_pmi_and_rog_are_zero_not_error() {
        // A single atom has a genuinely well-defined PMI/RoG of 0 -- not degenerate.
        // (Tolerance, not exact equality: the mass-weighted centroid divides by
        // the atom's own mass, so `p - com` is only ~0 to within float rounding,
        // not bit-exact zero -- squared in the tensor, that's a ~1e-30 residual.)
        let mol = mol_from_symbols(&[6]);
        let coords = coords_from(&[[1.0, 2.0, 3.0]]);
        assert!(rdkit_pmi1(&mol, &coords).unwrap().value < 1e-20);
        assert!(rdkit_pmi3(&mol, &coords).unwrap().value < 1e-20);
        assert!(rdkit_radius_of_gyration(&mol, &coords).unwrap().value < 1e-12);
    }

    #[test]
    fn single_atom_npr_is_degenerate_error() {
        // But NPR1 = PMI1/PMI3 = 0/0 for a single atom -- genuinely undefined.
        let mol = mol_from_symbols(&[6]);
        let coords = coords_from(&[[1.0, 2.0, 3.0]]);
        assert!(matches!(
            rdkit_npr1(&mol, &coords),
            Err(RdkitDescriptorError::DegenerateGeometry { .. })
        ));
        assert!(matches!(
            rdkit_asphericity(&mol, &coords),
            Err(RdkitDescriptorError::DegenerateGeometry { .. })
        ));
        assert!(matches!(
            rdkit_spherocity_index(&mol, &coords),
            Err(RdkitDescriptorError::DegenerateGeometry { .. })
        ));
    }

    #[test]
    fn all_coincident_atoms_is_degenerate_error() {
        let mol = mol_from_symbols(&[6, 6, 6, 6]);
        let coords = coords_from(&[[1.0, 1.0, 1.0]; 4]);
        assert!(matches!(
            rdkit_eccentricity(&mol, &coords),
            Err(RdkitDescriptorError::DegenerateGeometry { .. })
        ));
        assert!(matches!(
            rdkit_pbf(&mol, &coords),
            Err(RdkitDescriptorError::DegenerateGeometry { .. })
        ));
    }

    #[test]
    fn collinear_atoms_pbf_is_degenerate_error() {
        // 4 collinear points: a valid inertia/PMI shape (a rod) but PBF's
        // plane normal is not unique (chematic-added stricter guard).
        let mol = mol_from_symbols(&[6, 6, 6, 6]);
        let coords = coords_from(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ]);
        assert!(matches!(
            rdkit_pbf(&mol, &coords),
            Err(RdkitDescriptorError::DegenerateGeometry { .. })
        ));
        // But PMI/NPR are perfectly well-defined for a rod.
        assert!(rdkit_npr1(&mol, &coords).unwrap().value < 1e-9);
    }

    #[test]
    fn pbf_too_few_atoms_is_typed_error() {
        let mol = mol_from_symbols(&[6, 6, 6]);
        let coords = coords_from(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        assert_eq!(
            rdkit_pbf(&mol, &coords),
            Err(RdkitDescriptorError::TooFewAtoms {
                required: 4,
                found: 3
            })
        );
    }

    // --- Invariance -----------------------------------------------------------
    // Not all descriptors behave the same way -- see per-block comments.

    fn tetrahedron() -> (Molecule, Coords3D) {
        // 4 distinct, non-degenerate, non-symmetric points -- a real 3D shape.
        let mol = mol_from_symbols(&[6, 7, 8, 9]);
        let coords = coords_from(&[
            [0.0, 0.0, 0.0],
            [1.3, 0.0, 0.0],
            [0.2, 1.1, 0.3],
            [0.4, 0.3, 1.7],
        ]);
        (mol, coords)
    }

    fn translate(coords: &Coords3D, t: [f64; 3]) -> Coords3D {
        let mut out = Coords3D::new_zeroed(coords.atom_count());
        for i in 0..coords.atom_count() {
            let p = coords.get(AtomIdx(i as u32));
            out.set(
                AtomIdx(i as u32),
                Point3::new(p.x + t[0], p.y + t[1], p.z + t[2]),
            );
        }
        out
    }

    // Rotate all points by a fixed axis-angle rotation (Rodrigues' formula)
    // about an arbitrary axis, so the test isn't limited to axis-aligned
    // rotations (which could hide a bug that only shows up off-axis).
    fn rotate(coords: &Coords3D) -> Coords3D {
        let axis = {
            let n = (1.0_f64 + 4.0 + 9.0).sqrt();
            [1.0 / n, 2.0 / n, 3.0 / n]
        };
        let theta = 0.7_f64;
        let (s, c) = theta.sin_cos();
        let mut out = Coords3D::new_zeroed(coords.atom_count());
        for i in 0..coords.atom_count() {
            let p = coords.get(AtomIdx(i as u32));
            let v = [p.x, p.y, p.z];
            let dot = axis[0] * v[0] + axis[1] * v[1] + axis[2] * v[2];
            let cross = [
                axis[1] * v[2] - axis[2] * v[1],
                axis[2] * v[0] - axis[0] * v[2],
                axis[0] * v[1] - axis[1] * v[0],
            ];
            let r = [
                v[0] * c + cross[0] * s + axis[0] * dot * (1.0 - c),
                v[1] * c + cross[1] * s + axis[1] * dot * (1.0 - c),
                v[2] * c + cross[2] * s + axis[2] * dot * (1.0 - c),
            ];
            out.set(AtomIdx(i as u32), Point3::new(r[0], r[1], r[2]));
        }
        out
    }

    fn reflect_x(coords: &Coords3D) -> Coords3D {
        let mut out = Coords3D::new_zeroed(coords.atom_count());
        for i in 0..coords.atom_count() {
            let p = coords.get(AtomIdx(i as u32));
            out.set(AtomIdx(i as u32), Point3::new(-p.x, p.y, p.z));
        }
        out
    }

    fn permute(mol: &Molecule, coords: &Coords3D, order: &[usize]) -> (Molecule, Coords3D) {
        let mut b = MoleculeBuilder::new();
        let mut out_coords = Coords3D::new_zeroed(order.len());
        for (new_i, &old_i) in order.iter().enumerate() {
            b.add_atom(mol.atom(AtomIdx(old_i as u32)).clone());
            out_coords.set(AtomIdx(new_i as u32), coords.get(AtomIdx(old_i as u32)));
        }
        (b.build(), out_coords)
    }

    const TOL: f64 = 1e-9;

    #[test]
    fn translation_invariant_for_all_scalar_descriptors() {
        let (mol, c) = tetrahedron();
        let c2 = translate(&c, [37.0, -12.0, 5.0]);
        assert!(
            (rdkit_pmi1(&mol, &c).unwrap().value - rdkit_pmi1(&mol, &c2).unwrap().value).abs()
                < TOL
        );
        assert!(
            (rdkit_pmi3(&mol, &c).unwrap().value - rdkit_pmi3(&mol, &c2).unwrap().value).abs()
                < TOL
        );
        assert!(
            (rdkit_npr1(&mol, &c).unwrap().value - rdkit_npr1(&mol, &c2).unwrap().value).abs()
                < TOL
        );
        assert!(
            (rdkit_radius_of_gyration(&mol, &c).unwrap().value
                - rdkit_radius_of_gyration(&mol, &c2).unwrap().value)
                .abs()
                < TOL
        );
        assert!(
            (rdkit_asphericity(&mol, &c).unwrap().value
                - rdkit_asphericity(&mol, &c2).unwrap().value)
                .abs()
                < TOL
        );
        assert!(
            (rdkit_eccentricity(&mol, &c).unwrap().value
                - rdkit_eccentricity(&mol, &c2).unwrap().value)
                .abs()
                < TOL
        );
        assert!(
            (rdkit_inertial_shape_factor(&mol, &c).unwrap().value
                - rdkit_inertial_shape_factor(&mol, &c2).unwrap().value)
                .abs()
                < TOL
        );
        assert!(
            (rdkit_spherocity_index(&mol, &c).unwrap().value
                - rdkit_spherocity_index(&mol, &c2).unwrap().value)
                .abs()
                < TOL
        );
        assert!(
            (rdkit_pbf(&mol, &c).unwrap().value - rdkit_pbf(&mol, &c2).unwrap().value).abs() < TOL
        );
    }

    #[test]
    fn rotation_invariant_for_all_scalar_descriptors() {
        // Every G1 descriptor here is a function of the inertia/gyration
        // tensor's EIGENVALUES only (never raw eigenvector components), so
        // all of them -- including PBF's scalar mean-distance value -- are
        // rotation invariant.
        let (mol, c) = tetrahedron();
        let c2 = rotate(&c);
        assert!(
            (rdkit_pmi1(&mol, &c).unwrap().value - rdkit_pmi1(&mol, &c2).unwrap().value).abs()
                < 1e-6
        );
        assert!(
            (rdkit_pmi2(&mol, &c).unwrap().value - rdkit_pmi2(&mol, &c2).unwrap().value).abs()
                < 1e-6
        );
        assert!(
            (rdkit_pmi3(&mol, &c).unwrap().value - rdkit_pmi3(&mol, &c2).unwrap().value).abs()
                < 1e-6
        );
        assert!(
            (rdkit_radius_of_gyration(&mol, &c).unwrap().value
                - rdkit_radius_of_gyration(&mol, &c2).unwrap().value)
                .abs()
                < 1e-6
        );
        assert!(
            (rdkit_asphericity(&mol, &c).unwrap().value
                - rdkit_asphericity(&mol, &c2).unwrap().value)
                .abs()
                < 1e-6
        );
        assert!(
            (rdkit_pbf(&mol, &c).unwrap().value - rdkit_pbf(&mol, &c2).unwrap().value).abs() < 1e-6
        );
    }

    #[test]
    fn atom_permutation_invariant() {
        let (mol, c) = tetrahedron();
        let (mol2, c2) = permute(&mol, &c, &[3, 1, 0, 2]);
        assert!(
            (rdkit_pmi1(&mol, &c).unwrap().value - rdkit_pmi1(&mol2, &c2).unwrap().value).abs()
                < TOL
        );
        assert!(
            (rdkit_npr2(&mol, &c).unwrap().value - rdkit_npr2(&mol2, &c2).unwrap().value).abs()
                < TOL
        );
        assert!(
            (rdkit_pbf(&mol, &c).unwrap().value - rdkit_pbf(&mol2, &c2).unwrap().value).abs() < TOL
        );
        assert!(
            (rdkit_spherocity_index(&mol, &c).unwrap().value
                - rdkit_spherocity_index(&mol2, &c2).unwrap().value)
                .abs()
                < TOL
        );
    }

    #[test]
    fn reflection_does_not_change_any_g1_scalar() {
        // Every G1 descriptor here is built purely from tensor eigenVALUES
        // (never a signed/chiral combination of eigenvectors), so a mirror
        // reflection leaves every one of them unchanged -- unlike, say, a
        // signed dihedral or a CIP descriptor, none of these 11 values are
        // chirality-sensitive. This is a real property of these *specific*
        // formulas, not a blanket assumption: PBF in particular takes
        // |distance| (abs value), so a normal-vector sign flip under
        // reflection cancels out.
        let (mol, c) = tetrahedron();
        let c2 = reflect_x(&c);
        assert!(
            (rdkit_pmi1(&mol, &c).unwrap().value - rdkit_pmi1(&mol, &c2).unwrap().value).abs()
                < TOL
        );
        assert!(
            (rdkit_pmi2(&mol, &c).unwrap().value - rdkit_pmi2(&mol, &c2).unwrap().value).abs()
                < TOL
        );
        assert!(
            (rdkit_pmi3(&mol, &c).unwrap().value - rdkit_pmi3(&mol, &c2).unwrap().value).abs()
                < TOL
        );
        assert!(
            (rdkit_eccentricity(&mol, &c).unwrap().value
                - rdkit_eccentricity(&mol, &c2).unwrap().value)
                .abs()
                < TOL
        );
        assert!(
            (rdkit_pbf(&mol, &c).unwrap().value - rdkit_pbf(&mol, &c2).unwrap().value).abs() < TOL
        );
    }

    #[test]
    fn coordinate_precision_round_trip() {
        // Re-reading the same f64 coordinates back out and recomputing must
        // reproduce bit-identical results (no hidden internal rounding/state).
        let (mol, c) = tetrahedron();
        let r1 = rdkit_pmi1(&mol, &c).unwrap().value;
        let r2 = rdkit_pmi1(&mol, &c).unwrap().value;
        assert_eq!(r1.to_bits(), r2.to_bits());
        let a1 = rdkit_asphericity(&mol, &c).unwrap().value;
        let a2 = rdkit_asphericity(&mol, &c).unwrap().value;
        assert_eq!(a1.to_bits(), a2.to_bits());
    }

    // --- Macrocycle screening -----------------------------------------------

    #[test]
    fn macrocycle_ring_sets_warning_flag() {
        // 12-membered carbocycle (cyclododecane): real bonds needed for SSSR,
        // so parse from SMILES. Coordinates here are chematic's OWN
        // conformer generator -- fine for this qualitative flag-only check
        // (no RDKit numeric comparison happens in this test), but never used
        // for the RDKit-agreement fixture (see tests/rdkit_g1_parity.rs).
        let mol = chematic_smiles::parse("C1CCCCCCCCCCC1").unwrap();
        let coords = crate::dg::generate_coords(&mol);
        let status = detect_macrocycle_status(&mol);
        assert_eq!(
            status,
            MacrocycleStatus::Macrocyclic(MacrocycleWarning { max_ring_size: 12 })
        );

        let result = rdkit_radius_of_gyration(&mol, &coords).unwrap();
        assert!(
            matches!(result.macrocycle_status, MacrocycleStatus::Macrocyclic(_)),
            "successful descriptor result must carry the macrocycle status, not silently drop it"
        );
    }

    #[test]
    fn small_ring_has_no_macrocycle_warning() {
        let mol = chematic_smiles::parse("C1CCCCC1").unwrap(); // cyclohexane, 6-ring
        let coords = crate::dg::generate_coords(&mol);
        assert_eq!(
            detect_macrocycle_status(&mol),
            MacrocycleStatus::NotMacrocyclic
        );
        let result = rdkit_radius_of_gyration(&mol, &coords).unwrap();
        assert_eq!(result.macrocycle_status, MacrocycleStatus::NotMacrocyclic);
    }

    #[test]
    fn non_macrocycle_acyclic_molecule_has_no_warning() {
        let mol = chematic_smiles::parse("CCO").unwrap();
        let coords = crate::dg::generate_coords(&mol);
        assert_eq!(
            rdkit_pmi1(&mol, &coords).unwrap().macrocycle_status,
            MacrocycleStatus::NotMacrocyclic
        );
    }

    #[test]
    fn bond_free_molecule_is_unscreenable_not_silently_not_macrocyclic() {
        // A Molecule built with atoms but no bonds (e.g. our own
        // `mol_from_symbols` test helper, or any caller who constructs a
        // Molecule directly from raw coordinates without a bond-perception
        // step) cannot be ring-perceived at all. Reporting `NotMacrocyclic`
        // here would be exactly the silent, normal-looking-but-wrong success
        // value this module's failure policy forbids -- even for a 12-atom
        // ring's worth of atoms, laid out on an actual ring geometry.
        let mol = mol_from_symbols(&[6; 12]);
        let mut coords = Coords3D::new_zeroed(12);
        for i in 0..12u32 {
            let theta = std::f64::consts::TAU * (i as f64) / 12.0;
            coords.set(AtomIdx(i), Point3::new(theta.cos(), theta.sin(), 0.0));
        }
        assert_eq!(
            detect_macrocycle_status(&mol),
            MacrocycleStatus::Unscreenable
        );
        let result = rdkit_radius_of_gyration(&mol, &coords).unwrap();
        assert_eq!(result.macrocycle_status, MacrocycleStatus::Unscreenable);
    }
}
