//! MOL V2000 (Ctab) parser and writer.
//!
//! Reference format:
//!   Line 1  — molecule name (may be blank)
//!   Line 2  — program/date info (may be blank)
//!   Line 3  — comment (may be blank)
//!   Line 4  — counts line: fixed-width fields for atom count, bond count, version tag
//!   Lines 5..5+natoms — atom block (one line per atom)
//!   Lines 5+natoms..5+natoms+nbonds — bond block (one line per bond)
//!   "M  END" — molecule terminator

use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Chirality, Coords3D, Element, Molecule, MoleculeBuilder,
    Point3, STEREO_H_SENTINEL, SquarePlanarPermutation, StereoGeometry,
};
use chematic_perception::{
    EzDirectionDiagnostic, StereoDiagnostic, apply_ez_directions_from_2d_ex,
    apply_local_parity_from_wedges_with_diagnostics,
};

use crate::error::MolParseError;

/// Maximum number of atoms allowed in a MOL file (prevents memory exhaustion).
const MAX_ATOMS: usize = 100_000;

/// Maximum number of bonds allowed in a MOL file (prevents memory exhaustion).
const MAX_BONDS: usize = 200_000;

/// Metadata extracted from the three-line MOL header.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct MolMetadata {
    /// Molecule name from header line 1.
    pub name: String,
    /// Comment string from header line 3.
    pub comment: String,
}

impl MolMetadata {
    /// Set the molecule name and return `self` (builder style).
    pub fn with_name(mut self, name: &str) -> Self {
        self.name = name.to_owned();
        self
    }

    /// Set the comment and return `self` (builder style).
    pub fn with_comment(mut self, comment: &str) -> Self {
        self.comment = comment.to_owned();
        self
    }
}

/// The MOL header's own declared dimensionality -- the literal `2D`/`3D`
/// tag a writer stamps in columns 20..22 (0-indexed) of header line 2 (the
/// "program/date" line, MDL Ctfile spec), e.g. RDKit's
/// `"     RDKit          3D"`. Empirically confirmed against a live RDKit
/// `2026.03.3` oracle for both V2000 and V3000 (identical column in both).
///
/// This is the file's own CLAIM, independent of what the atom block's
/// actual `(x, y, z)` values look like -- see [`GeometryRank`] for the
/// observed side, and [`Stereo3DDiagnostic`] for whether the two agree.
/// `Unknown` is the common case: most real-world writers (including every
/// hand-written fixture already in this crate's own test suite) leave line
/// 2 blank or shorter than 22 columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CoordinateDimension {
    TwoD,
    ThreeD,
    #[default]
    Unknown,
}

/// What the atom block's actual `(x, y, z)` values look like, independent of
/// what the header declares. Deliberately kept separate from
/// [`CoordinateDimension`] (the file's claim) -- collapsing "does this file
/// have 3D data" to one boolean would conflate "header says 3D but is
/// mislabeled" with "header says 3D and the molecule is just, correctly,
/// flat" (e.g. benzene from a real conformer generator).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeometryRank {
    /// No atoms at all -- not enough information to say anything.
    Indeterminate,
    /// Every atom's z is within [`Z_EPS`] of exactly 0.0 -- the trivial flat
    /// case, indistinguishable from "no z data was ever written".
    FlatZero,
    /// Not all-zero-z, but every point still lies within [`COPLANAR_EPS`] of
    /// a single best-fit plane (which need not be the z=0 plane, and need
    /// not be axis-aligned) -- includes fewer-than-3-atom records, which
    /// are trivially coplanar. A genuinely flat molecule (e.g. benzene)
    /// legitimately lands here; on its own this is not evidence of a bug.
    Coplanar,
    /// At least one point lies measurably off the best-fit plane through
    /// the rest -- genuinely three-dimensional.
    ThreeD,
}

/// Coordinate magnitude below which a z value is treated as exactly zero.
const Z_EPS: f64 = 1e-4;
/// Maximum perpendicular deviation (Angstrom) from a best-fit plane still
/// classified as [`GeometryRank::Coplanar`]. Chosen well above float noise
/// and well below a real 3D ring pucker (~0.2-0.5 A for e.g. cyclohexane).
const COPLANAR_EPS: f64 = 1e-2;
/// Minimum cross-product norm (Angstrom^2) accepted as "not collinear" when
/// searching for a plane-defining pair of vectors.
const COLLINEAR_EPS: f64 = 1e-6;
/// Minimum |signed volume| (Angstrom^3) treated as a reliable, non-degenerate
/// tetrahedral sign in [`wedge_vs_3d_conflicts`].
const VOLUME_EPS: f64 = 1e-6;
/// Minimum bond-vector norm (Angstrom) accepted as non-degenerate by
/// [`classify_square_planar_geometry`] -- below this, a neighbor's position
/// coincides with (or is implausibly close to) the center's own position,
/// which can never be a real bond length for any element. Many orders of
/// magnitude below any real M-L bond (~1.5-2.5 Angstrom for a transition
/// metal): a pure corrupt-data guard, not a chemistry judgment call.
const SQUARE_PLANAR_MIN_BOND_NORM: f64 = 1e-3;
/// `cos(135 degrees)` -- the geometric midpoint between an ideal
/// square-planar cis angle (90 degrees) and trans angle (180 degrees), used
/// by [`classify_square_planar_geometry`] to decide whether a neighbor pair
/// (viewed from the center) is "trans-like". Self-justifying: over 45
/// degrees of distortion from either ideal is required to misclassify,
/// comfortably covering real (non-idealized) square-planar geometry while
/// still rejecting a tetrahedral center's ~109.5 degree angles outright
/// (already excluded earlier by the coplanarity check, since a real
/// tetrahedral center's 4 substituents are never coplanar with it -- this
/// threshold only has to disambiguate *which* pairing is trans for input
/// that already passed that gate). Chosen the same way `COPLANAR_EPS`/
/// `VOLUME_EPS` were: a deliberately-derived, stated-reasoning constant, not
/// a copied or guessed one. See `docs/rfcs/square_planar_mol_io_rfc.md` §7,
/// which also cites an RDKit 2026.03.4 oracle observation (empirical cutoff
/// between 30 and 40 degrees of distortion) in the same family as this
/// value, without copying RDKit's number.
const SQUARE_PLANAR_TRANS_COS_MAX: f64 = -std::f64::consts::FRAC_1_SQRT_2;

/// Classify a point cloud's [`GeometryRank`]. See that type's docs for the
/// four cases.
pub(crate) fn classify_geometry_rank(points: &[Point3]) -> GeometryRank {
    if points.is_empty() {
        return GeometryRank::Indeterminate;
    }
    if points.iter().all(|p| p.z.abs() < Z_EPS) {
        return GeometryRank::FlatZero;
    }
    if points.len() < 3 {
        return GeometryRank::Coplanar;
    }
    // Find a non-degenerate plane-defining pair of vectors from point 0.
    // ponytail: fixing the origin at point 0 makes this O(n) instead of
    // O(n^3); it misses only the vanishingly rare pathological case where
    // point 0 is collinear with every other point yet some other, unrelated
    // triple isn't -- not worth a full O(n^3) search for real molecular
    // input up to `MAX_ATOMS`.
    let origin = points[0];
    let mut normal: Option<Point3> = None;
    'outer: for (j, pj) in points.iter().enumerate().skip(1) {
        let v1 = pj.sub(&origin);
        for pk in points.iter().skip(j + 1) {
            let v2 = pk.sub(&origin);
            let n = v1.cross(&v2);
            if n.norm() > COLLINEAR_EPS {
                normal = Some(n);
                break 'outer;
            }
        }
    }
    let Some(n) = normal else {
        // Every point is collinear with point 0 -- trivially coplanar.
        return GeometryRank::Coplanar;
    };
    let n = n.normalize();
    let max_dev = points
        .iter()
        .map(|p| p.sub(&origin).dot(&n).abs())
        .fold(0.0_f64, f64::max);
    if max_dev < COPLANAR_EPS {
        GeometryRank::Coplanar
    } else {
        GeometryRank::ThreeD
    }
}

/// Parse the dimensional-code field (columns 20..22, 0-indexed) from a MOL
/// header's line 2. Returns [`CoordinateDimension::Unknown`] when the line
/// is shorter than 22 columns, blank in that field, or contains anything
/// other than the literal `2D`/`3D` tokens -- never an error (this field is
/// legacy/optional and most writers, including this crate's own 2D writer,
/// never populate it).
pub(crate) fn parse_dimension_code(line2: &str) -> CoordinateDimension {
    match line2.get(20..22).map(str::trim) {
        Some("2D") => CoordinateDimension::TwoD,
        Some("3D") => CoordinateDimension::ThreeD,
        _ => CoordinateDimension::Unknown,
    }
}

/// One 3D-geometry-related diagnostic surfaced while parsing a MOL/SDF
/// record -- see [`MolReadReport::stereo3d_diagnostics`].
///
/// An empty `stereo3d_diagnostics` vec does NOT mean "3D stereo was verified
/// correct" -- it means there was nothing to check (no wedge/hash bond was
/// present at any center, and the header/geometry dimensionality agreed),
/// mirroring how `stereo_diagnostics`/`ez_diagnostics` are only ever
/// populated when something was actually declared and rejected. A record
/// with genuine 3D coordinates and zero declared stereo is common and
/// entirely valid -- it produces no diagnostics, which must not be read as
/// "stereo check passed".
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Stereo3DDiagnostic {
    /// The header declares `2D`, but the atom block's actual z values are
    /// not all (near-)zero.
    DeclaredTwoDButNonzeroZ { observed: GeometryRank },
    /// The header declares `3D`, but the observed geometry is flat -- either
    /// every z is exactly 0 ([`GeometryRank::FlatZero`]) or all points lie
    /// on a single, possibly tilted, plane ([`GeometryRank::Coplanar`]). The
    /// two are reported distinctly via the payload: "every z is literally
    /// 0" is near-certainly a 2D file mislabeled 3D, while a real but
    /// planar/degenerate 3D geometry (e.g. benzene) is not on its own a bug.
    DeclaredThreeDButFlat { observed: GeometryRank },
    /// A wedge/hash bond was present at `atom`, and its 2D-perceived local
    /// parity (`Atom.chirality`, computed on this record's own `(x, y)` by
    /// the unmodified [`apply_local_parity_from_wedges_with_diagnostics`]
    /// pathway) disagrees with the parity read directly off the atom's real
    /// 3D geometry (`from_3d_geometry`). Named for what is actually being
    /// compared: for a genuinely-3D record, `(x, y)` alone is not a
    /// meaningful 2D depiction, so this is "the wedge disagrees with the
    /// shape", not "the declared stereo failed verification". When a wedge
    /// and nonzero-z geometry coexist and *agree*, nothing is emitted here
    /// -- 3D geometry is treated as the higher-fidelity signal, but it is
    /// never written back to `Atom.chirality` (which stays exactly what the
    /// 2D wedge pathway produced); only genuine disagreement is surfaced.
    WedgeVs3DParityConflict {
        atom: AtomIdx,
        wedge_2d: Chirality,
        from_3d_geometry: Chirality,
    },
}

/// Signed volume of the tetrahedron `(p1-p4, p2-p4, p3-p4)` from real 3D
/// coordinates. An independent, deliberately tiny (single determinant)
/// re-implementation of the same triple product
/// `chematic_perception::stereo2d_local` uses internally on a synthetic
/// wedge-derived z -- that helper is `pub(crate)` there and its z is not a
/// real coordinate, so it cannot be reused directly from this crate.
fn signed_volume3(p1: Point3, p2: Point3, p3: Point3, p4: Point3) -> f64 {
    let a = p1.sub(&p4);
    let b = p2.sub(&p4);
    let c = p3.sub(&p4);
    a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) + a.z * (b.x * c.y - b.y * c.x)
}

/// Compare every wedge/hash-perceived tetrahedral parity against the same
/// atom's real 3D geometry, using the exact sign convention
/// `chematic_perception::stereo2d_local::local_parity_from_wedges`
/// establishes (mirrored here from its doc comments, not re-derived): apex =
/// first neighbor in `stereo_neighbor_order`, viewed = the rest (in order),
/// for a 4-explicit-neighbor center; pivot = the center atom itself (no
/// synthetic H position needed -- the triple product of the 3 real bond
/// vectors from the center already carries full parity) for a 3-heavy +
/// implicit-H center. Only ever consulted when a wedge/hash bond was
/// actually present (`atom.chirality.is_tetrahedral()`), so this never
/// fires on the common case of a real 3D SDF record with no wedge notation
/// at all.
///
/// Gated on [`chematic_core::Chirality::is_tetrahedral`], not `!=
/// Chirality::None`: the computation below assumes a tetrahedral shape (it
/// only ever produces `Chirality::Clockwise`/`CounterClockwise` and compares
/// that against `atom.chirality`), so any non-tetrahedral-but-non-`None`
/// chirality (e.g. `Chirality::SquarePlanar`, produced by this crate's own
/// MOL/SDF readers as of this file's `perceive_square_planar_from_3d`) must
/// be excluded up front rather than coerced through a check that can only
/// ever disagree with it. `is_tetrahedral()` is an allowlist (`true` only
/// for the two known-tetrahedral variants), which makes this fix
/// exhaustive-match safe by construction: any future non-tetrahedral
/// geometry this crate adds later (trigonal-bipyramidal, octahedral --
/// sketched but unimplemented in `chematic_core::stereo_geometry`) is
/// automatically excluded too, without a new arm here. Same
/// equality-vs-exhaustive-match bug shape fixed twice before in this
/// project's history (`chematic-3d/src/stereo_constraints.rs`,
/// `chematic-chem/src/cip.rs`).
pub(crate) fn wedge_vs_3d_conflicts(
    mol: &Molecule,
    conformer: &Coords3D,
) -> Vec<Stereo3DDiagnostic> {
    let mut out = Vec::new();
    let get = |a: u32| conformer.points.get(a as usize).copied();

    for (idx, atom) in mol.atoms() {
        if !atom.chirality.is_tetrahedral() {
            continue;
        }
        let Some(order) = mol.stereo_neighbor_order(idx) else {
            continue;
        };

        let computed = if order.len() == 4 && order[3] == STEREO_H_SENTINEL {
            // 3 heavy neighbors + 1 implicit H: pivot = center.
            match (get(order[0]), get(order[1]), get(order[2]), get(idx.0)) {
                (Some(p0), Some(p1), Some(p2), Some(center)) => {
                    let v = signed_volume3(p0, p1, p2, center);
                    match v {
                        v if v.abs() < VOLUME_EPS => None,
                        v if v < 0.0 => Some(Chirality::Clockwise),
                        _ => Some(Chirality::CounterClockwise),
                    }
                }
                _ => None,
            }
        } else if order.len() == 4 {
            // 4 explicit neighbors: apex = order[0], viewed = order[1..4].
            match (get(order[0]), get(order[1]), get(order[2]), get(order[3])) {
                (Some(p0), Some(p1), Some(p2), Some(p3)) => {
                    let v = signed_volume3(p1, p2, p3, p0);
                    match v {
                        v if v.abs() < VOLUME_EPS => None,
                        v if v < 0.0 => Some(Chirality::CounterClockwise),
                        _ => Some(Chirality::Clockwise),
                    }
                }
                _ => None,
            }
        } else {
            None
        };

        if let Some(computed) = computed
            && computed != atom.chirality
        {
            out.push(Stereo3DDiagnostic::WedgeVs3DParityConflict {
                atom: idx,
                wedge_2d: atom.chirality,
                from_3d_geometry: computed,
            });
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Square-planar stereo: 3D-coordinate-derived reperception (read) and
// pre-write validation (write) -- see
// `docs/rfcs/square_planar_mol_io_rfc.md` for the full design and RDKit
// oracle provenance. MDL/CTfile has no symbolic field for a non-tetrahedral
// stereo tag (RFC §2); this is the *only* mechanism that can represent it.
// ---------------------------------------------------------------------------

/// Why a coplanar, undefined-valence, 4-coordinate center's neighbor
/// arrangement failed to resolve to exactly one of SP1/SP2/SP3 in
/// [`classify_square_planar_geometry`]. Used both as a read-side diagnostic
/// reason ([`SquarePlanarPerceptionDiagnostic`]) and, wrapped, as a write-side
/// error reason ([`UnsupportedStereoReason::GeometryRejected`]) -- one
/// classifier, one reason vocabulary, shared by both directions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SquarePlanarRejectionReason {
    /// A neighbor's position coincides with (or is implausibly close to) the
    /// center's own position -- see [`SQUARE_PLANAR_MIN_BOND_NORM`].
    DegenerateBondVector,
    /// The center and its 4 neighbors are not coplanar within tolerance --
    /// the common case for a genuinely tetrahedral/octahedral/other-shaped
    /// 4-coordinate center. Never surfaced as a read-side diagnostic (see
    /// [`perceive_square_planar_from_3d`]'s doc comment) since it is the
    /// expected, non-noteworthy outcome for the vast majority of real
    /// undefined-valence-element 4-coordinate centers.
    NotCoplanar,
    /// Coplanar, but no single pairing of neighbors had *both* its pairs at
    /// a trans-like angle (see [`SQUARE_PLANAR_TRANS_COS_MAX`]) -- or more
    /// than one pairing did, which can only happen for near-degenerate
    /// input. Either way, the arrangement does not unambiguously name one
    /// of SP1/SP2/SP3.
    AmbiguousTransPairing,
}

impl core::fmt::Display for SquarePlanarRejectionReason {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::DegenerateBondVector => write!(
                f,
                "a neighbor position coincides with the center (degenerate bond vector)"
            ),
            Self::NotCoplanar => write!(f, "center and neighbors are not coplanar"),
            Self::AmbiguousTransPairing => {
                write!(f, "no single trans-pairing of neighbors was unambiguous")
            }
        }
    }
}

/// Classify a candidate square-planar center's local geometry: `center` plus
/// its 4 neighbor positions, in [`Molecule::neighbors`]/`stereo_neighbor_order`
/// order (position `i` of the returned tag's
/// [`SquarePlanarPermutation::trans_pairs`] refers to `neighbors[i]`).
///
/// The single shared implementation used by both
/// [`perceive_square_planar_from_3d`] (read) and
/// [`validate_square_planar_for_write`] (write) -- see
/// `docs/rfcs/square_planar_mol_io_rfc.md` §7 for the full derivation of each
/// step's tolerance.
pub(crate) fn classify_square_planar_geometry(
    center: Point3,
    neighbors: [Point3; 4],
) -> Result<SquarePlanarPermutation, SquarePlanarRejectionReason> {
    let vecs: [Point3; 4] = [
        neighbors[0].sub(&center),
        neighbors[1].sub(&center),
        neighbors[2].sub(&center),
        neighbors[3].sub(&center),
    ];
    let norms: [f64; 4] = [
        vecs[0].norm(),
        vecs[1].norm(),
        vecs[2].norm(),
        vecs[3].norm(),
    ];
    if norms.iter().any(|&n| n < SQUARE_PLANAR_MIN_BOND_NORM) {
        return Err(SquarePlanarRejectionReason::DegenerateBondVector);
    }

    // Coplanarity of center + 4 neighbors, reusing `classify_geometry_rank`
    // directly rather than a second plane-fit implementation. `FlatZero` is
    // accepted here (it is the *local* 5-point subset that may legitimately
    // sit at z=0 even inside a molecule whose overall conformer is real 3D
    // data -- the whole-molecule flat/indeterminate case is rejected
    // separately, by the caller checking `conformer.is_some()` on read and
    // by `validate_square_planar_for_write` on write).
    let five = [
        center,
        neighbors[0],
        neighbors[1],
        neighbors[2],
        neighbors[3],
    ];
    match classify_geometry_rank(&five) {
        GeometryRank::Coplanar | GeometryRank::FlatZero => {}
        GeometryRank::ThreeD | GeometryRank::Indeterminate => {
            return Err(SquarePlanarRejectionReason::NotCoplanar);
        }
    }

    let cos_angle = |i: usize, j: usize| vecs[i].dot(&vecs[j]) / (norms[i] * norms[j]);
    let is_trans = |i: u8, j: u8| cos_angle(i as usize, j as usize) <= SQUARE_PLANAR_TRANS_COS_MAX;

    let mut matched: Option<SquarePlanarPermutation> = None;
    for tag in [
        SquarePlanarPermutation::SP1,
        SquarePlanarPermutation::SP2,
        SquarePlanarPermutation::SP3,
    ] {
        let [(a, b), (c, d)] = tag.trans_pairs();
        if is_trans(a, b) && is_trans(c, d) {
            if matched.is_some() {
                return Err(SquarePlanarRejectionReason::AmbiguousTransPairing);
            }
            matched = Some(tag);
        }
    }
    matched.ok_or(SquarePlanarRejectionReason::AmbiguousTransPairing)
}

/// One candidate square-planar center whose geometry looked plausible
/// (coplanar) but did not unambiguously resolve to SP1/SP2/SP3 -- surfaced
/// instead of silently leaving `Chirality::None`, matching this crate's
/// existing "explain why, don't guess" diagnostic discipline
/// (`StereoDiagnostic`/`Stereo3DDiagnostic`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SquarePlanarPerceptionDiagnostic {
    pub atom: AtomIdx,
    pub reason: SquarePlanarRejectionReason,
}

/// Scan `mol` for undefined-valence-element, 4-coordinate, not-already-tagged
/// atoms and, using `conformer`'s real 3D positions, reperceive
/// `Chirality::SquarePlanar` directly from geometry -- see
/// `docs/rfcs/square_planar_mol_io_rfc.md` for why this (not a symbolic MOL
/// field, which does not exist) is the mechanism.
///
/// Element eligibility reuses `Element::normal_valences().is_empty()`
/// (transition metals and other elements this crate already treats as
/// "valence undefined") -- reproduces an RDKit 2026.03.4 oracle observation
/// for every element tested, with one documented divergence (RFC §6).
///
/// Only [`SquarePlanarRejectionReason::DegenerateBondVector`]/
/// `AmbiguousTransPairing` are surfaced as diagnostics.
/// `NotCoplanar` is silently skipped -- it is the expected, non-noteworthy
/// outcome for the overwhelming majority of real 4-coordinate
/// undefined-valence-element centers (genuinely tetrahedral or octahedral,
/// not square-planar), mirroring `StereoDiagnostic`'s own "no wedge present
/// -> nothing to say" precedent rather than warning on every ordinary
/// tetrahedral transition-metal complex in a 3D SDF file.
pub(crate) fn perceive_square_planar_from_3d(
    mol: &mut Molecule,
    conformer: &Coords3D,
) -> Vec<SquarePlanarPerceptionDiagnostic> {
    let mut out = Vec::new();

    // Collect candidates first (immutable borrow of `mol`) before mutating.
    let candidates: Vec<(AtomIdx, [AtomIdx; 4])> = mol
        .atoms()
        .filter(|(_, atom)| {
            atom.chirality == Chirality::None && atom.element.normal_valences().is_empty()
        })
        .filter_map(|(idx, _)| {
            let nbs: Vec<AtomIdx> = mol.neighbors(idx).map(|(n, _)| n).collect();
            (nbs.len() == 4).then(|| (idx, [nbs[0], nbs[1], nbs[2], nbs[3]]))
        })
        .collect();

    for (idx, nbs) in candidates {
        let Some(center) = conformer.points.get(idx.0 as usize).copied() else {
            continue;
        };
        let mut pts = [Point3::zero(); 4];
        let mut all_present = true;
        for (slot, n) in nbs.iter().enumerate() {
            match conformer.points.get(n.0 as usize).copied() {
                Some(p) => pts[slot] = p,
                None => {
                    all_present = false;
                    break;
                }
            }
        }
        if !all_present {
            continue;
        }

        match classify_square_planar_geometry(center, pts) {
            Ok(tag) => {
                mol.set_chirality(idx, Chirality::SquarePlanar(tag));
                mol.set_stereo_neighbor_order(idx, nbs.iter().map(|a| a.0).collect());
            }
            Err(SquarePlanarRejectionReason::NotCoplanar) => {}
            Err(reason) => out.push(SquarePlanarPerceptionDiagnostic { atom: idx, reason }),
        }
    }

    out
}

/// Result of parsing a MOL/SDF record with stereo-perception diagnostics.
///
/// `stereo_diagnostics` is empty unless a wedge/hash bond was actually
/// present at some center and got rejected (contradictory, missing
/// coordinates, degenerate geometry, or an unsupported neighbor shape) --
/// see [`chematic_perception::StereoDiagnostic`]. It is never populated for
/// an atom with no wedge/hash bond at all.
///
/// `ez_diagnostics` is the E/Z (double-bond cis/trans) counterpart: empty
/// unless a stereogenic double bond (terminal alkenes, carbonyls, and
/// symmetric-substituent alkenes are never stereogenic in the first place)
/// had its direction rejected -- see
/// [`chematic_perception::EzDirectionDiagnostic`].
///
/// `conformer` is `Some` exactly when the atom block's real z values are not
/// all (near-)zero (i.e. `geometry_rank` is [`GeometryRank::Coplanar`] or
/// [`GeometryRank::ThreeD`]) -- this is "does this file actually have a 3D
/// conformer", independent of `coordinate_dimension` (the header's own,
/// possibly wrong or absent, claim). `coordinate_dimension` and
/// `geometry_rank` are always populated; `stereo3d_diagnostics` cross-checks
/// the two and checks any wedge/hash-declared stereo against real geometry
/// -- see [`Stereo3DDiagnostic`].
///
/// `square_planar_diagnostics` is populated only when `conformer` is `Some`
/// (square-planar reperception never runs against a flat/absent conformer --
/// see `docs/rfcs/square_planar_mol_io_rfc.md`), and even then only for a
/// candidate center whose geometry looked plausible but didn't
/// unambiguously resolve to SP1/SP2/SP3; a successfully-perceived center is
/// reflected directly in `mol`'s `Atom.chirality`/`stereo_neighbor_order`,
/// not listed here.
#[derive(Clone)]
pub struct MolReadReport {
    pub mol: Molecule,
    pub metadata: MolMetadata,
    pub coords: Vec<(f64, f64)>,
    pub stereo_diagnostics: Vec<StereoDiagnostic>,
    pub ez_diagnostics: Vec<EzDirectionDiagnostic>,
    pub conformer: Option<Coords3D>,
    pub coordinate_dimension: CoordinateDimension,
    pub geometry_rank: GeometryRank,
    pub stereo3d_diagnostics: Vec<Stereo3DDiagnostic>,
    pub square_planar_diagnostics: Vec<SquarePlanarPerceptionDiagnostic>,
}

// ---------------------------------------------------------------------------
// Charge encoding table (V2000 ccc field → formal charge)
// ---------------------------------------------------------------------------

/// Decode a V2000 charge code into a formal charge value.
fn decode_charge(code: i8) -> i8 {
    match code {
        1 => 3,
        2 => 2,
        3 => 1,
        4 => 0, // doublet radical — treated as neutral
        5 => -1,
        6 => -2,
        7 => -3,
        _ => 0,
    }
}

/// Encode a formal charge into a V2000 charge code.
fn encode_charge(charge: i8) -> u8 {
    match charge {
        3 => 1,
        2 => 2,
        1 => 3,
        -1 => 5,
        -2 => 6,
        -3 => 7,
        _ => 0,
    }
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a fixed-width 3-character integer field from a string slice.
///
/// Returns an error using `make_err` when the slice is missing or the text
/// cannot be parsed as an integer.
fn parse_field3(
    line: &str,
    start: usize,
    line_num: usize,
    make_err: impl Fn(usize, String) -> MolParseError,
) -> Result<usize, MolParseError> {
    let field = line
        .get(start..start + 3)
        .ok_or_else(|| make_err(line_num, format!("line too short at column {start}")))?;
    field
        .trim()
        .parse::<usize>()
        .map_err(|_| make_err(line_num, format!("cannot parse integer from '{field}'")))
}

/// Parse a MOL V2000 string, running stereo perception and returning every
/// rejected wedge/hash center as a structured [`StereoDiagnostic`].
///
/// The parser follows the MDL/CTfile fixed-width column layout. `coords[i]`
/// is the `(x, y)` position for atom `i` extracted from the atom block.
/// Z-coordinates are discarded. This is the one parsing core for V2000 MOL
/// text -- [`parse_mol_with_coords`]/[`parse_mol`] are thin wrappers that
/// discard `stereo_diagnostics`, and [`crate::sdf`]'s readers delegate here
/// per record.
///
/// Local tetrahedral parity (`Atom.chirality` + `Molecule::stereo_neighbor_order`)
/// is perceived unconditionally whenever a wedge/hash bond is present --
/// mirroring RDKit's own `assignChiralTypesFromBondDirs`, which runs
/// regardless of a `sanitize`-equivalent flag. It never touches `Atom.cip_code`
/// and never depends on CIP ranking.
pub fn read_mol_with_diagnostics(input: &str) -> Result<MolReadReport, MolParseError> {
    read_mol_internal(input, true)
}

/// Parse a MOL V2000 block without optional stereo/3D diagnostics.
///
/// This is used by the RDKit-compatible `SDMolSupplier` path, whose contract
/// is the molecule graph plus metadata and SD properties. Keeping the parser
/// and fixed-width validation shared with [`read_mol_with_diagnostics`] avoids
/// semantic drift while skipping perception and geometry work.
#[allow(clippy::type_complexity)]
pub(crate) fn parse_mol_fast(input: &str) -> Result<(Molecule, MolMetadata), MolParseError> {
    let report = read_mol_internal(input, false)?;
    Ok((report.mol, report.metadata))
}

fn read_mol_internal(
    input: &str,
    include_diagnostics: bool,
) -> Result<MolReadReport, MolParseError> {
    // Yields (1-based line number, line text); short-circuits on EOF.
    let mut lines = input.lines().enumerate().map(|(i, l)| (i + 1, l));
    let mut next_line = || lines.next().ok_or(MolParseError::UnexpectedEnd);

    // -- Header block: lines 1–3 -------------------------------------------

    let name = next_line()?.1.to_string();
    // Line 2 (program/date info) is mostly discarded, except for the
    // dimensional-code field (columns 20..22, "2D"/"3D") -- see
    // `parse_dimension_code`. Empirically confirmed against a live RDKit
    // 2026.03.3 oracle: `"     RDKit          3D"`/`"...2D"`.
    let (_, line2_raw) = next_line()?;
    let comment = next_line()?.1.to_string();

    let metadata = MolMetadata { name, comment };

    // -- Counts line (line 4) -----------------------------------------------

    let (counts_lineno, counts_line) = next_line()?;

    // Be lenient with shorter lines — just check the V2000 tag exists.
    if !counts_line.contains("V2000") {
        return Err(MolParseError::InvalidCountLine {
            line: counts_lineno,
            detail: "missing V2000 version tag".to_string(),
        });
    }

    let make_count_err = |ln: usize, d: String| MolParseError::InvalidCountLine {
        line: ln,
        detail: d,
    };

    let natoms = parse_field3(counts_line, 0, counts_lineno, make_count_err)?;
    let nbonds = parse_field3(counts_line, 3, counts_lineno, make_count_err)?;

    if natoms > MAX_ATOMS {
        return Err(MolParseError::InvalidCountLine {
            line: counts_lineno,
            detail: format!(
                "atom count {} exceeds maximum allowed {}",
                natoms, MAX_ATOMS
            ),
        });
    }

    if nbonds > MAX_BONDS {
        return Err(MolParseError::InvalidCountLine {
            line: counts_lineno,
            detail: format!(
                "bond count {} exceeds maximum allowed {}",
                nbonds, MAX_BONDS
            ),
        });
    }

    // -- Atom block ---------------------------------------------------------

    let mut builder = MoleculeBuilder::with_capacity(natoms, nbonds);
    let mut coords: Vec<(f64, f64)> = if include_diagnostics {
        Vec::with_capacity(natoms)
    } else {
        Vec::new()
    };
    let mut raw_z: Vec<f64> = if include_diagnostics {
        Vec::with_capacity(natoms)
    } else {
        Vec::new()
    };
    let make_atom_err = |ln: usize, d: String| MolParseError::InvalidAtomLine {
        line: ln,
        detail: d,
    };

    for atom_i in 0..natoms {
        let (raw_lineno, atom_line) = next_line()?;

        // Keep validating Z even on the graph-only supplier path. A malformed
        // value was already a typed error there, so skipping all coordinate
        // parsing would silently loosen strictParsing. Only X/Y conversion
        // and coordinate storage are diagnostic-path work.
        let z: f64 = match atom_line.get(20..30) {
            None => 0.0,
            Some(raw) => {
                let trimmed = raw.trim();
                if trimmed.is_empty() {
                    0.0
                } else {
                    let val: f64 = trimmed.parse().map_err(|_| {
                        make_atom_err(
                            raw_lineno,
                            format!("cannot parse z coordinate from '{trimmed}'"),
                        )
                    })?;
                    if !val.is_finite() {
                        return Err(make_atom_err(
                            raw_lineno,
                            format!("z coordinate is not finite (NaN/Infinite): '{trimmed}'"),
                        ));
                    }
                    val
                }
            }
        };

        if include_diagnostics {
            // X/Y are only needed by coordinate and stereo diagnostics. The
            // graph-only supplier does not expose either, so avoid converting
            // and storing values it immediately discards.
            let x: f64 = atom_line
                .get(0..10)
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0.0);
            let y: f64 = atom_line
                .get(10..20)
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0.0);
            coords.push((x, y));
            raw_z.push(z);
        }

        // Element symbol: bytes 31–33 (3 chars, left-padded with a space in
        // the spec, but writers vary; trim both ends).
        let sym = atom_line
            .get(31..34)
            .ok_or_else(|| {
                make_atom_err(
                    raw_lineno,
                    format!("atom line {atom_i} too short for element field"),
                )
            })?
            .trim();

        let element = Element::from_symbol(sym).ok_or_else(|| MolParseError::UnknownElement {
            symbol: sym.to_string(),
            line: raw_lineno,
        })?;

        // Charge code: bytes 36–38 (3 chars).
        let charge = atom_line
            .get(36..39)
            .map(|ccc| decode_charge(ccc.trim().parse().unwrap_or(0)))
            .unwrap_or(0);

        let mut atom = Atom::new(element);
        atom.charge = charge;
        builder.add_atom(atom);
    }

    // -- Bond block ---------------------------------------------------------

    let make_bond_err = |ln: usize, d: String| MolParseError::InvalidBondLine {
        line: ln,
        detail: d,
    };

    // Double bonds whose stereo field is 3 ("cis or trans -- either"), MDL's
    // crossed-bond convention for explicitly unspecified E/Z. This is a
    // THIRD, distinct state from a resolved direction, exactly analogous to
    // single-bond code 4 above -- confirmed against a live RDKit 2026.03.3
    // oracle (B0 diagnosis): RDKit sets `Bond::STEREOANY` and never prints
    // `/`/`\`, even when raw coordinates would otherwise resolve to a
    // definite E or Z. There is no per-bond Molecule-level field for this
    // (see `docs/rfcs/stereo2d_reader_integration_rfc.md` for why one wasn't
    // added); it is threaded directly into
    // `apply_ez_directions_from_2d_ex` below instead.
    let mut explicitly_unspecified_ez: std::collections::HashSet<BondIdx> =
        std::collections::HashSet::new();

    for bond_i in 0..nbonds {
        let (raw_lineno, bond_line) = next_line()?;

        let a1_raw = parse_field3(bond_line, 0, raw_lineno, make_bond_err)?;
        let a2_raw = parse_field3(bond_line, 3, raw_lineno, make_bond_err)?;
        let btype_raw = parse_field3(bond_line, 6, raw_lineno, make_bond_err)?;

        if a1_raw == 0 || a2_raw == 0 {
            return Err(MolParseError::InvalidBondLine {
                line: raw_lineno,
                detail: format!("bond {bond_i}: atom indices are 1-based; got {a1_raw}/{a2_raw}"),
            });
        }

        let a1 = AtomIdx((a1_raw - 1) as u32);
        let a2 = AtomIdx((a2_raw - 1) as u32);

        // Stereo field (columns 9-11, 0-indexed): only meaningful for single bonds.
        let stereo_raw: usize = if bond_line.len() >= 12 {
            parse_field3(bond_line, 9, raw_lineno, make_bond_err).unwrap_or(0)
        } else {
            0
        };

        let order = match btype_raw {
            // Code 4 ("either"/unspecified direction) is deliberately NOT
            // folded into `Up`: it is a third, distinct MDL state (RDKit
            // maps it to `Bond::BondDir::UNKNOWN`, never a definite wedge)
            // and once stereo perception reads `BondOrder::Up` as a
            // confident wedge, conflating the two would fabricate a
            // stereocenter from a bond whose direction the file explicitly
            // declares unknown. Falls through to `Single`, i.e. "no defined
            // direction" -- a documented, accepted round-trip-lossy case
            // (see `docs/rfcs/stereo2d_reader_integration_rfc.md`).
            1 => match stereo_raw {
                1 => BondOrder::Up,
                6 => BondOrder::Down,
                _ => BondOrder::Single,
            },
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            4 => BondOrder::Aromatic,
            5 => BondOrder::QuerySingleOrDouble,
            6 => BondOrder::QuerySingleOrAromatic,
            7 => BondOrder::QueryDoubleOrAromatic,
            8 => BondOrder::QueryAny,
            // 9 = dative/coordinate bond, a widely used (if formally
            // non-standard) MDL extension -- RDKit emits it for
            // `Bond::BondType::DATIVE` (V3000 only; see mol3000.rs), with
            // atom1/atom2 in the same donor/acceptor order as
            // `BondOrder::Dative` itself already documents. Previously fell
            // through to `Single`, silently discarding the coordination-bond
            // distinction (and, downstream, corrupting the donor atom's
            // implicit hydrogen count) on read -- see valence.rs.
            9 => BondOrder::Dative,
            _ => BondOrder::Single,
        };

        let bidx = builder
            .add_bond(a1, a2, order)
            .map_err(|e| MolParseError::InvalidBondLine {
                line: raw_lineno,
                detail: format!("bond {bond_i}: {e}"),
            })?;

        if btype_raw == 2 && stereo_raw == 3 {
            explicitly_unspecified_ez.insert(bidx);
        }
    }

    // Skip property lines until "M  END" (or EOF if absent).
    for (_, l) in lines.by_ref() {
        if l.trim_start().starts_with("M  END") {
            break;
        }
    }

    let mut mol = builder.build();

    if !include_diagnostics {
        return Ok(MolReadReport {
            mol,
            metadata,
            coords,
            stereo_diagnostics: Vec::new(),
            ez_diagnostics: Vec::new(),
            conformer: None,
            coordinate_dimension: CoordinateDimension::Unknown,
            geometry_rank: GeometryRank::Indeterminate,
            stereo3d_diagnostics: Vec::new(),
            square_planar_diagnostics: Vec::new(),
        });
    }

    // Tetrahedral parity first (raw wedge/hash still fully intact on
    // `bond.order`), THEN E/Z direction -- the E/Z stage only ever writes to
    // the separate `bond_direction` side channel, never to `bond.order`, so
    // running it after can never disturb the wedge/hash the tetrahedral
    // stage just read. See `docs/rfcs/stereo2d_reader_integration_rfc.md`.
    let stereo_diagnostics = apply_local_parity_from_wedges_with_diagnostics(&mut mol, &coords);
    let ez_diagnostics =
        apply_ez_directions_from_2d_ex(&mut mol, &coords, &explicitly_unspecified_ez);

    let coordinate_dimension = parse_dimension_code(line2_raw);
    let points: Vec<Point3> = coords
        .iter()
        .zip(raw_z.iter())
        .map(|(&(x, y), &z)| Point3::new(x, y, z))
        .collect();
    let geometry_rank = classify_geometry_rank(&points);
    let conformer = match geometry_rank {
        GeometryRank::Coplanar | GeometryRank::ThreeD => Some(Coords3D { points }),
        GeometryRank::FlatZero | GeometryRank::Indeterminate => None,
    };

    let mut stereo3d_diagnostics = Vec::new();
    match (coordinate_dimension, geometry_rank) {
        (CoordinateDimension::TwoD, GeometryRank::Coplanar | GeometryRank::ThreeD) => {
            stereo3d_diagnostics.push(Stereo3DDiagnostic::DeclaredTwoDButNonzeroZ {
                observed: geometry_rank,
            });
        }
        (CoordinateDimension::ThreeD, GeometryRank::FlatZero | GeometryRank::Coplanar) => {
            stereo3d_diagnostics.push(Stereo3DDiagnostic::DeclaredThreeDButFlat {
                observed: geometry_rank,
            });
        }
        _ => {}
    }
    // Square-planar reperception BEFORE `wedge_vs_3d_conflicts`: it may set
    // `Chirality::SquarePlanar` on some atoms, and that function's own
    // `is_tetrahedral()` gate (see its doc comment) must see the final
    // chirality to correctly skip them.
    let mut square_planar_diagnostics = Vec::new();
    if let Some(ref conf) = conformer {
        square_planar_diagnostics = perceive_square_planar_from_3d(&mut mol, conf);
        stereo3d_diagnostics.extend(wedge_vs_3d_conflicts(&mol, conf));
    }

    Ok(MolReadReport {
        mol,
        metadata,
        coords,
        stereo_diagnostics,
        ez_diagnostics,
        conformer,
        coordinate_dimension,
        geometry_rank,
        stereo3d_diagnostics,
        square_planar_diagnostics,
    })
}

/// Parse a MOL V2000 string into a `(Molecule, MolMetadata, coords)` triple.
///
/// Thin wrapper around [`read_mol_with_diagnostics`] that discards
/// `stereo_diagnostics` -- signature and behavior unchanged from before
/// stereo perception was wired in, except that a wedge/hash bond now
/// populates `Atom.chirality`/`Molecule::stereo_neighbor_order` where it
/// previously left them unset.
#[allow(clippy::type_complexity)]
pub fn parse_mol_with_coords(
    input: &str,
) -> Result<(Molecule, MolMetadata, Vec<(f64, f64)>), MolParseError> {
    read_mol_with_diagnostics(input).map(|r| (r.mol, r.metadata, r.coords))
}

/// Parse a MOL V2000 string into a `(Molecule, MolMetadata)` pair.
///
/// This is a convenience wrapper around [`parse_mol_with_coords`] that discards
/// the 2D coordinate data.
pub fn parse_mol(input: &str) -> Result<(Molecule, MolMetadata), MolParseError> {
    parse_mol_with_coords(input).map(|(mol, meta, _coords)| (mol, meta))
}

/// Parse all molecules from an SDF string, running stereo perception and
/// returning every rejected wedge/hash center per record.
///
/// Stops and returns an error on the first parse failure.
pub fn read_sdf_with_diagnostics(input: &str) -> Result<Vec<MolReadReport>, MolParseError> {
    let mut result = Vec::new();
    let mut remaining = input;
    loop {
        // A blank first line is a legal MOL name line (issue #171) -- do not
        // skip it. A genuinely empty gap between/after `$$$$` delimiters is
        // already handled below via `block.trim().is_empty()`.
        if remaining.is_empty() {
            break;
        }

        // Find the $$$$ delimiter (line-by-line to avoid false matches inside data).
        let mut byte_offset = 0usize;
        let (end_byte, after_delim) = loop {
            let rest = &remaining[byte_offset..];
            match rest.find('\n') {
                Some(nl) => {
                    let line = rest[..nl].trim_end_matches('\r');
                    if line == "$$$$" {
                        break (byte_offset, &remaining[byte_offset + nl + 1..]);
                    }
                    byte_offset += nl + 1;
                }
                None => {
                    if rest.trim_end_matches('\r') == "$$$$" {
                        break (byte_offset, "");
                    }
                    break (remaining.len(), "");
                }
            }
        };

        let block = &remaining[..end_byte];
        remaining = after_delim;
        if block.trim().is_empty() {
            continue;
        }

        result.push(read_mol_with_diagnostics(block)?);
    }
    Ok(result)
}

/// Parse all molecules from an SDF string, returning 2D coordinates.
///
/// Each entry contains the molecule, its metadata, and a `Vec<(x, y)>` of
/// 2D coordinates in atom-insertion order (the same order as `.atoms()`).
///
/// Thin wrapper around [`read_sdf_with_diagnostics`] that discards
/// `stereo_diagnostics`. Stops and returns an error on the first parse failure.
#[allow(clippy::type_complexity)]
pub fn parse_sdf_with_coords(
    input: &str,
) -> Result<Vec<(Molecule, MolMetadata, Vec<(f64, f64)>)>, MolParseError> {
    Ok(read_sdf_with_diagnostics(input)?
        .into_iter()
        .map(|r| (r.mol, r.metadata, r.coords))
        .collect())
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Write a `Molecule` to MOL V2000 format.
///
/// Coordinates are written as 0.0 because the core `Molecule` type does not
/// store 2D/3D coordinates.  All other atom and bond fields are derived from
/// the molecule graph.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo.** This is a 2D-only
/// writer with no z channel; MOL/CTfile has no other field for a
/// non-tetrahedral stereo tag either (see
/// `docs/rfcs/square_planar_mol_io_rfc.md`). A square-planar-tagged atom is
/// written with no indication anything was dropped. If `mol` may carry
/// `Chirality::SquarePlanar`, use [`write_mol_with_conformer_checked`]
/// instead, which fails closed with a typed error rather than silently
/// discarding the tag.
pub fn write_mol(mol: &Molecule, metadata: &MolMetadata) -> String {
    write_mol_with_coords(mol, metadata, &[])
}

/// Serialize `mol` to a V2000 MOL block, using `coords` for atom positions.
///
/// `coords[i]` is the `(x, y)` position in Ångström for atom index `i`.
/// Atoms beyond `coords.len()` receive `(0.0, 0.0, 0.0)`.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo** -- see
/// [`write_mol`]'s doc comment; the same 2D-only limitation applies here.
/// Use [`write_mol_with_conformer_checked`] if `mol` may carry it.
#[allow(clippy::write_with_newline)]
pub fn write_mol_with_coords(
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
) -> String {
    use std::fmt::Write as _;

    // A MOL block is dominated by fixed-width atom and bond rows.  Reserve
    // the common-size output up front so serialization does not repeatedly
    // grow and copy the String for every row (the old empty String was a
    // noticeable cost when writing many small records).
    let mut out = String::with_capacity(128 + mol.atom_count() * 80 + mol.bond_count() * 16);

    // Header lines 1–3
    out.push_str(&metadata.name);
    out.push('\n');
    out.push_str("  chematic\n");
    out.push_str(&metadata.comment);
    out.push('\n');

    // Counts line (line 4)
    let natoms = mol.atom_count();
    let nbonds = mol.bond_count();
    write!(
        &mut out,
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        natoms, nbonds
    )
    .expect("writing to String cannot fail");

    // Atom block
    for (idx, atom) in mol.atoms() {
        let sym = atom.element.symbol();
        let charge_code = encode_charge(atom.charge);
        if let Some(&(x, y)) = coords.get(idx.0 as usize) {
            write!(
                &mut out,
                "{:>10.4}{:>10.4}{:>10.4} {:<3} 0{:>3}  0  0  0  0  0  0  0  0  0\n",
                x, y, 0.0_f64, sym, charge_code,
            )
            .expect("writing to String cannot fail");
        } else {
            // Serialization-only SDF output overwhelmingly has no coordinate
            // array. Avoid running the float formatter three times per atom
            // for a byte sequence that is known at compile time.
            out.push_str("    0.0000    0.0000    0.0000 ");
            writeln!(
                &mut out,
                "{:<3} 0{:>3}  0  0  0  0  0  0  0  0  0",
                sym, charge_code,
            )
            .expect("writing to String cannot fail");
        }
    }

    // Bond block
    for (_idx, bond) in mol.bonds() {
        let a1 = bond.atom1.0 + 1; // convert to 1-based
        let a2 = bond.atom2.0 + 1;
        let btype = match bond.order {
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            BondOrder::QuerySingleOrDouble => 5,
            BondOrder::QuerySingleOrAromatic => 6,
            BondOrder::QueryDoubleOrAromatic => 7,
            BondOrder::QueryAny | BondOrder::Zero => 8,
            BondOrder::Quadruple => 4,
        };
        // Stereo field: preserve a wedge/hash bond so a re-parse recovers
        // the same local parity. This is the only channel MOL/SDF has for
        // recovered stereo -- `Atom.chirality` itself has no direct MOL
        // field, so what round-trips is the wedge bond it was derived from.
        let stereo = match bond.order {
            BondOrder::Up => 1,
            BondOrder::Down => 6,
            _ => 0,
        };
        writeln!(&mut out, "{:>3}{:>3}{:>3}{:>3}", a1, a2, btype, stereo)
            .expect("writing to String cannot fail");
    }

    // Terminator
    out.push_str("M  END\n");

    out
}

/// Serialize `mol` to MOL V2000 format using `conformer`'s real 3D
/// coordinates, stamping the header's line-2 dimensional code as `3D`.
///
/// This is the writer counterpart of [`read_mol_with_diagnostics`]'s new 3D
/// support -- round-tripping this output back through the reader reproduces
/// [`CoordinateDimension::ThreeD`] and never manufactures a fresh
/// [`Stereo3DDiagnostic::WedgeVs3DParityConflict`] on its own output: no
/// wedge/hash stereo field is ever emitted here (always `0`), matching the
/// fact that a real 3D geometry makes a 2D wedge symbol redundant at best,
/// contradictory at worst -- see [`write_mol_with_coords`] (unchanged, still
/// the 2D writer) for the wedge-preserving counterpart. Atoms beyond
/// `conformer.atom_count()` receive `(0.0, 0.0, 0.0)`.
///
/// **Does not validate `Chirality::SquarePlanar` stereo against
/// `conformer`** -- it writes whatever coordinates it is given, trusting the
/// caller (the same "trust the caller" posture this crate's writers have
/// always had for tetrahedral wedge bonds). If `conformer`'s geometry
/// doesn't actually match a declared square-planar tag, this will silently
/// write a self-inconsistent file. Use
/// [`write_mol_with_conformer_checked`] instead to fail closed on that
/// mismatch (or on a missing/flat conformer) rather than trusting it.
pub fn write_mol_with_conformer(
    mol: &Molecule,
    metadata: &MolMetadata,
    conformer: &Coords3D,
) -> String {
    let mut out = String::new();

    out.push_str(&metadata.name);
    out.push('\n');
    // Columns 0..10 "  chematic" match the 2D writer exactly; columns
    // 10..20 are the (blank) date field; columns 20..22 carry the "3D"
    // dimensional code -- the same column `parse_dimension_code` reads,
    // empirically confirmed against RDKit 2026.03.3's own `MolToMolBlock`.
    out.push_str("  chematic          3D\n");
    out.push_str(&metadata.comment);
    out.push('\n');

    let natoms = mol.atom_count();
    let nbonds = mol.bond_count();
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        natoms, nbonds
    ));

    for (idx, atom) in mol.atoms() {
        let sym = atom.element.symbol();
        let charge_code = encode_charge(atom.charge);
        let p = conformer
            .points
            .get(idx.0 as usize)
            .copied()
            .unwrap_or(Point3::zero());
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0{:>3}  0  0  0  0  0  0  0  0  0\n",
            p.x, p.y, p.z, sym, charge_code,
        ));
    }

    for (_idx, bond) in mol.bonds() {
        let a1 = bond.atom1.0 + 1;
        let a2 = bond.atom2.0 + 1;
        let btype = match bond.order {
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            BondOrder::QuerySingleOrDouble => 5,
            BondOrder::QuerySingleOrAromatic => 6,
            BondOrder::QueryDoubleOrAromatic => 7,
            BondOrder::QueryAny | BondOrder::Zero => 8,
            BondOrder::Quadruple => 4,
        };
        // Stereo field is always 0 -- see doc comment above.
        out.push_str(&format!("{:>3}{:>3}{:>3}{:>3}\n", a1, a2, btype, 0));
    }

    out.push_str("M  END\n");
    out
}

// ---------------------------------------------------------------------------
// Square-planar write-side validation -- see the read-side section above
// and `docs/rfcs/square_planar_mol_io_rfc.md` §8-9 for the design (no
// coordinate fabrication, ever; existing coordinates are validated against
// the declared tag, never trusted silently).
// ---------------------------------------------------------------------------

/// Which MOL/SDF variant a [`MolStereoWriteError`] was raised against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MolFormat {
    V2000,
    V3000,
}

impl core::fmt::Display for MolFormat {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::V2000 => write!(f, "V2000"),
            Self::V3000 => write!(f, "V3000"),
        }
    }
}

/// Why [`validate_square_planar_for_write`] refused to write a molecule
/// carrying `Chirality::SquarePlanar`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnsupportedStereoReason {
    /// No 3D conformer was supplied at all. MDL/CTfile has no symbolic field
    /// for a non-tetrahedral stereo tag (RFC §2) -- geometry is the only
    /// mechanism, and there is none to write here. The writer never
    /// fabricates one (RFC §8).
    NoConformerSupplied,
    /// A conformer was supplied, but this atom has no recorded
    /// `stereo_neighbor_order` to validate/write positions against.
    MissingNeighborOrder,
    /// This atom's `stereo_neighbor_order` contains the implicit-hydrogen
    /// sentinel -- there is no real spatial position for an implicit H to
    /// validate or write against (RFC §11).
    ImplicitHydrogenNeighbor,
    /// A conformer was supplied, but it has no real position for the center
    /// or one of its declared neighbors.
    MissingConformerPosition,
    /// The whole-molecule conformer is flat (`GeometryRank::FlatZero`) or
    /// indeterminate -- indistinguishable from an ordinary 2D depiction with
    /// no 3D data at all (RFC §10); this crate's own reader would never
    /// reperceive a tag from it, so writing one out would silently produce
    /// an unrecoverable file.
    WholeMoleculeConformerFlat,
    /// The conformer's local geometry around this atom does not
    /// unambiguously resolve to *any* square-planar permutation.
    GeometryRejected(SquarePlanarRejectionReason),
    /// The conformer's local geometry resolves to a *different* permutation
    /// than the one declared on the atom.
    GeometryTagMismatch { computed: SquarePlanarPermutation },
}

impl core::fmt::Display for UnsupportedStereoReason {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoConformerSupplied => write!(f, "no 3D conformer supplied"),
            Self::MissingNeighborOrder => write!(f, "atom has no recorded stereo_neighbor_order"),
            Self::ImplicitHydrogenNeighbor => write!(
                f,
                "stereo_neighbor_order contains an implicit-hydrogen slot"
            ),
            Self::MissingConformerPosition => {
                write!(
                    f,
                    "conformer is missing a position for this center or a neighbor"
                )
            }
            Self::WholeMoleculeConformerFlat => {
                write!(f, "whole-molecule conformer is flat/indeterminate")
            }
            Self::GeometryRejected(reason) => write!(f, "geometry rejected: {reason}"),
            Self::GeometryTagMismatch { computed } => write!(
                f,
                "conformer geometry encodes {computed:?}, which does not match the declared tag"
            ),
        }
    }
}

/// A MOL/SDF write was refused because `atom`'s declared
/// `chematic_core::StereoGeometry` could not be represented in `format` --
/// see [`UnsupportedStereoReason`] for why. Never raised for tetrahedral
/// chirality (which round-trips via wedge/hash bonds, already written
/// faithfully by every writer in this file).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolStereoWriteError {
    pub atom: AtomIdx,
    pub geometry: StereoGeometry,
    pub format: MolFormat,
    pub reason: UnsupportedStereoReason,
}

impl core::fmt::Display for MolStereoWriteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "cannot write atom {} ({:?} stereo) to MOL {}: {}",
            self.atom.0, self.geometry, self.format, self.reason
        )
    }
}

impl std::error::Error for MolStereoWriteError {}

/// Validate every `Chirality::SquarePlanar` atom in `mol` against `conformer`
/// before writing `format`. `Ok(())` when `mol` has no square-planar atom at
/// all (the common case -- this is a no-op for every existing caller), or
/// when every square-planar atom's declared tag matches what its
/// `conformer` positions actually encode. `pub`, not `pub(crate)`: callers
/// building their own write pipeline can use this as a pre-flight check
/// even without going through this crate's `*_checked` writer wrappers.
///
/// `conformer: None` means "the 2D-only writer path was used" -- always an
/// error (`NoConformerSupplied`) when `mol` has any square-planar atom, since
/// no 2D-only writer in this crate has a z channel to encode one (RFC §2).
pub fn validate_square_planar_for_write(
    mol: &Molecule,
    conformer: Option<&Coords3D>,
    format: MolFormat,
) -> Result<(), MolStereoWriteError> {
    // Whole-molecule flatness check once, up front: if ANY square-planar
    // atom exists and the conformer is flat/indeterminate, every one of them
    // fails the same way (RFC §10) -- report the first by atom order for a
    // deterministic, non-HashMap-dependent error.
    let whole_mol_flat_or_missing = match conformer {
        None => true,
        Some(conf) => matches!(
            classify_geometry_rank(&conf.points),
            GeometryRank::FlatZero | GeometryRank::Indeterminate
        ),
    };

    for (idx, atom) in mol.atoms() {
        let Chirality::SquarePlanar(declared) = atom.chirality else {
            continue;
        };
        let geometry = StereoGeometry::SquarePlanar;

        let Some(conf) = conformer else {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::NoConformerSupplied,
            });
        };
        if whole_mol_flat_or_missing {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::WholeMoleculeConformerFlat,
            });
        }

        let Some(order) = mol.stereo_neighbor_order(idx) else {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::MissingNeighborOrder,
            });
        };
        if order.len() != 4 {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::MissingNeighborOrder,
            });
        }
        if order.contains(&STEREO_H_SENTINEL) {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::ImplicitHydrogenNeighbor,
            });
        }

        let Some(center) = conf.points.get(idx.0 as usize).copied() else {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::MissingConformerPosition,
            });
        };
        let mut pts = [Point3::zero(); 4];
        for (slot, &n) in order.iter().enumerate() {
            match conf.points.get(n as usize).copied() {
                Some(p) => pts[slot] = p,
                None => {
                    return Err(MolStereoWriteError {
                        atom: idx,
                        geometry,
                        format,
                        reason: UnsupportedStereoReason::MissingConformerPosition,
                    });
                }
            }
        }

        match classify_square_planar_geometry(center, pts) {
            Ok(computed) if computed == declared => {}
            Ok(computed) => {
                return Err(MolStereoWriteError {
                    atom: idx,
                    geometry,
                    format,
                    reason: UnsupportedStereoReason::GeometryTagMismatch { computed },
                });
            }
            Err(reason) => {
                return Err(MolStereoWriteError {
                    atom: idx,
                    geometry,
                    format,
                    reason: UnsupportedStereoReason::GeometryRejected(reason),
                });
            }
        }
    }

    Ok(())
}

/// [`write_mol_with_conformer`], but fails closed with a typed
/// [`MolStereoWriteError`] instead of silently writing coordinates that
/// don't actually match a molecule's declared square-planar stereo (or
/// don't exist at all). See [`validate_square_planar_for_write`].
///
/// For a molecule with no `Chirality::SquarePlanar` atom, this is exactly
/// [`write_mol_with_conformer`] wrapped in `Ok`.
pub fn write_mol_with_conformer_checked(
    mol: &Molecule,
    metadata: &MolMetadata,
    conformer: &Coords3D,
) -> Result<String, MolStereoWriteError> {
    validate_square_planar_for_write(mol, Some(conformer), MolFormat::V2000)?;
    Ok(write_mol_with_conformer(mol, metadata, conformer))
}

// ---------------------------------------------------------------------------
// SDF writer
// ---------------------------------------------------------------------------

/// Serialise one or more molecules to SDF format.
///
/// `records` — slice of `(molecule, metadata, coords)` tuples.
/// `coords` is optional; pass an empty slice to write zero coordinates.
/// Each molecule block is terminated with `$$$$`.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo** -- see
/// [`write_mol`]'s doc comment; this delegates to [`write_mol_with_coords`]
/// per record. There is no `_checked` SDF-multi-record entry point today;
/// call [`write_sdf_record_with_conformer_checked`] once per record and
/// concatenate if `records` may carry it.
#[allow(clippy::type_complexity)]
pub fn write_sdf(records: &[(&Molecule, &MolMetadata, &[(f64, f64)])]) -> String {
    let mut out = String::new();
    for (mol, meta, coords) in records {
        out.push_str(&write_mol_with_coords(mol, meta, coords));
        out.push_str("$$$$\n");
    }
    out
}

/// Serialise one or more molecules to SDF format, appending per-atom partial
/// charges as an SD property `<PARTIAL_CHARGES>`.
///
/// `records` — slice of `(molecule, metadata, 2D-coords, charges)` tuples.
/// `charges[i]` is the partial charge for atom `i` (heavy atoms only).
/// Pass an empty charges slice to omit the property block.
///
/// Example SD block appended after `M  END`:
/// ```text
/// > <PARTIAL_CHARGES>
/// -0.2359 0.1076 -0.4500 0.1806
///
/// $$$$
/// ```
///
/// **Does not preserve `Chirality::SquarePlanar` stereo** -- see
/// [`write_mol`]'s doc comment; this delegates to [`write_mol_with_coords`]
/// per record.
#[allow(clippy::type_complexity)]
pub fn write_sdf_with_charges(
    records: &[(&Molecule, &MolMetadata, &[(f64, f64)], &[f64])],
) -> String {
    let mut out = String::new();
    for (mol, meta, coords, charges) in records {
        out.push_str(&write_mol_with_coords(mol, meta, coords));
        if !charges.is_empty() {
            out.push_str("> <PARTIAL_CHARGES>\n");
            let vals: Vec<String> = charges.iter().map(|q| format!("{q:.4}")).collect();
            out.push_str(&vals.join(" "));
            out.push_str("\n\n");
        }
        out.push_str("$$$$\n");
    }
    out
}

/// Serialise a single molecule to one SDF record with arbitrary SD data fields.
///
/// Keys starting with `_` are treated as internal/computed properties and are
/// omitted from the SD block (e.g. `_Name` is written into the MOL header, not
/// as an SD field).  The record is terminated with `$$$$`.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo** -- see
/// [`write_mol`]'s doc comment; this delegates to [`write_mol_with_coords`].
/// Use [`write_sdf_record_with_conformer_checked`] if `mol` may carry it.
pub fn write_sdf_record(
    mol: &Molecule,
    meta: &MolMetadata,
    coords: &[(f64, f64)],
    props: &std::collections::HashMap<String, String>,
) -> String {
    let mut out = write_mol_with_coords(mol, meta, coords);
    for (k, v) in props {
        if !k.starts_with('_') {
            out.push_str(&format!("> <{k}>\n{v}\n\n"));
        }
    }
    out.push_str("$$$$\n");
    out
}

/// Like [`write_sdf_record`] but emits a MOL V3000 (Extended Ctab) block
/// instead of V2000 — required for molecules with more than 999 atoms or
/// bonds, which don't fit V2000's fixed-width count fields.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo** -- see
/// [`write_mol`]'s doc comment; this delegates to [`crate::mol3000::write_mol_v3000`],
/// itself 2D-only.
pub fn write_sdf_record_v3000(
    mol: &Molecule,
    meta: &MolMetadata,
    coords: &[(f64, f64)],
    props: &std::collections::HashMap<String, String>,
) -> String {
    let mut out = crate::mol3000::write_mol_v3000(mol, meta, coords);
    for (k, v) in props {
        if !k.starts_with('_') {
            out.push_str(&format!("> <{k}>\n{v}\n\n"));
        }
    }
    out.push_str("$$$$\n");
    out
}

/// Like [`write_sdf_record`] but writes `conformer`'s real 3D coordinates
/// (V2000, via [`write_mol_with_conformer`]) instead of a 2D `(x, y)` slice
/// -- the 3D counterpart for a single SDF record. To write several
/// conformers of the same molecule as separate, repeated records (readable
/// back as one [`crate::sdf::ConformerEnsemble`] by
/// [`crate::sdf::read_sdf_conformer_ensembles`]), call this once per
/// conformer and concatenate the results.
///
/// **Does not validate `Chirality::SquarePlanar` stereo against
/// `conformer`** -- see [`write_mol_with_conformer`]'s doc comment for why.
/// Use [`write_sdf_record_with_conformer_checked`] instead to fail closed on
/// a geometry/tag mismatch (or a missing/flat conformer) rather than
/// trusting it.
pub fn write_sdf_record_with_conformer(
    mol: &Molecule,
    meta: &MolMetadata,
    conformer: &Coords3D,
    props: &std::collections::HashMap<String, String>,
) -> String {
    let mut out = write_mol_with_conformer(mol, meta, conformer);
    for (k, v) in props {
        if !k.starts_with('_') {
            out.push_str(&format!("> <{k}>\n{v}\n\n"));
        }
    }
    out.push_str("$$$$\n");
    out
}

/// [`write_sdf_record_with_conformer`], but fails closed with a typed
/// [`MolStereoWriteError`] instead of silently writing coordinates that
/// don't actually match a molecule's declared square-planar stereo (or
/// don't exist at all) -- the SDF-record counterpart of
/// [`write_mol_with_conformer_checked`]. See
/// [`validate_square_planar_for_write`].
///
/// For a molecule with no `Chirality::SquarePlanar` atom, this is exactly
/// [`write_sdf_record_with_conformer`] wrapped in `Ok`. To write a
/// multi-record SDF (e.g. several distinct molecules, or several conformers
/// of the same one), call this once per record and concatenate the
/// results -- there is no separate multi-record entry point, matching
/// [`write_sdf`]'s own per-record-call convention.
pub fn write_sdf_record_with_conformer_checked(
    mol: &Molecule,
    meta: &MolMetadata,
    conformer: &Coords3D,
    props: &std::collections::HashMap<String, String>,
) -> Result<String, MolStereoWriteError> {
    validate_square_planar_for_write(mol, Some(conformer), MolFormat::V2000)?;
    Ok(write_sdf_record_with_conformer(mol, meta, conformer, props))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal ethanol MOL V2000 block (CCO, 3 atoms, 2 bonds).
    const ETHANOL_MOL: &str = "\
ethanol
  chematic

  3  2  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
";

    #[test]
    fn test_parse_ethanol_counts() {
        let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse should succeed");
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        assert_eq!(meta.name, "ethanol");
    }

    #[test]
    fn test_parse_elements() {
        let (mol, _) = parse_mol(ETHANOL_MOL).expect("parse should succeed");
        let atoms: Vec<_> = mol.atoms().collect();
        assert_eq!(atoms[0].1.element, Element::C);
        assert_eq!(atoms[1].1.element, Element::C);
        assert_eq!(atoms[2].1.element, Element::O);
    }

    #[test]
    fn test_parse_bond_types() {
        // Two carbons: single, double, triple, aromatic bonds.
        let mol_str = "\
test
  chematic

  8  4  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  3  4  2  0
  5  6  3  0
  7  8  4  0
M  END
";
        let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
        let bonds: Vec<_> = mol.bonds().collect();
        assert_eq!(bonds[0].1.order, BondOrder::Single);
        assert_eq!(bonds[1].1.order, BondOrder::Double);
        assert_eq!(bonds[2].1.order, BondOrder::Triple);
        assert_eq!(bonds[3].1.order, BondOrder::Aromatic);
    }

    /// MDL bond type 9 (dative/coordinate) must read as `BondOrder::Dative`,
    /// not silently fall through to `Single` -- regression test for the
    /// platinum coordination-chemistry benchmark
    /// (`validation/platinum/FEASIBILITY.md`). See `mol3000.rs`'s
    /// `test_dative_bond_type_9_round_trips` for the RDKit-generated V3000
    /// version of the same finding (RDKit itself only ever writes bond type
    /// 9 in V3000, never V2000 -- this V2000 case defensively covers any
    /// other tool that does).
    #[test]
    fn test_parse_bond_type_9_is_dative() {
        let mol_str = "\
test
  chematic

  3  2  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 Pt  0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  9  0
  2  3  1  0
M  END
";
        let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
        let bonds: Vec<_> = mol.bonds().collect();
        assert_eq!(bonds[0].1.order, BondOrder::Dative);
        assert_eq!(mol.atom(bonds[0].1.atom1).element, Element::N);
        assert_eq!(mol.atom(bonds[0].1.atom2).element, Element::PT);
        assert_eq!(bonds[1].1.order, BondOrder::Single);
    }

    #[test]
    fn test_parse_query_bond_types_preserved() {
        let mol_str = "\
query_bonds
  chematic

  8  4  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  5  0
  3  4  6  0
  5  6  7  0
  7  8  8  0
M  END
";
        let (mol, meta) = parse_mol(mol_str).expect("parse query bonds");
        let bonds: Vec<_> = mol.bonds().collect();
        assert_eq!(bonds[0].1.order, BondOrder::QuerySingleOrDouble);
        assert_eq!(bonds[1].1.order, BondOrder::QuerySingleOrAromatic);
        assert_eq!(bonds[2].1.order, BondOrder::QueryDoubleOrAromatic);
        assert_eq!(bonds[3].1.order, BondOrder::QueryAny);

        let written = write_mol(&mol, &meta);
        assert!(written.contains("  1  2  5  0"), "{written}");
        assert!(written.contains("  7  8  8  0"), "{written}");
    }

    #[test]
    fn test_parse_charge() {
        // Nitrogen with charge code 3 (+1 formal charge).
        let mol_str = "\
charged
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
M  END
";
        let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
        assert_eq!(mol.atom(AtomIdx(0)).charge, 1);
    }

    #[test]
    fn test_parse_negative_charge() {
        // Oxygen with charge code 5 (-1 formal charge).
        let mol_str = "\
negcharge
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
M  END
";
        let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
        assert_eq!(mol.atom(AtomIdx(0)).charge, -1);
    }

    #[test]
    fn test_charged_atom_writes_implicit_h_in_smiles() {
        // MOL V2000 has no per-atom H-count field, so a charged atom read from
        // this format always has `hydrogen_count: None` (the format leaves
        // implicit H to be inferred, unlike bracket SMILES). Regression test
        // for a bug where the SMILES writer silently dropped these inferred
        // hydrogens on bracket atoms forced by charge/isotope/atom-map
        // (found via MRV oracle validation, see validation/mrv_io_parity_summary.json).
        let mol_str = "\
charged
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
M  END
";
        let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
        assert_eq!(mol.atom(AtomIdx(0)).hydrogen_count, None);
        assert_eq!(chematic_smiles::write(&mol), "[NH4+]");
        assert_eq!(chematic_smiles::canonical_smiles(&mol), "[NH4+]");
    }

    #[test]
    fn test_round_trip() {
        // Parse → write → parse again; atom and bond counts must match.
        let (mol1, meta1) = parse_mol(ETHANOL_MOL).expect("first parse");
        let written = write_mol(&mol1, &meta1);
        let (mol2, _meta2) = parse_mol(&written).expect("second parse");
        assert_eq!(mol1.atom_count(), mol2.atom_count());
        assert_eq!(mol1.bond_count(), mol2.bond_count());
    }

    #[test]
    fn test_round_trip_elements_preserved() {
        let (mol1, meta1) = parse_mol(ETHANOL_MOL).expect("first parse");
        let written = write_mol(&mol1, &meta1);
        let (mol2, _) = parse_mol(&written).expect("second parse");
        for ((_, a1), (_, a2)) in mol1.atoms().zip(mol2.atoms()) {
            assert_eq!(a1.element, a2.element);
        }
    }

    #[test]
    fn test_error_missing_v2000() {
        let bad = "\
bad
  prog

  3  2  0  0  0  0  0  0  0  0  0 V3000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
";
        assert!(matches!(
            parse_mol(bad),
            Err(MolParseError::InvalidCountLine { .. })
        ));
    }

    #[test]
    fn test_error_truncated_input() {
        // Counts line says 3 atoms but only 1 is provided.
        let bad = "\
trunc
  prog

  3  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
";
        assert!(matches!(parse_mol(bad), Err(MolParseError::UnexpectedEnd)));
    }

    #[test]
    fn test_error_invalid_counts_line() {
        // Counts line too short (no V2000 tag at all).
        let bad = "\
mol
  prog

  X  Y
M  END
";
        assert!(matches!(
            parse_mol(bad),
            Err(MolParseError::InvalidCountLine { .. })
        ));
    }

    #[test]
    fn test_write_contains_m_end() {
        let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse");
        let written = write_mol(&mol, &meta);
        assert!(written.contains("M  END"));
    }

    #[test]
    fn test_write_contains_v2000() {
        let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse");
        let written = write_mol(&mol, &meta);
        assert!(written.contains("V2000"));
    }

    #[test]
    fn test_parse_stereo_up_bond() {
        // MOL V2000 with a stereo=1 (Up) bond
        let mol_str = "\n\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  1  0  0  0\nM  END\n";
        let (mol, _) = crate::parse_mol(mol_str).unwrap();
        let bond = mol.bond(chematic_core::BondIdx(0));
        assert_eq!(bond.order, chematic_core::BondOrder::Up);
    }

    #[test]
    fn test_parse_stereo_down_bond() {
        let mol_str = "\n\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  6  0  0  0\nM  END\n";
        let (mol, _) = crate::parse_mol(mol_str).unwrap();
        let bond = mol.bond(chematic_core::BondIdx(0));
        assert_eq!(bond.order, chematic_core::BondOrder::Down);
    }

    #[test]
    fn test_parse_stereo_either_bond_code4_not_treated_as_wedge() {
        // MDL stereo code 4 ("either"/unspecified direction) is a third,
        // distinct state from a definite wedge (1) or hash (6) -- RDKit maps
        // it to `Bond::BondDir::UNKNOWN`, never a confident tag. Collapsing
        // it into `Up` (the pre-fix behavior) would fabricate a stereocenter
        // from a bond whose direction the file explicitly declares unknown,
        // now that stereo perception actually reads `BondOrder::Up`. This is
        // a documented, accepted lossy case: the "either" information itself
        // is not preserved (round-trips as a plain, unmarked single bond),
        // but it must never masquerade as a real wedge.
        let mol_str = "\n\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  4  0  0  0\nM  END\n";
        let (mol, _) = crate::parse_mol(mol_str).unwrap();
        let bond = mol.bond(chematic_core::BondIdx(0));
        assert_eq!(bond.order, chematic_core::BondOrder::Single);

        // Round-trip: writes as stereo field 0, not a false wedge/hash.
        let (mol, meta) = crate::parse_mol(mol_str).unwrap();
        let written = write_mol(&mol, &meta);
        assert!(written.contains("  1  2  1  0"), "{written}");
    }

    #[test]
    fn test_molmetadata_builder() {
        let meta = MolMetadata::default()
            .with_name("aspirin")
            .with_comment("test molecule");
        assert_eq!(meta.name, "aspirin");
        assert_eq!(meta.comment, "test molecule");
    }

    #[test]
    fn test_molmetadata_with_name_roundtrip() {
        // Build a two-atom molecule and write it → name appears on MOL header line 1.
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let mol = b.build();

        let meta = MolMetadata::default().with_name("acetic acid");
        let molblock = crate::write_mol(&mol, &meta);
        assert!(
            molblock.starts_with("acetic acid"),
            "MOL block must start with the molecule name"
        );
    }

    #[test]
    fn test_declared_max_atom_count_truncated_input_errors() {
        let bad = "\
max_atoms
  chematic

999  0  0  0  0  0  0  0  0  0  0 V2000
";
        assert!(matches!(parse_mol(bad), Err(MolParseError::UnexpectedEnd)));
    }

    #[test]
    fn test_declared_large_bond_count_truncated_input_errors() {
        let bad = "\
many_bonds
  chematic

  1 999  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
";
        assert!(matches!(parse_mol(bad), Err(MolParseError::UnexpectedEnd)));
    }

    #[test]
    fn test_write_sdf_record_props_roundtrip() {
        use crate::sdf::SdfRecordReader;
        use std::collections::HashMap;

        let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse");
        let coords = vec![(0.0, 0.0), (1.5, 0.0), (3.0, 0.0)];
        let mut props = HashMap::new();
        props.insert("Activity".to_string(), "7.2".to_string());
        props.insert("Source".to_string(), "test".to_string());
        props.insert("_Name".to_string(), "ethanol".to_string()); // internal, should be omitted

        let sdf = write_sdf_record(&mol, &meta, &coords, &props);

        // _Name must NOT appear as an SD field
        assert!(
            !sdf.contains("> <_Name>"),
            "internal prop leaked to SD block"
        );
        // Other props must appear
        assert!(sdf.contains("> <Activity>\n7.2"), "Activity missing");
        assert!(sdf.contains("> <Source>\ntest"), "Source missing");
        assert!(sdf.ends_with("$$$$\n"), "missing delimiter");

        // Parse back and verify
        let rec = SdfRecordReader::new(&sdf)
            .next()
            .expect("should have record")
            .expect("should parse");
        assert_eq!(rec.mol.atom_count(), mol.atom_count());
        assert_eq!(
            rec.properties.get("Activity").map(|s| s.as_str()),
            Some("7.2")
        );
        assert_eq!(
            rec.properties.get("Source").map(|s| s.as_str()),
            Some("test")
        );
        assert!(!rec.properties.contains_key("_Name"));
    }

    // ── Issue #171: blank MOL name line must not be eaten as inter-record
    // padding ──────────────────────────────────────────────────────────────
    // Fixture generated by a live RDKit 2026.03.3 oracle:
    //   AllChem.Compute2DCoords(Chem.MolFromSmiles("CC")); Chem.MolToMolBlock(mol)
    // RDKit's own MolToMolBlock leaves the name line blank for an unnamed
    // molecule -- this is the literal, spec-legal shape that previously broke
    // all 3 SDF-splitting entry points (see issue #171's repro).
    const RDKIT_BLANK_NAME_MOL: &str = "\
\n     RDKit          2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n   -0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n";

    fn crlf(s: &str) -> String {
        s.replace('\n', "\r\n")
    }

    #[test]
    fn test_read_sdf_with_diagnostics_blank_name_single_record() {
        let sdf = format!("{RDKIT_BLANK_NAME_MOL}$$$$\n");
        let reports = read_sdf_with_diagnostics(&sdf).expect("parse should succeed");
        assert_eq!(reports.len(), 1);
        assert_eq!(reports[0].metadata.name, "");
        assert_eq!(reports[0].mol.atom_count(), 2);
    }

    #[test]
    fn test_read_sdf_with_diagnostics_blank_name_first_record() {
        let sdf = format!("{RDKIT_BLANK_NAME_MOL}$$$$\n{ETHANOL_MOL}$$$$\n");
        let reports = read_sdf_with_diagnostics(&sdf).expect("parse should succeed");
        assert_eq!(reports.len(), 2);
        assert_eq!(reports[0].metadata.name, "");
        assert_eq!(reports[0].mol.atom_count(), 2);
        assert_eq!(reports[1].metadata.name, "ethanol");
        assert_eq!(reports[1].mol.atom_count(), 3);
    }

    #[test]
    fn test_read_sdf_with_diagnostics_blank_name_middle_record() {
        let sdf = format!("{ETHANOL_MOL}$$$$\n{RDKIT_BLANK_NAME_MOL}$$$$\n{ETHANOL_MOL}$$$$\n");
        let reports = read_sdf_with_diagnostics(&sdf).expect("parse should succeed");
        assert_eq!(reports.len(), 3);
        assert_eq!(reports[0].metadata.name, "ethanol");
        assert_eq!(reports[1].metadata.name, "");
        assert_eq!(reports[1].mol.atom_count(), 2);
        assert_eq!(reports[2].metadata.name, "ethanol");
    }

    #[test]
    fn test_read_sdf_with_diagnostics_blank_name_last_record() {
        let sdf = format!("{ETHANOL_MOL}$$$$\n{RDKIT_BLANK_NAME_MOL}$$$$\n");
        let reports = read_sdf_with_diagnostics(&sdf).expect("parse should succeed");
        assert_eq!(reports.len(), 2);
        assert_eq!(reports[0].metadata.name, "ethanol");
        assert_eq!(reports[1].metadata.name, "");
        assert_eq!(reports[1].mol.atom_count(), 2);
    }

    #[test]
    fn test_read_sdf_with_diagnostics_blank_name_crlf() {
        let sdf = crlf(&format!(
            "{ETHANOL_MOL}$$$$\n{RDKIT_BLANK_NAME_MOL}$$$$\n{ETHANOL_MOL}$$$$\n"
        ));
        let reports = read_sdf_with_diagnostics(&sdf).expect("parse should succeed (CRLF)");
        assert_eq!(reports.len(), 3);
        assert_eq!(reports[1].metadata.name, "");
        assert_eq!(reports[1].mol.atom_count(), 2);
    }

    #[test]
    fn test_parse_sdf_with_coords_blank_name_middle_record() {
        // parse_sdf_with_coords is a thin wrapper over read_sdf_with_diagnostics;
        // confirm the fix reaches it too.
        let sdf = format!("{ETHANOL_MOL}$$$$\n{RDKIT_BLANK_NAME_MOL}$$$$\n{ETHANOL_MOL}$$$$\n");
        let mols = parse_sdf_with_coords(&sdf).expect("parse should succeed");
        assert_eq!(mols.len(), 3);
        assert_eq!(mols[1].1.name, "");
        assert_eq!(mols[1].0.atom_count(), 2);
    }

    #[test]
    fn test_read_sdf_with_diagnostics_malformed_input_recovery_unaffected() {
        // Existing malformed-input behavior (bad counts line) must still error,
        // confirming the blank-line-skip removal didn't loosen error recovery.
        let bad_sdf = format!("{ETHANOL_MOL}$$$$\nbad\n  prog\n\n  X  Y\nM  END\n$$$$\n");
        assert!(read_sdf_with_diagnostics(&bad_sdf).is_err());
    }
}

/// Square-planar 3D-coordinate-derived stereo: classification, perception,
/// and write-side validation. See `docs/rfcs/square_planar_mol_io_rfc.md`
/// for the design and `tests/square_planar_mol_io.rs` for the
/// reader/writer/renumbering integration tests.
#[cfg(test)]
mod square_planar_tests {
    use super::*;

    /// Build the 4 neighbor positions (in trans_pairs() slot order) for an
    /// ideal, undistorted square-planar arrangement matching `tag` -- reuses
    /// `trans_pairs()` (not a hand-picked per-tag layout) so this generator
    /// and `classify_square_planar_geometry` agree by construction on what
    /// "matches the tag" means. All points share z = 1.5 (nonzero, to avoid
    /// the FlatZero ambiguity -- see RFC §10).
    fn ideal_square_planar_neighbors(tag: SquarePlanarPermutation) -> [Point3; 4] {
        let [(a, b), (c, d)] = tag.trans_pairs();
        let mut pts = [Point3::zero(); 4];
        let at = |deg: f64| -> Point3 {
            let r = deg.to_radians();
            Point3::new(1.5 * r.cos(), 1.5 * r.sin(), 1.5)
        };
        pts[a as usize] = at(45.0);
        pts[b as usize] = at(225.0); // trans to a
        pts[c as usize] = at(135.0);
        pts[d as usize] = at(315.0); // trans to c
        pts
    }

    #[test]
    fn classify_recovers_all_three_tags_from_ideal_geometry() {
        // Center must share the neighbors' z=1.5 plane (see
        // `ideal_square_planar_neighbors`) -- coplanarity is center-relative.
        let center = Point3::new(0.0, 0.0, 1.5);
        for tag in [
            SquarePlanarPermutation::SP1,
            SquarePlanarPermutation::SP2,
            SquarePlanarPermutation::SP3,
        ] {
            let neighbors = ideal_square_planar_neighbors(tag);
            assert_eq!(
                classify_square_planar_geometry(center, neighbors),
                Ok(tag),
                "failed to recover {tag:?}"
            );
        }
    }

    #[test]
    fn classify_rejects_genuinely_tetrahedral_geometry_as_not_coplanar() {
        // Ideal tetrahedral directions (unit vectors), scaled -- never
        // coplanar with the center.
        let center = Point3::zero();
        let neighbors = [
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(1.0, -1.0, -1.0),
            Point3::new(-1.0, 1.0, -1.0),
            Point3::new(-1.0, -1.0, 1.0),
        ];
        assert_eq!(
            classify_square_planar_geometry(center, neighbors),
            Err(SquarePlanarRejectionReason::NotCoplanar)
        );
    }

    #[test]
    fn classify_rejects_degenerate_coincident_neighbor() {
        let center = Point3::zero();
        let mut neighbors = ideal_square_planar_neighbors(SquarePlanarPermutation::SP1);
        neighbors[0] = center; // coincides with the center itself
        assert_eq!(
            classify_square_planar_geometry(center, neighbors),
            Err(SquarePlanarRejectionReason::DegenerateBondVector)
        );
    }

    #[test]
    fn classify_rejects_ambiguous_geometry_with_no_trans_pair() {
        // All 4 neighbors clustered within one 90-degree quadrant -- coplanar,
        // but no two are ever close to trans (180 degrees) apart.
        let center = Point3::new(0.0, 0.0, 1.5);
        let at = |deg: f64| -> Point3 {
            let r = deg.to_radians();
            Point3::new(1.5 * r.cos(), 1.5 * r.sin(), 1.5)
        };
        let neighbors = [at(0.0), at(20.0), at(40.0), at(60.0)];
        assert_eq!(
            classify_square_planar_geometry(center, neighbors),
            Err(SquarePlanarRejectionReason::AmbiguousTransPairing)
        );
    }

    /// Build a 5-atom (center + 4 ligand) molecule plus a matching `Coords3D`
    /// conformer encoding `tag`'s ideal geometry. `center_element` lets
    /// callers probe the element-eligibility gate (RFC §6).
    fn square_planar_fixture(
        center_element: Element,
        tag: SquarePlanarPermutation,
    ) -> (Molecule, Coords3D) {
        let mut b = MoleculeBuilder::new();
        let center = b.add_atom(Atom::new(center_element));
        let cl1 = b.add_atom(Atom::new(Element::CL));
        let cl2 = b.add_atom(Atom::new(Element::CL));
        let n1 = b.add_atom(Atom::new(Element::N));
        let n2 = b.add_atom(Atom::new(Element::N));
        b.add_bond(center, cl1, BondOrder::Single).unwrap();
        b.add_bond(center, cl2, BondOrder::Single).unwrap();
        b.add_bond(center, n1, BondOrder::Single).unwrap();
        b.add_bond(center, n2, BondOrder::Single).unwrap();
        let mol = b.build();

        let neighbor_pts = ideal_square_planar_neighbors(tag);
        let points = vec![
            Point3::new(0.0, 0.0, 1.5), // center
            neighbor_pts[0],
            neighbor_pts[1],
            neighbor_pts[2],
            neighbor_pts[3],
        ];
        (mol, Coords3D { points })
    }

    #[test]
    fn perceive_assigns_all_three_tags_for_an_undefined_valence_element() {
        for tag in [
            SquarePlanarPermutation::SP1,
            SquarePlanarPermutation::SP2,
            SquarePlanarPermutation::SP3,
        ] {
            let (mut mol, conformer) = square_planar_fixture(Element::PT, tag);
            let diagnostics = perceive_square_planar_from_3d(&mut mol, &conformer);
            assert!(
                diagnostics.is_empty(),
                "{tag:?}: unexpected {diagnostics:?}"
            );
            assert_eq!(mol.atom(AtomIdx(0)).chirality, Chirality::SquarePlanar(tag));
            assert_eq!(
                mol.stereo_neighbor_order(AtomIdx(0)),
                Some([1u32, 2, 3, 4].as_slice())
            );
        }
    }

    #[test]
    fn perceive_never_tags_an_element_with_a_defined_valence() {
        // Identical coplanar-square coordinates, but Carbon has a defined
        // (organic-subset) valence list -- RDKit itself does not perceive
        // non-tetrahedral chirality for such elements either (RFC §4/§6),
        // even from geometrically-square coordinates.
        let (mut mol, conformer) = square_planar_fixture(Element::C, SquarePlanarPermutation::SP1);
        let diagnostics = perceive_square_planar_from_3d(&mut mol, &conformer);
        assert!(diagnostics.is_empty());
        assert_eq!(mol.atom(AtomIdx(0)).chirality, Chirality::None);
    }

    #[test]
    fn perceive_reports_a_diagnostic_for_coplanar_but_ambiguous_geometry() {
        let mut b = MoleculeBuilder::new();
        let center = b.add_atom(Atom::new(Element::PT));
        let cl1 = b.add_atom(Atom::new(Element::CL));
        let cl2 = b.add_atom(Atom::new(Element::CL));
        let n1 = b.add_atom(Atom::new(Element::N));
        let n2 = b.add_atom(Atom::new(Element::N));
        b.add_bond(center, cl1, BondOrder::Single).unwrap();
        b.add_bond(center, cl2, BondOrder::Single).unwrap();
        b.add_bond(center, n1, BondOrder::Single).unwrap();
        b.add_bond(center, n2, BondOrder::Single).unwrap();
        let mut mol = b.build();

        let at = |deg: f64| -> Point3 {
            let r = deg.to_radians();
            Point3::new(1.5 * r.cos(), 1.5 * r.sin(), 1.5)
        };
        let points = vec![
            Point3::new(0.0, 0.0, 1.5),
            at(0.0),
            at(20.0),
            at(40.0),
            at(60.0),
        ];
        let conformer = Coords3D { points };

        let diagnostics = perceive_square_planar_from_3d(&mut mol, &conformer);
        assert_eq!(
            diagnostics,
            vec![SquarePlanarPerceptionDiagnostic {
                atom: AtomIdx(0),
                reason: SquarePlanarRejectionReason::AmbiguousTransPairing,
            }]
        );
        assert_eq!(mol.atom(AtomIdx(0)).chirality, Chirality::None);
    }

    #[test]
    fn validate_is_a_no_op_for_a_molecule_with_no_square_planar_atom() {
        let (mol, _conformer) = square_planar_fixture(Element::PT, SquarePlanarPermutation::SP1);
        // No chirality was set on this fixture yet (perceive was never
        // called) -- validate must be Ok regardless of conformer presence.
        assert!(validate_square_planar_for_write(&mol, None, MolFormat::V2000).is_ok());
    }

    #[test]
    fn validate_rejects_no_conformer_supplied() {
        let (mut mol, _conformer) =
            square_planar_fixture(Element::PT, SquarePlanarPermutation::SP1);
        mol.set_chirality(
            AtomIdx(0),
            Chirality::SquarePlanar(SquarePlanarPermutation::SP1),
        );
        mol.set_stereo_neighbor_order(AtomIdx(0), vec![1, 2, 3, 4]);
        let err = validate_square_planar_for_write(&mol, None, MolFormat::V2000).unwrap_err();
        assert_eq!(err.reason, UnsupportedStereoReason::NoConformerSupplied);
        assert_eq!(err.format, MolFormat::V2000);
        assert_eq!(err.geometry, StereoGeometry::SquarePlanar);
    }

    #[test]
    fn validate_rejects_whole_molecule_flat_conformer() {
        let (mut mol, mut conformer) =
            square_planar_fixture(Element::PT, SquarePlanarPermutation::SP1);
        mol.set_chirality(
            AtomIdx(0),
            Chirality::SquarePlanar(SquarePlanarPermutation::SP1),
        );
        mol.set_stereo_neighbor_order(AtomIdx(0), vec![1, 2, 3, 4]);
        // Flatten every z to exactly 0 -- the ordinary "2D depiction, no
        // real 3D data" case (RFC §10), indistinguishable from wedge-only
        // MOL input by this crate's own reader.
        for p in &mut conformer.points {
            p.z = 0.0;
        }
        let err =
            validate_square_planar_for_write(&mol, Some(&conformer), MolFormat::V3000).unwrap_err();
        assert_eq!(
            err.reason,
            UnsupportedStereoReason::WholeMoleculeConformerFlat
        );
    }

    #[test]
    fn validate_rejects_geometry_tag_mismatch() {
        // Declare SP1 but supply SP2-shaped coordinates.
        let (mut mol, conformer) = square_planar_fixture(Element::PT, SquarePlanarPermutation::SP2);
        mol.set_chirality(
            AtomIdx(0),
            Chirality::SquarePlanar(SquarePlanarPermutation::SP1),
        );
        mol.set_stereo_neighbor_order(AtomIdx(0), vec![1, 2, 3, 4]);
        let err =
            validate_square_planar_for_write(&mol, Some(&conformer), MolFormat::V2000).unwrap_err();
        assert_eq!(
            err.reason,
            UnsupportedStereoReason::GeometryTagMismatch {
                computed: SquarePlanarPermutation::SP2
            }
        );
    }

    #[test]
    fn validate_rejects_implicit_hydrogen_neighbor_slot() {
        let (mut mol, conformer) = square_planar_fixture(Element::PT, SquarePlanarPermutation::SP1);
        mol.set_chirality(
            AtomIdx(0),
            Chirality::SquarePlanar(SquarePlanarPermutation::SP1),
        );
        mol.set_stereo_neighbor_order(AtomIdx(0), vec![1, 2, 3, STEREO_H_SENTINEL]);
        let err =
            validate_square_planar_for_write(&mol, Some(&conformer), MolFormat::V2000).unwrap_err();
        assert_eq!(
            err.reason,
            UnsupportedStereoReason::ImplicitHydrogenNeighbor
        );
    }

    #[test]
    fn validate_accepts_matching_geometry_and_checked_writer_produces_output() {
        let (mut mol, conformer) = square_planar_fixture(Element::PT, SquarePlanarPermutation::SP1);
        mol.set_chirality(
            AtomIdx(0),
            Chirality::SquarePlanar(SquarePlanarPermutation::SP1),
        );
        mol.set_stereo_neighbor_order(AtomIdx(0), vec![1, 2, 3, 4]);
        assert!(validate_square_planar_for_write(&mol, Some(&conformer), MolFormat::V2000).is_ok());

        let block = write_mol_with_conformer_checked(&mol, &MolMetadata::default(), &conformer)
            .expect("matching geometry must write successfully");
        assert!(block.contains("V2000"));
        assert!(block.contains("Pt"));
    }

    /// The `wedge_vs_3d_conflicts` regression this PR fixes: before the
    /// `is_tetrahedral()` gate, ANY non-`None` chirality (not just
    /// `Chirality::None`) fell through into a computation that only ever
    /// produces `Chirality::Clockwise`/`CounterClockwise` and compares that
    /// against `atom.chirality` -- which can never equal
    /// `Chirality::SquarePlanar`, so every square-planar-tagged atom with a
    /// real conformer would have been flagged as a spurious
    /// `WedgeVs3DParityConflict`, regardless of whether its geometry was
    /// actually self-consistent.
    #[test]
    fn wedge_vs_3d_conflicts_skips_square_planar_chirality() {
        let (mut mol, conformer) = square_planar_fixture(Element::PT, SquarePlanarPermutation::SP1);
        mol.set_chirality(
            AtomIdx(0),
            Chirality::SquarePlanar(SquarePlanarPermutation::SP1),
        );
        mol.set_stereo_neighbor_order(AtomIdx(0), vec![1, 2, 3, 4]);
        let diagnostics = wedge_vs_3d_conflicts(&mol, &conformer);
        assert!(
            diagnostics.is_empty(),
            "square-planar chirality must never reach the tetrahedral-only conflict check: {diagnostics:?}"
        );
    }
}
