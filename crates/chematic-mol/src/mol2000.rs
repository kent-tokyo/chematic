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
    Point3, RGroupLabel, STEREO_H_SENTINEL, SquarePlanarPermutation, StereoGeometry,
};
use chematic_perception::{
    EzDirectionDiagnostic, StereoDiagnostic, apply_ctab_local_parity_from_wedges_with_diagnostics,
    apply_ez_directions_from_2d_ex,
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
    /// Opaque V3000 SGROUP logical lines preserved across a core V3000
    /// parse/write round trip. The molecule model does not interpret these
    /// polymer/query semantics; callers must not treat this as a typed SGroup
    /// representation.
    pub v3000_sgroups: Vec<String>,
    /// Opaque V3000 bond-line attributes such as `ENDPTS=` and `ATTACH=`,
    /// keyed by the writer's 1-based bond order. Rich query/coordination
    /// semantics are not interpreted by the core bond model, but are retained
    /// when a non-contiguous source V3000 id is normalized to writer order.
    pub v3000_bond_properties: Vec<(u32, String)>,
    /// Opaque V3000 atom-line attributes that the core atom model does not
    /// interpret. This currently preserves atom `CFG` parity emitted by
    /// third-party writers, so a V3000 parse/write round trip does not make a
    /// stereochemical CTAB invalid for its original writer. Keys use the
    /// writer's 1-based atom order, not the source V3000 atom id.
    pub v3000_atom_properties: Vec<(u32, String)>,
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
/// Atoms of an aromatic (type 4) bond are aromatic, as RDKit reads them: the
/// file records aromaticity on bonds only, and leaving the atoms aliphatic
/// gave `C:1:C:C:C:C:C:1` for benzene (SMARTS `c` then missed it).
pub(crate) fn flag_aromatic_bond_atoms(mol: &mut Molecule) {
    let atoms: Vec<AtomIdx> = mol
        .bonds()
        .filter(|(_, b)| b.order == BondOrder::Aromatic)
        .flat_map(|(_, b)| [b.atom1, b.atom2])
        .collect();
    for idx in atoms {
        if !mol.atom(idx).aromatic {
            mol.set_atom_aromatic(idx, true);
        }
    }
}

/// Bond orders to write for aromatic bonds. Aromatic bonds are written as
/// type 4 when a reader can recover every aromatic atom's H count from the
/// bonds alone; otherwise (pyrrole's `[nH]`, a charged aromatic atom) the
/// Kekulé orders are written instead, as RDKit always does, since type 4
/// would turn `c1cc[nH]c1` into `c1ccnc1`. `None`: write type 4.
pub(crate) fn kekule_orders_for_writing(
    mol: &Molecule,
) -> Option<std::collections::HashMap<chematic_core::BondIdx, BondOrder>> {
    if !mol.bonds().any(|(_, b)| b.order == BondOrder::Aromatic) {
        return None;
    }
    // An atom's inferred H count depends on its own bonds only, so each
    // stored count is compared with the inference directly.
    let inferable = mol.atoms().all(|(idx, atom)| match atom.hydrogen_count {
        Some(h) if atom.aromatic && !atom.wildcard => {
            chematic_core::valence::valence_inferred_hcount(mol, idx) == h
        }
        _ => true,
    });
    if inferable {
        return None;
    }
    deterministic_kekule(mol)
}

/// Kekulé bond orders that do not depend on `mol`'s atom numbering: the
/// molecule is kekulized with its atoms in canonical order and the orders
/// are mapped back, so two spellings of one molecule get the same Kekulé
/// structure (up to symmetry).
fn deterministic_kekule(
    mol: &Molecule,
) -> Option<std::collections::HashMap<chematic_core::BondIdx, BondOrder>> {
    let order = chematic_smiles::canonical_atom_order(mol);
    let mut new_of = vec![0u32; mol.atom_count()];
    let mut builder = chematic_core::MoleculeBuilder::new();
    for &old in &order {
        new_of[old] = builder
            .add_atom(mol.atom(chematic_core::AtomIdx(old as u32)).clone())
            .0;
    }
    let mut bonds: Vec<(u32, u32, chematic_core::BondIdx)> = mol
        .bonds()
        .map(|(b, e)| {
            let (x, y) = (new_of[e.atom1.0 as usize], new_of[e.atom2.0 as usize]);
            (x.min(y), x.max(y), b)
        })
        .collect();
    bonds.sort_unstable();
    let mut old_of_new_bond = Vec::with_capacity(bonds.len());
    for &(x, y, b) in &bonds {
        builder
            .add_bond(
                chematic_core::AtomIdx(x),
                chematic_core::AtomIdx(y),
                mol.bond(b).order,
            )
            .ok()?;
        old_of_new_bond.push(b);
    }
    let canon = builder.build();
    let kek = chematic_core::kekulize(&canon).ok()?;
    Some(
        kek.into_iter()
            .map(|(nb, order)| (old_of_new_bond[nb.0 as usize], order))
            .collect(),
    )
}

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

/// Decode the V2000 atom-line mass-difference field into chematic's isotope
/// mass number. MDL stores the difference from the element's most abundant
/// isotope; `Element::atomic_mass()` provides the corresponding rounded mass
/// number used by the rest of the crate.
fn decode_mass_difference(element: Element, field: &str) -> Option<u16> {
    if field == " 0" {
        return None;
    }
    let difference = field.trim().parse::<i16>().ok()?;
    if difference == 0 {
        return None;
    }
    let mass = (element.atomic_mass().round() as i16).checked_add(difference)?;
    (mass > 0).then_some(mass as u16)
}

/// Encode an isotope mass number for the fixed-width V2000 mass-difference
/// field. The field is three characters wide; values outside that representable
/// range cannot be faithfully written by V2000 and are left absent. Callers
/// requiring arbitrary isotope masses should use V3000 (`MASS=`).
fn encode_mass_difference(element: Element, isotope: Option<u16>) -> Option<i16> {
    let isotope = isotope? as i32;
    let base = element.atomic_mass().round() as i32;
    let difference = isotope.checked_sub(base)?;
    (i16::MIN as i32..=i16::MAX as i32)
        .contains(&difference)
        .then_some(difference as i16)
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
    parse_unsigned_ascii(field)
        .or_else(|| field.trim().parse::<usize>().ok())
        .ok_or_else(|| make_err(line_num, format!("cannot parse integer from '{field}'")))
}

/// Parse the common CTfile integer spelling without invoking the generic
/// `FromStr` machinery. Fixed-width count/bond fields are overwhelmingly
/// ASCII spaces followed by decimal digits; uncommon accepted spellings
/// (such as a leading `+`) fall back to `str::parse` at the call site.
#[inline]
fn parse_unsigned_ascii(field: &str) -> Option<usize> {
    // Common fixed-width fields: avoid a checked multiply per digit.
    // Every other spelling still takes the original parser/fallback path.
    if let [a, b, c] = *field.as_bytes()
        && c.is_ascii_digit()
    {
        let last = (c - b'0') as usize;
        if a == b' ' {
            if b == b' ' {
                return Some(last);
            }
            if b.is_ascii_digit() {
                return Some((b - b'0') as usize * 10 + last);
            }
        } else if a.is_ascii_digit() && b.is_ascii_digit() {
            return Some((a - b'0') as usize * 100 + (b - b'0') as usize * 10 + last);
        }
    }
    let mut value = 0usize;
    let mut saw_digit = false;
    for &byte in field.as_bytes() {
        match byte {
            b' ' if !saw_digit => {}
            b'0'..=b'9' => {
                saw_digit = true;
                value = value.checked_mul(10)?.checked_add((byte - b'0') as usize)?;
            }
            _ => return None,
        }
    }
    saw_digit.then_some(value)
}

mod reader;
pub use reader::*;
mod writer;
pub use writer::*;
#[cfg(test)]
use writer::{push_right_aligned_i16, push_right_aligned_u32};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "mol2000/tests.rs"]
mod tests;

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
        let result = classify_square_planar_geometry(center, neighbors);
        assert_eq!(result, Err(SquarePlanarRejectionReason::NotCoplanar));
        assert!(result.unwrap_err().to_string().contains("not coplanar"));
    }

    #[test]
    fn classify_rejects_degenerate_coincident_neighbor() {
        let center = Point3::zero();
        let mut neighbors = ideal_square_planar_neighbors(SquarePlanarPermutation::SP1);
        neighbors[0] = center; // coincides with the center itself
        let result = classify_square_planar_geometry(center, neighbors);
        assert_eq!(
            result,
            Err(SquarePlanarRejectionReason::DegenerateBondVector)
        );
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("degenerate bond vector")
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
        let result = classify_square_planar_geometry(center, neighbors);
        assert_eq!(
            result,
            Err(SquarePlanarRejectionReason::AmbiguousTransPairing)
        );
        assert!(result.unwrap_err().to_string().contains("trans-pairing"));
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
        assert!(
            err.to_string()
                .contains("atom 0 (SquarePlanar stereo) to MOL V2000")
        );
        assert!(err.to_string().contains("no 3D conformer supplied"));
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
        assert!(err.to_string().contains("MOL V3000"));
        assert!(err.to_string().contains("flat/indeterminate"));
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
        assert!(err.to_string().contains("geometry encodes SP2"));
        assert!(err.to_string().contains("does not match the declared tag"));
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
        assert!(err.to_string().contains("implicit-hydrogen slot"));
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

#[cfg(test)]
mod stereo_loss_message_boundary_tests {
    use super::*;

    #[test]
    fn actual_degenerate_layout_loss_reports_the_unexpressed_stereo_identity() {
        for (source, atom, bond) in [
            ("C[C@H](F)Cl", Some(AtomIdx(1)), None),
            ("F/C=C/F", None, Some(chematic_core::BondIdx(1))),
        ] {
            let mol = chematic_smiles::parse(source).unwrap();
            let coords: Vec<_> = (0..mol.atom_count())
                .map(|i| (1.5 * i as f64, 0.0))
                .collect();
            let (_, loss) = write_mol_with_stereo_report(&mol, &MolMetadata::default(), &coords);
            assert!(!loss.is_empty(), "{source}");
            let message = loss.to_string();
            if let Some(atom) = atom {
                assert!(loss.centres.contains(&atom));
                assert!(message.contains("tetrahedral centres without a wedge (atoms 1)"));
            }
            if let Some(bond) = bond {
                assert!(loss.double_bonds.contains(&bond));
                assert!(message.contains("E/Z written as either (bonds 1)"));
            }
        }
        let mol = chematic_smiles::parse("CCO").unwrap();
        let (_, loss) = write_mol_with_stereo_report(&mol, &MolMetadata::default(), &[]);
        assert!(loss.is_empty());
        assert_eq!(loss.to_string(), "no stereo lost");
    }
}
