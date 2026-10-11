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

    let metadata = MolMetadata {
        name,
        comment,
        v3000_sgroups: Vec::new(),
        v3000_bond_properties: Vec::new(),
        v3000_atom_properties: Vec::new(),
    };

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

    let mut valence_fields: Vec<(AtomIdx, u8)> = Vec::new();
    for atom_i in 0..natoms {
        let (raw_lineno, atom_line) = next_line()?;

        // Keep validating Z even on the graph-only supplier path. A malformed
        // value was already a typed error there, so skipping all coordinate
        // parsing would silently loosen strictParsing. Only X/Y conversion
        // and coordinate storage are diagnostic-path work.
        let z: f64 = match atom_line.get(20..30) {
            None => 0.0,
            Some("    0.0000") => 0.0,
            Some(raw) => {
                let trimmed = raw.trim();
                if trimmed.is_empty() {
                    0.0
                } else {
                    // Zero Z is by far the dominant V2000 spelling (especially
                    // in 2D SDF corpora). Avoid the generic float parser for the
                    // exact writer-produced value while retaining the existing
                    // strict parse/finite validation for every other spelling.
                    let val: f64 = if trimmed == "0.0000" || trimmed == "0" {
                        0.0
                    } else {
                        trimmed.parse().map_err(|_| {
                            make_atom_err(
                                raw_lineno,
                                format!("cannot parse z coordinate from '{trimmed}'"),
                            )
                        })?
                    };
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
        let sym = atom_line.get(31..34).ok_or_else(|| {
            make_atom_err(
                raw_lineno,
                format!("atom line {atom_i} too short for element field"),
            )
        })?;

        let trimmed_sym = sym.trim();
        let (mut atom, r_group) = if trimmed_sym == "*" {
            (Atom::wildcard(), None)
        } else if trimmed_sym == "R" || trimmed_sym == "R#" {
            (Atom::wildcard(), Some(RGroupLabel::unnumbered()))
        } else if let Some(number) = trimmed_sym.strip_prefix('R')
            && !number.is_empty()
            && number.bytes().all(|byte| byte.is_ascii_digit())
        {
            let label = number
                .parse::<u16>()
                .ok()
                .and_then(RGroupLabel::numbered)
                .ok_or_else(|| MolParseError::UnknownElement {
                    symbol: trimmed_sym.to_string(),
                    line: raw_lineno,
                })?;
            (Atom::wildcard(), Some(label))
        } else {
            let element = match sym {
                "C  " => Some(Element::C),
                "N  " => Some(Element::N),
                "O  " => Some(Element::O),
                "H  " => Some(Element::H),
                "S  " => Some(Element::S),
                "P  " => Some(Element::P),
                "F  " => Some(Element::F),
                "Cl " => Some(Element::CL),
                "Br " => Some(Element::BR),
                _ => Element::from_symbol(trimmed_sym),
            }
            .ok_or_else(|| MolParseError::UnknownElement {
                symbol: trimmed_sym.to_string(),
                line: raw_lineno,
            })?;
            (Atom::new(element), None)
        };

        // Charge code: bytes 36–38 (3 chars).
        let charge = atom_line
            .get(36..39)
            .map(|ccc| {
                if ccc == "  0" {
                    return 0;
                }
                let code = parse_unsigned_ascii(ccc)
                    .and_then(|value| i8::try_from(value).ok())
                    .unwrap_or_else(|| ccc.trim().parse().unwrap_or(0));
                decode_charge(code)
            })
            .unwrap_or(0);

        // V2000 mass difference: bytes 34–35 (two-character signed field).
        let isotope = atom_line
            .get(34..36)
            .and_then(|field| decode_mass_difference(atom.element, field));
        // V2000 atom-atom mapping number (columns 61-63, zero-based
        // 60..63). Zero means that the atom is not mapped.
        let atom_map = atom_line
            .get(60..63)
            .and_then(parse_unsigned_ascii)
            .and_then(|value| u16::try_from(value).ok())
            .filter(|&value| value != 0);

        // Valence field (vvv, columns 48-50): the atom's total valence, 15
        // meaning zero. Writers (RDKit, chematic) set it when the bonds alone
        // would give the atom a different hydrogen count.
        let valence_field = atom_line
            .get(48..51)
            .and_then(parse_unsigned_ascii)
            .filter(|&v| v != 0 && v <= 15)
            .map(|v| if v == 15 { 0u8 } else { v as u8 });

        atom.charge = charge;
        atom.isotope = isotope;
        atom.atom_map = atom_map;
        let idx = builder.add_atom(atom);
        if let Some(valence) = valence_field {
            valence_fields.push((idx, valence));
        }
        if let Some(label) = r_group {
            builder.set_r_group(idx, label);
        }
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

    // Read property lines until `M  END` (or EOF if absent). Even the graph-only
    // supplier must consume `M  CHG` and `M  ISO`: RDKit uses them for values
    // that do not fit, or are intentionally omitted from, the atom block's
    // fixed-width charge and mass-difference fields.
    let mut property_charges = Vec::new();
    let mut property_isotopes = Vec::new();
    let mut property_rgroups = Vec::new();
    for (line_number, line) in lines.by_ref() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 2 && fields[0] == "M" && fields[1] == "CHG" {
            let count = fields
                .get(2)
                .and_then(|value| value.parse::<usize>().ok())
                .ok_or_else(|| MolParseError::InvalidPropertyLine {
                    line: line_number,
                    detail: "M CHG lacks a valid entry count".to_string(),
                })?;
            let expected_fields = 3usize.checked_add(count.saturating_mul(2)).ok_or_else(|| {
                MolParseError::InvalidPropertyLine {
                    line: line_number,
                    detail: "M CHG entry count overflows".to_string(),
                }
            })?;
            if fields.len() != expected_fields {
                return Err(MolParseError::InvalidPropertyLine {
                    line: line_number,
                    detail: format!(
                        "M CHG declares {count} entries but has {} values",
                        fields.len().saturating_sub(3)
                    ),
                });
            }
            // The exact field-count check above guarantees pairs; use the
            // MSRV-compatible iterator rather than a newer slice API.
            for pair in fields[3..].chunks(2) {
                let atom_id = pair[0]
                    .parse::<usize>()
                    .ok()
                    .filter(|id| *id > 0 && *id <= natoms)
                    .ok_or_else(|| MolParseError::InvalidPropertyLine {
                        line: line_number,
                        detail: format!("M CHG atom id '{}' is outside 1..={natoms}", pair[0]),
                    })?;
                let charge =
                    pair[1]
                        .parse::<i8>()
                        .map_err(|_| MolParseError::InvalidPropertyLine {
                            line: line_number,
                            detail: format!("M CHG charge '{}' is not an i8", pair[1]),
                        })?;
                property_charges.push((AtomIdx((atom_id - 1) as u32), charge));
            }
        } else if fields.len() >= 2 && fields[0] == "M" && fields[1] == "ISO" {
            let count = fields
                .get(2)
                .and_then(|value| value.parse::<usize>().ok())
                .ok_or_else(|| MolParseError::InvalidPropertyLine {
                    line: line_number,
                    detail: "M ISO lacks a valid entry count".to_string(),
                })?;
            if fields.len() != 3 + count.saturating_mul(2) {
                return Err(MolParseError::InvalidPropertyLine {
                    line: line_number,
                    detail: format!(
                        "M ISO declares {count} entries but has {} values",
                        fields.len().saturating_sub(3)
                    ),
                });
            }
            for pair in fields[3..].chunks(2) {
                let atom_id = pair[0]
                    .parse::<usize>()
                    .ok()
                    .filter(|id| *id > 0 && *id <= natoms)
                    .ok_or_else(|| MolParseError::InvalidPropertyLine {
                        line: line_number,
                        detail: format!("M ISO atom id '{}' is outside 1..={natoms}", pair[0]),
                    })?;
                let isotope =
                    pair[1]
                        .parse::<u16>()
                        .map_err(|_| MolParseError::InvalidPropertyLine {
                            line: line_number,
                            detail: format!("M ISO mass '{}' is not a u16", pair[1]),
                        })?;
                if isotope == 0 {
                    return Err(MolParseError::InvalidPropertyLine {
                        line: line_number,
                        detail: "M ISO mass must be greater than zero".to_string(),
                    });
                }
                property_isotopes.push((AtomIdx((atom_id - 1) as u32), isotope));
            }
        } else if fields.len() >= 2 && fields[0] == "M" && fields[1] == "RGP" {
            let count = fields
                .get(2)
                .and_then(|value| value.parse::<usize>().ok())
                .ok_or_else(|| MolParseError::InvalidPropertyLine {
                    line: line_number,
                    detail: "M RGP lacks a valid entry count".to_string(),
                })?;
            if fields.len() != 3 + count.saturating_mul(2) {
                return Err(MolParseError::InvalidPropertyLine {
                    line: line_number,
                    detail: format!(
                        "M RGP declares {count} entries but has {} values",
                        fields.len().saturating_sub(3)
                    ),
                });
            }
            for pair in fields[3..].chunks(2) {
                let atom_id = pair[0]
                    .parse::<usize>()
                    .ok()
                    .filter(|id| *id > 0 && *id <= natoms)
                    .ok_or_else(|| MolParseError::InvalidPropertyLine {
                        line: line_number,
                        detail: format!("M RGP atom id '{}' is outside 1..={natoms}", pair[0]),
                    })?;
                let number = pair[1]
                    .parse::<u16>()
                    .ok()
                    .and_then(chematic_core::RGroupLabel::numbered)
                    .ok_or_else(|| MolParseError::InvalidPropertyLine {
                        line: line_number,
                        detail: format!("M RGP number '{}' must be in 1..=9999", pair[1]),
                    })?;
                property_rgroups.push((AtomIdx((atom_id - 1) as u32), number));
            }
        }
        if line.trim_start().starts_with("M  END") {
            break;
        }
    }

    let mut mol = builder.build();
    flag_aromatic_bond_atoms(&mut mol);
    for (atom, charge) in property_charges {
        mol.set_charge(atom, charge);
    }
    for (atom, isotope) in property_isotopes {
        mol.set_isotope(atom, Some(isotope));
    }
    // A valence field fixes the atom's hydrogens (as RDKit reads it); skipped
    // on aromatic-bond atoms, whose integer valence the bonds do not give.
    for (idx, valence) in valence_fields {
        let aromatic_bond = mol
            .neighbors(idx)
            .any(|(_, b)| mol.bond(b).order == BondOrder::Aromatic);
        let bonds = chematic_core::valence::bond_order_sum(&mol, idx);
        if !aromatic_bond && valence >= bonds {
            mol.set_hydrogen_count(idx, Some(valence - bonds));
        }
    }
    for (atom, label) in property_rgroups {
        mol.set_r_group(atom, label);
    }

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
    let stereo_diagnostics =
        apply_ctab_local_parity_from_wedges_with_diagnostics(&mut mol, &coords);
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
    let mut out = String::new();
    write_mol_with_coords_into(&mut out, mol, metadata, coords);
    out
}

/// Standard CTfile atom symbol for real atoms and the bounded pseudoatom
/// contract supported by this crate.
fn mdl_atom_symbol(mol: &Molecule, idx: AtomIdx, atom: &Atom) -> &'static str {
    if !atom.wildcard {
        return atom.element.symbol();
    }
    match mol.r_group_label(idx).and_then(|label| label.number()) {
        Some(_) => "R#",
        None if mol.r_group_label(idx).is_some() => "R",
        None => "*",
    }
}

fn write_v2000_rgroups(out: &mut String, mol: &Molecule) {
    use std::fmt::Write as _;

    let entries = mol
        .atoms()
        .filter_map(|(idx, _)| {
            mol.r_group_label(idx)
                .and_then(|label| label.number())
                .map(|number| (idx.0 + 1, number))
        })
        .collect::<Vec<_>>();
    // CTfile property records conventionally carry at most eight pairs.
    for chunk in entries.chunks(8) {
        write!(out, "M  RGP{:>3}", chunk.len()).expect("writing to String cannot fail");
        for &(atom_id, number) in chunk {
            write!(out, "{:>4}{:>4}", atom_id, number).expect("writing to String cannot fail");
        }
        out.push('\n');
    }
}

#[inline]
/// Append `x` as `format!("{x:>width$.4}")` would, without the general
/// float formatter (its exact fallback dominated writing laid-out
/// molecules). Values near a rounding tie, very large values and
/// non-finite values take the formatter, so the bytes are always the same.
pub(crate) fn push_fixed4(out: &mut String, x: f64, width: usize) {
    use std::fmt::Write as _;
    let scaled = x.abs() * 10_000.0;
    let frac = scaled - scaled.floor();
    if !x.is_finite() || scaled >= 1e15 || (frac - 0.5).abs() < 1e-6 {
        write!(out, "{x:>width$.4}").expect("writing to String cannot fail");
        return;
    }
    let units = scaled.round() as u64;
    let mut buf = [0u8; 24];
    let mut pos = buf.len();
    let mut push = |b: u8| {
        pos -= 1;
        buf[pos] = b;
    };
    let (mut int, mut dec) = (units / 10_000, units % 10_000);
    for _ in 0..4 {
        push(b'0' + (dec % 10) as u8);
        dec /= 10;
    }
    push(b'.');
    loop {
        push(b'0' + (int % 10) as u8);
        int /= 10;
        if int == 0 {
            break;
        }
    }
    if x.is_sign_negative() {
        push(b'-');
    }
    let digits = &buf[pos..];
    for _ in digits.len()..width {
        out.push(' ');
    }
    out.push_str(std::str::from_utf8(digits).expect("ASCII digits"));
}

fn push_right_aligned_u32(out: &mut String, mut value: u32, width: usize) {
    // V2000 counts and bond fields are almost always three columns wide.
    // Append the entire field in one operation, retaining the general path
    // for widths/values outside the fixed-width format.
    if width == 3 && value < 1000 {
        static FIELDS: [[u8; 3]; 1000] = {
            let mut fields = [[b' '; 3]; 1000];
            let mut n = 0;
            while n < fields.len() {
                if n >= 100 {
                    fields[n][0] = b'0' + (n / 100) as u8;
                }
                if n >= 10 {
                    fields[n][1] = b'0' + (n / 10 % 10) as u8;
                }
                fields[n][2] = b'0' + (n % 10) as u8;
                n += 1;
            }
            fields
        };
        // Validate ASCII at compile time, not at every emitted field.
        static TEXT: [&str; 1000] = {
            let mut text = [""; 1000];
            let mut n = 0;
            while n < text.len() {
                text[n] = match std::str::from_utf8(&FIELDS[n]) {
                    Ok(s) => s,
                    Err(_) => panic!("ASCII decimal field"),
                };
                n += 1;
            }
            text
        };
        out.push_str(TEXT[value as usize]);
        return;
    }
    let mut digits = [0_u8; 10];
    let mut start = digits.len();
    loop {
        start -= 1;
        digits[start] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }
    let len = digits.len() - start;
    for _ in len..width {
        out.push(' ');
    }
    for &digit in &digits[start..] {
        out.push(digit as char);
    }
}

#[inline]
fn push_right_aligned_i16(out: &mut String, value: i16, width: usize) {
    let negative = value < 0;
    let magnitude = value.unsigned_abs() as u32;
    let digits = if magnitude < 10 {
        1
    } else if magnitude < 100 {
        2
    } else if magnitude < 1_000 {
        3
    } else if magnitude < 10_000 {
        4
    } else {
        5
    };
    for _ in (digits + usize::from(negative))..width {
        out.push(' ');
    }
    if negative {
        out.push('-');
    }
    push_right_aligned_u32(out, magnitude, digits);
}

/// Serialize a MOL V2000 block into a reusable caller-owned buffer.
///
/// This is the allocation-amortized implementation used by streaming
/// writers. The buffer is cleared first; on return it contains exactly the
/// same bytes as [`write_mol_with_coords`].
#[doc(hidden)]
pub fn write_mol_with_coords_into(
    out: &mut String,
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
) {
    write_v2000_reporting(out, mol, metadata, coords, &[]);
}

/// [`write_mol_with_coords`] that also reports the stereo the block does not
/// carry. An empty `coords` lays the molecule out (as [`write_mol`] does).
///
/// ```
/// let mol = chematic_smiles::parse("C[C@H](N)O").unwrap();
/// let (block, loss) =
///     chematic_mol::write_mol_with_stereo_report(&mol, &Default::default(), &[]);
/// assert!(loss.is_empty());
/// assert!(block.contains("V2000"));
/// ```
pub fn write_mol_with_stereo_report(
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
) -> (String, MolStereoLoss) {
    let mut out = String::new();
    let loss = write_v2000_reporting(&mut out, mol, metadata, coords, &[]);
    (out, loss)
}

/// Stereo of a molecule that a written 2D MOL block does not carry. A reader
/// sees each listed centre as unspecified and each listed double bond as
/// "either" (V2000 stereo 3 / V3000 `CFG=2`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MolStereoLoss {
    /// Tetrahedral centres written without a wedge or hash: no wedge on the
    /// coordinates read back as the declared configuration (bridged and cage
    /// layouts, bonds drawn in nearly one direction).
    pub centres: Vec<AtomIdx>,
    /// Declared E/Z double bonds whose geometry the coordinates could not be
    /// given (ring bonds); written "either".
    pub double_bonds: Vec<chematic_core::BondIdx>,
    /// Square-planar centres: a 2D block has no field for them.
    pub non_tetrahedral_centres: Vec<AtomIdx>,
    /// Enhanced stereo groups other than absolute: V2000 has no field for
    /// them (V3000 writes them).
    pub stereo_groups_dropped: bool,
}

impl MolStereoLoss {
    /// Whether the block carries all of the molecule's stereo.
    pub fn is_empty(&self) -> bool {
        self.centres.is_empty()
            && self.double_bonds.is_empty()
            && self.non_tetrahedral_centres.is_empty()
            && !self.stereo_groups_dropped
    }

    pub(crate) fn from_depiction(
        mol: &Molecule,
        depiction: Option<&crate::stereo_depiction::StereoDepiction>,
        format: MolFormat,
    ) -> Self {
        let mut centres: Vec<AtomIdx> = depiction
            .map(|d| d.unexpressed_centres.clone())
            .unwrap_or_default();
        // A tetrahedral tag without a neighbour order has no configuration
        // to draw.
        centres.extend(mol.atoms().filter_map(|(idx, atom)| {
            (atom.chirality.is_tetrahedral() && mol.stereo_neighbor_order(idx).is_none())
                .then_some(idx)
        }));
        centres.sort_by_key(|a| a.0);
        centres.dedup();
        MolStereoLoss {
            centres,
            double_bonds: depiction
                .map(|d| d.unexpressed_double_bonds.clone())
                .unwrap_or_default(),
            non_tetrahedral_centres: mol
                .atoms()
                .filter(|(_, a)| {
                    a.chirality != chematic_core::Chirality::None && !a.chirality.is_tetrahedral()
                })
                .map(|(idx, _)| idx)
                .collect(),
            stereo_groups_dropped: format == MolFormat::V2000
                && mol
                    .stereo_groups()
                    .iter()
                    .any(|g| g.kind != chematic_core::StereoGroupKind::Absolute),
        }
    }
}

impl core::fmt::Display for MolStereoLoss {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let atoms = |v: &[AtomIdx]| v.iter().map(|a| a.0.to_string()).collect::<Vec<_>>();
        let mut parts = Vec::new();
        if !self.centres.is_empty() {
            parts.push(format!(
                "tetrahedral centres without a wedge (atoms {})",
                atoms(&self.centres).join(", ")
            ));
        }
        if !self.double_bonds.is_empty() {
            let bonds: Vec<String> = self.double_bonds.iter().map(|b| b.0.to_string()).collect();
            parts.push(format!(
                "E/Z written as either (bonds {})",
                bonds.join(", ")
            ));
        }
        if !self.non_tetrahedral_centres.is_empty() {
            parts.push(format!(
                "square-planar centres (atoms {})",
                atoms(&self.non_tetrahedral_centres).join(", ")
            ));
        }
        if self.stereo_groups_dropped {
            parts.push("enhanced stereo groups".to_string());
        }
        if parts.is_empty() {
            write!(f, "no stereo lost")
        } else {
            write!(f, "MOL block loses stereo: {}", parts.join("; "))
        }
    }
}

/// [`write_mol_with_stereo_report`] on chematic's own layout (the stereo
/// depiction when the molecule has stereo, the plain layout otherwise), as
/// the Python and WASM `to_mol_block` write it.
///
/// A layout draws every double bond cis or trans, and readers take drawn
/// geometry as E/Z. Double bonds that could carry E/Z but have none
/// declared are therefore written "either" (stereo 3), so the block never
/// gains a configuration the molecule does not have. Caller-supplied
/// coordinates ([`write_mol_with_coords`]) are data and stay unmarked.
///
/// ```
/// let mol = chematic_smiles::parse("CC=CC").unwrap();
/// let (block, loss) = chematic_mol::write_laid_out_mol(&mol, &Default::default());
/// assert!(loss.is_empty());
/// assert!(block.contains("  2  3  2  3\n"));
/// ```
pub fn write_laid_out_mol(mol: &Molecule, metadata: &MolMetadata) -> (String, MolStereoLoss) {
    let mut out = String::new();
    let loss = write_laid_out_mol_into(&mut out, mol, metadata);
    (out, loss)
}

/// [`write_laid_out_mol`] into a reusable buffer (cleared first).
#[doc(hidden)]
pub fn write_laid_out_mol_into(
    out: &mut String,
    mol: &Molecule,
    metadata: &MolMetadata,
) -> MolStereoLoss {
    let coords = crate::stereo_depiction::mol_block_coords(mol);
    let either = if crate::stereo_depiction::needs_stereo_depiction(mol) {
        // The stereo depiction marks its own unspecified double bonds.
        Vec::new()
    } else {
        crate::stereo_depiction::undeclared_stereo_double_bonds(mol)
    };
    write_v2000_reporting(out, mol, metadata, &coords, &either)
}

/// An SDF record written as [`write_laid_out_mol`] writes its MOL block.
pub fn write_laid_out_sdf_record_into(
    out: &mut String,
    mol: &Molecule,
    meta: &MolMetadata,
    props: &std::collections::HashMap<String, String>,
) {
    write_laid_out_mol_into(out, mol, meta);
    append_sd_fields_and_delimiter(out, props);
}

/// The V2000 valence field (vvv) for atom `idx`: 0 (unset) when a reader
/// infers the atom's hydrogens from its bonds, otherwise the total valence
/// (15 for zero), so `[I]`, `[S](=O)=O` or `[SH2]` keep their hydrogen
/// count (they used to read back as `[IH]`, `[SH]`, ...).
fn valence_field(mol: &Molecule, idx: AtomIdx, atom: &Atom) -> u8 {
    let Some(h) = atom.hydrogen_count else {
        return 0;
    };
    if atom.wildcard
        || chematic_core::valence::valence_inferred_hcount(mol, idx) == h
        || mol
            .neighbors(idx)
            .any(|(_, b)| mol.bond(b).order == BondOrder::Aromatic)
    {
        return 0;
    }
    match chematic_core::valence::bond_order_sum(mol, idx).saturating_add(h) {
        0 => 15,
        v if v < 15 => v,
        _ => 0,
    }
}

/// `either`: double bonds written with stereo 3 when there is no stereo
/// depiction (see [`write_laid_out_mol`]).
fn write_v2000_reporting(
    out: &mut String,
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
    either: &[chematic_core::BondIdx],
) -> MolStereoLoss {
    // A MOL block is dominated by fixed-width atom and bond rows.  Reserve
    // the common-size output up front so serialization does not repeatedly
    // grow and copy the String for every row (the old empty String was a
    // noticeable cost when writing many small records).
    out.clear();
    out.reserve(128 + mol.atom_count() * 80 + mol.bond_count() * 16);

    // Header lines 1–3
    out.push_str(&metadata.name);
    out.push('\n');
    out.push_str("  chematic\n");
    out.push_str(&metadata.comment);
    out.push('\n');

    // A molecule with stereo is written on coordinates that express it, with
    // one wedge or hash per tetrahedral centre (see `stereo_depiction`).
    let depiction = crate::stereo_depiction::needs_stereo_depiction(mol)
        .then(|| crate::stereo_depiction::stereo_depiction(mol, coords));
    let coords: &[(f64, f64)] = depiction.as_ref().map_or(coords, |d| &d.coords);

    // Counts line (line 4)
    let natoms = mol.atom_count();
    let nbonds = mol.bond_count();
    push_right_aligned_u32(out, natoms as u32, 3);
    push_right_aligned_u32(out, nbonds as u32, 3);
    out.push_str("  0  0  0  0  0  0  0  0999 V2000\n");

    // Atom block
    for (idx, atom) in mol.atoms() {
        let sym = mdl_atom_symbol(mol, idx, atom);
        let charge_code = encode_charge(atom.charge);
        let mass_difference = encode_mass_difference(atom.element, atom.isotope).unwrap_or(0);
        let atom_map = atom.atom_map.unwrap_or(0);
        let valence = valence_field(mol, idx, atom);
        if let Some(&(x, y)) = coords.get(idx.0 as usize) {
            push_fixed4(out, x, 10);
            push_fixed4(out, y, 10);
            // `"    0.0000 {:<3}{:>2}{:>3}  0  0  0  0  0  0  0{:>3}  0\n"`,
            // without the formatter (a tenth of laying out and writing).
            out.push_str("    0.0000 ");
            out.push_str(sym);
            for _ in sym.len()..3 {
                out.push(' ');
            }
            push_right_aligned_i16(out, mass_difference, 2);
            push_right_aligned_u32(out, charge_code as u32, 3);
            out.push_str("  0  0  0");
            push_right_aligned_u32(out, u32::from(valence), 3);
            out.push_str("  0  0  0");
            push_right_aligned_u32(out, atom_map as u32, 3);
            out.push_str("  0\n");
        } else {
            // Serialization-only SDF output overwhelmingly has no coordinate
            // array. Avoid running the float formatter three times per atom
            // for a byte sequence that is known at compile time.
            out.push_str("    0.0000    0.0000    0.0000 ");
            out.push_str(sym);
            for _ in sym.len()..3 {
                out.push(' ');
            }
            if mass_difference == 0 && charge_code == 0 {
                out.push_str(" 0  0");
            } else {
                push_right_aligned_i16(out, mass_difference, 2);
                push_right_aligned_u32(out, charge_code as u32, 3);
            }
            out.push_str("  0  0  0");
            push_right_aligned_u32(out, u32::from(valence), 3);
            out.push_str("  0  0  0");
            push_right_aligned_u32(out, atom_map as u32, 3);
            out.push_str("  0\n");
        }
    }

    // Bond block
    let kekule = kekule_orders_for_writing(mol);
    for (bond_idx, bond) in mol.bonds() {
        let a1 = bond.atom1.0 + 1; // convert to 1-based
        let a2 = bond.atom2.0 + 1;
        let order = kekule
            .as_ref()
            .and_then(|k| k.get(&bond_idx).copied())
            .unwrap_or(bond.order);
        let btype = match order {
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
        // Stereo field: with a stereo depiction, the chosen wedge/hash drawn
        // from its stereocentre (SMILES `/`/`\` markers are expressed by the
        // coordinates, not by a stereo code); otherwise a wedge/hash already
        // on the bond, so a re-parse recovers the same local parity.
        let (a1, a2, stereo) = match &depiction {
            Some(d) => match d.wedges.get(&bond_idx) {
                Some(w) => {
                    let (start, end) = if w.start == bond.atom1 {
                        (a1, a2)
                    } else {
                        (a2, a1)
                    };
                    (start, end, if w.order == BondOrder::Up { 1 } else { 6 })
                }
                // A stereo double bond these coordinates cannot draw (in a
                // ring) is written as "either" rather than with a wrong
                // geometry.
                None if d.unexpressed_double_bonds.contains(&bond_idx)
                    || d.unspecified_double_bonds.contains(&bond_idx) =>
                {
                    (a1, a2, 3)
                }
                None => (a1, a2, 0),
            },
            None => (
                a1,
                a2,
                match bond.order {
                    BondOrder::Up => 1,
                    BondOrder::Down => 6,
                    _ if either.contains(&bond_idx) => 3,
                    _ => 0,
                },
            ),
        };
        push_right_aligned_u32(out, a1, 3);
        push_right_aligned_u32(out, a2, 3);
        push_right_aligned_u32(out, btype, 3);
        push_right_aligned_u32(out, stereo, 3);
        out.push('\n');
    }

    write_v2000_rgroups(out, mol);

    // Terminator
    out.push_str("M  END\n");
    MolStereoLoss::from_depiction(mol, depiction.as_ref(), MolFormat::V2000)
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
        let sym = mdl_atom_symbol(mol, idx, atom);
        let charge_code = encode_charge(atom.charge);
        let mass_difference = encode_mass_difference(atom.element, atom.isotope).unwrap_or(0);
        let atom_map = atom.atom_map.unwrap_or(0);
        let p = conformer
            .points
            .get(idx.0 as usize)
            .copied()
            .unwrap_or(Point3::zero());
        push_fixed4(&mut out, p.x, 10);
        push_fixed4(&mut out, p.y, 10);
        push_fixed4(&mut out, p.z, 10);
        out.push_str(&format!(
            " {:<3}{:>2}{:>3}  0  0  0  0  0  0  0{:>3}  0\n",
            sym, mass_difference, charge_code, atom_map,
        ));
    }

    let kekule = kekule_orders_for_writing(mol);
    for (bond_idx, bond) in mol.bonds() {
        let a1 = bond.atom1.0 + 1;
        let a2 = bond.atom2.0 + 1;
        let order = kekule
            .as_ref()
            .and_then(|k| k.get(&bond_idx).copied())
            .unwrap_or(bond.order);
        let btype = match order {
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

    write_v2000_rgroups(&mut out, mol);
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
    let mut out = String::new();
    write_sdf_record_into(&mut out, mol, meta, coords, props);
    out
}

/// Serialize one SDF record into a reusable caller-owned buffer.
///
/// The buffer is cleared first. Reusing it across records amortizes the MOL
/// block allocation while preserving byte-for-byte output compatibility with
/// [`write_sdf_record`].
#[doc(hidden)]
pub fn write_sdf_record_into(
    out: &mut String,
    mol: &Molecule,
    meta: &MolMetadata,
    coords: &[(f64, f64)],
    props: &std::collections::HashMap<String, String>,
) {
    write_mol_with_coords_into(out, mol, meta, coords);
    append_sd_fields_and_delimiter(out, props);
}

/// Append visible SD fields and the record delimiter without allocating a
/// temporary formatted `String` for every property.  Reserving the aggregate
/// payload once also avoids repeated growth when records carry many fields.
fn append_sd_fields_and_delimiter(
    out: &mut String,
    props: &std::collections::HashMap<String, String>,
) {
    use std::fmt::Write as _;

    let additional = props
        .iter()
        .filter(|(key, _)| !key.starts_with('_'))
        .map(|(key, value)| key.len().saturating_add(value.len()).saturating_add(7))
        .sum::<usize>()
        .saturating_add(5);
    out.reserve(additional);
    for (key, value) in props {
        if !key.starts_with('_') {
            write!(out, "> <{key}>\n{value}\n\n").expect("writing to String cannot fail");
        }
    }
    out.push_str("$$$$\n");
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
    append_sd_fields_and_delimiter(&mut out, props);
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
    append_sd_fields_and_delimiter(&mut out, props);
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
