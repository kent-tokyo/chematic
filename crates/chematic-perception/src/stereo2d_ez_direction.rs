//! CIP-independent E/Z (double-bond cis/trans) *direction* from 2D
//! coordinates -- the double-bond counterpart to [`crate::stereo2d_local`]'s
//! tetrahedral parity stage.
//!
//! [`apply_ez_directions_from_2d`] computes writer-consumable SMILES `/`/`\`
//! direction for the substituent bonds flanking a stereogenic double bond,
//! straight from 2D geometry and topology -- mirroring RDKit's
//! `MolOps::detectBondStereochemistry`/`setDoubleBondNeighborDirections`
//! (see `docs/rfcs/stereo2d_reader_integration_rfc.md`'s RDKit source audit).
//! [`crate::stereo2d::assign_ez_from_2d`] already computes a CIP E/Z
//! *label* (`Atom.cip_code`) from the same kind of input, but a label is not
//! a direction: nothing in `chematic_smiles::write`/`canonical_smiles` reads
//! `cip_code`. This module produces the missing Stage-1-equivalent value a
//! SMILES writer can actually consume, and -- like the tetrahedral stage --
//! never reads `Atom.cip_code` and is never gated on the legacy CIP-label
//! engine succeeding.
//!
//! **Never mutates `BondOrder::Up`/`Down` on a bond's own `order` field.**
//! That field is the MOL/SDF reader's raw wedge/hash depiction (tetrahedral
//! stereocenter notation perceived by [`crate::stereo2d_local`] into
//! `Atom.chirality`) and must never be repurposed for E/Z. Every direction
//! this module writes goes through [`chematic_core::Molecule::set_bond_direction`]
//! -- the same side channel already used to stash an aromatic-bond-adjacent
//! E/Z direction, generalized here to a plain `Single`-order bond too -- a
//! strictly separate storage slot from `order`. A candidate substituent bond
//! whose *own* `order` is already `Up`/`Down` (an existing wedge) is never
//! usable as an E/Z carrier: see [`EzDirectionRejectionReason::CarrierConflict`].
//!
//! ## Scope limits (deliberate, see `docs/rfcs/stereo2d_reader_integration_rfc.md`)
//!
//! - Joint canonical carrier resolution across independently-stereogenic
//!   double bonds that share one physical candidate bond (issue #149) is out
//!   of scope. When two double bonds computed FROM RAW GEOMETRY ALONE (never
//!   from each other's output) happen to require the *same* literal
//!   direction on a shared bond, both succeed (this is the ordinary,
//!   expected shape of a conjugated diene, e.g. `(2E,4E)-hexa-2,4-diene`).
//!   When they disagree, *both* double bonds relying on that bond are
//!   rejected with [`EzDirectionRejectionReason::CarrierConflict`] rather
//!   than letting bond-index or processing order pick an arbitrary winner.
//! - No retry/search across double bonds is attempted beyond the natural
//!   sibling fallback within one alkene end (see [`resolve_end`]). A
//!   fixed-point joint resolver is exactly the issue #149 problem this PR
//!   does not solve.
//!
//! ## Aromatic (Kekulé) ring bonds are never candidates
//!
//! A MOL/SDF reader does not run aromaticity perception by default (per
//! `CLAUDE.md`: "Kekulé input requires explicit `apply_aromaticity`"), so a
//! benzene ring read from a Kekulized MOL file arrives as plain alternating
//! `Single`/`Double` bonds with `atom.aromatic == false` throughout --
//! structurally indistinguishable, to a naive per-bond check, from a real
//! open-chain conjugated diene. Found empirically (not anticipated in the
//! original fixture set): a broad-corpus run surfaced spurious
//! `CarrierConflict` cascades and, worse, wrongly-signed directions on
//! genuine adjacent alkenes, traced to this module attempting to assign E/Z
//! to Kekulé ring bonds that have no cis/trans isomerism at all (a ring's
//! geometry is fixed; RDKit's own pipeline never reaches this case because
//! its `sanitizeMol` step re-types these bonds as `AROMATIC` *before*
//! `detectBondStereochemistry` runs, so a literal `Double`-typed bond check
//! naturally excludes them there). Mirrored here via a one-time,
//! non-mutating [`crate::aromaticity::assign_aromaticity`] query (not
//! `apply_aromaticity`, which would return a new `Molecule` and silently
//! change the caller's aromatic flags as a side effect of adding E/Z
//! direction -- an unrelated, invasive change this module must not make):
//! any bond the Hückel model would classify as aromatic is excluded
//! up front, exactly like [`EzOutcome::NotRequested`] for a terminal alkene.

use std::collections::{HashMap, HashSet};

use crate::aromaticity::{AromaticityModel, assign_aromaticity};
use crate::cip_priority::compare_branches;
use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

/// Tolerance below which a 2D cross product (substituent vs. double-bond
/// axis) or a bond-length vector is treated as exactly zero.
const AXIS_EPS: f64 = 1e-6;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Why a candidate stereogenic double bond's direction was rejected by
/// [`apply_ez_directions_from_2d_with_diagnostics`] -- never emitted for a
/// double bond that isn't stereogenic in the first place (terminal alkene,
/// carbonyl/heteroatom terminus, or topologically-equivalent substituents at
/// either end): those are silently [`NotRequested`](self), matching how
/// [`crate::stereo2d_local::StereoDiagnostic`] only reports on a wedge/hash
/// bond that was actually drawn.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EzDirectionRejectionReason {
    /// The double bond's own, or a chosen substituent's, 2D coordinate is
    /// missing from `coords`.
    MissingCoordinate,
    /// A required coordinate is present but not finite (NaN or infinite).
    NonFiniteCoordinate,
    /// Zero-length double bond, or every candidate substituent at some end
    /// lies exactly on the double-bond axis (collinear geometry).
    DegenerateGeometry,
    /// The file explicitly marked this double bond's stereo as unspecified
    /// (MDL V2000 double-bond stereo code 3 / V3000 bond `CFG=2`, the
    /// "either"/crossed-bond convention) -- 2D coordinates are not consulted
    /// for direction in this case, matching RDKit (verified against a live
    /// RDKit 2026.03.3 oracle: an explicit "either" flag suppresses
    /// direction regardless of what the raw coordinates would otherwise
    /// imply).
    ExplicitlyUnspecified,
    /// A cumulated pi system (allene/cumulene: an endpoint participates in
    /// more than one double bond), or an endpoint with more than two
    /// non-double-bond substituents (not a valid trigonal alkene carbon).
    UnsupportedTopology,
    /// Every substituent bond usable as this double bond's carrier at some
    /// end is already claimed by something else that must not be
    /// overwritten or reinterpreted: an existing wedge/hash bond (a
    /// different physical meaning entirely), an unrelated pre-existing
    /// `bond_direction` stash, or a *different* double bond's own
    /// independently-computed, disagreeing requirement on the same shared
    /// bond (the issue #149 shape -- see module docs).
    CarrierConflict,
}

/// One rejected stereogenic double bond, as returned by
/// [`apply_ez_directions_from_2d_with_diagnostics`]. `bond` identifies the
/// double bond itself (not the substituent/carrier bond).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EzDirectionDiagnostic {
    pub bond: BondIdx,
    pub reason: EzDirectionRejectionReason,
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Apply [`apply_ez_directions_from_2d_with_diagnostics`], discarding
/// diagnostics.
pub fn apply_ez_directions_from_2d(mol: &mut Molecule, coords: &[(f64, f64)]) {
    apply_ez_directions_from_2d_with_diagnostics(mol, coords);
}

/// Compute and write E/Z bond directions for every stereogenic double bond
/// in `mol`, from `coords` alone (no CIP ranking, no `Atom.cip_code`).
///
/// `coords[i]` is the `(x, y)` position of atom `i`, in the same units/frame
/// used elsewhere in `chematic-perception`'s 2D stereo perception.
///
/// Returns one [`EzDirectionDiagnostic`] per double bond whose direction
/// could not be safely determined (see [`EzDirectionRejectionReason`]).
/// Never returns a diagnostic for a double bond that simply isn't
/// stereogenic (terminal alkene, carbonyl, or topologically-equivalent
/// substituents at either end) -- those are silently skipped, matching
/// `apply_local_parity_from_wedges_with_diagnostics`'s "nothing was
/// requested, nothing to reject" convention.
///
/// All-or-nothing per double bond: a bond either gets both of its carrier
/// directions written, or none at all.
pub fn apply_ez_directions_from_2d_with_diagnostics(
    mol: &mut Molecule,
    coords: &[(f64, f64)],
) -> Vec<EzDirectionDiagnostic> {
    apply_ez_directions_from_2d_ex(mol, coords, &HashSet::new())
}

/// Reader-integration entry point: same as
/// [`apply_ez_directions_from_2d_with_diagnostics`], but additionally takes
/// the set of double bonds a file explicitly marked as unspecified stereo
/// (MDL V2000 double-bond stereo code 3 / V3000 bond `CFG=2`) -- information
/// that exists only at parse time (there is no MOL/SDF-specific field on
/// `Molecule` to carry it, and this PR does not add one; see
/// `docs/rfcs/stereo2d_reader_integration_rfc.md` for why a Molecule-level
/// extension was avoided). A reader with no such bonds -- or a caller with
/// no reader-level information at all -- passes an empty set, exactly what
/// the two convenience entry points above do.
pub fn apply_ez_directions_from_2d_ex(
    mol: &mut Molecule,
    coords: &[(f64, f64)],
    explicitly_unspecified: &HashSet<BondIdx>,
) -> Vec<EzDirectionDiagnostic> {
    let double_bonds: Vec<BondIdx> = mol
        .bonds()
        .filter(|(_, b)| b.order == BondOrder::Double)
        .map(|(bidx, _)| bidx)
        .collect();

    // One-time, non-mutating aromaticity query (see module docs) so a
    // Kekulized aromatic ring bond -- read with `atom.aromatic == false`
    // since the reader never auto-perceives aromaticity -- is never treated
    // as an E/Z candidate.
    let aromaticity = assign_aromaticity(mol);

    // Pre-pass: a branch point (2-substituent alkene end) with a
    // substituent that is itself part of a DIFFERENT double bond poisons
    // BOTH double bonds, not just its own -- see
    // `poisoned_by_branch_ambiguity`'s doc comment for why rejecting only
    // the branch point's own bond is not sufficient (OpenSMILES `/`/`\`
    // markers are read by plain adjacency, not by which system "intended"
    // them, so the neighboring double bond's own, otherwise-legitimate
    // marker unavoidably leaks into the branch point's perceived stereo
    // too).
    let poisoned = poisoned_by_branch_ambiguity(mol, &aromaticity);

    // Phase 1: classify every double bond independently, from raw geometry
    // and topology alone -- never from another double bond's outcome, so
    // the result for one bond can't depend on which order this loop visits
    // bonds in.
    let mut outcomes: HashMap<BondIdx, EzOutcome> = HashMap::with_capacity(double_bonds.len());
    for &bidx in &double_bonds {
        let outcome = if poisoned.contains(&bidx) {
            EzOutcome::Rejected(EzDirectionRejectionReason::CarrierConflict)
        } else {
            classify_double_bond(mol, coords, bidx, explicitly_unspecified, &aromaticity)
        };
        outcomes.insert(bidx, outcome);
    }

    // Phase 2: a physical substituent bond can be the chosen carrier for TWO
    // different double bonds at once (the ordinary conjugated-diene shape,
    // see module docs). Group every Assigned outcome's claims by physical
    // bond and check they agree; any bond with disagreeing claims takes
    // every double bond relying on it down with it, symmetrically -- never
    // "first-claimed wins" (that would make the result depend on
    // `double_bonds` iteration order, i.e. bond-index/parse order, which the
    // spec explicitly rules out).
    let mut claims: HashMap<BondIdx, Vec<(BondIdx, BondOrder)>> = HashMap::new();
    for (&db, outcome) in &outcomes {
        if let EzOutcome::Assigned {
            carrier_a1,
            carrier_a2,
        } = outcome
        {
            claims
                .entry(carrier_a1.0)
                .or_default()
                .push((db, carrier_a1.1));
            claims
                .entry(carrier_a2.0)
                .or_default()
                .push((db, carrier_a2.1));
        }
    }
    let mut conflicted: HashSet<BondIdx> = HashSet::new();
    for claimants in claims.values() {
        if claimants.len() > 1 {
            let first = claimants[0].1;
            if claimants.iter().any(|&(_, v)| v != first) {
                for &(db, _) in claimants {
                    conflicted.insert(db);
                }
            }
        }
    }

    // Phase 3: commit. Only an Assigned outcome not caught by Phase 2 ever
    // writes to `mol`; everything else contributes at most one diagnostic.
    let mut diagnostics = Vec::with_capacity(double_bonds.len());
    for &bidx in &double_bonds {
        match outcomes.remove(&bidx).expect("classified in phase 1") {
            EzOutcome::NotRequested => {}
            EzOutcome::Rejected(reason) => {
                diagnostics.push(EzDirectionDiagnostic { bond: bidx, reason });
            }
            EzOutcome::Assigned {
                carrier_a1,
                carrier_a2,
            } => {
                if conflicted.contains(&bidx) {
                    diagnostics.push(EzDirectionDiagnostic {
                        bond: bidx,
                        reason: EzDirectionRejectionReason::CarrierConflict,
                    });
                } else {
                    mol.set_bond_direction(carrier_a1.0, carrier_a1.1);
                    mol.set_bond_direction(carrier_a2.0, carrier_a2.1);
                }
            }
        }
    }
    diagnostics
}

// ---------------------------------------------------------------------------
// Classification
// ---------------------------------------------------------------------------

/// Per-double-bond classification underlying [`apply_ez_directions_from_2d_ex`].
enum EzOutcome {
    /// Not a candidate for E/Z at all: terminal alkene, carbonyl/heteroatom
    /// terminus (an end with zero non-double-bond substituents), or
    /// topologically-equivalent substituents at either end (swapping them
    /// leaves the molecule unchanged, so there is no cis/trans to encode --
    /// matches RDKit's own behavior, empirically confirmed against a live
    /// RDKit 2026.03.3 oracle on a 2-methyl-2-butene-shaped fixture: no
    /// `BondDir`, no `Bond::Stereo`, no `/`/`\` in `MolToSmiles`).
    NotRequested,
    /// Both ends resolved to a specific carrier substituent bond, tagged
    /// with the literal `BondOrder` (`Up`/`Down`) value to stash there via
    /// [`chematic_core::Molecule::set_bond_direction`].
    Assigned {
        carrier_a1: (BondIdx, BondOrder),
        carrier_a2: (BondIdx, BondOrder),
    },
    Rejected(EzDirectionRejectionReason),
}

/// One alkene end's resolution, before the whole-molecule carrier-conflict
/// check in [`apply_ez_directions_from_2d_ex`] (which can still downgrade a
/// `Resolved` outcome to a rejection).
enum EndOutcome {
    /// This end's substituents are topologically equivalent -- the whole
    /// double bond is not stereogenic (see [`EzOutcome::NotRequested`]).
    NonStereogenic,
    /// `(carrier substituent bond, "up" per the axis-referenced cross
    /// product -- true means a positive sign)`.
    Resolved(BondIdx, bool),
    Failed(EzDirectionRejectionReason),
}

fn classify_double_bond(
    mol: &Molecule,
    coords: &[(f64, f64)],
    bond_idx: BondIdx,
    explicitly_unspecified: &HashSet<BondIdx>,
    aromaticity: &AromaticityModel,
) -> EzOutcome {
    let bond = mol.bond(bond_idx);
    debug_assert_eq!(bond.order, BondOrder::Double);
    let a1 = bond.atom1;
    let a2 = bond.atom2;

    // A Kekulé aromatic-ring bond has no cis/trans isomerism at all (see
    // module docs) -- checked first, ahead of every other classification,
    // since an aromatic bond is never a candidate regardless of what its
    // topology or geometry otherwise look like.
    if aromaticity.is_bond_aromatic(bond_idx) {
        return EzOutcome::NotRequested;
    }

    // Cumulated pi system (allene/cumulene): an endpoint with another double
    // bond besides this one. Checked before the terminal-alkene check so a
    // cumulated system is reported as UnsupportedTopology, never silently
    // folded into "not stereogenic" -- its real (orthogonal-plane) geometry
    // is a fundamentally different shape 2D coordinates in one plane can't
    // represent, not merely "no substituent here".
    if has_other_double_bond(mol, a1, bond_idx) || has_other_double_bond(mol, a2, bond_idx) {
        return EzOutcome::Rejected(EzDirectionRejectionReason::UnsupportedTopology);
    }

    if explicitly_unspecified.contains(&bond_idx) {
        return EzOutcome::Rejected(EzDirectionRejectionReason::ExplicitlyUnspecified);
    }

    let subs_a1 = substituents(mol, a1, a2);
    let subs_a2 = substituents(mol, a2, a1);
    if subs_a1.is_empty() || subs_a2.is_empty() {
        return EzOutcome::NotRequested; // terminal alkene / carbonyl / heteroatom terminus
    }
    if subs_a1.len() > 2 || subs_a2.len() > 2 {
        return EzOutcome::Rejected(EzDirectionRejectionReason::UnsupportedTopology);
    }

    let p1 = match coord(coords, a1) {
        Ok(p) => p,
        Err(e) => return EzOutcome::Rejected(e),
    };
    let p2 = match coord(coords, a2) {
        Ok(p) => p,
        Err(e) => return EzOutcome::Rejected(e),
    };
    let axis = (p2.0 - p1.0, p2.1 - p1.1);
    if axis.0.hypot(axis.1) < AXIS_EPS {
        return EzOutcome::Rejected(EzDirectionRejectionReason::DegenerateGeometry); // zero-length double bond
    }

    // Both ends' cross products are referenced from the SAME origin (`p1`,
    // a1's own position) -- valid because a2 also lies on the axis line
    // through p1, so "which side of the line" is well-defined regardless of
    // which endpoint the origin sits at. Mirrors
    // `crate::stereo2d::ez_from_coords`'s existing, already-tested
    // convention exactly, so a later CIP E/Z *label* computed on top of
    // these directions (`chematic_chem::cip::assign_ez`) reproduces the
    // same Z/E verdict `assign_ez_from_2d` would compute directly from the
    // same coordinates.
    let end1 = resolve_end(mol, coords, a1, &subs_a1, p1, axis, aromaticity);
    let end2 = resolve_end(mol, coords, a2, &subs_a2, p1, axis, aromaticity);

    match (end1, end2) {
        (EndOutcome::NonStereogenic, _) | (_, EndOutcome::NonStereogenic) => {
            EzOutcome::NotRequested
        }
        (EndOutcome::Failed(reason), _) => EzOutcome::Rejected(reason),
        (_, EndOutcome::Failed(reason)) => EzOutcome::Rejected(reason),
        (EndOutcome::Resolved(b1, up1), EndOutcome::Resolved(b2, up2)) => {
            let dir1 = direction_for_up(mol.bond(b1).atom1, a1, up1);
            let dir2 = direction_for_up(mol.bond(b2).atom1, a2, up2);
            EzOutcome::Assigned {
                carrier_a1: (b1, dir1),
                carrier_a2: (b2, dir2),
            }
        }
    }
}

/// Resolve one alkene end's carrier substituent bond and geometric side.
///
/// `subs` has already been checked non-empty and at most 2 elements by the
/// caller. When there are 2 substituents and they are topologically
/// equivalent (per [`compare_branches`], used here purely as a symmetric
/// equality oracle -- never to rank/choose a "winner", so this stays
/// CIP-*independent* in the sense that matters: no substituent is ever
/// preferred by priority), the end contributes no stereogenic signal at all.
/// Otherwise, candidates are tried in adjacency order (never CIP rank,
/// matching [`crate::stereo2d_local`]'s own precedent): a candidate already
/// claimed by an existing wedge/hash bond or an unrelated pre-existing
/// `bond_direction` stash is skipped in favor of its sibling (fixture: a
/// tetrahedral wedge on the OTHER substituent of a tri-substituted alkene
/// end must not block this end's direction); a candidate lying on the
/// double-bond axis (collinear) or with a missing/non-finite coordinate is
/// likewise skipped in favor of a sibling when one exists.
fn resolve_end(
    mol: &Molecule,
    coords: &[(f64, f64)],
    end: AtomIdx,
    subs: &[(AtomIdx, BondIdx)],
    axis_origin: (f64, f64),
    axis: (f64, f64),
    aromaticity: &AromaticityModel,
) -> EndOutcome {
    if subs.len() == 2
        && compare_branches(mol, end, subs[0].0, subs[1].0) == std::cmp::Ordering::Equal
    {
        return EndOutcome::NonStereogenic;
    }

    // A 2-substituted end where EITHER candidate is itself an endpoint of a
    // DIFFERENT, genuinely-classifiable (non-aromatic) double bond is a
    // branch-point-adjacent-to-conjugation shape this module deliberately
    // does not attempt (see `is_conjugated_to_another_double_bond`'s doc
    // comment for the full root-cause writeup: found empirically on the
    // broad-corpus run, not anticipated by the original design). Choosing
    // either candidate here risks corruption downstream in
    // `chematic_smiles::canonical`'s pre-existing `resolve_ez_markers`,
    // which cannot tell "a marker scoped to a different double bond's axis"
    // apart from "no marker at all" -- and reusing the conjugated
    // candidate as a shared carrier is NOT guaranteed to agree with the
    // neighboring system's own requirement the way a straight-chain
    // conjugated diene's shared bond is (that guarantee relies on neither
    // flanking atom having an extra branch; a branch point breaks it, and a
    // real disagreeing pair confirmed this directly). Reject rather than
    // guess, exactly like the ordinary conjugated-diene shared-carrier
    // agreement check, just detected one step earlier (before ever writing
    // a value) since a branch point can't be resolved without solving the
    // Issue #149 joint-carrier problem this module is out of scope for.
    if subs.len() == 2
        && subs
            .iter()
            .any(|&(a, b)| is_conjugated_to_another_double_bond(mol, a, b, aromaticity))
    {
        return EndOutcome::Failed(EzDirectionRejectionReason::CarrierConflict);
    }

    let mut first_geometry_failure: Option<EzDirectionRejectionReason> = None;
    for &(sub_atom, sub_bond) in subs {
        if matches!(mol.bond(sub_bond).order, BondOrder::Up | BondOrder::Down) {
            continue; // reserved for tetrahedral wedge/hash notation
        }
        if mol.bond_direction(sub_bond).is_some() {
            continue; // already occupied by an unrelated pre-existing stash
        }
        match coord(coords, sub_atom) {
            Err(reason) => {
                first_geometry_failure.get_or_insert(reason);
                continue;
            }
            Ok((sx, sy)) => {
                let side = cross2d(axis.0, axis.1, sx - axis_origin.0, sy - axis_origin.1);
                if side.abs() < AXIS_EPS {
                    first_geometry_failure
                        .get_or_insert(EzDirectionRejectionReason::DegenerateGeometry);
                    continue; // substituent lies on the double-bond axis
                }
                return EndOutcome::Resolved(sub_bond, side > 0.0);
            }
        }
    }
    match first_geometry_failure {
        Some(reason) => EndOutcome::Failed(reason),
        // Every candidate was skipped for occupancy reasons (wedge/hash or
        // an unrelated stash), never for a geometry reason.
        None => EndOutcome::Failed(EzDirectionRejectionReason::CarrierConflict),
    }
}

// ---------------------------------------------------------------------------
// Small geometry/topology helpers
// ---------------------------------------------------------------------------

/// True when `atom` has a `BondOrder::Double` neighbor other than `exclude`
/// -- i.e. `atom` sits in a cumulated pi system (allene/cumulene).
/// Find every double bond that must be rejected because a branch point (a
/// 2-substituent alkene end) somewhere in the molecule has a substituent
/// that is itself an endpoint of a DIFFERENT, genuinely-classifiable
/// (non-aromatic) double bond -- BOTH the branch point's own double bond
/// AND that other double bond are poisoned, not just the former.
///
/// Rejecting only the branch point's own double bond is NOT sufficient,
/// confirmed empirically against a live RDKit oracle (not assumed): a
/// directional marker in OpenSMILES is read by plain textual adjacency, not
/// by which system produced it, so the OTHER double bond's own,
/// individually-correct marker on the shared bond unavoidably becomes a
/// (possibly wrong) reference substituent for the branch point's double
/// bond too, once re-parsed by any standards-compliant consumer -- RDKit
/// re-parsing chematic's own SMILES output for a real corpus molecule of
/// this shape was directly observed inferring a definite (and wrong)
/// stereo for the "rejected" bond purely from the neighboring bond's
/// legitimate marker. This is exactly the Issue #149 joint-carrier problem,
/// just discovered one step earlier than the ordinary shared-carrier
/// agreement check in [`apply_ez_directions_from_2d_ex`] (which still
/// handles the ordinary, non-branched conjugated-diene case correctly --
/// this pre-pass only fires when a branch point is involved).
fn poisoned_by_branch_ambiguity(
    mol: &Molecule,
    aromaticity: &AromaticityModel,
) -> HashSet<BondIdx> {
    let mut poisoned = HashSet::new();
    for (bidx, bond) in mol.bonds() {
        if bond.order != BondOrder::Double || aromaticity.is_bond_aromatic(bidx) {
            continue;
        }
        for (end, other_end) in [(bond.atom1, bond.atom2), (bond.atom2, bond.atom1)] {
            let subs = substituents(mol, end, other_end);
            if subs.len() != 2 {
                continue;
            }
            for &(sub_atom, sub_bond) in &subs {
                let other_db = mol.neighbors(sub_atom).find(|&(_, nb_bidx)| {
                    nb_bidx != sub_bond
                        && mol.bond(nb_bidx).order == BondOrder::Double
                        && !aromaticity.is_bond_aromatic(nb_bidx)
                });
                if let Some((_, other_db)) = other_db {
                    poisoned.insert(bidx);
                    poisoned.insert(other_db);
                }
            }
        }
    }
    poisoned
}

fn has_other_double_bond(mol: &Molecule, atom: AtomIdx, exclude: BondIdx) -> bool {
    mol.neighbors(atom)
        .any(|(_, bidx)| bidx != exclude && mol.bond(bidx).order == BondOrder::Double)
}

/// True when `sub_atom` (reached from an alkene end via `sub_bond`) is
/// itself an endpoint of a DIFFERENT, genuinely-classifiable (non-aromatic)
/// double bond -- i.e. `sub_atom` is part of a longer conjugated system
/// (an azine/hydrazone chain, a polyene, etc.), not a terminal/unrelated
/// substituent.
///
/// Used by [`resolve_end`] to REJECT (not choose between) a 2-substituted
/// end when either candidate has this shape -- a branch point immediately
/// adjacent to a different double bond's own conjugated system. Found
/// empirically on the broad-corpus run (not anticipated by the original
/// design), via two failed attempts, both confirmed wrong by direct
/// atom-level RDKit comparison before landing on this one:
///
/// 1. An earlier version of this module let the branch point's OTHER
///    (unrelated) substituent carry its own, independently-computed
///    direction. That value was individually correct for THIS axis, but
///    `chematic_smiles::canonical`'s pre-existing `resolve_ez_markers`
///    carrier-selection (which predates this module and cannot distinguish
///    "a marker scoped to a different double bond's axis" from "no marker
///    at all") could then silently discard it in favor of the conjugated
///    substituent's OWN marker -- which is real, but scoped to the OTHER
///    double bond's axis, not this one -- corrupting the result.
/// 2. A second attempt tried reusing the conjugated substituent's bond as
///    THIS end's own carrier too (mirroring how an ordinary conjugated
///    diene's shared middle bond legitimately serves both flanking double
///    bonds at once). Measured directly: this shared-bond-agreement
///    guarantee holds for a straight chain (verified by hand for a real
///    diene fixture) but does NOT generally hold once one of the two
///    flanking atoms is a branch point with an extra substituent -- a real
///    corpus molecule of exactly this shape produced two independently-
///    computed, genuinely DISAGREEING requirements for the same bond.
///
/// Both failure modes are avoided by rejecting outright: this is exactly
/// the Issue #149 joint-carrier-resolution problem this module is out of
/// scope for, just detected one step earlier (before ever writing a value)
/// rather than via the whole-molecule agreement check in
/// [`apply_ez_directions_from_2d_ex`], which still catches the ordinary
/// (non-branched) shared-carrier case correctly.
fn is_conjugated_to_another_double_bond(
    mol: &Molecule,
    sub_atom: AtomIdx,
    sub_bond: BondIdx,
    aromaticity: &AromaticityModel,
) -> bool {
    mol.neighbors(sub_atom).any(|(_, bidx)| {
        bidx != sub_bond
            && mol.bond(bidx).order == BondOrder::Double
            && !aromaticity.is_bond_aromatic(bidx)
    })
}

/// Non-double-bond neighbors of `end`, excluding `other_end` (the double
/// bond's other atom) -- an alkene carbon's up-to-two sigma substituents.
fn substituents(mol: &Molecule, end: AtomIdx, other_end: AtomIdx) -> Vec<(AtomIdx, BondIdx)> {
    mol.neighbors(end)
        .filter(|&(nb, bidx)| nb != other_end && mol.bond(bidx).order != BondOrder::Double)
        .collect()
}

/// A finite `(x, y)` coordinate for `idx`, or the specific reason it isn't
/// usable.
fn coord(coords: &[(f64, f64)], idx: AtomIdx) -> Result<(f64, f64), EzDirectionRejectionReason> {
    let (x, y) = coords
        .get(idx.0 as usize)
        .copied()
        .ok_or(EzDirectionRejectionReason::MissingCoordinate)?;
    if !x.is_finite() || !y.is_finite() {
        return Err(EzDirectionRejectionReason::NonFiniteCoordinate);
    }
    Ok((x, y))
}

/// 2D cross product scalar: `vx*uy - vy*ux`. Same convention as
/// `crate::stereo2d`'s private helper of the same name/shape (not shared
/// directly -- that one is private to a different module -- but
/// deliberately identical so both engines agree on what "side" means).
fn cross2d(vx: f64, vy: f64, ux: f64, uy: f64) -> f64 {
    vx * uy - vy * ux
}

/// Which literal `BondOrder` (`Up`/`Down`) to store on a bond so that
/// reading it back "from `alkene_end`'s side" (the same
/// `bond.atom1 == alkene_end` convention `chematic_chem::cip::substituent_is_up`
/// uses) yields `want_up`. Reimplemented locally (chematic-perception cannot
/// depend on chematic-smiles/chematic-chem for production code -- see the
/// crate dependency graph in the workspace `CLAUDE.md`) but verified to
/// match `chematic_chem::cip::substituent_is_up`'s convention exactly by
/// direct algebraic inversion of its match arms, not by inspection alone.
fn direction_for_up(bond_atom1: AtomIdx, alkene_end: AtomIdx, want_up: bool) -> BondOrder {
    if (bond_atom1 == alkene_end) == want_up {
        BondOrder::Up
    } else {
        BondOrder::Down
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stereo2d::cip_ez_descriptor;
    use chematic_core::{Atom, CipCode, Element, MoleculeBuilder};

    /// Read a substituent's "up"-ness back the same way
    /// `chematic_chem::cip::substituent_is_up` does (literal `Up`/`Down`
    /// bond order, or the `bond_direction` stash), so tests can verify this
    /// module's convention round-trips to the same Z/E verdict the legacy
    /// CIP engine (`cip_ez_descriptor`) computes directly from the same
    /// coordinates -- without chematic-perception depending on
    /// chematic-chem.
    fn stored_is_up(mol: &Molecule, alkene_end: AtomIdx, sub_bond: BondIdx) -> Option<bool> {
        let bond = mol.bond(sub_bond);
        let effective = mol.bond_direction(sub_bond).unwrap_or(bond.order);
        match effective {
            BondOrder::Up => Some(bond.atom1 == alkene_end),
            BondOrder::Down => Some(bond.atom1 != alkene_end),
            _ => None,
        }
    }

    /// but-2-ene skeleton (`Me0-C1=C2-Me3`) with the given 2D layout,
    /// returning `(mol, coords, double_bond_idx)`.
    fn but2ene(
        c0: (f64, f64),
        c1: (f64, f64),
        c2: (f64, f64),
        c3: (f64, f64),
    ) -> (Molecule, Vec<(f64, f64)>, BondIdx) {
        let mut b = MoleculeBuilder::new();
        let m0 = b.add_atom(Atom::new(Element::C));
        let m1 = b.add_atom(Atom::new(Element::C));
        let m2 = b.add_atom(Atom::new(Element::C));
        let m3 = b.add_atom(Atom::new(Element::C));
        b.add_bond(m0, m1, BondOrder::Single).unwrap();
        let db = b.add_bond(m1, m2, BondOrder::Double).unwrap();
        b.add_bond(m2, m3, BondOrder::Single).unwrap();
        (b.build(), vec![c0, c1, c2, c3], db)
    }

    #[test]
    fn z_but2ene_assigns_and_matches_legacy_cip_engine() {
        // Same coordinates as `stereo2d::tests::test_ez_but2ene_z`.
        let (mut mol, coords, db) = but2ene((-0.866, 0.5), (0.0, 0.0), (1.5, 0.0), (2.366, 0.5));
        let legacy = cip_ez_descriptor(&mol, db, &coords);
        assert_eq!(legacy, Some(CipCode::Z));

        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert!(diagnostics.is_empty(), "{diagnostics:?}");

        let bond = mol.bond(db);
        let (c1, c2) = (bond.atom1, bond.atom2);
        let sub1 = mol.bond_between(c1, AtomIdx(0)).unwrap().0;
        let sub2 = mol.bond_between(c2, AtomIdx(3)).unwrap().0;
        assert!(mol.bond_direction(sub1).is_some());
        assert!(mol.bond_direction(sub2).is_some());
        let up1 = stored_is_up(&mol, c1, sub1).unwrap();
        let up2 = stored_is_up(&mol, c2, sub2).unwrap();
        let reconstructed = if up1 == up2 { CipCode::Z } else { CipCode::E };
        assert_eq!(
            reconstructed,
            legacy.unwrap(),
            "direction convention must reproduce the legacy CIP engine's own Z/E verdict"
        );
        // The double bond's own order is never touched.
        assert_eq!(mol.bond(db).order, BondOrder::Double);
    }

    #[test]
    fn e_but2ene_assigns_and_matches_legacy_cip_engine() {
        let (mut mol, coords, db) = but2ene((-0.866, 0.5), (0.0, 0.0), (1.5, 0.0), (2.366, -0.5));
        let legacy = cip_ez_descriptor(&mol, db, &coords);
        assert_eq!(legacy, Some(CipCode::E));

        apply_ez_directions_from_2d(&mut mol, &coords);
        let bond = mol.bond(db);
        let (c1, c2) = (bond.atom1, bond.atom2);
        let sub1 = mol.bond_between(c1, AtomIdx(0)).unwrap().0;
        let sub2 = mol.bond_between(c2, AtomIdx(3)).unwrap().0;
        let up1 = stored_is_up(&mol, c1, sub1).unwrap();
        let up2 = stored_is_up(&mol, c2, sub2).unwrap();
        let reconstructed = if up1 == up2 { CipCode::Z } else { CipCode::E };
        assert_eq!(reconstructed, legacy.unwrap());
    }

    #[test]
    fn terminal_alkene_not_requested() {
        // CH2=CH-CH3 (propene): one end (=CH2) has zero non-double-bond
        // heavy-atom substituents.
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::new(Element::C)); // =CH2
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C)); // methyl
        b.add_bond(c0, c1, BondOrder::Double).unwrap();
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let mut mol = b.build();
        let coords = vec![(0.0, 0.0), (1.5, 0.0), (2.366, 0.5)];
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert!(diagnostics.is_empty());
        assert!(mol.bond_direction(BondIdx(1)).is_none());
    }

    #[test]
    fn carbonyl_not_requested() {
        // CH3-CH=O: the oxygen end has zero non-double-bond substituents.
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::new(Element::C));
        let c1 = b.add_atom(Atom::new(Element::C));
        let o = b.add_atom(Atom::new(Element::O));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, o, BondOrder::Double).unwrap();
        let mut mol = b.build();
        let coords = vec![(-1.0, 0.0), (0.0, 0.0), (0.5, 1.0)];
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn equivalent_substituents_not_requested() {
        // (CH3)2C=CH-CH3, 2-methyl-2-butene: the disubstituted end carries
        // two plain, topologically-identical methyl groups -- confirmed
        // against a live RDKit 2026.03.3 oracle (B0 diagnosis) to set NO
        // BondDir/Stereo and print NO `/`/`\` at all, matching this
        // assertion.
        let mut b = MoleculeBuilder::new();
        let center = b.add_atom(Atom::new(Element::C));
        let me_a = b.add_atom(Atom::new(Element::C));
        let me_b = b.add_atom(Atom::new(Element::C));
        let ch = b.add_atom(Atom::new(Element::C));
        let me_c = b.add_atom(Atom::new(Element::C));
        b.add_bond(center, me_a, BondOrder::Single).unwrap();
        b.add_bond(center, me_b, BondOrder::Single).unwrap();
        let db = b.add_bond(center, ch, BondOrder::Double).unwrap();
        b.add_bond(ch, me_c, BondOrder::Single).unwrap();
        let mut mol = b.build();
        let coords = vec![
            (0.0, 0.0),
            (-0.866, 0.5),
            (-0.866, -0.5),
            (1.5, 0.0),
            (2.366, 0.5),
        ];
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert!(diagnostics.is_empty(), "{diagnostics:?}");
        for bidx in 0..mol.bond_count() {
            assert!(mol.bond_direction(BondIdx(bidx as u32)).is_none());
        }
        let _ = db;
    }

    #[test]
    fn trisubstituted_alkene_assigns() {
        // Cl(Br)C=CH-CH3: one end has two DIFFERENT substituents (Cl, Br) --
        // must resolve to whichever is eligible, not abstain merely because
        // there are two candidates.
        let mut b = MoleculeBuilder::new();
        let center = b.add_atom(Atom::new(Element::C));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let ch = b.add_atom(Atom::new(Element::C));
        let me = b.add_atom(Atom::new(Element::C));
        b.add_bond(center, cl, BondOrder::Single).unwrap();
        b.add_bond(center, br, BondOrder::Single).unwrap();
        b.add_bond(center, ch, BondOrder::Double).unwrap();
        b.add_bond(ch, me, BondOrder::Single).unwrap();
        let mut mol = b.build();
        let coords = vec![
            (0.0, 0.0),
            (-0.866, 0.5),
            (-0.866, -0.5),
            (1.5, 0.0),
            (2.366, 0.5),
        ];
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert!(diagnostics.is_empty(), "{diagnostics:?}");
        let has_direction =
            (0..mol.bond_count()).any(|i| mol.bond_direction(BondIdx(i as u32)).is_some());
        assert!(has_direction);
    }

    #[test]
    fn missing_coordinate_rejected() {
        let (mut mol, mut coords, _db) =
            but2ene((-0.866, 0.5), (0.0, 0.0), (1.5, 0.0), (2.366, 0.5));
        coords.truncate(3); // drop the last substituent's coordinate
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert_eq!(diagnostics.len(), 1);
        assert_eq!(
            diagnostics[0].reason,
            EzDirectionRejectionReason::MissingCoordinate
        );
        for bidx in 0..mol.bond_count() {
            assert!(mol.bond_direction(BondIdx(bidx as u32)).is_none());
        }
    }

    #[test]
    fn non_finite_coordinate_rejected() {
        let (mut mol, coords, _db) =
            but2ene((-0.866, 0.5), (0.0, 0.0), (1.5, 0.0), (f64::NAN, 0.5));
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert_eq!(diagnostics.len(), 1);
        assert_eq!(
            diagnostics[0].reason,
            EzDirectionRejectionReason::NonFiniteCoordinate
        );
    }

    #[test]
    fn zero_length_double_bond_rejected() {
        // Both alkene atoms placed at the exact same point.
        let (mut mol, coords, _db) = but2ene((-0.866, 0.5), (0.0, 0.0), (0.0, 0.0), (2.366, 0.5));
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert_eq!(diagnostics.len(), 1);
        assert_eq!(
            diagnostics[0].reason,
            EzDirectionRejectionReason::DegenerateGeometry
        );
    }

    #[test]
    fn collinear_substituent_and_all_or_nothing() {
        // Mono-substituted end (Me0) placed exactly on the double-bond axis
        // (same y as both alkene carbons): geometrically ambiguous, and
        // with no sibling substituent to fall back to. The OTHER end (Me3)
        // has perfectly good geometry -- it must NOT get a direction either
        // (all-or-nothing per double bond).
        let (mut mol, coords, db) = but2ene((-1.0, 0.0), (0.0, 0.0), (1.5, 0.0), (2.366, 0.5));
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert_eq!(diagnostics.len(), 1);
        assert_eq!(
            diagnostics[0].reason,
            EzDirectionRejectionReason::DegenerateGeometry
        );
        let bond = mol.bond(db);
        let sub2 = mol.bond_between(bond.atom2, AtomIdx(3)).unwrap().0;
        assert!(
            mol.bond_direction(sub2).is_none(),
            "the good end must not get a direction when the other end fails"
        );
    }

    #[test]
    fn explicitly_unspecified_rejected() {
        let (mut mol, coords, db) = but2ene((-0.866, 0.5), (0.0, 0.0), (1.5, 0.0), (2.366, 0.5));
        let mut unspecified = HashSet::new();
        unspecified.insert(db);
        let diagnostics = apply_ez_directions_from_2d_ex(&mut mol, &coords, &unspecified);
        assert_eq!(diagnostics.len(), 1);
        assert_eq!(
            diagnostics[0].reason,
            EzDirectionRejectionReason::ExplicitlyUnspecified
        );
        for bidx in 0..mol.bond_count() {
            assert!(mol.bond_direction(BondIdx(bidx as u32)).is_none());
        }
    }

    #[test]
    fn cumulene_allene_rejected() {
        // H2C=C=CH2 (propa-1,2-diene / allene): central atom has two double
        // bonds. Both double bonds must be rejected as UnsupportedTopology,
        // never silently treated as ordinary (and definitely never
        // assigned a direction from naive 2D coordinates, which cannot
        // represent an allene's orthogonal substituent planes anyway).
        let mut b = MoleculeBuilder::new();
        let t1 = b.add_atom(Atom::new(Element::C));
        let central = b.add_atom(Atom::new(Element::C));
        let t2 = b.add_atom(Atom::new(Element::C));
        let s1 = b.add_atom(Atom::new(Element::C));
        let s2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(t1, central, BondOrder::Double).unwrap();
        b.add_bond(central, t2, BondOrder::Double).unwrap();
        b.add_bond(t1, s1, BondOrder::Single).unwrap();
        b.add_bond(t2, s2, BondOrder::Single).unwrap();
        let mut mol = b.build();
        let coords = vec![(-1.0, 0.0), (0.0, 0.0), (1.0, 0.0), (-1.5, 1.0), (1.5, 1.0)];
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert_eq!(diagnostics.len(), 2);
        assert!(
            diagnostics
                .iter()
                .all(|d| d.reason == EzDirectionRejectionReason::UnsupportedTopology)
        );
    }

    #[test]
    fn existing_wedge_on_only_candidate_is_carrier_conflict() {
        // The ONE substituent bond at an end already carries a literal
        // wedge (as if a tetrahedral-parity stage ran first and this same
        // physical bond was drawn as a wedge from the OTHER endpoint's
        // perspective) -- must never be reinterpreted as an E/Z marker, and
        // with no sibling to fall back to, the whole double bond rejects.
        let (mol, coords, db) = but2ene((-0.866, 0.5), (0.0, 0.0), (1.5, 0.0), (2.366, 0.5));
        let bond = mol.bond(db);
        let sub1 = mol.bond_between(bond.atom1, AtomIdx(0)).unwrap().0;
        // There is no in-place bond-order setter, so rebuild through the
        // builder to force sub1's order to a literal wedge (`Up`), keeping
        // every other atom/bond identical.
        let mut b = MoleculeBuilder::new();
        for (_, atom) in mol.atoms() {
            b.add_atom(atom.clone());
        }
        for (bidx, bond) in mol.bonds() {
            let order = if bidx == sub1 {
                BondOrder::Up
            } else {
                bond.order
            };
            b.add_bond(bond.atom1, bond.atom2, order).unwrap();
        }
        let mut mol = b.build();

        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert_eq!(diagnostics.len(), 1);
        assert_eq!(
            diagnostics[0].reason,
            EzDirectionRejectionReason::CarrierConflict
        );
        // The wedge itself must survive untouched.
        assert_eq!(mol.bond(sub1).order, BondOrder::Up);
        assert!(mol.bond_direction(sub1).is_none());
    }

    #[test]
    fn existing_wedge_on_sibling_falls_back_to_other_substituent() {
        // Trisubstituted end (Cl, Br); Cl's bond is already a literal wedge
        // (unrelated tetrahedral notation elsewhere) -- the E/Z stage must
        // fall back to Br rather than abstaining, and must never touch the
        // Cl bond's order or write anything into its direction stash.
        let mut b = MoleculeBuilder::new();
        let center = b.add_atom(Atom::new(Element::C));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let ch = b.add_atom(Atom::new(Element::C));
        let me = b.add_atom(Atom::new(Element::C));
        let cl_bond = b.add_bond(center, cl, BondOrder::Up).unwrap(); // pre-existing wedge
        let br_bond = b.add_bond(center, br, BondOrder::Single).unwrap();
        b.add_bond(center, ch, BondOrder::Double).unwrap();
        b.add_bond(ch, me, BondOrder::Single).unwrap();
        let mut mol = b.build();
        let coords = vec![
            (0.0, 0.0),
            (-0.866, 0.5),
            (-0.866, -0.5),
            (1.5, 0.0),
            (2.366, 0.5),
        ];
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert!(diagnostics.is_empty(), "{diagnostics:?}");
        assert_eq!(mol.bond(cl_bond).order, BondOrder::Up, "wedge untouched");
        assert!(
            mol.bond_direction(cl_bond).is_none(),
            "must not stash a direction on the wedge bond"
        );
        assert!(
            mol.bond_direction(br_bond).is_some(),
            "must fall back to the sibling substituent"
        );
    }

    #[test]
    fn conjugated_diene_shared_bond_agrees() {
        // (2E,4E)-hexa-2,4-diene, standard all-anti zigzag layout: the
        // middle single bond (Cb-Cc) is simultaneously the carrier
        // candidate for BOTH double bonds. Hand-verified: both bonds'
        // independent geometric computations require `Down` on the shared
        // bond -- they must agree, not conflict, and both must be Assigned.
        let mut b = MoleculeBuilder::new();
        let me1 = b.add_atom(Atom::new(Element::C));
        let ca = b.add_atom(Atom::new(Element::C));
        let cb = b.add_atom(Atom::new(Element::C));
        let cc = b.add_atom(Atom::new(Element::C));
        let cd = b.add_atom(Atom::new(Element::C));
        let me2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(me1, ca, BondOrder::Single).unwrap();
        let db1 = b.add_bond(ca, cb, BondOrder::Double).unwrap();
        let shared = b.add_bond(cb, cc, BondOrder::Single).unwrap();
        let db2 = b.add_bond(cc, cd, BondOrder::Double).unwrap();
        b.add_bond(cd, me2, BondOrder::Single).unwrap();
        let mut mol = b.build();
        let coords = vec![
            (-2.0, 0.5), // Me1
            (-1.0, 0.0), // Ca
            (0.0, 0.5),  // Cb
            (1.0, 0.0),  // Cc
            (2.0, 0.5),  // Cd
            (3.0, 0.0),  // Me2
        ];
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert!(diagnostics.is_empty(), "{diagnostics:?}");
        assert_eq!(mol.bond_direction(shared), Some(BondOrder::Down));
        let _ = (db1, db2);
    }

    #[test]
    fn conjugated_diene_shared_bond_conflict() {
        // Deliberately non-planar-consistent layout (hand-verified): bond1
        // requires `Up` on the shared Cb-Cc bond while bond2 independently
        // requires `Down` on the SAME physical bond. Per module scope, both
        // double bonds relying on it must reject with CarrierConflict --
        // neither wins by bond-index/parse order, and nothing is written.
        let mut b = MoleculeBuilder::new();
        let me1 = b.add_atom(Atom::new(Element::C));
        let ca = b.add_atom(Atom::new(Element::C));
        let cb = b.add_atom(Atom::new(Element::C));
        let cc = b.add_atom(Atom::new(Element::C));
        let cd = b.add_atom(Atom::new(Element::C));
        let me2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(me1, ca, BondOrder::Single).unwrap();
        b.add_bond(ca, cb, BondOrder::Double).unwrap();
        let shared = b.add_bond(cb, cc, BondOrder::Single).unwrap();
        b.add_bond(cc, cd, BondOrder::Double).unwrap();
        b.add_bond(cd, me2, BondOrder::Single).unwrap();
        let mut mol = b.build();
        let coords = vec![
            (-1.0, 1.0), // Me1
            (0.0, 0.0),  // Ca
            (1.0, 0.0),  // Cb
            (2.0, 1.0),  // Cc
            (2.0, 2.0),  // Cd
            (3.0, 3.0),  // Me2
        ];
        let diagnostics = apply_ez_directions_from_2d_with_diagnostics(&mut mol, &coords);
        assert_eq!(diagnostics.len(), 2, "{diagnostics:?}");
        assert!(
            diagnostics
                .iter()
                .all(|d| d.reason == EzDirectionRejectionReason::CarrierConflict)
        );
        assert!(
            mol.bond_direction(shared).is_none(),
            "a conflicting shared carrier must end up with NO direction written"
        );
    }
}
