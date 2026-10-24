//! Typed indexed-graph relation comparison and diagnostics.

use super::*;

// ─── Typed indexed-graph-relation diagnostics ──────────────────────────────

/// Typed reason [`compare_indexed_graph_relation`] could not reach a
/// conclusive [`GraphRelation`] -- never a free-form string, so callers can
/// match on it instead of parsing text. Scoped to exactly what this API
/// produces (unlike an earlier draft, which defined a single, larger
/// `IdentityDiagnostic` enum shared with a sibling, still-unmerged PR's
/// accurate-CIP dedup preflight (issue #161, `compare_molecules_with_accurate_cip_preflight`
/// et al. -- see #178). That sibling PR's diagnostics (CIP-tie/budget/engine-error
/// reasons) don't apply to this API at all, so this version keeps its own,
/// smaller enum rather than forward-declaring unused variants. If both PRs
/// land, reconciling the two enums (or keeping them separate, since they
/// serve genuinely different purposes -- CIP-ranking tie-breaking vs. graph
/// structural comparison) is a merge-order decision for a human, not
/// resolved here since it would require depending on unmerged code.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IndexedGraphRelationDiagnostic {
    /// Native InChI is unavailable: the `native-inchi` feature is disabled,
    /// or the C library failed for a reason unrelated to CIP ranking
    /// (kekulization limit, library error).
    NativeInchiUnavailable,
    /// The molecule has no valid InChI at all (e.g. zero heavy atoms).
    InvalidMolecule,
    /// [`compare_indexed_graph_relation`]'s bounded structural comparator
    /// only establishes atom correspondence by matching raw input index (see
    /// that function's module docs) -- the two molecules' bond structure did
    /// not line up completely under that direct index correspondence, so no
    /// conclusive relation was reached. A real, literal graph-isomorphism
    /// search (atom orders free to differ) is out of scope for this API's
    /// current version -- see the module docs for the full reasoning.
    AtomOrderCorrespondenceNotEstablished,
    /// A specified tetrahedral stereocentre's raw SMILES/MOL neighbor order
    /// (`Molecule::stereo_neighbor_order`) was missing, or not exactly 4
    /// entries, on at least one side -- its parity could not be expressed in
    /// the canonical (index-sorted) reference frame [`compare_indexed_graph_relation`]
    /// compares in, so that specific centre's contribution to the relation is
    /// unverified.
    TetrahedralParityNotRepresented,
    /// [`GraphStrictness::ChemicalGraphExact`] could not Kekulize one (or
    /// both) of the input molecules' aromatic bonds (the same limitation
    /// `crate::native::standard_inchi` already surfaces as
    /// `InchiError::KekulizationFailed` for some fused/bridgehead-heteroatom
    /// ring systems) -- comparing an un-Kekulizable aromatic structure
    /// literally against the other input would risk a false `Distinct` (the
    /// aromatic-order bonds simply would not match a Kekulized structure's
    /// Single/Double spelling, for a reason that has nothing to do with
    /// whether the molecules are actually the same), so the comparison
    /// fails closed to "inconclusive" instead.
    ChemicalGraphKekulizationFailed,
}

// ─── Indexed graph relation (IndexedGraphRelationMode / IndexedGraphRelationEvidence) ──
//
// A separate, standalone identity API -- not a fifth `IdentityPolicy`
// variant (that enum, and `identity_verify`/`verified_identity_key`, are
// matched exhaustively by callers with no wildcard arm; adding a variant
// would be a breaking change -- see `ALL_POLICIES` in
// `tests/dedup_fixtures.rs`). This is the follow-up the module docs above
// flagged as out of scope for the original four InChI-normalization-based
// policies: no tautomer/mobile-H normalization at all, literal
// element/isotope/charge/aromaticity/bond-order/tetrahedral-parity/E-Z/
// fragment-composition equality.
//
// # Naming: "indexed graph relation", deliberately not "graph identity"
//
// An earlier draft of this API called itself `compare_graph_identity` /
// `GraphIdentityMode` and described its result as "exact graph identity."
// That overclaimed what the comparator actually does. It is **not** a
// general graph-isomorphism test: [`compare_indexed_graph_relation`] only
// recognizes two molecules as the same graph when they already share the
// same atom indexing (`Molecule` index `i` in one input corresponds to
// index `i` in the other) -- see "Scope boundary: atom ordering" below. Two
// molecules that represent the identical graph but happen to have their
// atoms numbered/ordered differently are **not** detected as identical here
// -- the comparison reports `None` (inconclusive), never a guess, in that
// case. The name `compare_indexed_graph_relation` (and
// `IndexedGraphRelationMode`/`IndexedGraphRelationEvidence`) states that
// precondition directly instead of implying a capability (full isomorphism
// search) this version does not have. This is a real, motivated use case on
// its own -- e.g. conformer/ensemble grouping over records that already
// share one connection table -- just not "graph identity" in the general
// sense.
//
// **Canonical SMILES is never used as the proof here** (same standing rule
// as the rest of this module): [`compare_indexed_graph_relation`] establishes
// atom correspondence via direct input-index matching plus a bond-set check,
// and verifies every attribute in-place -- never via a canonical string.
//
// # Scope boundary: atom ordering
//
// This first version's structural comparator only recognizes two molecules
// as the same graph when they share **the same atom indexing** (see above)
// -- it is not a general-purpose graph-isomorphism search over arbitrary
// atom permutations. This is not a corner deliberately cut for convenience:
// a real isomorphism search needs either an isomorphism-library dependency
// (`chematic-smarts`'s VF2 implementation is not on `chematic-inchi`'s
// dependency graph today) or a from-scratch bounded backtracking matcher
// with its own, non-trivial correctness-testing burden -- both a larger
// undertaking than fits this change, and risky to get subtly wrong under
// time pressure. Renaming the API to state its real precondition (this
// section) was judged the more honest and feasible choice for this version;
// implementing real isomorphism search is flagged as a legitimate, separate
// follow-up, not attempted here. The motivating real-world consumer named in
// this program (conformer/ensemble grouping over records that share one
// connection table) already satisfies same-indexing by construction, so this
// scope boundary does not block that use case. When atom ordering doesn't
// line up, [`IndexedGraphRelationEvidence::graph_relation`] is `None` (never
// a guess) with [`IndexedGraphRelationDiagnostic::AtomOrderCorrespondenceNotEstablished`],
// and the native-InChI corroboration fields are still populated when
// possible so the caller has *some* signal. A full isomorphism search is a
// natural follow-up, not attempted here.
//
// # Scope boundary: enhanced stereo
//
// The task brief for this API asked for "enhanced stereo (if representable)"
// to be part of what's compared. `chematic-core`/`chematic-perception` have
// no enhanced-stereochemistry representation today (no `StereoGroup`
// AND/OR/ABS-style relative/absolute grouping over multiple centres beyond
// plain per-atom `Chirality` -- confirmed by grep across both crates before
// writing this module). Building a new enhanced-stereo model as a side
// effect of this task would be exactly the kind of unrequested, speculative
// abstraction this codebase's own conventions warn against; this is
// documented here as a real scope boundary, not silently ignored.
//
// # Scope boundary: `GraphStrictness::ChemicalGraphExact` and resonance forms
//
// See that variant's own doc comment: Kekulization normalizes "aromatic-
// flagged" vs. "one particular already-Kekulized" spelling of the same ring,
// but does not recognize two independently-written, already-Kekulized
// alternative resonance structures of the same aromatic system as
// equivalent (that needs Hückel aromaticity perception, which is not on this
// crate's dependency graph either).

/// How strictly [`compare_indexed_graph_relation`] treats bond-order/aromaticity
/// spelling. Orthogonal to [`AtomMapPolicy`] -- see [`IndexedGraphRelationMode`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphStrictness {
    /// Literal bond-order and `aromatic`-flag equality. Two molecules that
    /// differ only in Kekulé vs. aromatic-order spelling for the same ring
    /// system (e.g. one writer emits `Aromatic`-order bonds, another emits
    /// the same ring pre-Kekulized as alternating Single/Double) are **not**
    /// identical here -- the strictest, most literal "same graph as written."
    RawGraphExact,
    /// Chemically-meaningful equality: any input molecule with `Aromatic`-
    /// order bonds is Kekulized first, via the same deterministic algorithm
    /// (`chematic_core::kekulize`/`apply_kekule`) this crate's native InChI
    /// conversion already relies on -- **not** canonical-SMILES string
    /// equality (see this module's standing rule). The `aromatic` flag
    /// itself is then excluded from atom comparison: Kekulization does not
    /// clear it (see `apply_kekule`'s own doc comment), so comparing it
    /// post-Kekulization would silently reintroduce the exact spelling-
    /// sensitivity this mode exists to remove.
    ///
    /// Known limitation (see the module docs above): two independently-
    /// written, already-Kekulized structures representing *different* valid
    /// resonance forms of the same aromatic system (e.g. two of
    /// naphthalene's several valid Kekulé structures) are not recognized as
    /// equivalent here -- that needs Hückel aromaticity perception, not on
    /// this crate's dependency graph. Never a false
    /// [`GraphRelation::Identical`] either way: worst case this degrades to
    /// an honest mismatch or, if Kekulization itself fails, an inconclusive
    /// result (see [`IndexedGraphRelationDiagnostic::ChemicalGraphKekulizationFailed`]).
    ChemicalGraphExact,
}

/// Whether [`compare_indexed_graph_relation`] treats reaction atom-mapping numbers
/// (`Atom::atom_map`, the `:1` in `[CH4:1]`) as part of the molecule's
/// identity. Orthogonal to [`GraphStrictness`] -- any combination of the two
/// axes is valid, see [`IndexedGraphRelationMode`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomMapPolicy {
    /// Atom-mapping numbers must correspond exactly (`None` only matches
    /// `None`; `Some(n)` only matches the literal same `n`) at every atom
    /// under the discovered correspondence.
    Include,
    /// Atom-mapping numbers are not compared at all -- appropriate for
    /// almost every caller outside reaction-atom-mapping bookkeeping (e.g.
    /// conformer/ensemble grouping, MOL/SDF round-trip checks, ordinary
    /// dedup: none of these contexts carry or care about reaction atom
    /// maps).
    Ignore,
}

/// Mode for [`compare_indexed_graph_relation`]: two independent, orthogonal axes
/// (see [`GraphStrictness`] / [`AtomMapPolicy`]) rather than one flat,
/// 3-or-4-variant enum -- a caller can freely combine either strictness with
/// either atom-map policy (e.g. `ChemicalGraphExact` + ignore atom maps, or
/// `RawGraphExact` + include atom maps) instead of being limited to whatever
/// preset combinations a flat enum happened to name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexedGraphRelationMode {
    pub graph_strictness: GraphStrictness,
    pub atom_maps: AtomMapPolicy,
}

/// A conclusive indexed-graph-relation result between two molecules (see
/// [`compare_indexed_graph_relation`]). Two values only -- both mean "the search
/// examined every dimension `graph_strictness`/`atom_maps` requires and
/// reached a decisive answer." When the search cannot reach one (atom
/// ordering doesn't line up, a stereocentre's parity can't be represented,
/// Kekulization failed, ...), [`IndexedGraphRelationEvidence::graph_relation`] is `None`
/// instead -- never a third "maybe" value smuggled into this enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphRelation {
    Identical,
    Distinct,
}

/// Evidence produced by [`compare_indexed_graph_relation`] for one pair of
/// molecules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexedGraphRelationEvidence {
    /// The conclusive structural relation, if the bounded comparator reached
    /// one. `None` means "inconclusive" -- never "distinct" or "identical"
    /// by omission. Always check `diagnostics` for why.
    pub graph_relation: Option<GraphRelation>,
    /// The shared native standard InChI string for both molecules, IF native
    /// InChI generation succeeded independently for both AND their full
    /// (stereo-inclusive) strings agree. Supplementary corroboration only:
    /// per this module's standing rule, `graph_relation` is never derived
    /// from this field alone -- native InChI performs its own IUPAC-defined
    /// tautomer/mobile-H normalization, a different (weaker) question than
    /// literal graph identity (see the module-level docs' `IdentityPolicy`
    /// discussion of the same distinction).
    pub standard_inchi: Option<String>,
    /// The shared standard InChIKey, under the same "both succeeded and
    /// agreed" condition as `standard_inchi`.
    pub inchikey: Option<String>,
    /// Typed reasons for anything short of a fully conclusive, fully
    /// corroborated result.
    pub diagnostics: Vec<IndexedGraphRelationDiagnostic>,
}

/// Structural bond-order class used by [`compare_indexed_graph_relation`]'s
/// comparisons: `Up`/`Down`/`Dative` normalize to the same class as `Single`
/// (mirroring `crate::native::convert::mol_to_inchi_atoms`'s own
/// `bond_type` mapping) since their directional/donor-acceptor meaning is
/// verified separately (E/Z via [`directed_bond_marker`], not via this raw
/// class); everything else keeps its own distinct class.
fn structural_bond_order(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 5,
        BondOrder::Zero => 6,
        BondOrder::QueryAny => 7,
        BondOrder::QuerySingleOrDouble => 8,
        BondOrder::QuerySingleOrAromatic => 9,
        BondOrder::QueryDoubleOrAromatic => 10,
    }
}

/// Cheap, atom-order-independent composition signature: fragment sizes,
/// atom multiset (element/isotope/charge/aromatic-per-`strictness`), bond
/// multiset (structural order + endpoint elements). A mismatch here is
/// conclusive (`GraphRelation::Distinct`) regardless of atom ordering --
/// this is a real graph invariant, not canonical-SMILES string equality.
#[derive(PartialEq, Eq)]
struct CompositionSignature {
    fragment_sizes: Vec<usize>,
    atom_multiset: Vec<(u8, Option<u16>, i8, bool)>,
    bond_multiset: Vec<(u8, u8, u8)>,
}

fn composition_signature(mol: &Molecule, strictness: GraphStrictness) -> CompositionSignature {
    let mut fragment_sizes: Vec<usize> = mol.fragments().iter().map(|f| f.atom_count()).collect();
    fragment_sizes.sort_unstable();

    let mut atom_multiset: Vec<(u8, Option<u16>, i8, bool)> = mol
        .atoms()
        .map(|(_, a)| {
            let aromatic = match strictness {
                GraphStrictness::RawGraphExact => a.aromatic,
                GraphStrictness::ChemicalGraphExact => false,
            };
            (a.element.atomic_number(), a.isotope, a.charge, aromatic)
        })
        .collect();
    atom_multiset.sort_unstable();

    let mut bond_multiset: Vec<(u8, u8, u8)> = mol
        .bonds()
        .map(|(_, b)| {
            let e1 = mol.atom(b.atom1).element.atomic_number();
            let e2 = mol.atom(b.atom2).element.atomic_number();
            (structural_bond_order(b.order), e1.min(e2), e1.max(e2))
        })
        .collect();
    bond_multiset.sort_unstable();

    CompositionSignature {
        fragment_sizes,
        atom_multiset,
        bond_multiset,
    }
}

/// Attempt to Kekulize `mol` for [`GraphStrictness::ChemicalGraphExact`].
/// `Ok(None)` = no aromatic bonds present (nothing to do, compare `mol`
/// as-is); `Ok(Some(k))` = Kekulized working copy; `Err(())` = Kekulization
/// failed (caller must not fall back to comparing the un-Kekulized aromatic
/// form literally against the other input -- see
/// [`IndexedGraphRelationDiagnostic::ChemicalGraphKekulizationFailed`]'s doc comment for
/// why that would risk a false `Distinct`).
fn kekulize_for_comparison(mol: &Molecule) -> Result<Option<Molecule>, ()> {
    if !mol.bonds().any(|(_, b)| b.order == BondOrder::Aromatic) {
        return Ok(None);
    }
    match kekulize(mol) {
        Ok(k) => Ok(Some(apply_kekule(mol, &k))),
        Err(_) => Err(()),
    }
}

/// Permutation parity (`true` = even) needed to sort `list` into ascending
/// order, via inversion counting. `list` always has exactly 4, pairwise-
/// distinct entries here (`STEREO_H_SENTINEL` counts as its own, always-
/// largest value) -- see [`canonical_tetrahedral_parity`].
fn is_even_permutation(list: &[u32]) -> bool {
    let mut inversions = 0usize;
    for i in 0..list.len() {
        for j in (i + 1)..list.len() {
            if list[i] > list[j] {
                inversions += 1;
            }
        }
    }
    inversions.is_multiple_of(2)
}

/// Canonical signed tetrahedral parity for atom `idx`, expressed in a common
/// reference frame (ascending raw `AtomIdx` order) so it is directly
/// comparable across two molecules that share the same atom indexing --
/// pure combinatorial permutation parity, no CIP engine involved (a
/// deliberately different, lower-level tool than any CIP-ranking-based
/// comparison: this comparator needs plain combinatorial parity, not
/// chemical priority ranking). `None` if the atom isn't a well-formed
/// tetrahedral centre in chematic's own representation here (no chirality
/// specified, no recorded neighbor order, or not exactly 4 neighbors).
fn canonical_tetrahedral_parity(mol: &Molecule, idx: AtomIdx) -> Option<bool> {
    let atom = mol.atom(idx);
    let chirality_is_even = match atom.chirality {
        Chirality::Clockwise => true,
        Chirality::CounterClockwise => false,
        _ => return None,
    };
    let order = mol.stereo_neighbor_order(idx)?;
    if order.len() != 4 {
        return None;
    }
    let sorted_is_even = is_even_permutation(order);
    Some(chirality_is_even == sorted_is_even)
}

/// Directed wedge/hash marker for the bond from `from` to `to`: `Some(true)`
/// = "up" traversed in that direction, `Some(false)` = "down", `None` = no
/// `Up`/`Down` direction stored for this bond. Mirrors
/// `crate::native::convert::mol_to_inchi_atoms`'s own `is_up` convention
/// exactly (same `atom1`-relative sign rule) rather than inventing a new
/// one.
fn directed_bond_marker(mol: &Molecule, from: AtomIdx, to: AtomIdx) -> Option<bool> {
    let (bidx, bond) = mol.bond_between(from, to)?;
    let effective = mol.bond_direction(bidx).unwrap_or(bond.order);
    match effective {
        BondOrder::Up => Some(bond.atom1 == from),
        BondOrder::Down => Some(bond.atom1 == to),
        _ => None,
    }
}

/// Outcome of [`same_structure_under_identity`].
enum StructuralComparison {
    Identical,
    Distinct,
    Inconclusive(IndexedGraphRelationDiagnostic),
}

/// The identity-atom-index-correspondence structural comparator (see this
/// module's "scope boundary: atom ordering" docs above). `a`/`b` are assumed
/// to already be [`GraphStrictness`]-normalized (Kekulized if
/// `ChemicalGraphExact`) by the caller.
fn same_structure_under_identity(
    a: &Molecule,
    b: &Molecule,
    mode: IndexedGraphRelationMode,
) -> StructuralComparison {
    if a.atom_count() != b.atom_count() {
        return StructuralComparison::Distinct;
    }
    let n = a.atom_count();

    let bond_key = |bd: &chematic_core::BondEntry| -> (u32, u32, u8) {
        let (x, y) = (bd.atom1.0, bd.atom2.0);
        (x.min(y), x.max(y), structural_bond_order(bd.order))
    };
    let bonds_a: std::collections::HashSet<(u32, u32, u8)> =
        a.bonds().map(|(_, bd)| bond_key(bd)).collect();
    let bonds_b: std::collections::HashSet<(u32, u32, u8)> =
        b.bonds().map(|(_, bd)| bond_key(bd)).collect();
    if bonds_a != bonds_b {
        return StructuralComparison::Inconclusive(
            IndexedGraphRelationDiagnostic::AtomOrderCorrespondenceNotEstablished,
        );
    }

    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let (aa, ab) = (a.atom(idx), b.atom(idx));
        let (aromatic_a, aromatic_b) = match mode.graph_strictness {
            GraphStrictness::RawGraphExact => (aa.aromatic, ab.aromatic),
            GraphStrictness::ChemicalGraphExact => (false, false),
        };
        if aa.element != ab.element
            || aa.isotope != ab.isotope
            || aa.charge != ab.charge
            || aromatic_a != aromatic_b
        {
            return StructuralComparison::Distinct;
        }
        if mode.atom_maps == AtomMapPolicy::Include && aa.atom_map != ab.atom_map {
            return StructuralComparison::Distinct;
        }
    }

    let mut inconclusive = false;
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        match (
            canonical_tetrahedral_parity(a, idx),
            canonical_tetrahedral_parity(b, idx),
        ) {
            (None, None) => {}
            (Some(pa), Some(pb)) => {
                if pa != pb {
                    return StructuralComparison::Distinct;
                }
            }
            _ => inconclusive = true,
        }
    }

    for (_, bd) in a.bonds() {
        if bd.order != BondOrder::Double {
            continue;
        }
        for from in [bd.atom1, bd.atom2] {
            let to_end = if from == bd.atom1 { bd.atom2 } else { bd.atom1 };
            for (nb, _) in a.neighbors(from) {
                if nb == to_end {
                    continue;
                }
                let da = directed_bond_marker(a, from, nb);
                let db = directed_bond_marker(b, from, nb);
                if da != db {
                    return StructuralComparison::Distinct;
                }
            }
        }
    }

    if inconclusive {
        StructuralComparison::Inconclusive(
            IndexedGraphRelationDiagnostic::TetrahedralParityNotRepresented,
        )
    } else {
        StructuralComparison::Identical
    }
}

/// Independent native-InChI corroboration for [`IndexedGraphRelationEvidence`] (see that
/// struct's docs): each molecule's own full-stereo standard InChI
/// ([`identity_verify`] under [`IdentityPolicy::StandardInchiString`],
/// unchanged/legacy-only -- this corroboration signal deliberately does not
/// pull in the accurate-CIP preflight, keeping the two features decoupled),
/// populated only when both succeed AND agree.
fn shared_native_inchi_evidence(
    mol_a: &Molecule,
    mol_b: &Molecule,
    diagnostics: &mut Vec<IndexedGraphRelationDiagnostic>,
) -> (Option<String>, Option<String>) {
    let (ra, rb) = match (
        identity_verify(mol_a, IdentityPolicy::StandardInchiString),
        identity_verify(mol_b, IdentityPolicy::StandardInchiString),
    ) {
        (Verify::Ok(ra), Verify::Ok(rb)) => (ra, rb),
        (Verify::Invalid, _) | (_, Verify::Invalid) => {
            diagnostics.push(IndexedGraphRelationDiagnostic::InvalidMolecule);
            return (None, None);
        }
        _ => {
            diagnostics.push(IndexedGraphRelationDiagnostic::NativeInchiUnavailable);
            return (None, None);
        }
    };
    if ra != rb {
        return (None, None);
    }
    let key = inchi_key(&ra).ok();
    (Some(ra), key)
}

/// Compare two molecules for indexed graph relation under `mode` --
/// **requires the two molecules to already share the same atom indexing**
/// (see the module docs above for the full rationale, the naming decision,
/// scope boundaries, and why canonical SMILES is never the proof).
pub fn compare_indexed_graph_relation(
    mol_a: &Molecule,
    mol_b: &Molecule,
    mode: IndexedGraphRelationMode,
) -> IndexedGraphRelationEvidence {
    let mut diagnostics: Vec<IndexedGraphRelationDiagnostic> = Vec::new();

    let (working_a, working_b) = match mode.graph_strictness {
        GraphStrictness::RawGraphExact => (None, None),
        GraphStrictness::ChemicalGraphExact => {
            match (
                kekulize_for_comparison(mol_a),
                kekulize_for_comparison(mol_b),
            ) {
                (Ok(a), Ok(b)) => (a, b),
                _ => {
                    diagnostics
                        .push(IndexedGraphRelationDiagnostic::ChemicalGraphKekulizationFailed);
                    let (inchi, key) = shared_native_inchi_evidence(mol_a, mol_b, &mut diagnostics);
                    return IndexedGraphRelationEvidence {
                        graph_relation: None,
                        standard_inchi: inchi,
                        inchikey: key,
                        diagnostics,
                    };
                }
            }
        }
    };
    let ref_a = working_a.as_ref().unwrap_or(mol_a);
    let ref_b = working_b.as_ref().unwrap_or(mol_b);

    if composition_signature(ref_a, mode.graph_strictness)
        != composition_signature(ref_b, mode.graph_strictness)
    {
        return IndexedGraphRelationEvidence {
            graph_relation: Some(GraphRelation::Distinct),
            standard_inchi: None,
            inchikey: None,
            diagnostics,
        };
    }

    match same_structure_under_identity(ref_a, ref_b, mode) {
        StructuralComparison::Distinct => IndexedGraphRelationEvidence {
            graph_relation: Some(GraphRelation::Distinct),
            standard_inchi: None,
            inchikey: None,
            diagnostics,
        },
        StructuralComparison::Inconclusive(reason) => {
            diagnostics.push(reason);
            let (inchi, key) = shared_native_inchi_evidence(mol_a, mol_b, &mut diagnostics);
            IndexedGraphRelationEvidence {
                graph_relation: None,
                standard_inchi: inchi,
                inchikey: key,
                diagnostics,
            }
        }
        StructuralComparison::Identical => {
            let (inchi, key) = shared_native_inchi_evidence(mol_a, mol_b, &mut diagnostics);
            IndexedGraphRelationEvidence {
                graph_relation: Some(GraphRelation::Identical),
                standard_inchi: inchi,
                inchikey: key,
                diagnostics,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    include!("../dedup_boundary_tests.rs");
}
