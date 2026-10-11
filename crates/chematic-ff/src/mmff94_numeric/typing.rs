// ── Atom type assignment ─────────────────────────────────────────────────────

/// Assign MMFF94 numeric atom types (1–99) to all atoms in the molecule.
///
/// This implements the core atom type perception rules for organic chemistry.
/// For atoms not handled, returns `Err`.
///
/// Thin wrapper over [`assign_mmff94_numeric_types_with_view`] that discards
/// the re-perceived molecule -- correct for callers that only need numeric
/// types (charges, vdW, atom-level reporting). Any caller that ALSO does
/// bond-order-dependent classification (`bond_type_for`/`angle_type_for`/
/// `torsion_type_for`/`stretch_bend_type_for`, all of which read
/// [`chematic_core::BondOrder`] directly) must call
/// [`assign_mmff94_numeric_types_with_view`] instead and use its returned
/// molecule for that purpose, not its own input `mol` -- see that function's
/// doc for why (issue #227 Phase 1, torsion parameter gap root cause).
pub fn assign_mmff94_numeric_types(mol: &Molecule) -> Result<Vec<u8>, NumericTypeError> {
    assign_mmff94_numeric_types_with_view(mol).map(|(types, _)| types)
}

/// Like [`assign_mmff94_numeric_types`], but also returns the MMFF-specific
/// re-perceived molecule ([`compute_mmff94_aromatic_view`]'s output, same
/// atom count/bond topology as `mol`, only [`chematic_core::BondOrder`]
/// values on ring bonds can differ) used to derive those types.
///
/// Root cause this exists to fix (issue #227 Phase 1, torsion parameter gap,
/// 2026-08-15): `bond_type_for`'s "`BondOrder::Aromatic` forces
/// `bond_type=0`" rule is itself correct (confirmed against a live RDKit
/// oracle for benzene's own ring bond), but every one of chematic-ff's
/// bond-order-dependent classification call sites
/// (`bond_type_for`/`angle_type_for`/`torsion_type_for`/
/// `stretch_bend_type_for`, both in production energy/gradient code and in
/// coverage-gate diagnostics) was being fed the CALLER'S original `mol`, not
/// this re-perceived view -- even though the numeric TYPES those same call
/// sites use were already correctly derived from it. For a ring system where
/// chematic's general/SMILES aromaticity perception and MMFF94's own
/// stricter, Kekule-based perception (`setMMFFAromaticity`) disagree (e.g.
/// caffeine's pyrimidinedione ring: chematic's general model treats the
/// whole fused bicyclic system as one delocalized aromatic ring, RDKit's
/// real sanitizer Kekulizes that ring to alternating single/double bonds
/// while leaving only the fused imidazole ring aromatic -- oracle-confirmed
/// via `MolFromSmiles(...).GetBondBetweenAtoms(5,6).GetIsAromatic() ==
/// False`, and independently confirmed against chematic's own
/// pre-existing, already-oracle-validated
/// `validation/results/mmff94_aromaticity_bond_parity_227_oracle.json` dump,
/// which already recorded `bond_aromatic["5-6"]: false` for caffeine before
/// this fix), this meant the classification code was computed from the
/// WRONG bond order, landing on a torsion/bond/angle-type code with no table
/// row even though chematic's own, unmodified parameter table already
/// carries the correct row at the code RDKit's real (Kekulized) bond type
/// resolves to. Measured on the 265-molecule Wave 1 corpus: 254/254 of the
/// `torsions_missing` instances with `present_at_different_classification =
/// Some` (`crates/chematic-3d/examples/mmff94_term_coverage_audit.rs`) have
/// this exact shape -- oracle-validated (all 254, not a sample) via a live
/// `GetMMFFTorsionParams` cross-check: the oracle's returned value matches
/// EXACTLY one of chematic's own pre-existing rows at a different
/// classification code in 254/254 cases, and RDKit's real bond object is
/// non-aromatic (`GetIsAromatic() == False`) in 254/254 cases. No Halgren
/// empirical-rule implementation was needed for any of them (see
/// `scripts/mmff94_provenance/PROVENANCE.md`'s Torsion entry for the two
/// falsified alternative hypotheses this ruled out first).
pub fn assign_mmff94_numeric_types_with_view(
    mol: &Molecule,
) -> Result<(Vec<u8>, Molecule), NumericTypeError> {
    let n = mol.atom_count();
    let mut types = vec![0u8; n];
    // MMFF94's aromaticity loop follows RDKit's ordered RingInfo semantics,
    // including the large-cycle representatives that control when charged
    // fused rings become resolved. `compute_mmff94_aromatic_view` restores
    // the compatibility-specific processing order without changing the
    // general perception API's canonical ring ordering.
    // RDKit's own ring list (its symmetrized SSSR, in its order) where the
    // port has it; RDKit types from that `RingInfo` throughout.
    let rings = chematic_perception::rdkit_sssr_ring_order(mol).unwrap_or_else(|| {
        chematic_perception::find_symmetrized_sssr(mol)
            .rings()
            .to_vec()
    });
    // MMFF94 has its own, stricter, Kekule-based aromaticity perception
    // (RDKit's `setMMFFAromaticity`), distinct from chematic's own general
    // Huckel model -- most visibly for "mancude" ring systems where a ring
    // atom carries a genuine exocyclic multiple bond (e.g. caffeine's
    // carbonyl carbons, which chematic's own model still treats as
    // aromatic). All typing below reads aromaticity/bond-order state from
    // this re-perceived view, never from `mol` directly, so every
    // `assign_*_type`/`aromatic_*_type` function's existing `atom.aromatic`
    // and `BondOrder::Aromatic` checks stay correct for MMFF94 purposes
    // without needing their own signatures changed.
    let mmff_mol = compute_mmff94_aromatic_view(mol, &rings)?;

    for (i, ty) in types.iter_mut().enumerate().take(n) {
        let idx = AtomIdx(i as u32);
        *ty = assign_heavy_atom_type(&mmff_mol, &rings, idx)?;
    }
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        if mmff_mol.atom(idx).element == Element::H {
            types[i] = assign_h_type(&mmff_mol, idx, &types)?;
        }
    }

    // Construction-time semantic-compatibility invariant (issue #227,
    // Phase 1B-0): a numeric type assigned to an atom must represent the
    // same element the registry says it does. This can only fire if a
    // `assign_*_type` function above has a bug (assigns a type belonging to
    // the wrong element's registry entry) -- it is a defense-in-depth
    // backstop, not the primary correctness mechanism, but it converts any
    // such bug from a silent wrong-parameter collision (issue #227's
    // `furan` finding: a real-but-semantically-wrong table row silently
    // accepted because the lookup only ever checked `Some`/`None`) into a
    // typed, fail-closed `Err` here instead.
    validate_numeric_type_elements(mol, &types)?;
    Ok((types, mmff_mol))
}

fn assign_heavy_atom_type(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
) -> Result<u8, NumericTypeError> {
    let atom = mol.atom(idx);
    let degree = total_degree(mol, idx);
    let is_ion = |charge| degree == 0 && atom.charge == charge;
    match atom.element {
        Element::C => assign_c_type(mol, rings, idx),
        Element::N => assign_n_type(mol, rings, idx),
        Element::O => assign_o_type(mol, rings, idx),
        Element::S => assign_s_type(mol, rings, idx),
        Element::P => assign_p_type(mol, idx),
        Element::SI => Ok(19),
        Element::F if is_ion(-1) => Ok(89),
        Element::F => Ok(11),
        Element::CL if degree == 4 => Ok(77),
        Element::CL if is_ion(-1) => Ok(90),
        Element::CL => Ok(12),
        Element::BR if is_ion(-1) => Ok(91),
        Element::BR => Ok(13),
        Element::I => Ok(14),
        Element::LI if is_ion(1) => Ok(92),
        Element::NA if is_ion(1) => Ok(93),
        Element::K if is_ion(1) => Ok(94),
        Element::MG if is_ion(2) => Ok(99),
        Element::CA if is_ion(2) => Ok(96),
        Element::ZN if is_ion(2) => Ok(95),
        Element::FE if degree == 0 && atom.charge == 2 => Ok(87),
        Element::FE if degree == 0 && atom.charge == 3 => Ok(88),
        Element::CU if degree == 0 && atom.charge == 1 => Ok(97),
        Element::CU if degree == 0 && atom.charge == 2 => Ok(98),
        Element::H => Ok(0),
        _ => Err(NumericTypeError(format!(
            "unsupported element {:?} at atom {}",
            atom.element, idx.0
        ))),
    }
}

fn validate_numeric_type_elements(mol: &Molecule, types: &[u8]) -> Result<(), NumericTypeError> {
    for (i, &assigned) in types.iter().enumerate() {
        let idx = AtomIdx(i as u32);
        let actual_element = mol.atom(idx).element;
        match chematic_ff_numeric_type_registry_lookup(assigned) {
            Some(info) if info.element == actual_element => {}
            Some(info) => {
                return Err(NumericTypeError(format!(
                    "internal error: atom {i} ({actual_element:?}) was assigned MMFF94 \
                     numeric type {assigned} ({}), whose registry element is {:?} -- \
                     semantic-compatibility invariant violated",
                    info.symbol, info.element
                )));
            }
            None => {
                return Err(NumericTypeError(format!(
                    "internal error: atom {i} ({actual_element:?}) was assigned MMFF94 \
                     numeric type {assigned}, which does not exist in the numeric type \
                     registry"
                )));
            }
        }
    }
    Ok(())
}

fn chematic_ff_numeric_type_registry_lookup(
    id: u8,
) -> Option<&'static crate::mmff94_numeric_type_registry::Mmff94NumericTypeInfo> {
    crate::mmff94_numeric_type_registry::mmff94_numeric_type_info(id)
}

// ── Helper: bond iteration ───────────────────────────────────────────────────

struct BondInfo {
    neighbor: AtomIdx,
    order: BondOrder,
}

/// The atom's bonds in bond-index order (the order a scan of `mol.bonds()`
/// gives), from its adjacency list rather than a scan of every bond.
fn bonds_of(mol: &Molecule, idx: AtomIdx) -> Vec<BondInfo> {
    let mut nbrs: Vec<(u32, BondInfo)> = mol
        .neighbors(idx)
        .map(|(nb, b)| {
            (
                b.0,
                BondInfo {
                    neighbor: nb,
                    order: mol.bond(b).order,
                },
            )
        })
        .collect();
    nbrs.sort_unstable_by_key(|&(b, _)| b);
    nbrs.into_iter().map(|(_, info)| info).collect()
}

fn count_bond_order(mol: &Molecule, idx: AtomIdx, order: BondOrder) -> usize {
    bonds_of(mol, idx)
        .iter()
        .filter(|b| b.order == order)
        .count()
}

// ── Aromatic ring helpers ────────────────────────────────────────────────────
//
// Faithful port of RDKit's `RingMembershipSize`/`isAtomNOxide`
// (`Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp`, tag
// `Release_2026_03_3`, commit `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f` --
// see `scripts/mmff94_provenance/PROVENANCE.md`), not a re-derivation. An
// "aromatic ring of size N" requires every consecutive bond in the ring to
// literally be `BondOrder::Aromatic` (chematic's SMILES parser already
// represents aromatic-SMILES ring bonds this way, matching RDKit's
// `Bond::AROMATIC` requirement in `isRingAromatic`), not just that the
// endpoint atoms carry `aromatic = true`.

fn ring_is_fully_aromatic(mol: &Molecule, ring: &[AtomIdx]) -> bool {
    let n = ring.len();
    (0..n).all(|i| {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        matches!(
            mol.bond_between(a, b).map(|(_, e)| e.order),
            Some(BondOrder::Aromatic)
        )
    })
}

fn atom_in_aromatic_ring_of_size(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    atom: AtomIdx,
    size: usize,
) -> bool {
    rings
        .iter()
        .any(|r| r.len() == size && r.contains(&atom) && ring_is_fully_aromatic(mol, r))
}

/// RDKit's `RingInfo::isAtomInRingOfSize` -- plain ring membership, no
/// aromaticity requirement (used for CR3R/CR4R strained-aliphatic-ring
/// typing, unlike `atom_in_aromatic_ring_of_size` above).
fn atom_in_ring_of_size(rings: &[Vec<AtomIdx>], atom: AtomIdx, size: usize) -> bool {
    rings.iter().any(|r| r.len() == size && r.contains(&atom))
}

fn atoms_share_aromatic_ring_of_size(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    a: AtomIdx,
    b: AtomIdx,
    size: usize,
) -> bool {
    rings.iter().any(|r| {
        r.len() == size && r.contains(&a) && r.contains(&b) && ring_is_fully_aromatic(mol, r)
    })
}

/// RDKit's `getTotalDegree()`: heavy-atom bond count plus implicit/explicit
/// hydrogens.
fn total_degree(mol: &Molecule, idx: AtomIdx) -> usize {
    bonds_of(mol, idx).len() + implicit_hcount(mol, idx) as usize
}

/// Sum of Kekule bond orders plus implicit H count -- RDKit's
/// `getValence(Atom::ValenceType::EXPLICIT) + getNumImplicitHs()`. Requires
/// real (non-`Aromatic`) bond orders on `mol`.
fn total_valence(mol: &Molecule, idx: AtomIdx) -> u32 {
    let explicit: u32 = bonds_of(mol, idx)
        .iter()
        .map(|b| b.order.order_int() as u32)
        .sum();
    explicit + implicit_hcount(mol, idx) as u32
}
