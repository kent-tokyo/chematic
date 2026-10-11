// ---------------------------------------------------------------------------
// Production entry points (A1-1b-1): fallible opt-in API
// ---------------------------------------------------------------------------

/// Error from the RDKit-parity experimental aromaticity API.
///
/// Unlike [`assign_aromaticity_ex`](crate::assign_aromaticity_ex)/
/// [`apply_aromaticity_ex`](crate::apply_aromaticity_ex) (infallible, and
/// unchanged by this addition), this engine requires an explicit
/// kekulization step it does not control the success of, so its entry
/// points return `Result` rather than silently falling back to another
/// algorithm or panicking.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AromaticityError {
    /// The input could not be reduced to a Kekulé form (no `BondOrder::Aromatic`
    /// bonds), which this engine requires as a precondition -- mirrors RDKit's
    /// own pipeline, where `Kekulize` always runs before `setAromaticity`.
    KekulizationFailed {
        /// Human-readable detail from the underlying `chematic_core::KekuleError`.
        reason: String,
    },
    /// A post-computation sanity check failed (e.g. an aromatic bond with a
    /// non-aromatic endpoint atom) -- should never happen for chemically
    /// valid input; surfaced as an error rather than a panic or a silently
    /// wrong result.
    InternalInvariantViolation {
        /// Human-readable detail of which invariant failed.
        reason: String,
    },
}

impl std::fmt::Display for AromaticityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AromaticityError::KekulizationFailed { reason } => {
                write!(f, "rdkit-parity aromaticity: kekulization failed: {reason}")
            }
            AromaticityError::InternalInvariantViolation { reason } => {
                write!(
                    f,
                    "rdkit-parity aromaticity: internal invariant violation: {reason}"
                )
            }
        }
    }
}

impl std::error::Error for AromaticityError {}

/// Clone `mol` with every atom's `aromatic` flag reset to `false`.
///
/// Bond orders (including any `BondOrder::Aromatic`) are copied unchanged --
/// only the atom-level annotation is cleared. This engine derives
/// aromaticity purely from element/charge/bond-order structure (never reads
/// `atom.aromatic`), so clearing stale flags here has no effect on the
/// computation itself; it only ensures the *output* molecule's flags come
/// entirely from this engine's own verdict, never from whatever annotation
/// the caller's input happened to carry in.
fn clear_aromatic_flags(mol: &Molecule) -> Molecule {
    use chematic_core::MoleculeBuilder;
    let mut builder = MoleculeBuilder::new();
    for (_, atom) in mol.atoms() {
        let mut a = atom.clone();
        a.aromatic = false;
        builder.add_atom(a);
    }
    for (_, bond) in mol.bonds() {
        let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
    }
    builder.copy_r_groups_from(mol);
    builder.copy_atom_tags_from(mol);
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.build()
}

/// Preserve a parser-supplied aromatic representation when it is internally
/// self-consistent and contains no locally impossible neutral group-16 donor.
/// RDKit accepts some fused heteroaromatic SMILES whose
/// aromatic graph has no single Kekulé assignment under chematic's matching
/// model; discarding that representation would reject a valid RDKit input.
/// This path is deliberately limited to the literal aromatic bond/endpoint
/// annotation and never invents aromaticity for an aliphatic graph.
/// The representation is preserved (cloned) exactly when this holds: the
/// molecule has at least one aromatic bond and every aromatic bond's
/// endpoints are flagged aromatic. Allocation-free.
fn explicit_aromaticity_is_consistent(mol: &Molecule) -> bool {
    // Parser-supplied lowercase flags are annotations, not an aromaticity
    // oracle.  Reject the specific locally impossible donor state from
    // #767 so it is Kekulized and re-perceived instead of taking the identity
    // shortcut.  Other valid explicit aromatic representations retain the
    // compatibility fallback below.
    if mol.atoms().any(|(idx, atom)| {
        atom.aromatic
            && atom.charge == 0
            && matches!(atom.element.atomic_number(), 8 | 16 | 34 | 52)
            && {
                let h = chematic_core::implicit_hcount(mol, idx) as usize;
                mol.degree(idx) + h > 2
            }
    }) {
        return false;
    }

    let mut any = false;
    for (_, bond) in mol.bonds() {
        if bond.order == BondOrder::Aromatic {
            any = true;
            if !mol.atom(bond.atom1).aromatic || !mol.atom(bond.atom2).aromatic {
                return false;
            }
        }
    }
    any
}

fn model_from_explicit_aromaticity(mol: &Molecule) -> AromaticityModel {
    AromaticityModel::from_atom_bond_sets(
        mol.atoms()
            .filter_map(|(idx, atom)| atom.aromatic.then_some(idx))
            .collect(),
        mol.bonds()
            .filter_map(|(idx, bond)| (bond.order == BondOrder::Aromatic).then_some(idx))
            .collect(),
    )
}

/// Return the already-complete aromatic representation when that can be
/// proven without constructing an SSSR/cycle basis.
///
/// The common SMILES cases are either acyclic, saturated rings, or complete
/// lowercase aromatic ring systems. A linear-time bridge pass is sufficient
/// to prove that none of those can gain another aromatic bond. Mixed
/// aromatic/Kekulé ring systems and non-aromatic cyclic multiple bonds still
/// take the full RDKit-parity perception path; this is deliberately a
/// conservative fast path, not a new aromaticity heuristic.
#[cfg(test)]
fn complete_preperceived_aromaticity(mol: &Molecule) -> Option<Molecule> {
    match preperceived_kind(mol)? {
        ParityShortcut::ClearFlags => Some(clear_aromatic_flags(mol)),
        _ => Some(mol.clone()),
    }
}

/// How the RDKit-parity view of a molecule can be produced.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ParityShortcut {
    /// The view is the molecule itself (an exact copy).
    Identity,
    /// The view is the molecule with every atom's aromatic flag cleared.
    ClearFlags,
    /// Kekulization and full re-perception are needed.
    Full,
}

/// Memoized decision between the shortcuts of
/// [`apply_aromaticity_rdkit_parity_uncached`]: the complete-preperceived
/// test ([`preperceived_kind_with`]), then the pre-kekulization test that the
/// explicit representation is kept ([`unchanged_without_kekulization_with`]).
fn parity_shortcut(mol: &Molecule) -> ParityShortcut {
    *mol.derived(chematic_core::DerivedSlot::RdkitParityShortcut, || {
        let model_may_disagree = crate::rdkit_model_hook::rdkit_model_may_disagree(mol);
        // Molecules whose view the RDKit-model hook may correct take the
        // full path, where the hook runs.
        if model_may_disagree && crate::rdkit_model_hook::hook_for(mol).is_some() {
            return ParityShortcut::Full;
        }
        // RDKit's sanitization rewrites some graphs (perchlorate,
        // organometallic dative bonds); the view is then never `mol`. Every
        // cleanup candidate is already included in the broader model gate,
        // so ordinary molecules do not need a second full atom scan.
        if model_may_disagree && crate::rdkit_cleanup::rdkit_sanitize_cleanup_needed(mol) {
            return ParityShortcut::Full;
        }
        let info = ComponentInfo::new(mol);
        match preperceived_kind_with(mol, &info) {
            Some(kind) => kind,
            None if unchanged_without_kekulization_with(mol, &info) => ParityShortcut::Identity,
            None if explicit_kept_by_kekule_verdict(mol, &info) => ParityShortcut::Identity,
            None => ParityShortcut::Full,
        }
    })
}

/// Whether the RDKit-parity aromatic view of `mol`
/// ([`apply_aromaticity_rdkit_parity_shared`]) is an exact copy of `mol`
/// (same atoms, flags, bonds and orders). Linear time, memoized; lets
/// callers that only read the view (SMARTS matching) use `mol` directly.
pub fn rdkit_parity_view_is_identity(mol: &Molecule) -> bool {
    parity_shortcut(mol) == ParityShortcut::Identity
}

/// Which copy of `mol` is already the complete RDKit-parity representation,
/// when that can be proven without an SSSR (`None` otherwise).
///
/// The common SMILES cases are either acyclic, saturated rings, or complete
/// lowercase aromatic ring systems. A linear-time bridge pass is sufficient
/// to prove that none of those can gain another aromatic bond. Mixed
/// aromatic/Kekulé ring systems and non-aromatic cyclic multiple bonds still
/// take the full RDKit-parity perception path; this is deliberately a
/// conservative fast path, not a new aromaticity heuristic.
///
/// Per connected component of the cyclic subgraph: a component with both
/// aromatic and non-aromatic cyclic bonds, or a non-aromatic component with
/// a cyclic multiple bond, needs the full path. A formally saturated
/// component can still become aromatic through exocyclic C=N/C=O/C=S bonds
/// (for example O=C1NNC(=O)N1); for a ring to qualify, its atoms that pass
/// RDKit's donor/candidate rules must themselves contain a cycle, so a
/// non-aromatic component with an incident multiple bond also needs the full
/// path exactly when its candidate-induced cyclic subgraph has a cycle.
#[cfg(test)]
fn preperceived_kind(mol: &Molecule) -> Option<ParityShortcut> {
    preperceived_kind_with(mol, &ComponentInfo::new(mol))
}

/// Per cyclic component facts shared by [`preperceived_kind_with`] and
/// [`unchanged_without_kekulization_with`].
struct ComponentInfo {
    cyclic: std::sync::Arc<Vec<bool>>,
    /// Component label of every atom with a cyclic bond (`u32::MAX` otherwise).
    root: Vec<u32>,
    /// Per label: `COMP_*` bits.
    kind: Vec<u8>,
    /// Per label: cyclic bonds and atoms.
    edges: Vec<u32>,
    atoms: Vec<u32>,
    has_aromatic_bond: bool,
    /// Every aromatic bond joins two flagged atoms (and one exists), and no
    /// flagged neutral group-16 atom is hypercoordinate.
    explicit: bool,
    non_ring_aromatic_bond: bool,
    any_flagged: bool,
}

/// Some cyclic bond of the component is aromatic.
const COMP_AROMATIC: u8 = 1;
/// Some cyclic bond is not aromatic.
const COMP_NON_AROMATIC: u8 = 2;
/// Some cyclic bond is neither aromatic nor single (up/down included).
const COMP_NON_AROMATIC_MULTIPLE: u8 = 4;
/// Some atom of the component has a double, triple or quadruple bond.
const COMP_INCIDENT_MULTIPLE: u8 = 8;
/// Some atom of the component is flagged aromatic.
const COMP_FLAGGED: u8 = 16;

impl ComponentInfo {
    fn new(mol: &Molecule) -> Self {
        let n = mol.atom_count();
        let cyclic = crate::sssr::ring_bond_flags_shared(mol);
        let components = crate::sssr::ring_components_shared(mol);
        let c = components.count as usize;
        let mut root = vec![u32::MAX; n];
        let mut kind = vec![0u8; c];
        let mut edges = vec![0u32; c];
        let mut atoms = vec![0u32; c];
        let mut has_aromatic_bond = false;
        let mut consistent = true;
        let mut non_ring_aromatic_bond = false;
        for (bond_idx, bond) in mol.bonds() {
            if bond.order == BondOrder::Aromatic {
                has_aromatic_bond = true;
                if !mol.atom(bond.atom1).aromatic || !mol.atom(bond.atom2).aromatic {
                    consistent = false;
                }
            }
            if cyclic[bond_idx.0 as usize] {
                let r = components.label[bond.atom1.0 as usize];
                root[bond.atom1.0 as usize] = r;
                root[bond.atom2.0 as usize] = r;
                let r = r as usize;
                edges[r] += 1;
                kind[r] |= match bond.order {
                    BondOrder::Aromatic => COMP_AROMATIC,
                    BondOrder::Single | BondOrder::Up | BondOrder::Down => COMP_NON_AROMATIC,
                    _ => COMP_NON_AROMATIC | COMP_NON_AROMATIC_MULTIPLE,
                };
            } else if bond.order == BondOrder::Aromatic {
                non_ring_aromatic_bond = true;
            }
        }
        let mut any_flagged = false;
        for (idx, atom) in mol.atoms() {
            any_flagged |= atom.aromatic;
            if atom.aromatic
                && atom.charge == 0
                && matches!(atom.element.atomic_number(), 8 | 16 | 34 | 52)
                && mol.degree(idx) + chematic_core::implicit_hcount(mol, idx) as usize > 2
            {
                consistent = false;
            }
            let r = root[idx.0 as usize];
            if r != u32::MAX {
                atoms[r as usize] += 1;
                if atom.aromatic {
                    kind[r as usize] |= COMP_FLAGGED;
                }
            }
        }
        for (_, bond) in mol.bonds() {
            if matches!(
                bond.order,
                BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
            ) {
                for end in [bond.atom1, bond.atom2] {
                    let r = root[end.0 as usize];
                    if r != u32::MAX {
                        kind[r as usize] |= COMP_INCIDENT_MULTIPLE;
                    }
                }
            }
        }
        Self {
            cyclic,
            root,
            kind,
            edges,
            atoms,
            has_aromatic_bond,
            explicit: has_aromatic_bond && consistent,
            non_ring_aromatic_bond,
            any_flagged,
        }
    }

    /// Labels of the cyclic components (components with a cyclic bond).
    fn roots(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.edges.len()).filter(|&c| self.edges[c] > 0)
    }
}

fn preperceived_kind_with(mol: &Molecule, info: &ComponentInfo) -> Option<ParityShortcut> {
    if info.has_aromatic_bond && !info.explicit {
        return None;
    }
    if info.non_ring_aromatic_bond {
        return None;
    }

    let mut needs_candidate_check = false;
    for r in info.roots() {
        let f = info.kind[r];
        let aromatic = f & COMP_AROMATIC != 0;
        if (aromatic && f & COMP_NON_AROMATIC != 0)
            || (!aromatic && f & COMP_NON_AROMATIC_MULTIPLE != 0)
        {
            return None;
        }
        if !aromatic && f & COMP_INCIDENT_MULTIPLE != 0 {
            needs_candidate_check = true;
        }
    }

    if needs_candidate_check {
        let n = mol.atom_count();
        let checked = |atom: usize| -> bool {
            let r = info.root[atom];
            r != u32::MAX && {
                let f = info.kind[r as usize];
                f & COMP_AROMATIC == 0 && f & COMP_INCIDENT_MULTIPLE != 0
            }
        };
        let mut candidate = vec![false; n];
        for (atom, slot) in candidate.iter_mut().enumerate() {
            if checked(atom) {
                let idx = AtomIdx(atom as u32);
                let donor = get_atom_electron_donor_type(mol, idx, info.cyclic.as_slice());
                *slot = is_atom_candidate_for_aromaticity(mol, idx, donor);
            }
        }
        // A cyclic bond between two candidates that joins atoms already
        // connected through candidate cyclic bonds closes a candidate cycle.
        let mut candidate_uf = UnionFind::new(n);
        for (bond_idx, bond) in mol.bonds() {
            let (a, b) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
            if info.cyclic[bond_idx.0 as usize]
                && candidate[a]
                && candidate[b]
                && !candidate_uf.union(a as u32, b as u32)
            {
                return None;
            }
        }
    }

    if info.explicit {
        Some(ParityShortcut::Identity)
    } else if info.any_flagged {
        Some(ParityShortcut::ClearFlags)
    } else {
        Some(ParityShortcut::Identity)
    }
}

/// Minimal union-find over atom indices.
struct UnionFind {
    parent: Vec<u32>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n as u32).collect(),
        }
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let p = self.parent[x as usize];
            self.parent[x as usize] = self.parent[p as usize];
            x = p;
        }
        x
    }

    /// Join the sets of `a` and `b`; `false` if they were already joined.
    fn union(&mut self, a: u32, b: u32) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        self.parent[ra as usize] = rb;
        true
    }
}

/// Whether `rdkit_parity_aromaticity(kekulized)` could mark any atom or bond
/// that the explicit representation `mol` does not already mark aromatic.
///
/// That verdict only marks atoms and bonds of SSSR rings made entirely of
/// aromaticity candidates (donor/candidate rules evaluated on `kekulized`
/// against the global ring-bond set, which equals the cyclic bonds). Every
/// such ring is a cycle of the subgraph induced by candidate atoms over cyclic
/// bonds, so each bond it marks is a non-bridge edge of that subgraph. If all
/// of those edges are already explicit aromatic bonds between flagged atoms,
/// the verdict is a subset of the explicit sets and cannot strictly extend
/// them — decided without an SSSR. `kekulized` must be index-aligned with
/// `mol` (same atoms and bonds).
fn re_perception_can_extend(mol: &Molecule, kekulized: &Molecule) -> bool {
    candidate_cycles_can_extend(mol, kekulized)
}

/// Whether the full RDKit-parity path would return `mol` unchanged, decided
/// without kekulizing or building an SSSR (`false` means "not proven").
///
/// Applies to molecules whose aromaticity is either explicit and consistent,
/// or absent altogether (no aromatic atom or bond). The full path then
/// returns an exact copy of `mol` whenever the RDKit-parity verdict marks
/// nothing outside the explicit aromatic atoms and bonds: with explicit
/// input it keeps the explicit form unless the verdict strictly extends it
/// (or kekulization fails, which also keeps it); with no aromaticity the
/// molecule is its own Kekulé form and an empty verdict changes nothing.
///
/// The verdict only marks outer-perimeter bonds (and their atoms) of SSSR
/// rings made entirely of aromaticity candidates, i.e. edges on a cycle of
/// the candidate-induced cyclic subgraph, and the cycles of distinct cyclic
/// components are disjoint. Per component:
///
/// * every cyclic bond aromatic: all such edges are already explicit;
/// * no aromatic atom or bond: kekulization leaves every bond of the
///   component as it is, so candidacy and donor types computed on `mol` are
///   exact. With no candidate cycle nothing is marked; a component that is a
///   single ring made entirely of candidates is its own SSSR ring and its
///   own fused group, so it is marked exactly when [`apply_huckel`] accepts
///   it; anything else is left to the full path;
/// * otherwise (mixed): only atoms flagged aromatic change bonds under
///   kekulization, and the donor/candidate rules read an atom's own flag,
///   element, charge, H count and incident bond orders plus its neighbours'
///   elements, so a non-aromatic atom's candidacy is the same on `mol` and on
///   any Kekulé form. Flagged atoms are taken as candidates unconditionally,
///   giving a supergraph of the real candidate subgraph; an edge on a cycle
///   of the real subgraph lies on a cycle of the supergraph, so if no such
///   edge is non-explicit here, none is there.
fn unchanged_without_kekulization_with(mol: &Molecule, info: &ComponentInfo) -> bool {
    if info.has_aromatic_bond {
        if !info.explicit {
            return false;
        }
    } else if info.any_flagged {
        return false;
    }

    let n = mol.atom_count();
    let cyclic: &[bool] = info.cyclic.as_slice();
    // Exact components: no aromatic cyclic bond and no flagged atom.
    let exact =
        |r: u32| r != u32::MAX && info.kind[r as usize] & (COMP_AROMATIC | COMP_FLAGGED) == 0;

    // Candidacy: exact for atoms not flagged aromatic, assumed for flagged ones.
    let mut candidate = vec![false; n];
    let mut donor: Vec<ElectronDonorType> = vec![ElectronDonorType::None; n];
    for (idx, atom) in mol.atoms() {
        let a = idx.0 as usize;
        if info.root[a] == u32::MAX {
            continue;
        }
        if atom.aromatic {
            candidate[a] = true;
        } else {
            let d = get_atom_electron_donor_type(mol, idx, cyclic);
            donor[a] = d;
            candidate[a] = is_atom_candidate_for_aromaticity(mol, idx, d);
        }
    }

    // Exact components: find candidate cycles.
    let mut candidate_uf = UnionFind::new(n);
    let mut has_candidate_cycle = vec![false; n];
    for (idx, bond) in mol.bonds() {
        let (a, b) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        if !cyclic[idx.0 as usize] || !candidate[a] || !candidate[b] {
            continue;
        }
        let r = info.root[a];
        if exact(r) && !candidate_uf.union(a as u32, b as u32) {
            has_candidate_cycle[r as usize] = true;
        }
    }
    for r in info.roots() {
        if !exact(r as u32) || !has_candidate_cycle[r] {
            continue;
        }
        if info.edges[r] != info.atoms[r] {
            return false; // a candidate cycle in a fused system: full path
        }
        // A single ring made entirely of candidates.
        let members: Vec<AtomIdx> = (0..n)
            .filter(|&a| info.root[a] as usize == r)
            .map(|a| AtomIdx(a as u32))
            .collect();
        if apply_huckel(mol, &members, &donor) {
            return false;
        }
    }
    // Remaining components: exact ones contribute nothing (above), so only
    // mixed components can hold a non-explicit candidate cycle edge.
    for (slot, &root) in candidate.iter_mut().zip(&info.root) {
        if exact(root) {
            *slot = false;
        }
    }
    !non_explicit_candidate_cycle_edge(mol, cyclic, &candidate)
}

/// Whether the full path keeps an explicit aromatic representation, decided
/// from the RDKit-parity verdict of its Kekulé form computed without
/// kekulizing (`false` means "not proven").
///
/// With consistent explicit aromaticity the full path returns `mol` unless
/// kekulization succeeds and the verdict on the Kekulé form contains every
/// explicit aromatic atom and bond and more. Kekulization only reassigns
/// aromatic bonds (between flagged atoms): when it succeeds, every
/// must-match atom (`atom_must_be_matched`) gets exactly one double bond on
/// an aromatic bond to another must-match atom and every other bond of the
/// molecule is single or unchanged. The one exception is the bridgehead-N
/// fallback, which can leave a must-match N with three or more such bonds
/// unmatched; molecules where it could apply are left to the full path.
/// When every aromatic bond of a must-match atom is cyclic, the typing rules
/// see the same thing whichever of them is the double one, so each flagged
/// atom's donor type and candidacy follow from its matched status alone
/// ([`KekuleAtom`]); unflagged atoms keep all their bonds. The SSSR ring set
/// is used only when it is decidable from the graph alone (the Kekulé form
/// has the same graph and adjacency). Kekulization failure keeps `mol` as
/// well, so the conclusion "not strictly extended" holds either way.
fn explicit_kept_by_kekule_verdict(mol: &Molecule, info: &ComponentInfo) -> bool {
    if !info.explicit {
        return false;
    }
    let n = mol.atom_count();
    let cyclic: &[bool] = info.cyclic.as_slice();
    // Matched status and the bond carrying the double bond, per atom.
    let mut must = vec![false; n];
    let mut n_flagged = 0usize;
    for (idx, atom) in mol.atoms() {
        if atom.aromatic {
            n_flagged += 1;
            must[idx.0 as usize] = chematic_core::kekulization::atom_must_be_matched(mol, idx);
        }
    }
    let mut double: Vec<Option<BondIdx>> = vec![None; n];
    for (idx, atom) in mol.atoms() {
        if !atom.aromatic {
            continue;
        }
        let a = idx.0 as usize;
        let mut has_aromatic_bond = false;
        let mut must_neighbors = 0usize;
        let mut first = None;
        for &(nb, b) in mol.neighbor_slice(idx) {
            if mol.bond(b).order != BondOrder::Aromatic {
                continue;
            }
            has_aromatic_bond = true;
            if must[a] && !cyclic[b.0 as usize] {
                return false; // the double bond might be exocyclic
            }
            if must[nb.0 as usize] {
                must_neighbors += 1;
                first.get_or_insert(b);
            }
        }
        if !has_aromatic_bond {
            return false; // a flag the Kekulé form types differently
        }
        if must[a] {
            if must_neighbors == 0 {
                return true; // kekulization fails: `mol` is kept
            }
            if atom.element.atomic_number() == 7 && must_neighbors >= 3 {
                return false; // bridgehead-N fallback possible
            }
            double[a] = first;
        }
    }
    let n_aromatic_bonds = mol
        .bonds()
        .filter(|(_, b)| b.order == BondOrder::Aromatic)
        .count();

    let Some((atoms, bonds)) = verdict_with_typing(
        mol,
        6,
        |idx, cyclic| {
            if mol.atom(idx).aromatic {
                let view = KekuleAtom {
                    double: double[idx.0 as usize],
                };
                let d = get_atom_electron_donor_type_in(mol, idx, cyclic, &view);
                (d, is_atom_candidate_for_aromaticity_in(mol, idx, d, &view))
            } else {
                let d = get_atom_electron_donor_type(mol, idx, cyclic);
                (d, is_atom_candidate_for_aromaticity(mol, idx, d))
            }
        },
        false,
    ) else {
        return false;
    };
    let preserves = mol
        .atoms()
        .all(|(idx, a)| !a.aromatic || atoms.contains(&idx))
        && mol
            .bonds()
            .all(|(idx, b)| b.order != BondOrder::Aromatic || bonds.contains(&idx));
    let strictly_extends = atoms.len() > n_flagged || bonds.len() > n_aromatic_bonds;
    !(preserves && strictly_extends)
}

/// Body of [`re_perception_can_extend`]: candidacy of every ring atom read
/// from `typing` (index-aligned with `mol`) by the donor/candidate rules.
fn candidate_cycles_can_extend(mol: &Molecule, typing: &Molecule) -> bool {
    let n = mol.atom_count();
    let cyclic = crate::sssr::ring_bond_flags_shared(mol);
    let cyclic: &[bool] = cyclic.as_slice();
    let mut ring_atom = vec![false; n];
    for (idx, bond) in mol.bonds() {
        if cyclic[idx.0 as usize] {
            ring_atom[bond.atom1.0 as usize] = true;
            ring_atom[bond.atom2.0 as usize] = true;
        }
    }
    let mut candidate = vec![false; n];
    for i in 0..n {
        if ring_atom[i] {
            let idx = AtomIdx(i as u32);
            let donor = get_atom_electron_donor_type(typing, idx, cyclic);
            candidate[i] = is_atom_candidate_for_aromaticity(typing, idx, donor);
        }
    }
    non_explicit_candidate_cycle_edge(mol, cyclic, &candidate)
}

/// Whether some cyclic bond between two candidates lies on a cycle of the
/// candidate-induced cyclic subgraph and is not already an explicit aromatic
/// bond between two flagged atoms.
fn non_explicit_candidate_cycle_edge(mol: &Molecule, cyclic: &[bool], candidate: &[bool]) -> bool {
    // An edge between two candidates can only lie on a candidate cycle if
    // it is itself not explicit; with no such edge there is nothing to test.
    let explicit = |b: BondIdx, u: AtomIdx, v: AtomIdx| {
        mol.bond(b).order == BondOrder::Aromatic && mol.atom(u).aromatic && mol.atom(v).aromatic
    };
    let is_edge = |b: BondIdx, u: AtomIdx, v: AtomIdx| {
        cyclic[b.0 as usize] && candidate[u.0 as usize] && candidate[v.0 as usize]
    };
    if !mol
        .bonds()
        .any(|(b, bond)| is_edge(b, bond.atom1, bond.atom2) && !explicit(b, bond.atom1, bond.atom2))
    {
        return false;
    }
    // Bridges of the candidate-induced cyclic subgraph (iterative low-link
    // DFS over the molecule's adjacency, restricted to candidate edges): an
    // edge lies on a cycle exactly when it is a back edge or a tree edge
    // whose child subtree reaches its parent or above.
    let n = mol.atom_count();
    const UNSEEN: u32 = u32::MAX;
    const NO_BOND: u32 = u32::MAX;
    let mut dl = vec![(UNSEEN, UNSEEN); n];
    let mut on_cycle = vec![false; mol.bond_count()];
    let mut time = 0u32;
    let mut stack: Vec<(u32, u32, u32)> = Vec::new();
    for root in 0..n {
        if dl[root].0 != UNSEEN || !candidate[root] {
            continue;
        }
        dl[root] = (time, time);
        time += 1;
        stack.push((root as u32, NO_BOND, 0));
        while let Some(&mut (atom, parent_bond, ref mut next)) = stack.last_mut() {
            let a = atom as usize;
            let nbrs = mol.neighbor_slice(AtomIdx(atom));
            if (*next as usize) < nbrs.len() {
                let (nb, bond) = nbrs[*next as usize];
                *next += 1;
                if bond.0 == parent_bond || !is_edge(bond, AtomIdx(atom), nb) {
                    continue;
                }
                let b = nb.0 as usize;
                if dl[b].0 == UNSEEN {
                    dl[b] = (time, time);
                    time += 1;
                    stack.push((nb.0, bond.0, 0));
                } else {
                    if dl[b].0 < dl[a].1 {
                        dl[a].1 = dl[b].0;
                    }
                    on_cycle[bond.0 as usize] = true; // back edge
                }
                continue;
            }
            stack.pop();
            if parent_bond != NO_BOND {
                let parent = stack
                    .last()
                    .expect("non-root DFS frame must retain its parent")
                    .0 as usize;
                let low_a = dl[a].1;
                if low_a < dl[parent].1 {
                    dl[parent].1 = low_a;
                }
                if low_a <= dl[parent].0 {
                    on_cycle[parent_bond as usize] = true;
                }
            }
        }
    }
    mol.bonds().any(|(b, bond)| {
        on_cycle[b.0 as usize]
            && is_edge(b, bond.atom1, bond.atom2)
            && !explicit(b, bond.atom1, bond.atom2)
    })
}

/// `apply_kekule(&clear_aromatic_flags(mol), k)`, built as one modified copy
/// of `mol` when that is provably the same molecule.
///
/// `clear_aromatic_flags` rebuilds `mol` (same atoms, bonds and side
/// channels) with every aromatic flag cleared; `apply_kekule` returns that
/// unchanged for an empty `k` and otherwise rebuilds it with the new orders,
/// copying every side channel except the R-group labels. A rebuild re-adds
/// bonds in index order, so it equals a copy exactly when `mol`'s adjacency
/// lists are already in that order.
fn kekulized_copy(mol: &Molecule, k: &chematic_core::KekuleResult) -> Molecule {
    if !mol.adjacency_in_bond_order() {
        return chematic_core::apply_kekule(&clear_aromatic_flags(mol), k);
    }
    let mut out = mol.clone();
    for (idx, atom) in mol.atoms() {
        if atom.aromatic {
            out.set_atom_aromatic(idx, false);
        }
    }
    if !k.is_empty() {
        for (&b, &order) in k {
            out.set_bond_order(b, order);
        }
        out.clear_r_group_labels();
    }
    out
}

/// The aromatic atom/bond sets `rdkit_parity_aromaticity` computes on a
/// kekulized molecule.
type AromaticVerdict = (FxHashSet<AtomIdx>, FxHashSet<BondIdx>);

/// Normalize `mol` for RDKit-parity aromaticity. A mixed aromatic/Kekulé
/// input can contain additional atoms that RDKit's sanitizer promotes into
/// the final aromatic system. Accept that re-perception only when it strictly
/// extends the supplied aromatic atom/bond sets without removing anything;
/// otherwise preserve the explicit representation. The latter matters for
/// large fused cages where a valid alternate Kekulé assignment can change
/// the aromatic partition. A non-representable graph without an explicit
/// fallback returns an error; no partial rewrite is exposed.
///
/// When the re-perception is accepted, its verdict on the returned Kekulé
/// molecule is returned too, so callers need not compute it a second time.
fn kekulize_for_rdkit_parity_with_verdict(
    mol: &Molecule,
) -> Result<(Molecule, Option<AromaticVerdict>), AromaticityError> {
    let explicit_ok = explicit_aromaticity_is_consistent(mol);
    // Kekulization never reads aromatic flags, so kekulizing `mol` is the
    // same as kekulizing its flag-cleared copy.
    match chematic_core::kekulize(mol) {
        Ok(k) => {
            let kekulized = kekulized_copy(mol, &k);
            // Same graph; kekulization only turns aromatic bonds into single
            // or double ones, all ring-eligible, so the cyclic bonds agree.
            crate::sssr::seed_ring_data_from(&kekulized, mol);
            if !explicit_ok {
                return Ok((kekulized, None));
            }
            if !re_perception_can_extend(mol, &kekulized) {
                // The verdict could only be a subset of the explicit sets, so
                // the explicit representation is kept (see below).
                return Ok((mol.clone(), None));
            }

            let (candidate_atoms, candidate_bonds) = rdkit_parity_aromaticity(&kekulized);
            let mut explicit_atoms = 0usize;
            let mut preserves_explicit = true;
            for (idx, atom) in mol.atoms() {
                if atom.aromatic {
                    explicit_atoms += 1;
                    preserves_explicit &= candidate_atoms.contains(&idx);
                }
            }
            let mut explicit_bonds = 0usize;
            for (idx, bond) in mol.bonds() {
                if bond.order == BondOrder::Aromatic {
                    explicit_bonds += 1;
                    preserves_explicit &= candidate_bonds.contains(&idx);
                }
            }
            let strictly_extends =
                candidate_atoms.len() > explicit_atoms || candidate_bonds.len() > explicit_bonds;

            if preserves_explicit && strictly_extends {
                Ok((kekulized, Some((candidate_atoms, candidate_bonds))))
            } else {
                Ok((mol.clone(), None))
            }
        }
        Err(e) => {
            if explicit_ok {
                Ok((mol.clone(), None))
            } else {
                Err(AromaticityError::KekulizationFailed { reason: e.detail })
            }
        }
    }
}

/// Every aromatic bond's two endpoint atoms must themselves be in the
/// aromatic atom set -- a basic well-formedness property of any Hückel
/// verdict. Cheap to check, catches a class of bug that would otherwise
/// surface downstream as a confusing SMILES/valence inconsistency instead
/// of a clear error at the source.
fn validate_aromaticity_invariants(
    mol: &Molecule,
    atoms: &FxHashSet<AtomIdx>,
    bonds: &FxHashSet<BondIdx>,
) -> Result<(), AromaticityError> {
    for &bidx in bonds {
        let bond = mol.bond(bidx);
        if !atoms.contains(&bond.atom1) || !atoms.contains(&bond.atom2) {
            return Err(AromaticityError::InternalInvariantViolation {
                reason: format!(
                    "aromatic bond {bidx:?} ({:?}-{:?}) has a non-aromatic endpoint atom",
                    bond.atom1, bond.atom2
                ),
            });
        }
    }
    Ok(())
}

fn assign_from_kekulized(kekulized: &Molecule) -> Result<AromaticityModel, AromaticityError> {
    assign_from_verdict(kekulized, rdkit_parity_aromaticity(kekulized))
}

fn assign_from_verdict(
    kekulized: &Molecule,
    (atoms, bonds): AromaticVerdict,
) -> Result<AromaticityModel, AromaticityError> {
    validate_aromaticity_invariants(kekulized, &atoms, &bonds)?;
    Ok(AromaticityModel::from_atom_bond_sets(atoms, bonds))
}

/// Assign aromaticity using the RDKit-parity reference engine
/// (`rdkit_parity_aromaticity`, see the module doc comment).
///
/// Explicitly opt-in and separate from [`assign_aromaticity_ex`]/
/// [`AromaticityAlgorithm`] -- those remain infallible and unchanged. This
/// function is fallible because it performs its own kekulization
/// internally (this engine requires pre-kekulized input); on failure,
/// `mol` is never touched and no partial result is produced.
///
/// The returned model's atom/bond indices correspond 1:1 with `mol`'s own
/// indices (kekulization here is index-preserving: it only clears stale
/// aromatic flags and normalizes bond orders, never adds/removes/reorders
/// atoms or bonds).
///
/// [`ring_classifications`](AromaticityModel::ring_classifications) and
/// [`antiaromatic_rings`](AromaticityModel::antiaromatic_rings) are always
/// empty on the returned model -- this engine (like RDKit's own) determines
/// only the aromatic atom/bond sets, not a per-ring classification or
/// antiaromaticity verdict.
///
/// [`assign_aromaticity_ex`]: crate::assign_aromaticity_ex
/// [`AromaticityAlgorithm`]: crate::AromaticityAlgorithm
pub fn assign_aromaticity_rdkit_parity_experimental(
    mol: &Molecule,
) -> Result<AromaticityModel, AromaticityError> {
    match parity_shortcut(mol) {
        ParityShortcut::Identity => return Ok(model_from_explicit_aromaticity(mol)),
        ParityShortcut::ClearFlags => {
            return Ok(model_from_explicit_aromaticity(&clear_aromatic_flags(mol)));
        }
        ParityShortcut::Full => {}
    }
    let (kekulized, verdict) = kekulize_for_rdkit_parity_with_verdict(mol)?;
    if kekulized
        .bonds()
        .any(|(_, bond)| bond.order == BondOrder::Aromatic)
    {
        return Ok(model_from_explicit_aromaticity(&kekulized));
    }
    match verdict {
        Some(v) => assign_from_verdict(&kekulized, v),
        None => assign_from_kekulized(&kekulized),
    }
}

/// Apply aromaticity using the RDKit-parity reference engine, returning a
/// new [`Molecule`] with atom/bond flags set according to the computed
/// model.
///
/// See [`assign_aromaticity_rdkit_parity_experimental`] for the fallibility
/// contract (kekulization failure is reported, never silently substituted
/// or partially applied) and the index-correspondence guarantee.
pub fn apply_aromaticity_rdkit_parity_experimental(
    mol: &Molecule,
) -> Result<Molecule, AromaticityError> {
    // Memoized on `mol`: descriptors and the RDKit-compatible fingerprints
    // all start from this perceived copy.
    let view = (*apply_aromaticity_rdkit_parity_shared(mol)).clone()?;
    // Clones start with an empty cache: re-mark the copy as already holding
    // RDKit's dative bonds.
    crate::rdkit_sssr_order::seed_no_further_dative_bonds(&view);
    Ok(view)
}

/// Shared, memoized [`apply_aromaticity_rdkit_parity_experimental`] result
/// (no copy of the perceived molecule). The perceived molecule keeps its own
/// derived caches (rings, ring flags) across calls.
pub fn apply_aromaticity_rdkit_parity_shared(
    mol: &Molecule,
) -> std::sync::Arc<Result<Molecule, AromaticityError>> {
    mol.derived(chematic_core::DerivedSlot::RdkitParityAromatic, || {
        apply_aromaticity_rdkit_parity_uncached(mol)
    })
}

/// Run `f` on the RDKit-parity aromatic view of `mol`
/// ([`apply_aromaticity_rdkit_parity_shared`]) without copying it: `mol`
/// itself when the view would be an exact copy of it
/// ([`rdkit_parity_view_is_identity`]), the memoized view otherwise.
pub fn with_rdkit_parity_view<R>(
    mol: &Molecule,
    f: impl FnOnce(Result<&Molecule, &AromaticityError>) -> R,
) -> R {
    if rdkit_parity_view_is_identity(mol) {
        return f(Ok(mol));
    }
    let view = apply_aromaticity_rdkit_parity_shared(mol);
    f(view.as_ref().as_ref())
}

fn apply_aromaticity_rdkit_parity_uncached(mol: &Molecule) -> Result<Molecule, AromaticityError> {
    // The gate flags only molecules whose (memoized) shortcut is `Full`.
    let hook = if parity_shortcut(mol) == ParityShortcut::Full {
        crate::rdkit_model_hook::hook_for(mol)
    } else {
        None
    };
    let Some(hook) = hook else {
        return apply_aromaticity_rdkit_parity_parity_only(mol);
    };
    let view = apply_aromaticity_rdkit_parity_parity_only(mol);
    let base = match &view {
        Ok(v) => v,
        Err(_) => mol,
    };
    match hook(mol, base) {
        Some(corrected) => {
            crate::rdkit_sssr_order::seed_no_further_dative_bonds(&corrected);
            Ok(corrected)
        }
        None => view,
    }
}

/// The RDKit-parity view without the RDKit-model hook (diagnostics).
#[doc(hidden)]
pub fn aromaticity_rdkit_parity_only(mol: &Molecule) -> Result<Molecule, AromaticityError> {
    apply_aromaticity_rdkit_parity_parity_only(mol)
}

fn apply_aromaticity_rdkit_parity_parity_only(
    mol: &Molecule,
) -> Result<Molecule, AromaticityError> {
    if let Some(cleaned) = crate::rdkit_cleanup::rdkit_sanitize_cleanup(mol) {
        // The view of RDKit's sanitized graph; its ring data is its own.
        let view = apply_aromaticity_rdkit_parity_unseeded(&cleaned)?;
        crate::rdkit_sssr_order::seed_no_further_dative_bonds(&view);
        return Ok(view);
    }
    let view = apply_aromaticity_rdkit_parity_unseeded(mol)?;
    crate::rdkit_sssr_order::seed_no_further_dative_bonds(&view);
    // Same graph with the same ring-eligible bonds (the view only changes
    // aromatic flags and turns bonds aromatic, single or double), so the
    // cyclic-bond flags carry over.
    crate::sssr::seed_ring_data_from(&view, mol);
    Ok(view)
}

fn apply_aromaticity_rdkit_parity_unseeded(mol: &Molecule) -> Result<Molecule, AromaticityError> {
    match parity_shortcut(mol) {
        ParityShortcut::Identity => return Ok(mol.clone()),
        ParityShortcut::ClearFlags => return Ok(clear_aromatic_flags(mol)),
        ParityShortcut::Full => {}
    }
    let (kekulized, verdict) = kekulize_for_rdkit_parity_with_verdict(mol)?;
    if kekulized
        .bonds()
        .any(|(_, bond)| bond.order == BondOrder::Aromatic)
    {
        return Ok(kekulized);
    }
    let model = match verdict {
        Some(v) => assign_from_verdict(&kekulized, v)?,
        None => assign_from_kekulized(&kekulized)?,
    };
    Ok(crate::aromaticity::build_molecule_from_model(
        &kekulized, &model,
    ))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

