//! ECFP (Extended Connectivity Fingerprints) based on the Morgan algorithm.
//!
//! Uses FNV-1a 64-bit hashing for reproducibility and WASM-compatibility.

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};
use chematic_perception::find_sssr;
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

use crate::bitvec::BitVec2048;

const FNV_OFFSET: u64 = 14695981039346656037;
const FNV_PRIME: u64 = 1099511628211;

/// Compute the FNV-1a 64-bit hash of `bytes`.
pub(crate) fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h = FNV_OFFSET;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

/// Hash one atom's neighbourhood at iteration `r` of the Morgan expansion.
///
/// Byte layout: `[r as u8, self_id (8 bytes), (bond_type (1) ++ nb_id (8))*]`
/// Neighbours are sorted before hashing to make the result order-independent.
pub(crate) fn expand_atom_id(mol: &Molecule, i: usize, r: u32, ids: &[u64]) -> u64 {
    let idx = AtomIdx(i as u32);
    // SmallVec<6>: typical atoms have ≤4 heavy neighbors; avoids heap alloc for ~95% of calls.
    let mut neighbor_info: SmallVec<[(u8, u64); 6]> = mol
        .neighbors(idx)
        .map(|(nb_idx, bond_idx)| {
            (
                bond_type_int(mol.bond(bond_idx).order),
                ids[nb_idx.0 as usize],
            )
        })
        .collect();
    neighbor_info.sort_unstable();

    // 1 (radius) + 8 (self id) + up to 6 × 9 (bond_type + nb_id) = 63 bytes max on stack.
    let mut bytes: SmallVec<[u8; 64]> = SmallVec::new();
    bytes.push(r as u8);
    bytes.extend_from_slice(&ids[i].to_le_bytes());
    for (btype, nb_id) in &neighbor_info {
        bytes.push(*btype);
        bytes.extend_from_slice(&nb_id.to_le_bytes());
    }
    fnv1a(&bytes)
}

/// Which atom-invariant definition [`initial_atom_id`] (and everything built on
/// it) uses for iteration 0 of the Morgan expansion.
///
/// Only the *atom invariant* is affected — hashing (FNV-1a), environment
/// deduplication, and bit-folding are unchanged in both modes, so neither mode
/// is bit-compatible with RDKit's own Morgan fingerprint bit positions. See
/// [`crate::ecfp::EcfpInvariantMode::RdkitMorgan`] for what it does and does
/// not claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EcfpInvariantMode {
    /// chematic's original invariant: atomic number, degree (counts explicit
    /// H neighbors too), implicit H count, formal charge, ring membership,
    /// aromaticity. Unchanged from every prior release — this is what
    /// [`ecfp`]/[`ecfp4`]/[`ecfp6`]/[`ecfp_with_bitinfo`]/[`morgan_fp_counts`]
    /// have always computed.
    #[default]
    Chematic,
    /// RDKit's default `GetConnectivityInvariants` atom invariant, matching
    /// RDKit's own component set: atomic number, total degree (heavy
    /// neighbors + H, implicit or explicit — RDKit's `getTotalDegree()`),
    /// total H count as a *separate* component, formal charge (full `i8`
    /// range), ring membership, isotope mass delta from the element's
    /// average (IUPAC standard) atomic weight, truncated toward zero — no
    /// aromaticity component. Empirically verified against RDKit
    /// 2026.03.3's `GetConnectivityInvariants` (partition agreement, not raw
    /// hash-value agreement — the two implementations use different hash
    /// functions by construction). This is *atom-invariant* parity only, not
    /// a claim of RDKit fingerprint bit-compatibility — hence "RdkitMorgan",
    /// not "RdkitCompatible": that name is reserved for a possible future
    /// mode that also matches RDKit's hash, environment deduplication, and
    /// folding.
    RdkitMorgan,
}

/// RDKit's `getTotalDegree()`-equivalent: heavy-atom neighbors plus H count
/// (implicit or explicit — RDKit's connectivity invariant uses total degree
/// including Hs as one component, and total H count as a *separate* second
/// component; see `rdkit_total_h_count`). Explicit H atoms are already
/// counted once via `mol.neighbors`, so only `implicit_hcount` is added on
/// top (verified: `[H]C([H])([H])[H]`'s carbon gets the same invariant as
/// plain `C`'s — both total 4).
pub(crate) fn rdkit_total_degree(mol: &Molecule, idx: AtomIdx) -> u16 {
    let graph_degree = mol.neighbors(idx).count() as u16;
    let inferred_or_bracket_h = implicit_hcount(mol, idx) as u16;
    graph_degree + inferred_or_bracket_h
}

/// RDKit's total H count: implicit H plus any explicit H *atoms* in the
/// graph, combined into one count regardless of how they're spelled
/// (verified: a mixed implicit/explicit-H atom gets the same invariant as
/// its all-implicit or all-explicit spelling of the same total H count).
pub(crate) fn rdkit_total_h_count(mol: &Molecule, idx: AtomIdx) -> u8 {
    let explicit_h = mol
        .neighbors(idx)
        .filter(|&(nb, _)| mol.atom(nb).element.atomic_number() == 1)
        .count() as u8;
    implicit_hcount(mol, idx).saturating_add(explicit_h)
}

/// Exact `deltaMass` for a given `(atomic_number, mass_number)` isotope, via
/// binary search over [`crate::rdkit_isotope_delta_table::RDKIT_ISOTOPE_DELTA_TABLE`] --
/// a table generated directly from RDKit's own `PeriodicTable`
/// (`scripts/gen_rdkit_isotope_delta_table.py`), covering every isotope
/// RDKit's `GetMassForIsotope` recognizes for every element, not a
/// hand-picked subset. `None` means `(atomic_number, mass_number)` isn't a
/// real isotope RDKit itself recognizes either.
fn rdkit_isotope_delta_for(atomic_number: u8, mass_number: u16) -> Option<i16> {
    use crate::rdkit_isotope_delta_table::RDKIT_ISOTOPE_DELTA_TABLE;
    RDKIT_ISOTOPE_DELTA_TABLE
        .binary_search_by_key(&(atomic_number, mass_number), |&(z, a, _)| (z, a))
        .ok()
        .map(|i| RDKIT_ISOTOPE_DELTA_TABLE[i].2)
}

/// Isotope-mass delta, matching RDKit's actual `getConnectivityInvariants`
/// source (`ichi`/Morgan fingerprint invariant computation): the explicit
/// isotope's exact mass minus the element's average (standard) atomic
/// weight, truncated toward zero — **not** an integer mass-number
/// difference from the most common isotope, which was this function's
/// earlier (incorrect) definition and does not reproduce RDKit's actual
/// rounding behavior (e.g. it would treat carbon-13 as different from
/// unspecified/carbon-12, when RDKit's real mass-based truncation puts all
/// three at delta 0: `13.00335 - 12.011 = 0.992`, which truncates to `0`,
/// not `1`).
///
/// `atom.isotope == None` always gives delta 0, since `Atom::getMass()`
/// itself returns the average atomic weight for an unspecified isotope
/// (confirmed directly, not inferred: RDKit's `C`'s `GetMass()` is `12.011`,
/// not the monoisotopic `12.000`) — so both sides of the subtraction are the
/// same value. An explicit isotope RDKit itself doesn't recognize
/// (`rdkit_isotope_delta_for` returns `None`) also has defined RDKit
/// behavior, confirmed the same way: `Atom::getMass()` falls back to the
/// raw mass number itself (`[500CH4]`'s `GetMass()` is exactly `500.0`, not
/// 0 and not a nearby real isotope's mass), so the delta is
/// `mass_number - average_atomic_weight` — using
/// [`RDKIT_ATOMIC_WEIGHTS`](crate::rdkit_isotope_delta_table::RDKIT_ATOMIC_WEIGHTS)
/// for `average_atomic_weight` here, **not** chematic-core's
/// `Element::atomic_mass()` (a different, monoisotopic quantity — using it
/// would silently reproduce the exact carbon-11-shaped bug this table was
/// built to close, just for a different isotope).
pub(crate) fn rdkit_isotope_delta(mol: &Molecule, idx: AtomIdx) -> i32 {
    use crate::rdkit_isotope_delta_table::RDKIT_ATOMIC_WEIGHTS;

    let atom = mol.atom(idx);
    match atom.isotope {
        None => 0,
        Some(mass_number) => {
            match rdkit_isotope_delta_for(atom.element.atomic_number(), mass_number) {
                Some(delta) => delta as i32,
                None => {
                    let average = RDKIT_ATOMIC_WEIGHTS[atom.element.atomic_number() as usize];
                    (mass_number as f64 - average) as i32
                }
            }
        }
    }
}

/// Build the pre-chirality invariant byte sequence for `idx` under `mode`.
///
/// `Chematic` mode returns the exact same 6 bytes the pre-mode-aware
/// implementation always produced (byte order, byte count, and values
/// unchanged) — this is what makes [`ecfp`]/[`ecfp4`]/[`ecfp6`]/
/// [`ecfp_with_bitinfo`]/[`morgan_fp_counts`] fully non-regressing.
fn rdkit_connectivity_invariant_bytes(
    mol: &Molecule,
    idx: AtomIdx,
    ring_set: &chematic_perception::RingSet,
) -> SmallVec<[u8; 16]> {
    let atom = mol.atom(idx);
    let mut bytes: SmallVec<[u8; 16]> = SmallVec::new();
    bytes.push(atom.element.atomic_number());
    bytes.extend_from_slice(&rdkit_total_degree(mol, idx).to_le_bytes());
    bytes.push(rdkit_total_h_count(mol, idx));
    // Full i8 range, injective (unlike Chematic mode's `+8` clamp, which
    // collides charges outside roughly -8..+247) — a new public
    // "RdkitMorgan" invariant has no reason to inherit that collision.
    bytes.push(atom.charge as u8);
    bytes.push(ring_set.contains_atom(idx) as u8);
    bytes.extend_from_slice(&rdkit_isotope_delta(mol, idx).to_le_bytes());
    bytes
}

fn invariant_bytes(
    mol: &Molecule,
    idx: AtomIdx,
    ring_set: &chematic_perception::RingSet,
    mode: EcfpInvariantMode,
) -> SmallVec<[u8; 16]> {
    match mode {
        EcfpInvariantMode::Chematic => {
            let atom = mol.atom(idx);
            let charge_adjusted = (atom.charge as i16 + 8).clamp(0, 255) as u8;
            SmallVec::from_slice(&[
                atom.element.atomic_number(),
                mol.neighbors(idx).count().min(255) as u8,
                implicit_hcount(mol, idx),
                charge_adjusted,
                ring_set.contains_atom(idx) as u8,
                atom.aromatic as u8,
            ])
        }
        EcfpInvariantMode::RdkitMorgan => rdkit_connectivity_invariant_bytes(mol, idx, ring_set),
    }
}

/// Compute the FNV-1a atom identifier for iteration 0 of the Morgan algorithm.
///
/// See [`EcfpInvariantMode`] for what `mode` covers. When `use_chirality` is
/// true an extra chirality byte is appended; this preserves bit-compatibility
/// with the default (`use_chirality=false`) fingerprints.
pub(crate) fn initial_atom_id(
    mol: &Molecule,
    idx: AtomIdx,
    ring_set: &chematic_perception::RingSet,
    use_chirality: bool,
    mode: EcfpInvariantMode,
) -> u64 {
    let mut bytes = invariant_bytes(mol, idx, ring_set, mode);
    if use_chirality {
        use chematic_core::Chirality;
        // Reads the raw stored tag, not any permutation-corrected form -- same
        // reordering-sensitivity tetrahedral chirality already has here, not a new
        // inconsistency introduced by adding SquarePlanar.
        let chirality_byte = match mol.atom(idx).chirality {
            Chirality::None => 0u8,
            Chirality::CounterClockwise => 1u8,
            Chirality::Clockwise => 2u8,
            Chirality::SquarePlanar(p) => 3 + p as u8,
        };
        bytes.push(chirality_byte);
    }
    fnv1a(&bytes)
}

/// Configuration for ECFP computation.
#[derive(Debug, Clone)]
pub struct EcfpConfig {
    /// Number of iterations (radius). ECFP4 = 2, ECFP6 = 3.
    pub radius: u32,
    /// Output bitvector size (default 2048).
    pub nbits: usize,
    /// When `true`, include tetrahedral chirality in the initial atom hash so
    /// that R and S enantiomers produce different fingerprints.
    ///
    /// Defaults to `false` (chirality ignored, matching RDKit's `useChirality=False`).
    pub use_chirality: bool,
    /// When `true`, each hash sets two bit positions (using single and double-folded hash)
    /// to reduce bitvector collisions. This reduces collision probability but changes
    /// fingerprint values — **not backwards-compatible** with stored fingerprints.
    ///
    /// Defaults to `false` (single-bit folding, current behavior preserved).
    pub use_double_fold: bool,
}

impl Default for EcfpConfig {
    fn default() -> Self {
        Self {
            radius: 2,
            nbits: 2048,
            use_chirality: false,
            use_double_fold: false,
        }
    }
}

/// Map a `BondOrder` to the integer code used in the ECFP hash.
///
/// - Single / Up / Down → 1
/// - Double            → 2
/// - Triple            → 3
/// - Aromatic          → 4
/// - Quadruple         → 5  (not in standard ECFP; assigned a distinct value)
#[inline]
pub(crate) fn bond_type_int(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
        BondOrder::Quadruple => 5,
        BondOrder::Zero => 0,
        BondOrder::QueryAny => 6,
        BondOrder::QuerySingleOrDouble => 7,
        BondOrder::QuerySingleOrAromatic => 8,
        BondOrder::QueryDoubleOrAromatic => 9,
    }
}

/// Compute an ECFP fingerprint for `mol` using the given configuration.
///
/// # Algorithm overview
/// 1. Compute initial atom identifiers from atomic properties.
/// 2. Iteratively expand each identifier by incorporating neighbour identifiers
///    (with their bond types) for `config.radius` rounds.
/// 3. After each iteration (including iteration 0), map every identifier to a
///    bit in the output bitvector.
///
/// Maximum supported radius for `ecfp`.  Matches the cap in `morgan_fp_counts`.
/// Beyond this, `r as u8` would silently truncate, producing hash collisions.
pub const MAX_ECFP_RADIUS: u32 = 20;

pub fn ecfp(mol: &Molecule, config: &EcfpConfig) -> BitVec2048 {
    ecfp_with_invariant_mode(mol, config, EcfpInvariantMode::Chematic)
}

/// Like [`ecfp`], but with an explicit [`EcfpInvariantMode`] choice for the
/// iteration-0 atom invariant. [`ecfp`] is exactly
/// `ecfp_with_invariant_mode(mol, config, EcfpInvariantMode::Chematic)`.
pub fn ecfp_with_invariant_mode(
    mol: &Molecule,
    config: &EcfpConfig,
    mode: EcfpInvariantMode,
) -> BitVec2048 {
    let n = mol.atom_count();
    let nbits = config.nbits;
    // Cap radius to prevent `r as u8` truncation at r > 255 (hash collision bug).
    let config = &EcfpConfig {
        radius: config.radius.min(MAX_ECFP_RADIUS),
        ..*config
    };
    let mut fp = BitVec2048::new();

    if n == 0 {
        return fp;
    }

    let ring_set = find_sssr(mol);

    // Step 1: initial atom identifiers (iteration 0).
    let mut ids: Vec<u64> = Vec::with_capacity(n);
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let id = initial_atom_id(mol, idx, &ring_set, config.use_chirality, mode);
        fp.set((id % nbits as u64) as usize);
        if config.use_double_fold {
            fp.set(((id >> 11) % nbits as u64) as usize);
        }
        ids.push(id);
    }

    // Step 2: iterative expansion.
    let mut new_ids: Vec<u64> = vec![0u64; n];
    for r in 1..=config.radius {
        for (i, slot) in new_ids.iter_mut().enumerate() {
            let new_id = expand_atom_id(mol, i, r, &ids);
            *slot = new_id;
            fp.set((new_id % nbits as u64) as usize);
            if config.use_double_fold {
                fp.set(((new_id >> 11) % nbits as u64) as usize);
            }
        }
        core::mem::swap(&mut ids, &mut new_ids);
    }

    fp
}

/// Like [`ecfp`] but also returns, for each set bit, the list of
/// `(atom_idx, radius)` environments that produced it — the data behind
/// RDKit's `bitInfo` map.
///
/// The fingerprint bits are identical to [`ecfp`] with the same config
/// (same hash, same fold), so the recorded environments are the true origin
/// of each bit. Bit positions still differ from RDKit (FNV-1a vs MurmurHash),
/// so this is shape-compatible and internally consistent, not bit-identical.
pub fn ecfp_with_bitinfo(
    mol: &Molecule,
    config: &EcfpConfig,
) -> (BitVec2048, FxHashMap<usize, Vec<(u32, u32)>>) {
    ecfp_with_bitinfo_and_mode(mol, config, EcfpInvariantMode::Chematic)
}

/// Like [`ecfp_with_bitinfo`], but with an explicit [`EcfpInvariantMode`]
/// choice for the iteration-0 atom invariant — shares the same invariant
/// computation as [`ecfp_with_invariant_mode`], so the two never diverge.
pub fn ecfp_with_bitinfo_and_mode(
    mol: &Molecule,
    config: &EcfpConfig,
    mode: EcfpInvariantMode,
) -> (BitVec2048, FxHashMap<usize, Vec<(u32, u32)>>) {
    let n = mol.atom_count();
    let nbits = config.nbits;
    let config = &EcfpConfig {
        radius: config.radius.min(MAX_ECFP_RADIUS),
        ..*config
    };
    let mut fp = BitVec2048::new();
    let mut info: FxHashMap<usize, Vec<(u32, u32)>> = FxHashMap::default();

    if n == 0 {
        return (fp, info);
    }

    let ring_set = find_sssr(mol);

    // Iteration 0: initial atom identifiers.
    let mut ids: Vec<u64> = Vec::with_capacity(n);
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let id = initial_atom_id(mol, idx, &ring_set, config.use_chirality, mode);
        record_bit(
            &mut fp,
            &mut info,
            id,
            i as u32,
            0,
            nbits,
            config.use_double_fold,
        );
        ids.push(id);
    }

    // Iterations 1..=radius: expansion.
    let mut new_ids: Vec<u64> = vec![0u64; n];
    for r in 1..=config.radius {
        for (i, slot) in new_ids.iter_mut().enumerate() {
            let new_id = expand_atom_id(mol, i, r, &ids);
            *slot = new_id;
            record_bit(
                &mut fp,
                &mut info,
                new_id,
                i as u32,
                r,
                nbits,
                config.use_double_fold,
            );
        }
        core::mem::swap(&mut ids, &mut new_ids);
    }

    (fp, info)
}

/// Set the bit(s) for hash `id` and record the `(atom, radius)` origin.
///
/// Shared by [`ecfp_with_bitinfo`] and [`crate::fcfp::fcfp_with_bitinfo`] — the
/// bit-recording step is identical regardless of how `id` was derived.
pub(crate) fn record_bit(
    fp: &mut BitVec2048,
    info: &mut FxHashMap<usize, Vec<(u32, u32)>>,
    id: u64,
    atom: u32,
    radius: u32,
    nbits: usize,
    use_double_fold: bool,
) {
    let bit = (id % nbits as u64) as usize;
    fp.set(bit);
    info.entry(bit).or_default().push((atom, radius));
    if use_double_fold {
        let bit2 = ((id >> 11) % nbits as u64) as usize;
        fp.set(bit2);
        info.entry(bit2).or_default().push((atom, radius));
    }
}

/// Count-based Morgan fingerprint: returns a map of `hash → count` for all
/// atom environments up to `radius` iterations.
///
/// Each (atom, iteration) pair contributes its hash to the map.  Unlike the
/// default RDKit behavior, redundant (duplicate) environments are **not**
/// suppressed — every atom contributes at every iteration level.
///
/// This corresponds to `GetMorganFingerprint(mol, radius,
/// useFeatures=False, includeRedundantEnvironments=True)` in RDKit.
pub fn morgan_fp_counts(mol: &Molecule, radius: u32) -> FxHashMap<u64, u32> {
    const MAX_RADIUS: u32 = 20;
    let radius = radius.min(MAX_RADIUS);

    let n = mol.atom_count();
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();

    if n == 0 {
        return counts;
    }

    let ring_set = find_sssr(mol);

    // Radius-0: initial atom identifiers.
    let mut ids: Vec<u64> = (0..n)
        .map(|i| {
            initial_atom_id(
                mol,
                AtomIdx(i as u32),
                &ring_set,
                false,
                EcfpInvariantMode::Chematic,
            )
        })
        .collect();

    for &id in &ids {
        *counts.entry(id).or_insert(0) += 1;
    }

    // Radius 1..=radius: iterative expansion (same hash scheme as ecfp).
    let mut new_ids = vec![0u64; n];
    for r in 1..=radius {
        for (i, slot) in new_ids.iter_mut().enumerate() {
            let new_id = expand_atom_id(mol, i, r, &ids);
            *slot = new_id;
            *counts.entry(new_id).or_insert(0) += 1;
        }
        core::mem::swap(&mut ids, &mut new_ids);
    }

    counts
}

/// ECFP4 fingerprint (radius = 2, 2048 bits).
pub fn ecfp4(mol: &Molecule) -> BitVec2048 {
    ecfp(mol, &EcfpConfig::default())
}

/// ECFP6 fingerprint (radius = 3, 2048 bits).
pub fn ecfp6(mol: &Molecule) -> BitVec2048 {
    ecfp(
        mol,
        &EcfpConfig {
            radius: 3,
            ..EcfpConfig::default()
        },
    )
}

/// Tanimoto similarity between two molecules using ECFP4.
pub fn tanimoto_ecfp4(a: &Molecule, b: &Molecule) -> f64 {
    ecfp4(a).tanimoto(&ecfp4(b))
}

/// ECFP4 fingerprint (radius = 2, 2048 bits) using RDKit's default Morgan
/// atom invariant instead of chematic's own. See [`EcfpInvariantMode::RdkitMorgan`].
pub fn ecfp4_rdkit_invariants(mol: &Molecule) -> BitVec2048 {
    ecfp_with_invariant_mode(mol, &EcfpConfig::default(), EcfpInvariantMode::RdkitMorgan)
}

/// ECFP6 fingerprint (radius = 3, 2048 bits) using RDKit's default Morgan
/// atom invariant instead of chematic's own. See [`EcfpInvariantMode::RdkitMorgan`].
pub fn ecfp6_rdkit_invariants(mol: &Molecule) -> BitVec2048 {
    ecfp_with_invariant_mode(
        mol,
        &EcfpConfig {
            radius: 3,
            ..EcfpConfig::default()
        },
        EcfpInvariantMode::RdkitMorgan,
    )
}

/// **Experimental.** ECFP4 fingerprint (radius = 2, 2048 bits) additionally
/// applying RDKit's redundant-environment suppression on top of
/// [`EcfpInvariantMode::RdkitMorgan`] — see
/// [`ecfp_with_bitinfo_rdkit_environment_experimental`] for what this mode
/// does and does not claim to match.
pub fn ecfp4_rdkit_environment_experimental(mol: &Molecule) -> BitVec2048 {
    ecfp_with_bitinfo_rdkit_environment_experimental(mol, &EcfpConfig::default()).0
}

/// **Experimental.** ECFP6 fingerprint (radius = 3, 2048 bits) additionally
/// applying RDKit's redundant-environment suppression. See
/// [`ecfp_with_bitinfo_rdkit_environment_experimental`].
pub fn ecfp6_rdkit_environment_experimental(mol: &Molecule) -> BitVec2048 {
    ecfp_with_bitinfo_rdkit_environment_experimental(
        mol,
        &EcfpConfig {
            radius: 3,
            ..EcfpConfig::default()
        },
    )
    .0
}

/// **Experimental.** Like [`ecfp_with_bitinfo_and_mode`], but also applies
/// RDKit's redundant-environment suppression (an atom is not re-emitted at a
/// higher radius once its cumulative bond-environment duplicates one already
/// emitted) on top of [`EcfpInvariantMode::RdkitMorgan`] — see
/// `crate::morgan_environment` for the algorithm, verified directly against
/// RDKit's own Morgan generator source.
///
/// **Not** bit-compatible with RDKit's fingerprint: FNV-1a still doesn't
/// match RDKit's hash, so representative selection when multiple atoms
/// produce an identical bond-environment in the same round (RDKit breaks
/// this tie by hash value) is not guaranteed to match RDKit's specific pick
/// when the tied atoms are chemically non-equivalent. Emitted-`(atom,
/// radius)`-pair-set agreement with RDKit is high (measured against RDKit's
/// own default Morgan generator; see `scripts/ecfp_rdkit_suppression_parity.py`),
/// not the raw bit values.
pub fn ecfp_with_bitinfo_rdkit_environment_experimental(
    mol: &Molecule,
    config: &EcfpConfig,
) -> (BitVec2048, FxHashMap<usize, Vec<(u32, u32)>>) {
    crate::morgan_environment::ecfp_environments(
        mol,
        config,
        EcfpInvariantMode::RdkitMorgan,
        crate::morgan_environment::EnvironmentEmissionMode::SuppressRdkitRedundant,
    )
}

/// Raw (unfolded) iteration-0 atom invariant for every atom in `mol`, under
/// `mode`, with `use_chirality=false`. One entry per atom, in `AtomIdx`
/// order. Exposed for atom-invariant-partition comparison against an
/// external oracle (e.g. RDKit's `GetConnectivityInvariants`) — the raw `u64`
/// values themselves are not meant to be compared directly across
/// implementations (different hash functions), only whether two atoms get
/// the *same* value as each other (the equivalence partition).
pub fn atom_invariants(mol: &Molecule, mode: EcfpInvariantMode) -> Vec<u64> {
    let ring_set = find_sssr(mol);
    (0..mol.atom_count())
        .map(|i| initial_atom_id(mol, AtomIdx(i as u32), &ring_set, false, mode))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn benzene() -> Molecule {
        parse("c1ccccc1").unwrap()
    }

    fn ethane() -> Molecule {
        parse("CC").unwrap()
    }

    fn toluene() -> Molecule {
        parse("Cc1ccccc1").unwrap()
    }

    fn aspirin() -> Molecule {
        // Acetylsalicylic acid
        parse("CC(=O)Oc1ccccc1C(=O)O").unwrap()
    }

    fn methane() -> Molecule {
        parse("C").unwrap()
    }

    fn water() -> Molecule {
        parse("O").unwrap()
    }

    #[test]
    fn benzene_ecfp4_nonzero() {
        let fp = ecfp4(&benzene());
        assert!(fp.popcount() > 0, "benzene ECFP4 must be non-zero");
    }

    #[test]
    fn benzene_ecfp4_deterministic() {
        let fp1 = ecfp4(&benzene());
        let fp2 = ecfp4(&benzene());
        assert_eq!(fp1, fp2, "ECFP4 must be deterministic");
    }

    #[test]
    fn ethane_vs_benzene_tanimoto_lt1() {
        let t = tanimoto_ecfp4(&ethane(), &benzene());
        assert!(t < 1.0, "ethane and benzene must differ (tanimoto={t})");
    }

    #[test]
    fn benzene_vs_benzene_tanimoto_eq1() {
        let t = tanimoto_ecfp4(&benzene(), &benzene());
        assert_eq!(t, 1.0, "identical molecules must have tanimoto == 1.0");
    }

    #[test]
    fn bitinfo_fp_matches_ecfp4() {
        // The fp returned alongside bitInfo must equal the plain ecfp4 output.
        let mol = aspirin();
        let (fp, _info) = ecfp_with_bitinfo(&mol, &EcfpConfig::default());
        assert_eq!(fp, ecfp4(&mol), "bitInfo fp must match ecfp4");
    }

    #[test]
    fn bitinfo_keys_are_set_bits_and_valid() {
        let mol = aspirin();
        let n = mol.atom_count() as u32;
        let (fp, info) = ecfp_with_bitinfo(&mol, &EcfpConfig::default());
        assert!(!info.is_empty());
        for (&bit, envs) in &info {
            assert!(fp.get(bit), "every bitInfo key must be a set bit");
            for &(atom, radius) in envs {
                assert!(atom < n, "atom idx in range");
                assert!(radius <= 2, "radius within ECFP4 (<=2)");
            }
        }
    }

    #[test]
    fn benzene_vs_toluene_tanimoto_between() {
        let t = tanimoto_ecfp4(&benzene(), &toluene());
        assert!(t > 0.0, "benzene and toluene share bits (tanimoto={t})");
        assert!(
            t < 1.0,
            "benzene and toluene are not identical (tanimoto={t})"
        );
    }

    #[test]
    fn aspirin_ecfp4_many_bits() {
        let fp = ecfp4(&aspirin());
        assert!(
            fp.popcount() > 5,
            "aspirin ECFP4 must have more than 5 bits set (got {})",
            fp.popcount()
        );
    }

    #[test]
    fn ecfp6_vs_ecfp4_benzene_differ() {
        let fp4 = ecfp4(&benzene());
        let fp6 = ecfp6(&benzene());
        // Larger radius explores more environment — the bit counts should differ
        // because radius-3 adds new hash values not present at radius-2.
        assert_ne!(
            fp4.popcount(),
            fp6.popcount(),
            "ECFP6 and ECFP4 should produce different bit counts for benzene"
        );
    }

    #[test]
    fn methane_ecfp4_nonzero() {
        let fp = ecfp4(&methane());
        assert!(fp.popcount() > 0, "methane ECFP4 must be non-zero");
    }

    #[test]
    fn water_ecfp4_nonzero() {
        let fp = ecfp4(&water());
        assert!(fp.popcount() > 0, "water ECFP4 must be non-zero");
    }

    #[test]
    fn tanimoto_ecfp4_benzene_self_is_one() {
        let t = tanimoto_ecfp4(&benzene(), &benzene());
        assert_eq!(t, 1.0, "tanimoto_ecfp4 of identical molecules must be 1.0");
    }

    #[test]
    fn tanimoto_ecfp4_methane_vs_benzene_lt_half() {
        let t = tanimoto_ecfp4(&methane(), &benzene());
        assert!(
            t < 0.5,
            "methane and benzene should be very dissimilar (tanimoto={t})"
        );
    }

    // ── morgan_fp_counts ─────────────────────────────────────────────────────

    #[test]
    fn morgan_counts_radius0_atom_count() {
        // At radius 0, one hash per atom.
        let m = benzene();
        let counts = morgan_fp_counts(&m, 0);
        let total: u32 = counts.values().sum();
        assert_eq!(
            total,
            m.atom_count() as u32,
            "radius-0 total count == atom_count"
        );
    }

    #[test]
    fn morgan_counts_radius2_total_grows() {
        // Each additional radius adds one hash per atom → total = n * (radius+1).
        let m = methane();
        let n = m.atom_count() as u32;
        let r = 2u32;
        let counts = morgan_fp_counts(&m, r);
        let total: u32 = counts.values().sum();
        assert_eq!(
            total,
            n * (r + 1),
            "methane total = atom_count * (radius+1)"
        );
    }

    #[test]
    fn morgan_counts_benzene_symmetry() {
        // All 6 benzene C atoms are equivalent → radius-0 yields 1 unique hash.
        let m = benzene();
        let counts = morgan_fp_counts(&m, 0);
        assert_eq!(
            counts.len(),
            1,
            "benzene has one unique radius-0 environment"
        );
        assert_eq!(
            *counts.values().next().unwrap(),
            6,
            "that environment appears 6 times"
        );
    }

    #[test]
    fn morgan_counts_empty_mol_is_empty() {
        use chematic_core::MoleculeBuilder;
        let m = MoleculeBuilder::new().build();
        let counts = morgan_fp_counts(&m, 2);
        assert!(counts.is_empty(), "empty molecule yields empty count map");
    }

    #[test]
    fn morgan_counts_deterministic() {
        let m = aspirin();
        let c1 = morgan_fp_counts(&m, 2);
        let c2 = morgan_fp_counts(&m, 2);
        assert_eq!(c1, c2, "morgan_fp_counts must be deterministic");
    }

    #[test]
    fn morgan_counts_consistent_with_ecfp_bits() {
        // Every hash in the count map should be reachable from the ecfp bit set
        // (after folding to 2048 bits).  This checks the same hash scheme.
        let m = toluene();
        let fp = ecfp(
            &m,
            &EcfpConfig {
                radius: 2,
                nbits: 2048,
                use_chirality: false,
                use_double_fold: false,
            },
        );
        let counts = morgan_fp_counts(&m, 2);
        for &hash in counts.keys() {
            let bit = (hash % 2048) as usize;
            assert!(
                fp.get(bit),
                "bit {bit} from count map not set in ECFP bitvec"
            );
        }
    }

    // -- Chirality tests ------------------------------------------------------

    #[test]
    fn ecfp4_ignores_chirality_by_default() {
        // L-alanine and D-alanine should produce the same ECFP4 when
        // use_chirality=false (default), since chirality is not in the hash.
        let l_ala = parse("N[C@@H](C)C(=O)O").unwrap();
        let d_ala = parse("N[C@H](C)C(=O)O").unwrap();
        let fp_l = ecfp4(&l_ala);
        let fp_d = ecfp4(&d_ala);
        assert_eq!(
            fp_l, fp_d,
            "L/D-alanine ECFP4 should be identical when use_chirality=false"
        );
    }

    #[test]
    fn ecfp4_distinguishes_enantiomers_with_chirality() {
        // With use_chirality=true, L-alanine and D-alanine must have different FPs.
        let l_ala = parse("N[C@@H](C)C(=O)O").unwrap();
        let d_ala = parse("N[C@H](C)C(=O)O").unwrap();
        let config = EcfpConfig {
            radius: 2,
            nbits: 2048,
            use_chirality: true,
            use_double_fold: false,
        };
        let fp_l = ecfp(&l_ala, &config);
        let fp_d = ecfp(&d_ala, &config);
        assert_ne!(
            fp_l, fp_d,
            "L/D-alanine ECFP4 must differ when use_chirality=true"
        );
        // Tanimoto < 1.0 confirms they are not identical.
        assert!(
            fp_l.tanimoto(&fp_d) < 1.0,
            "Tanimoto of L/D-alanine must be < 1.0 with use_chirality"
        );
    }

    #[test]
    fn ecfp4_non_chiral_generates_with_chirality_flag() {
        // Non-chiral molecules like benzene should generate valid ECFP with use_chirality flag.
        // This test verifies that use_chirality=true doesn't break on achiral molecules.
        let mol = parse("c1ccccc1").unwrap(); // Benzene — no stereo centers.
        let config = EcfpConfig {
            radius: 2,
            nbits: 2048,
            use_chirality: true,
            use_double_fold: false,
        };
        let fp = ecfp(&mol, &config);
        assert!(
            fp.popcount() > 0,
            "Benzene should generate non-empty ECFP4 with use_chirality=true"
        );
    }

    // -----------------------------------------------------------------------
    // Implicit vs. explicit hydrogen representation (CDK issue #1084 pattern)
    // -----------------------------------------------------------------------
    //
    // chematic's molecule model stores only heavy atoms; implicit H counts are
    // computed on demand via `implicit_hcount()`.  When a SMILES contains
    // explicit H atoms (e.g. `[H]O` or `[OH2]`), those atoms ARE stored in
    // the molecular graph and WILL change the ECFP invariant for the heavy
    // atom they're bonded to (its degree increases and its implicit_hcount
    // decreases).  This is the expected behaviour for a graph-based fingerprint
    // and mirrors CDK / RDKit behaviour.  These tests document it explicitly.

    #[test]
    fn ecfp4_implicit_h_water_vs_no_atoms() {
        // "O" has 1 heavy atom (O) with 2 implicit H.
        // "[OH2]" is the same molecule: O with explicit H count = 2 but still
        // only 1 heavy atom (implicit H also = 0 because OH2 bracket sets it).
        // Both should produce the same fingerprint.
        let implicit = parse("O").unwrap();
        let bracketed = parse("[OH2]").unwrap();
        assert_eq!(
            ecfp4(&implicit),
            ecfp4(&bracketed),
            "[OH2] and O should give the same ECFP4 (same heavy-atom graph)"
        );
    }

    #[test]
    fn ecfp4_explicit_h_atom_changes_fingerprint() {
        // "[H]O[H]" parses as 3 atoms: H-O-H.  The O atom now has degree=2
        // and implicit_hcount=0, so its Morgan invariant differs from the
        // single-atom "O" (degree=0, implicit_hcount=2).  The fingerprint
        // must differ — this documents the expected behaviour when explicit
        // H atoms are present in the molecular graph.
        let implicit = parse("O").unwrap();
        let explicit_h = parse("[H]O[H]").unwrap();
        assert_ne!(
            ecfp4(&implicit),
            ecfp4(&explicit_h),
            "explicit H atoms in the graph change the ECFP4"
        );
    }

    #[test]
    fn ecfp4_implicit_vs_explicit_h_in_organic_molecule() {
        // Methanol "CO" vs "C([H])([H])([H])O" — the second form has 3
        // explicit H atoms on C, which changes C's degree and implicit_hcount.
        // These are different molecular graphs → different ECFP4.
        let implicit = parse("CO").unwrap();
        let explicit_h = parse("C([H])([H])([H])O").unwrap();
        assert_ne!(
            ecfp4(&implicit),
            ecfp4(&explicit_h),
            "methanol with explicit H atoms has a different heavy-atom neighbourhood"
        );
    }

    // ── EcfpInvariantMode::RdkitMorgan edge-case fixtures ──────────────────
    //
    // Every equivalence/non-equivalence asserted here was independently
    // verified against RDKit 2026.03.3's `GetConnectivityInvariants` (see
    // scripts/ecfp_rdkit_invariant_parity.py and
    // scripts/ecfp_rdkit_edge_fixtures.csv) -- 100% partition agreement on
    // this fixture set (including cross-atom isotope/charge stress cases)
    // and the full 5,000-molecule ChEMBL corpus, no residual.

    fn rdkit_inv(mol: &Molecule) -> Vec<u64> {
        atom_invariants(mol, EcfpInvariantMode::RdkitMorgan)
    }

    #[test]
    fn rdkit_mode_ignores_explicit_h_in_degree() {
        // Unlike Chematic mode (see `ecfp4_implicit_vs_explicit_h_in_organic_molecule`
        // above), RDKit's total-degree component counts H the same way
        // regardless of spelling -- verified: methane's C gets the same
        // RDKit invariant whether all 4 H are implicit or all 4 are
        // explicit atoms.
        let implicit = parse("C").unwrap();
        let explicit_h = parse("[H]C([H])([H])[H]").unwrap();
        assert_eq!(
            rdkit_inv(&implicit)[0],
            rdkit_inv(&explicit_h)[1],
            "RdkitMorgan mode must not distinguish implicit vs explicit H representation"
        );
    }

    #[test]
    fn rdkit_mode_ignores_aromaticity() {
        // The whole point of this mode: benzene's aromatic ring carbon and
        // cyclohexadiene's non-aromatic sp2 ring carbon (same degree/H/ring
        // membership, only aromaticity differs) must be invariant-equivalent.
        let benzene = parse("c1ccccc1").unwrap();
        let kekule_diene = parse("C1=CC=CC=C1").unwrap();
        assert_eq!(
            rdkit_inv(&benzene)[0],
            rdkit_inv(&kekule_diene)[0],
            "RdkitMorgan mode must not distinguish aromatic from Kekule sp2 ring atoms"
        );
    }

    #[test]
    fn rdkit_mode_distinguishes_ring_from_acyclic() {
        let cyclohexane = parse("C1CCCCC1").unwrap();
        let hexane = parse("CCCCCC").unwrap();
        assert_ne!(
            rdkit_inv(&cyclohexane)[0],
            rdkit_inv(&hexane)[1], // hexane's atom 1 is a chain CH2, same degree/H as the ring atom
            "RdkitMorgan mode must still distinguish ring membership"
        );
    }

    #[test]
    fn rdkit_mode_distinguishes_charge() {
        let pyridine_n = parse("c1ccncc1").unwrap();
        let pyridinium_n = parse("c1cc[nH+]cc1").unwrap();
        // Ring N in each case is atom index 3.
        assert_ne!(
            rdkit_inv(&pyridine_n)[3],
            rdkit_inv(&pyridinium_n)[3],
            "RdkitMorgan mode must still distinguish formal charge (and the H it brings)"
        );
    }

    #[test]
    fn rdkit_mode_charge_stress_full_i8_range_no_collision() {
        // The Chematic-mode `+8` clamp would collide e.g. charge -10 and
        // charge -11 (both clamp to the same byte). RdkitMorgan mode must
        // not inherit that: every formal charge in a wide, deliberately
        // out-of-clamp-range sweep must map to a distinct invariant (same
        // degree/H/ring/isotope otherwise, only charge varies).
        let mut invs = Vec::new();
        for charge in [-20i8, -10, -9, -8, -1, 0, 1, 8, 9, 10, 20, 100, 127] {
            let mut mol = parse("C").unwrap();
            let idx = AtomIdx(0);
            let mut atom = mol.atom(idx).clone();
            atom.charge = charge;
            let mut builder = chematic_core::MoleculeBuilder::new();
            builder.add_atom(atom);
            mol = builder.build();
            invs.push(rdkit_inv(&mol)[0]);
        }
        let unique: std::collections::HashSet<_> = invs.iter().collect();
        assert_eq!(
            unique.len(),
            invs.len(),
            "every charge in the sweep must produce a distinct invariant, got {invs:?}"
        );
    }

    #[test]
    fn rdkit_mode_isotope_unspecified_equals_most_common() {
        // Verified against RDKit: an atom with no isotope specified and one
        // with its isotope explicitly set to the element's most common value
        // are invariant-equivalent (both are the delta=0 baseline).
        let unspecified = parse("C").unwrap();
        let explicit_12 = parse("[12CH4]").unwrap();
        assert_eq!(
            rdkit_inv(&unspecified)[0],
            rdkit_inv(&explicit_12)[0],
            "unspecified isotope and the most-common isotope must both be delta=0"
        );
    }

    #[test]
    fn rdkit_mode_isotope_carbon_12_13_same_class() {
        // Fixed: RDKit's real deltaMass truncates `13.00335 - 12.011 =
        // 0.992` toward zero to `0`, the same as carbon-12's `12.00000 -
        // 12.011 = -0.011 -> 0` -- verified against RDKit, not assumed.
        // Cross-atom (not single-atom) so the partition comparison is
        // actually exercised: a 1-heavy-atom molecule's partition is
        // trivially "1 class" regardless of the invariant.
        let mol = parse("[12CH3][13CH3]").unwrap();
        let inv = rdkit_inv(&mol);
        assert_eq!(
            inv[0], inv[1],
            "carbon-12 and carbon-13 must be the same invariant class, matching RDKit"
        );
    }

    #[test]
    fn rdkit_mode_isotope_carbon_11_12_same_class() {
        // The counterexample that caught the previous partial-table
        // fallback: carbon-11 (11.0114336, not in the original hand-picked
        // 12/13/14 table) truncated via the `mass_number as f64`
        // approximation to `(11.0 - 12.011) as i32 = -1`, disagreeing with
        // RDKit's real `(11.0114336 - 12.011) as i32 = 0` -- same as
        // carbon-12. Now backed by the full generated table
        // (RDKIT_ISOTOPE_DELTA_TABLE), not an approximation.
        let mol = parse("[11CH3][12CH3]").unwrap();
        let inv = rdkit_inv(&mol);
        assert_eq!(
            inv[0], inv[1],
            "carbon-11 and carbon-12 must be the same invariant class, matching RDKit"
        );
    }

    #[test]
    fn rdkit_mode_isotope_carbon_10_11_different_class() {
        // `10.0168532 - 12.011 = -1.994` truncates to `-1`, different from
        // carbon-11's `0`.
        let mol = parse("[10CH3][11CH3]").unwrap();
        let inv = rdkit_inv(&mol);
        assert_ne!(
            inv[0], inv[1],
            "carbon-10 and carbon-11 must be different invariant classes, matching RDKit"
        );
    }

    #[test]
    fn rdkit_mode_isotope_carbon_10_12_different_class() {
        let mol = parse("[10CH3][12CH3]").unwrap();
        let inv = rdkit_inv(&mol);
        assert_ne!(
            inv[0], inv[1],
            "carbon-10 and carbon-12 must be different invariant classes, matching RDKit"
        );
    }

    #[test]
    fn rdkit_isotope_delta_table_is_sorted_for_binary_search() {
        // rdkit_isotope_delta_for's binary_search_by_key requires the table
        // to be sorted by (atomic_number, mass_number) -- a corrupted or
        // wrongly-regenerated table would silently give wrong lookups
        // rather than an error, so check the precondition explicitly rather
        // than trusting the generator forever.
        let table = crate::rdkit_isotope_delta_table::RDKIT_ISOTOPE_DELTA_TABLE;
        for w in table.windows(2) {
            let (z0, a0, _) = w[0];
            let (z1, a1, _) = w[1];
            assert!(
                (z0, a0) < (z1, a1),
                "table must be strictly sorted by (atomic_number, mass_number): \
                 ({z0}, {a0}) is not before ({z1}, {a1})"
            );
        }
    }

    #[test]
    fn rdkit_isotope_delta_table_exhaustive_lookup_round_trip() {
        // Every one of the table's 3,111 (atomic_number, mass_number,
        // delta) rows -- generated directly from RDKit 2026.03.3's
        // PeriodicTable, covering every isotope RDKit itself recognizes for
        // every element -- must round-trip through rdkit_isotope_delta_for
        // exactly. This is the actual "full RDKit-supported isotope delta
        // table: mismatch = 0" gate: not a hand-picked sample, all 3,111
        // rows, checked against chematic's own lookup path (the table's
        // *values* are RDKit's by construction; this exhaustively verifies
        // the Rust binary-search lookup never returns the wrong one).
        let table = crate::rdkit_isotope_delta_table::RDKIT_ISOTOPE_DELTA_TABLE;
        let mut mismatches = 0usize;
        for &(z, a, expected_delta) in table.iter() {
            let got = rdkit_isotope_delta_for(z, a);
            if got != Some(expected_delta) {
                mismatches += 1;
                if mismatches <= 5 {
                    eprintln!("MISMATCH: Z={z} A={a} expected={expected_delta} got={got:?}");
                }
            }
        }
        assert_eq!(
            mismatches,
            0,
            "{mismatches}/{} table entries did not round-trip through rdkit_isotope_delta_for",
            table.len()
        );
    }

    /// Parse a single-atom SMILES and return atom 0's raw RdkitMorgan
    /// isotope delta directly (not just whether it matches another atom's
    /// class) -- needed to catch a uniform off-by-N shift that a
    /// partition-only comparison can't distinguish from "correct".
    fn isotope_delta_for_smiles(smi: &str) -> i32 {
        let mol = parse(smi).unwrap();
        rdkit_isotope_delta(&mol, AtomIdx(0))
    }

    #[test]
    fn rdkit_mode_isotope_unrecognized_uses_average_atomic_weight_fallback() {
        // Mass numbers 499/500 aren't real carbon isotopes RDKit's
        // PeriodicTable recognizes (confirmed: GetMassForIsotope(6, 500) ==
        // 0.0), so they fall back to RDKit's own defined behavior for an
        // unrecognized explicit isotope: Atom::getMass() returns the raw
        // mass number itself (confirmed directly via RDKit's
        // Atom.GetMass(), not inferred -- [500CH4]'s GetMass() is exactly
        // 500.0). deltaMass = mass_number - average_atomic_weight,
        // truncated toward zero: `500.0 - 12.011 = 487.989 -> 487`,
        // `499.0 - 12.011 = 486.989 -> 486`. Using chematic-core's
        // monoisotopic Element::atomic_mass() (12.000) here instead would
        // give 488/487 -- silently wrong by exactly 1, a uniform shift a
        // partition-only comparison (same molecule's atoms all shifted
        // together) would never catch.
        assert_eq!(isotope_delta_for_smiles("[500CH4]"), 487);
        assert_eq!(isotope_delta_for_smiles("[499CH4]"), 486);
    }

    #[test]
    fn rdkit_mode_isotope_recognized_raw_delta_values() {
        // Exact raw deltas (not just partition membership) for a few
        // in-table isotopes, cross-checked against
        // RDKIT_ISOTOPE_DELTA_TABLE directly.
        assert_eq!(isotope_delta_for_smiles("[12CH4]"), 0);
        assert_eq!(isotope_delta_for_smiles("[13CH4]"), 0);
        assert_eq!(isotope_delta_for_smiles("[14CH4]"), 1);
        assert_eq!(isotope_delta_for_smiles("[11CH4]"), 0);
        assert_eq!(isotope_delta_for_smiles("[10CH4]"), -1);
    }

    #[test]
    fn rdkit_mode_isotope_carbon_12_14_different_class() {
        // `14.00324 - 12.011 = 1.992` truncates to `1`, different from
        // carbon-12/13's `0`. `[12CH4]`/`[14CH4]` are each a single bracket
        // atom (H4 is the bracket's compact hydrogen-count field, not
        // separate graph atoms) -- two disconnected fragments = 2 atoms
        // total, indices 0 and 1.
        let mol = parse("[12CH4].[14CH4]").unwrap();
        let inv = rdkit_inv(&mol);
        assert_ne!(
            inv[0], inv[1],
            "carbon-12 and carbon-14 must be different invariant classes, matching RDKit"
        );
    }

    #[test]
    fn rdkit_mode_isotope_oxygen_15_16_same_class() {
        // `15.00011 - 15.999 = -0.999` and `15.99491 - 15.999 = -0.004`
        // both truncate to `0`.
        let mol = parse("[15OH2].[16OH2]").unwrap();
        let inv = rdkit_inv(&mol);
        assert_eq!(
            inv[0], inv[1],
            "oxygen-15 and oxygen-16 must be the same invariant class, matching RDKit"
        );
    }

    #[test]
    fn rdkit_mode_isotope_oxygen_16_17_different_class() {
        // `16.99913 - 15.999 = 1.000` truncates to `1`, different from
        // oxygen-16's `0`.
        let mol = parse("[16OH2].[17OH2]").unwrap();
        let inv = rdkit_inv(&mol);
        assert_ne!(
            inv[0], inv[1],
            "oxygen-16 and oxygen-17 must be different invariant classes, matching RDKit"
        );
    }

    #[test]
    fn rdkit_mode_isotope_oxygen_16_18_different_class() {
        // `17.99916 - 15.999 = 2.000` truncates to `2`, different from both
        // oxygen-16 (`0`) and oxygen-17 (`1`).
        let mol = parse("[16OH2].[18OH2]").unwrap();
        let inv = rdkit_inv(&mol);
        assert_ne!(
            inv[0], inv[1],
            "oxygen-16 and oxygen-18 must be different invariant classes, matching RDKit"
        );
    }

    #[test]
    fn rdkit_mode_isotope_multi_element_stress() {
        // A single connected molecule carrying isotope labels on three
        // different elements at once (carbon-13, nitrogen-15, oxygen-18) --
        // each atom's isotope delta must be independently correct, not just
        // in isolation. `[13CH2]([15NH2])[18OH]` -- aminomethanol-shaped.
        let mol = parse("[13CH2]([15NH2])[18OH]").unwrap();
        let unlabeled = parse("C(N)O").unwrap();
        let inv = rdkit_inv(&mol);
        let inv_ref = rdkit_inv(&unlabeled);
        // C: 13.00335-12.011=0.992->0 (same as unlabeled C, atom 0)
        assert_eq!(
            inv[0], inv_ref[0],
            "carbon-13 in a multi-isotope molecule must still be delta=0"
        );
        // N: 15.00011-14.007=0.993->0 (same as unlabeled N, atom 1)
        assert_eq!(
            inv[1], inv_ref[1],
            "nitrogen-15 in a multi-isotope molecule must still be delta=0"
        );
        // O: 17.99916-15.999=2.000->2 (different from unlabeled O, atom 2)
        assert_ne!(
            inv[2], inv_ref[2],
            "oxygen-18 in a multi-isotope molecule must be delta=2, not delta=0"
        );
    }
}
