// ─── Type classification helpers ──────────────────────────────────────────────

/// Atom types with "extended" multiple-bond character (sp2/aromatic/sp
/// centers) — the MMFF94 property table's `mltb ≠ 0` set, i.e. atoms that
/// can be the single-bond end of an MMFF94 "single bond between
/// multiply-bonded atoms" (sbmb) pair (e.g. biphenyl's inter-ring bond, or
/// buta-1,3-diene's central C-C bond).
///
/// Exactly the 16 types that appear in `MMFF94_BOND_ENERGY`'s own
/// `(1, ti, tj, ...)` rows: {2,3,4,9,30,37,39,54,57,58,63,64,67,78,80,81}.
/// Measured, not assumed: a broader candidate (`oop_energy`'s trigonal-only
/// `SP2_TYPES` plus type 4) was tried first and made the 58-molecule corpus
/// harness's bond-miss count *worse* (38 vs. 16) — e.g. amide C(=O)-N
/// (types 3+10) and furan's aromatic C-O (types 37/38+43) both have no `(1,
/// ...)` row for that specific pair, so classifying them bt=1 under the
/// broader set was itself a fresh miss the tighter, evidence-derived set
/// avoids by routing them to their real `(0, ...)` row instead.
pub const MLTB_TYPES: &[u8] = &[2, 3, 4, 9, 30, 37, 39, 54, 57, 58, 63, 64, 67, 78, 80, 81];

/// Real bond order between two bonded atoms, defaulting to `Single` if no
/// bond exists between them (defensive only — every call site here passes
/// atom pairs already known to be bonded via `mol.neighbors`/`mol.bonds`).
/// The bond's order, with a directional (`/`, `\`) single bond read as
/// single: RDKit's MMFF typing sees `Bond::SINGLE` there.
fn bond_order_between(mol: &Molecule, a: usize, b: usize) -> BondOrder {
    match mol
        .bond_between(AtomIdx(a as u32), AtomIdx(b as u32))
        .map(|(_, bond)| bond.order)
    {
        Some(BondOrder::Up | BondOrder::Down) | None => BondOrder::Single,
        Some(order) => order,
    }
}

/// True if ALL of `atoms` belong to a single common SSSR ring of exactly
/// `size` atoms (not just individually each in *some* ring of that size —
/// e.g. an exocyclic substituent atom on a ring must not count).
fn atoms_share_ring_of_size(rings: &[Vec<AtomIdx>], atoms: &[usize], size: usize) -> bool {
    rings
        .iter()
        .any(|ring| ring.len() == size && atoms.iter().all(|&a| ring.contains(&AtomIdx(a as u32))))
}

/// `isAngleInRingOfSize3or4` (`AtomTyper.cpp`, pinned commit): LOCAL bond
/// adjacency, NOT SSSR-based -- deliberately distinct from
/// [`atoms_share_ring_of_size`] above, which IS SSSR-based and feeds angle
/// *type* classification. RDKit's real empirical-rule ring gate and its real
/// angle-type ring gate are two different, independently-defined mechanisms
/// (confirmed by direct source read, issue #227 Stage C). Returns 3 if `i`
/// and `k` are directly bonded (i-j-k-i, a 3-ring), 4 if `i` and `k` share a
/// common neighbor other than `j` (a 4-ring through that neighbor), else 0.
/// `i`/`k` are assumed already bonded to `j` (true for every call site,
/// which only ever iterates real `mol.neighbors(j)` pairs).
///
/// `pub` so downstream coverage-checkers (e.g. `chematic-3d`'s independent
/// MMFF94 coverage measurement) can pass the same `ring_size` into
/// [`crate::mmff94_energy::mmff94_angle_energy_resolved`] this module's own
/// `angle_energy`/`stretch_bend_energy` use, instead of hand-copying this
/// logic and risking drift (same rationale as [`OOP_SP2_TYPES`]).
pub fn is_angle_in_ring_of_size_3_or_4(mol: &Molecule, i: usize, j: usize, k: usize) -> u8 {
    if mol
        .bond_between(AtomIdx(i as u32), AtomIdx(k as u32))
        .is_some()
    {
        return 3;
    }
    let j = j as u32;
    let i_neighbors: std::collections::HashSet<u32> = mol
        .neighbors(AtomIdx(i as u32))
        .map(|(nb, _)| nb.0)
        .filter(|&n| n != j)
        .collect();
    let shares_neighbor = mol
        .neighbors(AtomIdx(k as u32))
        .any(|(nb, _)| nb.0 != j && i_neighbors.contains(&nb.0));
    if shares_neighbor { 4 } else { 0 }
}

/// Determine the MMFF94 bond-type index (Halgren 1996 bond-type index BT).
///
/// BT=1 only for a formally SINGLE, *non-aromatic* bond directly linking two
/// atoms that are BOTH independently "conjugation-capable" ([`MLTB_TYPES`]) —
/// the sbmb special case (e.g. biphenyl's inter-ring bond, r0=1.436 Å).
/// `BondOrder::Aromatic` gets BT=0 unconditionally, same as a real
/// double/triple bond — confirmed against a live RDKit oracle (issue #227
/// Phase 1B-0): benzene's own ring bond (types 37-37, `BondOrder::Aromatic`)
/// resolves to RDKit's `(bondType=0, kb=5.573, r0=1.374)`, not `(1, 5.178,
/// 1.436)` — i.e. the aromatic ring bond and biphenyl's non-aromatic
/// inter-ring bond between the *same two atom types* get *different* BT
/// values, which only the bond order (not the MLTB_TYPES membership alone)
/// can distinguish. A prior version of this function treated `Aromatic`
/// bonds the same as `Single` here and was never checked against a real
/// oracle for that specific claim — the `(1, 63, 63, ...)` row it was
/// "confirmed" against belonged to the *pre-fix*, wrong (63, not 37) atom
/// typing, so the match was coincidental, not a correctness check. A real
/// double/triple bond always gets BT=0 too (its atom types already encode
/// the multiply-bonded hybridization; the `(1, ...)` table rows are the
/// single-bond exception, not an alternate double-bond row), and a bond
/// where either atom is a plain non-conjugating type (e.g. sp3 CR) also
/// always gets BT=0 — confirmed empirically: zero `(1, 1, x, ...)` rows
/// exist in the 493-row bond table for sp3 carbon type 1 with any partner.
pub fn bond_type_for(ti: u8, tj: u8, order: BondOrder) -> u8 {
    if matches!(
        order,
        BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple | BondOrder::Aromatic
    ) {
        return 0;
    }
    if MLTB_TYPES.binary_search(&ti).is_ok() && MLTB_TYPES.binary_search(&tj).is_ok() {
        1
    } else {
        0
    }
}

/// Determine the MMFF94 angle-type index (Halgren 1996, types 0-8).
///
/// `bt_sum` = sum of the two flanking bonds' [`bond_type_for`] indices (0-2).
/// Ring-embedded angles (i,j,k all in one common 3- or 4-membered SSSR ring)
/// get dedicated types; everything else uses `bt_sum` directly:
///
/// | ring   | bt_sum=0 | bt_sum=1 | bt_sum=2 |
/// |--------|----------|----------|----------|
/// | none   | 0        | 1        | 2        |
/// | 3-ring | 3        | 5        | 6        |
/// | 4-ring | 4        | 7        | 8        |
///
/// Matches RDKit's real `getMMFFAngleType` formula (`AtomTyper.cpp:2412-`
/// `2447`, pinned commit — see `scripts/mmff94_provenance/PROVENANCE.md`'s
/// "Stretch-bend" row, Priority 2C addendum): `angleType = ring_size; if
/// bond_type_sum != 0 { angleType += bond_type_sum + ring_size - 2 }`. A
/// prior version of this table gave 3-ring bt_sum=2 -> 8 and 4-ring
/// bt_sum=1 -> 6 / bt_sum=2 -> 7, which disagreed with RDKit's formula —
/// fixed here (issue #227). Measured LATENT on the 265-molecule Wave 1
/// corpus (0/113 reachable ring-embedded angle triples hit these branches),
/// but it is a real, independently-provable formula bug and a needed
/// correct input to [`crate::mmff94_energy::mmff94_stbn`]'s stretch-bend
/// type classification (a wrong angle_type feeds directly into
/// `getMMFFStretchBendType`'s first argument).
///
/// Confirmed empirically against the angle table: `(3, 22, 22, 22)`
/// (all-CR3R, cyclopropane) has θ0≈60°; `(4, 6, 20, 20)` (4-ring) has
/// θ0≈93° — i.e. 3-ring/4-ring are not swapped.
pub fn angle_type_for(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    i: usize,
    j: usize,
    k: usize,
    types: &[u8],
) -> u8 {
    let bt_ij = bond_type_for(types[i], types[j], bond_order_between(mol, i, j));
    let bt_jk = bond_type_for(types[j], types[k], bond_order_between(mol, j, k));
    let bt_sum = bt_ij + bt_jk;

    if atoms_share_ring_of_size(rings, &[i, j, k], 3) {
        return match bt_sum {
            0 => 3,
            1 => 5,
            _ => 6,
        };
    }
    if atoms_share_ring_of_size(rings, &[i, j, k], 4) {
        return match bt_sum {
            0 => 4,
            1 => 7,
            _ => 8,
        };
    }
    bt_sum
}

/// Determine the MMFF94 stretch-bend-type index (Halgren MMFF.V, types
/// 0-11) for angle i-j-k.
///
/// Ported verbatim from RDKit's `getMMFFStretchBendType`
/// (`AtomTyper.cpp:2480-2508`, pinned commit — see
/// `scripts/mmff94_provenance/PROVENANCE.md`'s "Stretch-bend" row,
/// Priority 2C addendum) and the diagnostic's own `resolve_rdkit` (issue
/// #227, `mmff94_stbn_equivalence_diagnostic_227.rs`). Angle types 1, 5, 7
/// each split into two distinct stretch-bend types depending on whether
/// either flanking bond individually has MMFF bond type 1 — [`angle_type_for`]'s
/// `bt_sum` only records the *sum*, discarding exactly the information this
/// split needs, which is why stretch-bend cannot reuse `angle_type` directly
/// as its own [`crate::mmff94_energy::MMFF94_STBN`] table key (issue #227
/// root cause).
///
/// `ta`/`tc` are the outer atom MMFF types (i and k); `bond_type_ij`/
/// `bond_type_jk` are the individual (unswapped) [`bond_type_for`] results
/// for the i-j and j-k bonds. The two bond-type arguments fed into the
/// classification match are canonicalized on `ta <= tc` — this is **not**
/// the same swap rule as [`mmff94_stbn_type_only`](crate::mmff94_energy::mmff94_stbn_type_only)'s
/// own i/k table-lookup canonicalization (that one additionally tie-breaks
/// on bond type when `ta == tc`; see `AtomTyper.cpp:3598-3600`). Do not
/// conflate the two swap rules.
pub fn stretch_bend_type_for(
    angle_type: u8,
    ta: u8,
    tc: u8,
    bond_type_ij: u8,
    bond_type_jk: u8,
) -> u8 {
    let (arg1, arg2) = if ta <= tc {
        (
            bond_type_ij,
            if ta < tc { bond_type_jk } else { bond_type_ij },
        )
    } else {
        (bond_type_jk, bond_type_ij)
    };
    match angle_type {
        1 => {
            if arg1 != 0 || arg1 == arg2 {
                1
            } else {
                2
            }
        }
        2 => 3,
        4 => 4,
        3 => 5,
        5 => {
            if arg1 != 0 || arg1 == arg2 {
                6
            } else {
                7
            }
        }
        6 => 8,
        7 => {
            if arg1 != 0 || arg1 == arg2 {
                9
            } else {
                10
            }
        }
        8 => 11,
        _ => 0,
    }
}

/// Port of RDKit's `isTorsionInRingOfSize4or5` (`AtomTyper.cpp:403-447`,
/// pinned commit — see `scripts/mmff94_provenance/PROVENANCE.md`'s
/// "Torsion" row) and the diagnostic's own `rdkit_ring_size_4_or_5` (issue
/// #227, `mmff94_torsion_equivalence_diagnostic_227.rs`). Purely LOCAL
/// bond-adjacency, NOT SSSR-based: 4-ring iff i-l are directly bonded;
/// 5-ring iff i and l, excluding their ring neighbours j and k
/// respectively, share a common neighbour. Returns 0 if neither.
fn ring_size_4_or_5(mol: &Molecule, i: AtomIdx, j: AtomIdx, k: AtomIdx, l: AtomIdx) -> u8 {
    if mol.bond_between(i, l).is_some() {
        return 4;
    }
    let nbrs_i: Vec<AtomIdx> = mol
        .neighbors(i)
        .map(|(n, _)| n)
        .filter(|&n| n != j)
        .collect();
    let has_common = mol
        .neighbors(l)
        .map(|(n, _)| n)
        .filter(|&n| n != k)
        .any(|n| nbrs_i.contains(&n));
    if has_common { 5 } else { 0 }
}

/// True when RDKit's real torsion resolution generates NO term at all for
/// this i-j-k-l torsion, by design -- not a missing-parameter gap.
///
/// Halgren's empirical-rule cascade (`getMMFFTorsionEmpiricalRuleParams`,
/// `AtomTyper.cpp`, rule (a) per the public transcription cited in
/// `scripts/mmff94_provenance/PROVENANCE.md`'s Torsion entry) omits the
/// torsion term entirely whenever either central atom (`type_j`/`type_k`,
/// the j-k bond this torsion rotates around) has MMFF's `lin` flag
/// ([`crate::mmff94_numeric_type_registry::Mmff94NumericTypeInfo::linear`],
/// e.g. type 4 CSP/nitrile-or-acetylenic carbon, type 53 `=N=`/cumulated
/// azide nitrogen): rotating around a bond whose other end is a linear
/// (180°) center changes no real geometry, so there is nothing to
/// parameterize. Issue #227 Phase 1: measured as the exact, complete
/// explanation for the only 3 genuine `table_gap` Torsion instances left in
/// the 265-molecule Wave 1 corpus after the classification fix below —
/// oracle-confirmed (`GetMMFFTorsionParams` returns `None` for all 3, and
/// each central atom's registered `linear` flag matches RDKit's own MMFF
/// atom type exactly at that atom). A caller that already omits the term for
/// an unresolved lookup (chematic-ff's own `torsion_energy`, which just adds
/// nothing on `mmff94_torsion_energy(..) == None`) needs no code change for
/// physics; this exists so coverage-reporting/diagnostic callers (the
/// `Mmff94BondAngleStrict` gate under `include_torsion_oop_in_gate`, the
/// Phase 1A audit) can tell "correctly no term" apart from "genuinely
/// missing" instead of counting both as the same failure.
pub fn torsion_no_term_by_design(type_j: u8, type_k: u8) -> bool {
    let linear = |t: u8| {
        crate::mmff94_numeric_type_registry::mmff94_numeric_type_info(t)
            .is_some_and(|info| info.linear)
    };
    linear(type_j) || linear(type_k)
}

/// Determine the MMFF94 torsion-type index (Halgren 1996, types 0-8).
///
/// Ported verbatim from RDKit's `getMMFFTorsionType` (`AtomTyper.cpp:2528-`
/// `2571`, pinned commit — see `scripts/mmff94_provenance/PROVENANCE.md`'s
/// "Torsion" row) and the diagnostic's own `rdkit_torsion_type` (issue #227,
/// `mmff94_torsion_equivalence_diagnostic_227.rs`).
///
/// The base (non-ring) code is the j-k bond's OWN [`bond_type_for`] result —
/// NOT atom-type ([`MLTB_TYPES`]) membership, which is what a prior version
/// of this function used and which disagreed with RDKit's real
/// classification on 100% of the diagnostic's 1,107 routing-bug candidates
/// and 76.3% of ALL 13,530 torsion instances in the 265-molecule corpus (a
/// double/triple/aromatic j-k bond always gets `bond_type_jk=0` from
/// `bond_type_for` regardless of `tj`/`tk`'s own MLTB membership — the prior
/// atom-type-membership formula couldn't see that at all). An
/// empirically-required override bumps the base code to 2 when
/// `bond_type_jk==0 && order_jk==Single && (bond_type_ij==1 ||
/// bond_type_kl==1)` — MMFF.IV page 609's simple condition fails CYGUAN01 in
/// RDKit's own test suite; this corrected condition is RDKit's real code,
/// not derivable from Halgren's paper alone.
///
/// Ring-embedded torsions get dedicated types 4 (4-ring) / 5 (5-ring),
/// determined by [`ring_size_4_or_5`] — a purely LOCAL bond-adjacency check,
/// NOT SSSR-based, replacing this function's prior
/// [`atoms_share_ring_of_size`]-based check entirely (the diagnostic's whole
/// point is that the SSSR-based check was the wrong mechanism here). The
/// 4-ring override additionally requires i-k and j-l to NOT be directly
/// bonded (excludes degenerate/bridged cases); the 5-ring override
/// additionally requires at least one of the four atoms to be MMFF numeric
/// type 1 (`ti==1 || tj==1 || tk==1 || tl==1`) — a condition the prior
/// SSSR-based check had no equivalent of at all.
#[allow(clippy::too_many_arguments)]
pub fn torsion_type_for(
    mol: &Molecule,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    ti: u8,
    tj: u8,
    tk: u8,
    tl: u8,
) -> u8 {
    torsion_type_pair_for(mol, i, j, k, l, ti, tj, tk, tl).0
}

/// RDKit's `getMMFFTorsionType`: the torsion type and, when a 4- or
/// 5-membered ring sets it to 4 or 5, the type the bonds alone give
/// (RDKit's second lookup pass uses it; 0 otherwise).
#[allow(clippy::too_many_arguments)]
pub fn torsion_type_pair_for(
    mol: &Molecule,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    ti: u8,
    tj: u8,
    tk: u8,
    tl: u8,
) -> (u8, u8) {
    let order_ij = bond_order_between(mol, i, j);
    let order_jk = bond_order_between(mol, j, k);
    let order_kl = bond_order_between(mol, k, l);
    let bond_type_ij = bond_type_for(ti, tj, order_ij);
    let bond_type_jk = bond_type_for(tj, tk, order_jk);
    let bond_type_kl = bond_type_for(tk, tl, order_kl);

    let mut torsion_type = bond_type_jk;
    if bond_type_jk == 0
        && order_jk == BondOrder::Single
        && (bond_type_ij == 1 || bond_type_kl == 1)
    {
        torsion_type = 2;
    }

    let (ai, aj, ak, al) = (
        AtomIdx(i as u32),
        AtomIdx(j as u32),
        AtomIdx(k as u32),
        AtomIdx(l as u32),
    );
    let ring_size = ring_size_4_or_5(mol, ai, aj, ak, al);
    let mut second = 0;
    if ring_size == 4 && mol.bond_between(ai, ak).is_none() && mol.bond_between(aj, al).is_none() {
        second = torsion_type;
        torsion_type = 4;
    } else if ring_size == 5 && (ti == 1 || tj == 1 || tk == 1 || tl == 1) {
        second = torsion_type;
        torsion_type = 5;
    }

    (torsion_type, second)
}

/// MMFF94 torsion parameters of `i-j-k-l` as RDKit's force field takes
/// them: the MMFFTOR table through RDKit's equivalence-level stages
/// (`MMFFTorCollection::getMMFFTorParams`), else Halgren's empirical rule
/// (`getMMFFTorsionEmpiricalRuleParams`). `None` when every barrier is zero
/// (RDKit adds no term).
pub fn mmff94_torsion_term_params(
    mol: &Molecule,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    types: &[u8],
) -> Option<TorsionEnergyParams> {
    mmff_torsion_term_params_variant(mol, i, j, k, l, types, false)
}

/// [`mmff94_torsion_term_params`] with the MMFF94s torsion table when
/// `mmffs` is set (RDKit `mmffVariant="MMFF94s"`).
pub(crate) fn mmff_torsion_term_params_variant(
    mol: &Molecule,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    types: &[u8],
    mmffs: bool,
) -> Option<TorsionEnergyParams> {
    let (ti, tj, tk, tl) = (types[i], types[j], types[k], types[l]);
    let pair = torsion_type_pair_for(mol, i, j, k, l, ti, tj, tk, tl);
    let table = if mmffs {
        crate::mmff94_energy::MMFF94S_TORSION_ENERGY
    } else {
        crate::mmff94_energy::MMFF94_TORSION_ENERGY
    };
    let p = torsion_table_lookup_rdkit(table, pair, ti, tj, tk, tl)
        .unwrap_or_else(|| torsion_empirical_rule(mol, j, k, tj, tk));
    (p.v1 != 0.0 || p.v2 != 0.0 || p.v3 != 0.0).then_some(p)
}

type TorsionTable = [(u8, u8, u8, u8, u8, f64, f64, f64)];

fn torsion_table_search(
    table: &TorsionTable,
    tt: u8,
    i: u8,
    j: u8,
    k: u8,
    l: u8,
) -> Option<TorsionEnergyParams> {
    table
        .binary_search_by_key(&(tt, i, j, k, l), |&(t0, t1, t2, t3, t4, _, _, _)| {
            (t0, t1, t2, t3, t4)
        })
        .ok()
        .map(|idx| {
            let (_, _, _, _, _, v1, v2, v3) = table[idx];
            TorsionEnergyParams { v1, v2, v3 }
        })
}

/// RDKit's five-stage lookup: equivalence levels 1-1-1-1, 2-2-2-2, 3-2-2-5,
/// 5-2-2-3, 5-2-2-5 on the torsion type, then again on the second type.
fn torsion_table_lookup_rdkit(
    table: &TorsionTable,
    (first, second): (u8, u8),
    ti: u8,
    tj: u8,
    tk: u8,
    tl: u8,
) -> Option<TorsionEnergyParams> {
    let eq = |t: u8| {
        crate::mmff94_numeric_type_registry::mmff94_numeric_type_info(t)
            .map(|info| info.equivalence_levels)
    };
    let (eq_i, eq_l) = (eq(ti)?, eq(tl)?);
    let mut found = None;
    let mut iter = 0usize;
    let mut max_iter = 5usize;
    let mut tt = first;
    while (iter < max_iter && (found.is_none() || max_iter == 4))
        || (iter == 4 && first == 5 && second != 0)
    {
        if max_iter == 5 && iter == 4 {
            max_iter = 4;
            iter = 0;
            tt = second;
        }
        let (wi, wl) = match iter {
            1 => (1, 3),
            2 => (3, 1),
            n => (n, n),
        };
        // Level index 4 does not exist (RDKit reads past its array there).
        let (Some(&ci), Some(&cl)) = (eq_i.get(wi), eq_l.get(wl)) else {
            break;
        };
        let (mut ci, mut cj, mut ck, mut cl) = (ci, tj, tk, cl);
        if cj > ck {
            std::mem::swap(&mut cj, &mut ck);
            std::mem::swap(&mut ci, &mut cl);
        } else if cj == ck && ci > cl {
            std::mem::swap(&mut ci, &mut cl);
        }
        if let Some(p) = torsion_table_search(table, tt, ci, cj, ck, cl) {
            found = Some(p);
            if max_iter == 4 {
                break;
            }
        }
        iter += 1;
    }
    found
}

/// Halgren's empirical torsion rule (MMFF.IV), as RDKit's
/// `getMMFFTorsionEmpiricalRuleParams` writes it, for the central bond j-k.
fn torsion_empirical_rule(
    mol: &Molecule,
    j: usize,
    k: usize,
    tj: u8,
    tk: u8,
) -> TorsionEnergyParams {
    use crate::mmff94_numeric_type_registry::mmff94_numeric_type_info;
    let zero = TorsionEnergyParams {
        v1: 0.0,
        v2: 0.0,
        v3: 0.0,
    };
    let (Some(pj), Some(pk)) = (mmff94_numeric_type_info(tj), mmff94_numeric_type_info(tk)) else {
        return zero;
    };
    let order = bond_order_between(mol, j, k);
    let z = [
        mol.atom(AtomIdx(j as u32)).element.atomic_number(),
        mol.atom(AtomIdx(k as u32)).element.atomic_number(),
    ];
    let constants = z.map(torsion_element_constants);
    let u = [constants[0].0, constants[1].0];
    let v = [constants[0].1, constants[1].1];
    let w = [constants[0].2, constants[1].2];
    let n_jk = f64::from((pj.coordination - 1) * (pk.coordination - 1));
    let v3 = TorsionEnergyParams {
        v1: 0.0,
        v2: 0.0,
        v3: (v[0] * v[1]).sqrt() / n_jk,
    };
    if pj.linear || pk.linear {
        return zero;
    }
    if let Some(barrier) = multiple_bond_torsion_barrier(pj, pk, order, u) {
        return barrier;
    }
    if let Some(barrier) = coordination_torsion_barrier(pj, pk, v3, zero) {
        return barrier;
    }
    if let Some(barrier) = conjugated_torsion_barrier(pj, pk, order, z, u) {
        return barrier;
    }
    if matches!(z[0], 8 | 16) && matches!(z[1], 8 | 16) {
        TorsionEnergyParams {
            v1: 0.0,
            v2: -(w[0] * w[1]).sqrt(),
            v3: 0.0,
        }
    } else {
        v3
    }
}

fn v2_torsion_barrier(u: [f64; 2], beta: f64, pi: f64) -> TorsionEnergyParams {
    TorsionEnergyParams {
        v1: 0.0,
        v2: beta * pi * (u[0] * u[1]).sqrt(),
        v3: 0.0,
    }
}

fn multiple_bond_torsion_barrier(
    left: &crate::mmff94_numeric_type_registry::Mmff94NumericTypeInfo,
    right: &crate::mmff94_numeric_type_registry::Mmff94NumericTypeInfo,
    order: BondOrder,
    u: [f64; 2],
) -> Option<TorsionEnergyParams> {
    if left.aromatic && right.aromatic && order == BondOrder::Aromatic {
        let beta = if matches!((left.valence, right.valence), (3, 4) | (4, 3)) {
            3.0
        } else {
            6.0
        };
        let pi = if !left.has_pi_lone_pair && !right.has_pi_lone_pair {
            0.5
        } else {
            0.3
        };
        Some(v2_torsion_barrier(u, beta, pi))
    } else if order == BondOrder::Double {
        let pi = if left.multiple_bond_count == 2 && right.multiple_bond_count == 2 {
            1.0
        } else {
            0.4
        };
        Some(v2_torsion_barrier(u, 6.0, pi))
    } else {
        None
    }
}

fn coordination_torsion_barrier(
    left: &crate::mmff94_numeric_type_registry::Mmff94NumericTypeInfo,
    right: &crate::mmff94_numeric_type_registry::Mmff94NumericTypeInfo,
    v3: TorsionEnergyParams,
    zero: TorsionEnergyParams,
) -> Option<TorsionEnergyParams> {
    if left.coordination == 4 && right.coordination == 4 {
        Some(v3)
    } else if left.coordination == 4 {
        Some(if torsion_end_conjugates(right) { zero } else { v3 })
    } else if right.coordination == 4 {
        Some(if torsion_end_conjugates(left) { zero } else { v3 })
    } else {
        None
    }
}

fn conjugated_torsion_barrier(
    left: &crate::mmff94_numeric_type_registry::Mmff94NumericTypeInfo,
    right: &crate::mmff94_numeric_type_registry::Mmff94NumericTypeInfo,
    order: BondOrder,
    atomic_numbers: [u8; 2],
    u: [f64; 2],
) -> Option<TorsionEnergyParams> {
    let (left_multiple, right_multiple) =
        (left.multiple_bond_count, right.multiple_bond_count);
    let (left_lone_pair, right_lone_pair) =
        (left.has_pi_lone_pair, right.has_pi_lone_pair);
    let conjugated = (order == BondOrder::Single && left_multiple != 0 && right_multiple != 0)
        || (left_multiple != 0 && right_lone_pair)
        || (left_lone_pair && right_multiple != 0);
    if !conjugated {
        return None;
    }
    let pi = if left_lone_pair && right_lone_pair {
        return Some(TorsionEnergyParams {
            v1: 0.0,
            v2: 0.0,
            v3: 0.0,
        });
    } else if left_lone_pair && right_multiple != 0 {
        pi_lone_pair_factor(left_multiple, atomic_numbers)
    } else if right_lone_pair && left_multiple != 0 {
        pi_lone_pair_factor(right_multiple, atomic_numbers)
    } else if (left_multiple == 1 || right_multiple == 1)
        && (atomic_numbers[0] != 6 || atomic_numbers[1] != 6)
    {
        0.4
    } else {
        0.15
    };
    Some(v2_torsion_barrier(u, 6.0, pi))
}

fn torsion_element_constants(atomic_number: u8) -> (f64, f64, f64) {
    match atomic_number {
        6 => (2.0, 2.12, 0.0),
        7 => (2.0, 1.5, 0.0),
        8 => (2.0, 0.2, 2.0),
        14 => (1.25, 1.22, 0.0),
        15 => (1.25, 2.40, 0.0),
        16 => (1.25, 0.49, 8.0),
        _ => (0.0, 0.0, 0.0),
    }
}

fn torsion_end_conjugates(
    info: &crate::mmff94_numeric_type_registry::Mmff94NumericTypeInfo,
) -> bool {
    (info.coordination == 3
        && matches!(info.valence, 4 | 34)
            || info.coordination == 3 && info.multiple_bond_count != 0)
        || (info.coordination == 2 && (info.valence == 3 || info.multiple_bond_count != 0))
}

fn pi_lone_pair_factor(multiple_bonds: u8, atomic_numbers: [u8; 2]) -> f64 {
    if multiple_bonds == 1 {
        0.5
    } else if atomic_numbers.map(periodic_table_row) == [2, 2] {
        0.3
    } else {
        0.15
    }
}

fn periodic_table_row(atomic_number: u8) -> u8 {
    match atomic_number {
        3..=10 => 1,
        11..=18 => 2,
        19..=36 => 3,
        37..=54 => 4,
        _ => 0,
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────
