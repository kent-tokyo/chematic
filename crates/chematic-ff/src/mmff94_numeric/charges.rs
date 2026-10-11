// ── Partial charge calculation ───────────────────────────────────────────────

/// RDKit's real `computeMMFFCharges` does NOT feed the molecule's raw/literal
/// formal charge (`atom->getFormalCharge()`) into equation 15 -- it first
/// computes a separate, MMFF-atom-TYPE-derived formal charge ("MMFFFormalCharge"
/// in RDKit's own naming) via a dedicated per-type switch statement that runs
/// BEFORE the main charge loop (`AtomTyper.cpp` lines ~3095-3350, "We need to
/// set formal charges upfront"), and uses THAT value everywhere a formal
/// charge is needed: as this atom's own q0 in `(1-M*v)*q0`, and as the
/// NEIGHBOR formal charge source for both the `v*sumFormalCharge` term and
/// the anionic-neighbor-leak adjustment. For most types (every type not
/// named in the switch) the derived charge defaults to 0.0 -- RDKit's own
/// pre-switch initial value -- even when the atom's raw SMILES formal charge
/// is nonzero (e.g. nitro nitrogen type 45, azide types 47/53, sulfoxide
/// type 17: none of these are switch cases, so RDKit treats them as
/// formal-charge-neutral for equation 15's purposes despite how the input
/// structure happened to place the literal charge).
///
/// Issue #227 Phase 2 Step 6: this was the root cause of the 67/6693-atom
/// residual left after the BCI bond-type fix (PR #331) -- see
/// `scripts/mmff94_provenance/PROVENANCE.md`'s Charges/BCI section for the
/// full molecule-by-molecule evidence.
///
/// **Faithful but intentionally partial port.** Implemented and
/// independently verified against a live RDKit 2026.03.4 oracle query (not
/// merely re-derived from this port's own output): the unconditional
/// +1/+2/+3/-1 "simple" type groups; the O2CM/SM (types 32/72)
/// neighbor-counting redistribution for its carbon-neighbor
/// (carboxylate/thiocarboxylate), nitro/nitrate-nitrogen-neighbor (type 45,
/// 3-terminal-O case only), and sulfone/sulfonate/sulfonamide-sulfur-neighbor
/// (type 18) branches (synthetic fixtures: acetate, nitrate ion,
/// methanesulfonate, dimethyl sulfone; see the `o2cm_sm_*`/`nitrate_ion_*`/
/// `sulfone_and_sulfonate_*` tests below). Implemented but **not**
/// independently oracle-verified (zero corpus exposure either way, so
/// nothing to falsify against -- see below): type 62's (NM) extra
/// "subtract half of each positively-charged neighbor's derived charge"
/// adjustment (`AtomTyper.cpp` lines ~3378-3383, applied in
/// `mmff94_charges_numeric`'s Step 1, not here, since -- unlike the
/// `isDoubleZero(v)` leak -- this adjustment happens BEFORE the `(1-M*v)`
/// multiplication and therefore cannot be expressed as a separate additive
/// term when `v != 0`, which it does for type 62, `fcadj=0.25`).
///
/// **NOT ported** (falls through to the 0.0 default, same as RDKit's own
/// fallback when no switch condition matches):
/// - O2CM/SM's phosphorus-neighbor (type 25, phosphate/phosphonate O),
///   thiosulfinate-sulfur-neighbor (type 73), and perchlorate-chlorine-
///   neighbor (type 77) branches.
/// - Type 76 (N5M): formal charge shared across N5M nitrogens co-membered in
///   a 5-ring (needs ring perception).
/// - Types 55/56/81 (NIM+/N5A+/N5B+): formal charge averaged across a
///   conjugated cationic-nitrogen network reached via alternating type-57/80
///   carbons (needs BFS over a subset of the molecule graph).
/// - Type 61's diazonium special case (+1 bump when bonded to a type-42
///   diazonium nitrogen).
///
/// None of these gaps are guesses papered over with an untested formula --
/// zero atoms of types 62/76/55/56/61/81, and zero O2CM/SM atoms with a
/// phosphorus/type-73/type-77 neighbor, appear anywhere in the 264-molecule
/// Wave 1 corpus this port was measured against (confirmed by a dedicated,
/// committed, independently re-runnable survey,
/// `crates/chematic-3d/examples/mmff94_fchg_type_exposure_survey_227.rs` --
/// issue #227 Phase 2 Step 6), so none of this can be silently masking a
/// corpus-visible bug. Flagged as a follow-up for whoever next touches
/// MMFF94 formal-charge handling on a corpus that does exercise these
/// types.
///
/// **Blast radius, stated explicitly in both directions**: within the
/// 264-molecule corpus, this fix changes the computed charge for exactly
/// 5/6,693 atoms (see the per-atom-join measurement in
/// `scripts/mmff94_provenance/PROVENANCE.md`) -- small, and the reason no
/// downstream 3D-pipeline re-measurement was run for this step (contrast
/// the prior BCI bond-type fix, which moved 1,620 atoms and produced one
/// genuine new stereo violation). Outside this corpus, the behavioral
/// change is broader than that number suggests: it applies to *any*
/// molecule containing a carboxylate, sulfonate/sulfamate, nitrate,
/// nitro, azide, sulfoxide, or quaternary-ammonium group -- the corpus
/// simply happens to contain none of the anionic O2CM/SM forms (every
/// type-32 atom in it has a sulfone/nitro/sulfoxide neighbor, never a
/// carbon or phosphorus neighbor) and only 3 molecules combining nitro/
/// azide/sulfoxide with the absent-from-switch types this fix's Step
/// 1/Step 3 change directly affects.
fn mmff_derived_formal_charge(
    mol: &Molecule,
    view: &Molecule,
    rings: &[Vec<AtomIdx>],
    types: &[u8],
    idx: AtomIdx,
) -> f64 {
    match types[idx.0 as usize] {
        // "Non-complicated" +1/+2/+3/-1 atom types (`AtomTyper.cpp`'s
        // `computeMMFFCharges` switch, cases with a single hardcoded `fChg`
        // assignment) -- independent of the atom's own raw formal charge.
        34 | 49 | 51 | 54 | 58 | 92 | 93 | 94 | 97 => 1.0,
        87 | 95 | 96 | 98 | 99 => 2.0,
        88 => 3.0,
        35 | 62 | 89 | 90 | 91 => -1.0,
        // O2CM (32) / SM (72): formal charge shared/localized across
        // terminal O/S atoms bonded to a common neighbor.
        32 | 72 => o2cm_sm_formal_charge(view, types, idx),
        // N5M: -1 shared by the N5M nitrogens of the (first) ring holding it.
        76 => rings
            .iter()
            .find(|r| r.contains(&idx))
            .map(|r| r.iter().filter(|a| types[a.0 as usize] == 76).count())
            .filter(|&k| k > 0)
            .map_or(0.0, |k| -1.0 / k as f64),
        // NIM+ / N5A+ / N5B+ / N5+: the formal charges of the nitrogens of
        // these types conjugated through C5A/CIM+-type carbons (57, 80),
        // shared evenly among them.
        55 | 56 | 81 => conjugated_cation_formal_charge(mol, types, idx),
        // N=N+ next to a diazonium N (42): +1.
        61 if mol.neighbors(idx).any(|(nb, _)| types[nb.0 as usize] == 42) => 1.0,
        _ => 0.0,
    }
}

/// RDKit's `case 55: case 56: case 81:` in `computeMMFFCharges`.
fn conjugated_cation_formal_charge(mol: &Molecule, types: &[u8], idx: AtomIdx) -> f64 {
    let n = mol.atom_count();
    let mut conj = vec![false; n];
    conj[idx.0 as usize] = true;
    let mut total = f64::from(mol.atom(idx).charge);
    let mut n_conj = 1usize;
    let mut old = 0usize;
    while n_conj > old {
        old = n_conj;
        for i in 0..n {
            if !conj[i] {
                continue;
            }
            for (nb, _) in mol.neighbors(AtomIdx(i as u32)) {
                if !matches!(types[nb.0 as usize], 57 | 80) {
                    continue;
                }
                for (nb2, _) in mol.neighbors(nb) {
                    let j = nb2.0 as usize;
                    if matches!(types[j], 55 | 56 | 81) && !conj[j] {
                        conj[j] = true;
                        total += f64::from(mol.atom(nb2).charge);
                        n_conj += 1;
                    }
                }
            }
        }
    }
    total / n_conj as f64
}

/// O2CM/SM formal-charge redistribution, RDKit's `case 32: case 72:`: the
/// first neighbour that matches a branch decides, from its terminal O/S
/// (degree 1) and secondary nitrogen (degree 2, not aromatic) neighbours.
fn o2cm_sm_formal_charge(view: &Molecule, types: &[u8], idx: AtomIdx) -> f64 {
    for (nbr, _) in view.neighbors(idx) {
        let nbr_elem = view.atom(nbr).element;
        let nbr_type = types[nbr.0 as usize];
        let mut n_sec_n = 0usize;
        let mut n_term_os = 0usize;
        for (nb2, _) in view.neighbors(nbr) {
            let a = view.atom(nb2);
            if a.element == Element::N && view.degree(nb2) == 2 && !a.aromatic {
                n_sec_n += 1;
            }
            if matches!(a.element, Element::O | Element::S) && view.degree(nb2) == 1 {
                n_term_os += 1;
            }
        }
        // Deprotonated-sulfonamide fixup: a sulfur with 2 terminal O/S and 1
        // secondary N is not treated as having a "replaceable" secondary N.
        if nbr_elem == Element::S && n_term_os == 2 && n_sec_n == 1 {
            n_sec_n = 0;
        }
        let shared = |k: usize, zero_at_one: bool| {
            if k == 1 {
                if zero_at_one { 0.0 } else { -1.0 }
            } else {
                -((k - 1) as f64) / (k as f64)
            }
        };
        if nbr_elem == Element::C && n_term_os > 0 {
            return shared(n_term_os, false);
        }
        if nbr_type == 45 && n_term_os == 3 {
            return -1.0 / 3.0;
        }
        if nbr_type == 25 && n_term_os > 0 {
            return shared(n_term_os, true);
        }
        if nbr_type == 18 && n_term_os > 0 {
            let total = n_sec_n + n_term_os;
            return if total == 2 {
                0.0
            } else {
                -((total as f64) - 2.0) / (n_term_os as f64)
            };
        }
        if nbr_type == 73 && n_term_os > 0 {
            return shared(n_term_os, true);
        }
        if nbr_type == 77 && n_term_os > 0 {
            return -1.0 / n_term_os as f64;
        }
    }
    0.0
}

/// Compute MMFF94 partial charges using the full PBCI+CHG tables (Halgren 1996).
///
/// Implements equation 15 from MMFF.V paper. For most neutral organic atoms
/// (fcadj=0, no formal charge), this reduces to:
///   q_i = Σ_{j bonded} bci(j→i)
///
/// Returns per-atom partial charges in units of elementary charge.
///
/// Issue #227 Phase 2 (BCI investigation): the BCI (bond-charge-increment)
/// step below reads bond order from [`assign_mmff94_numeric_types_with_view`]'s
/// re-perceived molecule, the same MMFF-specific Kekulized view Phase 1
/// already threads through `bond_type_for`/`angle_type_for`/
/// `torsion_type_for`/`stretch_bend_type_for` -- not this function's own
/// `mol` argument. RDKit's real `computeMMFFCharges` calls the identical
/// `getMMFFBondType(bond)` its bond-stretch code calls, on the identical
/// (sanitized/Kekulized) `mol` object (`AtomTyper.cpp:3071-3488`, pinned
/// commit) -- there is no separate "charge bond order" in RDKit's own
/// algorithm, so there must not be one here either.
pub fn mmff94_charges_numeric(mol: &Molecule) -> Result<Vec<f64>, NumericTypeError> {
    let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
    let n = mol.atom_count();
    let mut charges = vec![0.0f64; n];

    // Derived MMFF formal charges ("MMFFFormalCharge" in RDKit's own naming)
    // -- NOT the molecule's raw/literal `atom.charge`. See
    // `mmff_derived_formal_charge`'s doc for why this distinction is load-
    // bearing (issue #227 Phase 2 Step 6 BCI residual fix).
    let rings = chematic_perception::rdkit_sssr_ring_order(mol).unwrap_or_else(|| {
        chematic_perception::find_symmetrized_sssr(mol)
            .rings()
            .to_vec()
    });
    let fchg: Vec<f64> = (0..n)
        .map(|i| mmff_derived_formal_charge(mol, &mmff_mol, &rings, &types, AtomIdx(i as u32)))
        .collect();

    // Step 1: formal charge contribution (scaled by fcadj)
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let (_, fcadj) = pbci_for(types[i]);
        let mut q0 = fchg[i];
        // Type 62 (NM, anionic divalent N) special case (`AtomTyper.cpp`
        // lines ~3378-3383): subtract half of each positively-charged
        // neighbor's derived formal charge from q0, purely locally -- this
        // adjustment is never written back to `fchg`, so it is invisible to
        // any OTHER atom reading this atom's charge as a neighbor value
        // (matching RDKit's own `getMMFFFormalCharge`, which always returns
        // the switch-only stored value). Must happen here, before the
        // `(1-M*v)` multiplication below, not as a separate additive term
        // in Step 3 -- unlike the `isDoubleZero(v)` leak, this fires
        // whenever `v != 0` (fcadj(62) = 0.25), so it cannot be expressed
        // as a no-op-when-v-is-zero additive term the way that leak can.
        if types[i] == 62 {
            for b in bonds_of(mol, idx) {
                let nbr_fchg = fchg[b.neighbor.0 as usize];
                if nbr_fchg > 0.0 {
                    q0 -= nbr_fchg / 2.0;
                }
            }
        }
        // (1 - M*v)*q0 simplified for fcadj=0 (most atoms): charge[i] = 0
        // For charged atoms with fcadj > 0:
        let m = bonds_of(mol, idx).len() as f64;
        charges[i] = (1.0 - m * fcadj) * q0;
    }

    // Step 2: BCI contributions from each bond. `mmff_mol` (not `mol`) is
    // the bond-order source -- same reperceived view as Step 1's atom
    // types, and the same view Phase 1 already requires for every other
    // bond-order-dependent MMFF94 classification. `mmff_mol` has the same
    // atom count/bond topology as `mol` (only ring-bond `BondOrder` values
    // can differ), so atom-index-based `types[i]`/`types[j]` stay valid.
    for (_, bond) in mmff_mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        let ti = types[i];
        let tj = types[j];
        let bt = crate::mmff94_minimizer::bond_type_for(ti, tj, bond.order);

        // Contribution to atom i
        let ci =
            lookup_chg_contribution(bt, ti, tj).unwrap_or_else(|| pbci_for(ti).0 - pbci_for(tj).0);

        // Contribution to atom j
        let cj =
            lookup_chg_contribution(bt, tj, ti).unwrap_or_else(|| pbci_for(tj).0 - pbci_for(ti).0);

        charges[i] += ci;
        charges[j] += cj;
    }

    // Step 3: formal charge redistribution -- equation 15's `v*sumFormalCharge`
    // term and RDKit's `isDoubleZero(v)` anionic-neighbor-leak adjustment
    // (`AtomTyper.cpp` lines ~3384-3399). These are RDKit's own two
    // MUTUALLY EXCLUSIVE branches on whether the ipso atom's own fcadj is
    // (numerically) zero, not two independent unconditional additions: the
    // leak only fires when fcadj_i == 0, and `v*sumFormalCharge` is
    // multiplied by v so it is a no-op exactly when fcadj_i == 0 anyway --
    // the two forms are equivalent, but computing them as an if/else makes
    // the "leak only applies when v==0" condition explicit rather than
    // relying on a multiply-by-zero coincidence. Both branches read the
    // derived `fchg` (not raw `atom.charge`) for the same reason Step 1
    // does. `ISDOUBLEZERO_EPS` mirrors RDKit's own `isDoubleZero` helper
    // (`Params.h`: `x < 1e-10 && x > -1e-10`), not an arbitrary tolerance.
    const ISDOUBLEZERO_EPS: f64 = 1.0e-10;
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let (_, fcadj_i) = pbci_for(types[i]);
        if fcadj_i.abs() < ISDOUBLEZERO_EPS {
            for b in bonds_of(mol, idx) {
                let nbr_fchg = fchg[b.neighbor.0 as usize];
                if nbr_fchg < 0.0 {
                    let deg = bonds_of(mol, b.neighbor).len() as f64;
                    charges[i] += nbr_fchg / (2.0 * deg);
                }
            }
        } else {
            let sum_fc: f64 = bonds_of(mol, idx)
                .iter()
                .map(|b| fchg[b.neighbor.0 as usize])
                .sum();
            charges[i] += fcadj_i * sum_fc;
        }
    }

    Ok(charges)
}

/// MMFF94 partial charges evaluated in RDKit's exact floating-point order
/// (`MMFFMolProperties::computeMMFFCharges`): per atom, neighbours are
/// visited in `adj` order (RDKit's bond-creation adjacency order) and the
/// charge is `(1 - M*v)*q0 + v*sumFormalCharge + sumPartialCharge`.
pub(crate) fn mmff94_charges_rdkit_order(
    mol: &Molecule,
    adj: &[Vec<usize>],
) -> Result<Vec<f64>, NumericTypeError> {
    let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
    let n = mol.atom_count();
    let rings = chematic_perception::rdkit_sssr_ring_order(mol).unwrap_or_else(|| {
        chematic_perception::find_symmetrized_sssr(mol)
            .rings()
            .to_vec()
    });
    let fchg: Vec<f64> = (0..n)
        .map(|i| mmff_derived_formal_charge(mol, &mmff_mol, &rings, &types, AtomIdx(i as u32)))
        .collect();
    let mut charges = vec![0.0f64; n];
    for i in 0..n {
        let ti = types[i];
        let (pbci_i, v) = pbci_for(ti);
        let m = adj[i].len() as f64;
        let mut q0 = fchg[i];
        if v < 1.0e-10 && v > -1.0e-10 {
            for &j in &adj[i] {
                let nf = fchg[j];
                if nf < 0.0 {
                    q0 += nf / (2.0 * adj[j].len() as f64);
                }
            }
        }
        if ti == 62 {
            for &j in &adj[i] {
                let nf = fchg[j];
                if nf > 0.0 {
                    q0 -= nf / 2.0;
                }
            }
        }
        let mut sum_formal = 0.0f64;
        let mut sum_partial = 0.0f64;
        for &j in &adj[i] {
            let tj = types[j];
            let order = mmff_mol
                .bond_between(AtomIdx(i as u32), AtomIdx(j as u32))
                .map(|(_, b)| b.order)
                .unwrap_or(BondOrder::Single);
            let bt = crate::mmff94_minimizer::bond_type_for(ti, tj, order);
            sum_partial +=
                lookup_chg_contribution(bt, ti, tj).unwrap_or_else(|| pbci_i - pbci_for(tj).0);
            sum_formal += fchg[j];
        }
        charges[i] = (1.0 - m * v) * q0 + v * sum_formal + sum_partial;
    }
    Ok(charges)
}

// ── Tests ────────────────────────────────────────────────────────────────────

