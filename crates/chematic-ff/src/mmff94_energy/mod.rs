//! MMFF94 full energy parameters (Halgren 1996, Tables IV–VII).
//
//! Data extracted verbatim from RDKit `Code/ForceField/MMFF/Params.cpp` (BSD license).
//! Original parameters: Copyright (c) Merck and Co., Inc., 1994, 1995, 1996.
//
//! Units: bond kb in md/Å, angle ka in md·Å/rad², theta0 in degrees,
//! torsion v1/v2/v3 in kcal/mol, vdW alpha_i in Å³.

#![allow(
    clippy::approx_constant,
    clippy::items_after_test_module,
    clippy::type_complexity
)]

/// Bond stretching parameters (Halgren 1996, MMFF.II eq. 1).
/// Energy = (143.9325 × kb / 2) × (ΔR)²  [kcal/mol]
#[derive(Debug, Clone, Copy)]
pub struct BondEnergyParams {
    /// Force constant (md/Å = millidyne/Å); multiply by 143.9325 for kcal/(mol·Å²)
    pub kb: f64,
    /// Equilibrium bond length (Å)
    pub r0: f64,
}

/// Angle bending parameters (Halgren 1996, MMFF.III eq. 2).
/// Energy = (0.043844 × ka / 2) × (Δθ)²  [kcal/mol, Δθ in degrees]
#[derive(Debug, Clone, Copy)]
pub struct AngleEnergyParams {
    /// Force constant (md·Å/rad²); 0.043844 conversion to kcal/mol per deg²
    pub ka: f64,
    /// Equilibrium angle (degrees)
    pub theta0: f64,
}

/// Torsion dihedral parameters (Halgren 1996, MMFF.IV).
/// Energy = (v1/2)(1+cosφ) + (v2/2)(1-cos2φ) + (v3/2)(1+cos3φ)  [kcal/mol]
#[derive(Debug, Clone, Copy)]
pub struct TorsionEnergyParams {
    /// 1-fold Fourier barrier (kcal/mol)
    pub v1: f64,
    /// 2-fold Fourier barrier (kcal/mol)
    pub v2: f64,
    /// 3-fold Fourier barrier (kcal/mol)
    pub v3: f64,
}

/// Van der Waals Slater-Kirkwood parameters (Halgren 1996, MMFF.I Table VII).
#[derive(Debug, Clone, Copy)]
pub struct VdwEnergyParams {
    /// Atomic polarizability (Å³)
    pub alpha_i: f64,
    /// Effective number of electrons (Slater-Kirkwood)
    pub n_i: f64,
    /// Scale factor A_i (for r*_ii = A_i × alpha_i^(1/4))
    pub a_i: f64,
    /// Scale factor G_i (for eps combining rule)
    pub g_i: f64,
    /// Donor/acceptor flag: 0=standard, 1=H-bond donor, 2=H-bond acceptor
    pub da: u8,
}

/// Which mechanism resolved a Bond-stretch or Angle-bend parameter lookup.
///
/// Exposed for diagnostics/tests, not required by production physics: both
/// `DirectTable`/`EquivalentType` and the `Empirical*` variants hand back a
/// real, usable [`BondEnergyParams`]/[`AngleEnergyParams`]. It exists
/// because a table hit alone does not prove correctness -- the Angle
/// `eqLevel` equivalence ladder (issue #227 Stage B) can substitute atom
/// types and land on a real row that is nonetheless the WRONG parameter for
/// the original triple's chemistry, the same failure class as the #236
/// furan collision. Only `mmff94_*_energy_resolved` (issue #227 Stage C)
/// return this; the plain `mmff94_*_energy` functions keep their original
/// `Option<Params>` signature unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mmff94Resolution {
    /// Exact `(type, ti, tj[, tk])` table row, no substitution.
    DirectTable,
    /// Angle only: RDKit's real `eqLevel` canonical-type-substitution ladder
    /// (`Code/ForceField/MMFF/Params.h`), `level` in `{3, 4, 5}` (MMFF.I note
    /// 68's Level 2 is always identity, already covered by `DirectTable`).
    EquivalentType { level: u8 },
    /// Angle only: chematic-specific safety net with no RDKit equivalent --
    /// the same triple re-searched with `angle_type` forced to `0` after the
    /// real eqLevel ladder is exhausted. Predates Stage B; kept so a
    /// correctly-typed triple the specialized angle-type table doesn't cover
    /// doesn't silently drop the term.
    GenericAngleTypeFallback,
    /// Halgren MMFF.V eq. 18-19 empirical bond-stretch rule -- no table row
    /// found at any stage (Bond has no eqLevel ladder at all).
    EmpiricalBond,
    /// Halgren MMFF.V eq. 20 empirical angle-bend rule. Covers both of
    /// RDKit's own sub-cases: no table row found anywhere (`ka`/`theta0`
    /// both derived from scratch), and a table row found with `ka == 0.0`
    /// (RDKit's `isDoubleZero` placeholder -- only `ka` is derived, `theta0`
    /// is reused verbatim from that row, no ring-size override applied;
    /// `getMMFFAngleBendEmpiricalRuleParams`, `AtomTyper.cpp`). Both are the
    /// same RDKit function and the same "no usable ka in the table" case, so
    /// they share this one variant.
    EmpiricalAngle,
    /// Torsion only (issue #227 Phase 1): RDKit's real empirical-rule
    /// cascade generates NO torsion term at all when either central (j-k)
    /// atom has MMFF's `lin` flag -- see
    /// `crate::mmff94_minimizer::torsion_no_term_by_design`'s doc. A missing
    /// lookup with this cause is correct, matching RDKit exactly, not a
    /// coverage gap; distinct from a genuine unresolved miss so callers
    /// (the `include_torsion_oop_in_gate` strict-policy gate, coverage
    /// audits) don't misclassify it as one.
    NoTermByDesign,
}

mod angle;
mod bond;
mod oop_stbn;
mod torsion;
mod vdw;

pub use angle::{MMFF94_ANGLE_ENERGY, mmff94_angle_energy, mmff94_angle_energy_resolved};
pub use bond::{MMFF94_BOND_ENERGY, mmff94_bond_energy, mmff94_bond_energy_resolved};
pub use oop_stbn::{MMFF94_OOP, MMFF94_STBN, mmff94_oop, mmff94_stbn, mmff94_stbn_type_only};
pub use torsion::{MMFF94_TORSION_ENERGY, mmff94_torsion_energy};
pub use vdw::{MMFF94_VDW_ENERGY, mmff94_vdw_combined, mmff94_vdw_energy};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_sizes() {
        assert_eq!(MMFF94_BOND_ENERGY.len(), 493);
        assert_eq!(MMFF94_ANGLE_ENERGY.len(), 2342);
        assert_eq!(MMFF94_TORSION_ENERGY.len(), 926);
        assert_eq!(MMFF94_VDW_ENERGY.len(), 95);
    }

    #[test]
    fn bond_cc_sp3() {
        // C(sp3)-C(sp3): type 1-1, bond_type 0
        let p = mmff94_bond_energy(0, 1, 1).expect("C-C sp3 bond");
        assert!((p.r0 - 1.508).abs() < 0.001, "r0={}", p.r0);
        assert!((p.kb - 4.258).abs() < 0.001, "kb={}", p.kb);
    }

    #[test]
    fn bond_ch_sp3() {
        // C(sp3)-H: type 1-5, bond_type 0
        let p = mmff94_bond_energy(0, 1, 5).expect("C-H sp3 bond");
        assert!((p.r0 - 1.093).abs() < 0.001, "r0={}", p.r0);
        assert!((p.kb - 4.766).abs() < 0.001, "kb={}", p.kb);
    }

    #[test]
    fn bond_symmetric() {
        // Order should not matter
        assert_eq!(
            mmff94_bond_energy(0, 1, 2).map(|p| p.r0),
            mmff94_bond_energy(0, 2, 1).map(|p| p.r0),
        );
    }

    #[test]
    fn angle_ccc_sp3() {
        // C(sp3)-C(sp3)-C(sp3): types 1-1-1, angle_type 0
        let p = mmff94_angle_energy(0, 1, 1, 1).expect("C-C-C sp3 angle");
        assert!((p.theta0 - 109.608).abs() < 0.1, "theta0={}", p.theta0);
        assert!(p.ka > 0.5, "ka={}", p.ka);
    }

    #[test]
    fn angle_symmetric() {
        // (1,1,2) and (2,1,1) should give same params
        let a = mmff94_angle_energy(0, 1, 1, 2).map(|p| p.theta0);
        let b = mmff94_angle_energy(0, 2, 1, 1).map(|p| p.theta0);
        assert_eq!(a, b, "angle lookup not symmetric: {:?} vs {:?}", a, b);
    }

    #[test]
    fn torsion_cccc() {
        // C-C-C-C (butane): tors_type=0, types 1-1-1-1
        // Expected from RDKit API: v1=0.103, v2=0.681, v3=0.332
        let p = mmff94_torsion_energy(0, 1, 1, 1, 1).expect("C-C-C-C torsion");
        assert!((p.v1 - 0.103).abs() < 0.001, "v1={}", p.v1);
        assert!((p.v2 - 0.681).abs() < 0.001, "v2={}", p.v2);
        assert!((p.v3 - 0.332).abs() < 0.001, "v3={}", p.v3);
    }

    #[test]
    fn torsion_hcch() {
        // H-C-C-H: types 5-1-1-5
        let p = mmff94_torsion_energy(0, 5, 1, 1, 5);
        assert!(p.is_some(), "H-C-C-H torsion should be found");
    }

    #[test]
    fn torsion_wildcard_fallback() {
        // An unusual type combo should fall back to wildcard
        // Use types that likely only have wildcard coverage
        // Any result (even zero-barrier) is acceptable
        let _ = mmff94_torsion_energy(0, 99, 1, 1, 99);
        // Just verify it doesn't panic
    }

    #[test]
    fn vdw_carbon_sp3() {
        // Type 1 = CR (sp3 carbon)
        let p = mmff94_vdw_energy(1).expect("sp3 C vdW");
        assert!(p.alpha_i > 0.0, "alpha_i={}", p.alpha_i);
        assert!(p.n_i > 0.0, "n_i={}", p.n_i);
    }

    #[test]
    fn vdw_combined_cc() {
        let (r_star, eps) = mmff94_vdw_combined(1, 1).expect("C-C vdW combined");
        // MMFF94 C(sp3)-C(sp3): r* ≈ 3.9 Å, eps ≈ 0.04 kcal/mol
        assert!(r_star > 2.0 && r_star < 6.0, "r_star={}", r_star);
        assert!(eps > 0.0, "eps={}", eps);
    }

    // --- mmff94_stbn's Dfsb periodic-row fallback (Priority 2B, issue #227) ---
    // Parity fixture: each expected value below is traceable to a specific row
    // in `scripts/mmff94_provenance/rdkit_defaultMMFFDfsb.txt` (the pinned
    // RDKit commit's `defaultMMFFDfsb` table, programmatically extracted, not
    // hand-transcribed -- see PROVENANCE.md's "Stretch-bend" row).

    #[test]
    fn stbn_type_table_hit_takes_priority_over_dfsb() {
        // C(sp3)-C(sp3)-C(sp3): type table has (0,1,1,1)=(0.2060,0.2060)
        // (see mmff94_stbn_type_only's own MMFF94_STBN row). Dfsb's own
        // (row=1,row=1,row=1) row is (0.30,0.30) -- a DIFFERENT value, so
        // this test would catch the Dfsb tier firing before/instead of a
        // real type-table hit.
        let carbon = 6u8;
        let via_type_only = mmff94_stbn_type_only(0, 1, 1, 1).expect("C-C-C stbn (type table)");
        let via_full = mmff94_stbn(0, 1, 1, 1, carbon, carbon, carbon).expect("C-C-C stbn (full)");
        assert_eq!(via_full, via_type_only, "type-table hit must win over Dfsb");
        assert!(
            (via_full.0 - 0.2060).abs() < 1e-9 && (via_full.1 - 0.2060).abs() < 1e-9,
            "got {:?}, expected the real (0.2060, 0.2060) type-table row, not Dfsb's (0.30, 0.30)",
            via_full
        );
    }

    #[test]
    fn stbn_dfsb_fallback_resolves_when_type_table_misses() {
        // F-C-Cl: no exact/generic MMFF94_STBN row exists for this triple
        // (guaranteed by using MMFF type 200, which appears in no real table
        // row, to isolate the Dfsb tier from the type-table entirely).
        // Periodic rows: F=9->1, C=6->1(center), Cl=17->2 -> canonical
        // (1,1,2) -> rdkit_defaultMMFFDfsb.txt row "1  1  2  0.30  0.50".
        assert!(
            mmff94_stbn_type_only(0, 200, 200, 200).is_none(),
            "type 200 must not accidentally exist in MMFF94_STBN"
        );
        let (fluorine, carbon, chlorine) = (9u8, 6u8, 17u8);
        let got = mmff94_stbn(0, 200, 200, 200, fluorine, carbon, chlorine)
            .expect("F-C-Cl should resolve via Dfsb once the type table misses");
        assert!(
            (got.0 - 0.30).abs() < 1e-9 && (got.1 - 0.50).abs() < 1e-9,
            "got {:?}, expected (0.30, 0.50) from rdkit_defaultMMFFDfsb.txt row (1,1,2)",
            got
        );
    }

    #[test]
    fn stbn_dfsb_canonicalizes_row_order_and_swaps_the_result() {
        // Same F-C-Cl chemistry as above but with i/k reversed at the call
        // site (Cl first, F last) -- Dfsb must canonicalize row_i<=row_k
        // internally (matching RDKit's own `MMFFDfsbCollection` swap logic)
        // AND swap the returned (kba_ijk, kba_kji) pair to match the
        // caller's actual atom order, not the table's internal one.
        let (fluorine, carbon, chlorine) = (9u8, 6u8, 17u8);
        let forward = mmff94_stbn(0, 200, 200, 200, fluorine, carbon, chlorine).unwrap();
        let reversed = mmff94_stbn(0, 200, 200, 200, chlorine, carbon, fluorine).unwrap();
        assert_eq!(
            (forward.0, forward.1),
            (reversed.1, reversed.0),
            "reversing the outer atoms must swap kba_ijk/kba_kji, forward={:?} reversed={:?}",
            forward,
            reversed
        );
    }

    #[test]
    fn stbn_dfsb_all_zero_row_is_not_resolved() {
        // H-S-H (or H-P-H): periodic rows H=1->0, S=16->2(center), H=1->0
        // -> canonical (0,2,0) -- rdkit_defaultMMFFDfsb.txt's ONLY all-zero
        // row ("0  2  0  0.00  0.00"). RDKit's own isDoubleZero(kbaIJK) &&
        // isDoubleZero(kbaKJI) check treats this as unresolved, not a real
        // (0.0, 0.0) hit -- replicated here.
        let (hydrogen, sulfur) = (1u8, 16u8);
        assert!(
            mmff94_stbn_type_only(0, 200, 200, 200).is_none(),
            "type 200 must not accidentally exist in MMFF94_STBN"
        );
        assert_eq!(
            mmff94_stbn(0, 200, 200, 200, hydrogen, sulfur, hydrogen),
            None,
            "the Dfsb table's one all-zero row must not be reported as resolved"
        );
    }

    #[test]
    fn stbn_dfsb_out_of_table_combination_is_not_resolved() {
        // No real angle center has periodic row 0 (H/He) in MMFF94 chemistry
        // in practice, but as a structural check: row_j=0 never appears as
        // the *second* column in any Dfsb row, so any triple canonicalizing
        // to (_, 0, _) must stay unresolved rather than silently matching an
        // unrelated row.
        let hydrogen = 1u8;
        assert_eq!(
            mmff94_stbn(0, 200, 200, 200, hydrogen, hydrogen, hydrogen),
            None,
            "row_j=0 never appears in MMFF94_DFSB -- H-H-H must stay unresolved"
        );
    }
}
