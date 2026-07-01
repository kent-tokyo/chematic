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

mod angle;
mod bond;
mod oop_stbn;
mod torsion;
mod vdw;

pub use angle::{MMFF94_ANGLE_ENERGY, mmff94_angle_energy};
pub use bond::{MMFF94_BOND_ENERGY, mmff94_bond_energy};
pub use oop_stbn::{MMFF94_OOP, MMFF94_STBN, mmff94_oop, mmff94_stbn};
pub use torsion::{MMFF94_TORSION_ENERGY, mmff94_torsion_energy};
pub use vdw::{MMFF94_VDW_ENERGY, mmff94_vdw_combined, mmff94_vdw_energy};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_sizes() {
        assert_eq!(MMFF94_BOND_ENERGY.len(), 493);
        assert_eq!(MMFF94_ANGLE_ENERGY.len(), 2245);
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
}
