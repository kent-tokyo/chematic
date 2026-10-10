use crate::MMFF94Type::{self, *};
use crate::mmff94_params::*;
use chematic_core::BondOrder;

const TYPES: &[MMFF94Type] = &[
    C_sp3,
    C_sp2_Alkene,
    C_sp_Alkyne,
    C_Aromatic,
    C_Carbonyl,
    C_Carboxylic,
    C_Carbamate,
    C_Ester,
    C_Amide,
    C_Imide,
    C_CarbamideN,
    N_sp3_Amine,
    N_sp3_AmineAromatic,
    N_sp2_Imine,
    N_sp2_Aromatic,
    N_sp2_Carbonyl,
    N_sp_Nitrile,
    N_Amide,
    N_Carbamate,
    N_Ester,
    N_Imide,
    N_Aromatic_5ring,
    N_Aromatic_6ring,
    N_Aromatic_Pyridine,
    N_Aromatic_Pyrrole,
    N_Aromatic_Imidazole,
    N_Aromatic_Triazole,
    N_Aromatic_Tetrazole,
    N_Aromatic_Pyrimidine,
    N_Aromatic_Pyrazine,
    O_Alcohol,
    O_Phenol,
    O_Ether,
    O_Carbonyl,
    O_Carboxylic,
    O_Carbamate,
    O_Ester,
    O_Amide,
    O_Imide,
    O_CarbamideN,
    O_Sulfoxide,
    O_Sulfone,
    S_Thiol,
    S_Thioether,
    S_Disulfide,
    S_Sulfoxide,
    S_Sulfone,
    S_Aromatic,
    P_sp3,
    P_Oxide,
    Si_sp3,
    Si_sp2,
    F,
    Cl,
    Br,
    I,
    H_Carbon,
    H_Nitrogen,
    H_Oxygen,
    H_Sulfur,
    H_Halogen,
    H_Aromatic,
    Generic,
];

#[test]
fn parameter_tables_have_physical_domains_and_symmetric_bond_lookups() {
    let mut supported_bonds = 0;
    let mut supported_angles = 0;
    for &a in TYPES {
        let vdw = mmff94_vdw_params(a);
        assert!(vdw.r_star.is_finite() && vdw.r_star > 0., "{a:?}");
        assert!(vdw.epsilon.is_finite() && vdw.epsilon > 0., "{a:?}");
        assert!(mmff94_charge_params(a).charge.is_finite());
        for &b in TYPES {
            let forward = mmff94_bond_params(a, b, BondOrder::Single);
            let backward = mmff94_bond_params(b, a, BondOrder::Single);
            assert_eq!(forward.is_some(), backward.is_some());
            if let Some(p) = forward {
                supported_bonds += 1;
                let reverse = backward.unwrap();
                assert_eq!((p.r0, p.kb), (reverse.r0, reverse.kb));
                assert!(p.r0.is_finite() && p.r0 > 0.);
                assert!(p.kb.is_finite() && p.kb > 0.);
            }
            let dipole = mmff94_bond_dipole(a, b, 1.5);
            let reverse = mmff94_bond_dipole(b, a, 1.5);
            assert_eq!(dipole.dipole_moment, reverse.dipole_moment);
            assert_eq!(
                dipole.electronegativity_diff,
                -reverse.electronegativity_diff
            );
            assert!((0. ..=4.).contains(&dipole.dipole_moment));
            assert!(mmff94_bond_dipole(a, b, 10.).dipole_moment >= dipole.dipole_moment);
            let scaling = mmff94_electrostatic_scaling_1_4(a, b);
            assert!((0. ..=1.).contains(&scaling.scale_1_4));
            assert!(scaling.scale_dielectric >= 1.);
            for &c in TYPES {
                if let Some(p) = mmff94_angle_params(a, b, c) {
                    supported_angles += 1;
                    assert!(
                        p.theta0.is_finite() && p.theta0 > 0. && p.theta0 <= std::f64::consts::PI
                    );
                    assert!(p.ka.is_finite() && p.ka > 0.);
                }
            }
        }
    }
    assert!(supported_bonds > 20);
    assert!(supported_angles > 10);
}
