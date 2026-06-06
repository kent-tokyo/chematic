use chematic_core::BondOrder;
use std::f64::consts::PI;

use crate::dreiding::DREIDINGType;

pub fn dreiding_vdw(t: DREIDINGType) -> (f64, f64) {
    let (r0, depth) = match t {
        DREIDINGType::H_ => (2.646, 0.044),
        DREIDINGType::C_3 => (3.473, 0.066),
        DREIDINGType::C_2 => (3.473, 0.066),
        DREIDINGType::C_1 => (3.473, 0.066),
        DREIDINGType::C_R => (3.473, 0.066),
        DREIDINGType::N_3 => (3.243, 0.068),
        DREIDINGType::N_2 => (3.243, 0.068),
        DREIDINGType::N_1 => (3.243, 0.068),
        DREIDINGType::N_R => (3.243, 0.068),
        DREIDINGType::O_3 => (3.144, 0.052),
        DREIDINGType::O_2 => (3.144, 0.052),
        DREIDINGType::O_R => (3.144, 0.052),
        DREIDINGType::S_3 => (3.952, 0.274),
        DREIDINGType::S_R => (3.952, 0.274),
        DREIDINGType::P_3 => (3.741, 0.210),
        DREIDINGType::F_ => (3.048, 0.050),
        DREIDINGType::Cl => (3.516, 0.227),
        DREIDINGType::Br => (3.641, 0.389),
        DREIDINGType::I_ => (3.984, 0.680),
        DREIDINGType::X_ => (3.473, 0.066),
    };
    (r0, depth)
}

pub fn dreiding_bond_len(t1: DREIDINGType, t2: DREIDINGType, order: BondOrder) -> f64 {
    let key = (
        dreiding_type_order(t1),
        dreiding_type_order(t2),
        bond_order_value(order),
    );

    match key {
        (1, 1, 1) => 1.090,
        (1, 6, 1) => 1.090,
        (1, 7, 1) => 1.010,
        (1, 8, 1) => 0.960,
        (1, 16, 1) => 1.336,
        (1, 17, 1) => 0.980,
        (1, 35, 1) => 1.084,
        (1, 53, 1) => 1.008,

        (6, 6, 1) => 1.540,
        (6, 6, 2) => 1.340,
        (6, 6, 3) => 1.204,
        (6, 6, 4) => 1.395,
        (6, 7, 1) => 1.469,
        (6, 7, 2) => 1.279,
        (6, 7, 3) => 1.158,
        (6, 7, 4) => 1.340,
        (6, 8, 1) => 1.427,
        (6, 8, 2) => 1.217,
        (6, 8, 4) => 1.370,
        (6, 16, 1) => 1.820,
        (6, 17, 1) => 1.770,
        (6, 35, 1) => 1.940,
        (6, 53, 1) => 2.140,

        (7, 7, 1) => 1.450,
        (7, 7, 2) => 1.225,
        (7, 7, 3) => 1.098,
        (7, 8, 1) => 1.400,
        (7, 8, 2) => 1.218,
        (7, 16, 1) => 1.754,
        (7, 17, 1) => 1.690,

        (8, 8, 1) => 1.480,
        (8, 8, 2) => 1.208,
        (8, 16, 1) => 1.680,
        (8, 17, 1) => 1.611,

        (16, 16, 1) => 2.048,
        (16, 16, 2) => 1.549,

        _ => match order {
            BondOrder::Triple => 1.20,
            BondOrder::Double => 1.34,
            BondOrder::Aromatic => 1.40,
            _ => 1.54,
        },
    }
}

pub fn dreiding_angle(t: DREIDINGType) -> f64 {
    let deg = match t {
        DREIDINGType::C_3 => 109.47,
        DREIDINGType::C_2 => 120.0,
        DREIDINGType::C_1 => 180.0,
        DREIDINGType::C_R => 120.0,
        DREIDINGType::N_3 => 106.7,
        DREIDINGType::N_2 => 120.0,
        DREIDINGType::N_1 => 180.0,
        DREIDINGType::N_R => 120.0,
        DREIDINGType::O_3 => 104.51,
        DREIDINGType::O_2 => 120.0,
        DREIDINGType::O_R => 120.0,
        DREIDINGType::S_3 => 99.0,
        DREIDINGType::S_R => 120.0,
        DREIDINGType::P_3 => 93.0,
        DREIDINGType::H_ => 109.47,
        DREIDINGType::F_ => 109.47,
        DREIDINGType::Cl => 109.47,
        DREIDINGType::Br => 109.47,
        DREIDINGType::I_ => 109.47,
        DREIDINGType::X_ => 109.47,
    };
    deg * PI / 180.0
}

pub fn dreiding_torsion_barrier(_t1: DREIDINGType, _t2: DREIDINGType) -> f64 {
    2.0
}

fn dreiding_type_order(t: DREIDINGType) -> u8 {
    match t {
        DREIDINGType::H_ => 1,
        DREIDINGType::C_3 | DREIDINGType::C_2 | DREIDINGType::C_1 | DREIDINGType::C_R => 6,
        DREIDINGType::N_3 | DREIDINGType::N_2 | DREIDINGType::N_1 | DREIDINGType::N_R => 7,
        DREIDINGType::O_3 | DREIDINGType::O_2 | DREIDINGType::O_R => 8,
        DREIDINGType::S_3 | DREIDINGType::S_R => 16,
        DREIDINGType::P_3 => 15,
        DREIDINGType::F_ => 9,
        DREIDINGType::Cl => 17,
        DREIDINGType::Br => 35,
        DREIDINGType::I_ => 53,
        DREIDINGType::X_ => 0,
    }
}

fn bond_order_value(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 4,
        _ => 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vdw_carbon() {
        let (r0, depth) = dreiding_vdw(DREIDINGType::C_3);
        assert!((r0 - 3.473).abs() < 0.001);
        assert!((depth - 0.066).abs() < 0.001);
    }

    #[test]
    fn test_bond_len_cc_single() {
        let r = dreiding_bond_len(DREIDINGType::C_3, DREIDINGType::C_3, BondOrder::Single);
        assert!((r - 1.540).abs() < 0.001);
    }

    #[test]
    fn test_bond_len_cc_double() {
        let r = dreiding_bond_len(DREIDINGType::C_2, DREIDINGType::C_2, BondOrder::Double);
        assert!((r - 1.340).abs() < 0.001);
    }

    #[test]
    fn test_angle_sp3() {
        let theta = dreiding_angle(DREIDINGType::C_3);
        let expected = 109.47 * PI / 180.0;
        assert!((theta - expected).abs() < 0.01);
    }

    #[test]
    fn test_angle_sp2() {
        let theta = dreiding_angle(DREIDINGType::C_2);
        let expected = 120.0 * PI / 180.0;
        assert!((theta - expected).abs() < 0.01);
    }
}
