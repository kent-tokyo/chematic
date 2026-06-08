use chematic_core::{AtomIdx, BondOrder, Molecule};

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DREIDINGType {
    C_3,
    C_2,
    C_1,
    C_R,
    N_3,
    N_2,
    N_1,
    N_R,
    O_3,
    O_2,
    O_R,
    S_3,
    S_R,
    P_3,
    H_,
    F_,
    Cl,
    Br,
    I_,
    X_,
}

impl DREIDINGType {
    pub fn symbol(&self) -> &'static str {
        match self {
            Self::C_3 => "C.3",
            Self::C_2 => "C.2",
            Self::C_1 => "C.1",
            Self::C_R => "C.r",
            Self::N_3 => "N.3",
            Self::N_2 => "N.2",
            Self::N_1 => "N.1",
            Self::N_R => "N.r",
            Self::O_3 => "O.3",
            Self::O_2 => "O.2",
            Self::O_R => "O.r",
            Self::S_3 => "S.3",
            Self::S_R => "S.r",
            Self::P_3 => "P.3",
            Self::H_ => "H",
            Self::F_ => "F",
            Self::Cl => "Cl",
            Self::Br => "Br",
            Self::I_ => "I",
            Self::X_ => "X",
        }
    }
}

pub fn assign_dreiding_types(mol: &Molecule) -> Vec<DREIDINGType> {
    let mut types = vec![DREIDINGType::X_; mol.atom_count()];

    for (idx, atom) in mol.atoms() {
        let sym = atom.element.symbol();

        types[idx.0 as usize] = match sym {
            "H" => DREIDINGType::H_,
            "C" => {
                if atom.aromatic {
                    DREIDINGType::C_R
                } else {
                    match highest_bond_order(mol, idx) {
                        BondOrder::Triple => DREIDINGType::C_1,
                        BondOrder::Double => DREIDINGType::C_2,
                        _ => DREIDINGType::C_3,
                    }
                }
            }
            "N" => {
                if atom.aromatic {
                    DREIDINGType::N_R
                } else {
                    match highest_bond_order(mol, idx) {
                        BondOrder::Triple => DREIDINGType::N_1,
                        BondOrder::Double => DREIDINGType::N_2,
                        _ => DREIDINGType::N_3,
                    }
                }
            }
            "O" => {
                if atom.aromatic {
                    DREIDINGType::O_R
                } else {
                    match highest_bond_order(mol, idx) {
                        BondOrder::Double => DREIDINGType::O_2,
                        _ => DREIDINGType::O_3,
                    }
                }
            }
            "S" => {
                if atom.aromatic {
                    DREIDINGType::S_R
                } else {
                    DREIDINGType::S_3
                }
            }
            "P" => DREIDINGType::P_3,
            "F" => DREIDINGType::F_,
            "Cl" => DREIDINGType::Cl,
            "Br" => DREIDINGType::Br,
            "I" => DREIDINGType::I_,
            _ => DREIDINGType::X_,
        };
    }

    types
}

fn highest_bond_order(mol: &Molecule, idx: AtomIdx) -> BondOrder {
    let mut max_order = BondOrder::Single;
    for (_, bond) in mol.neighbors(idx) {
        let bond_entry = &mol.bond(bond);
        match bond_entry.order {
            BondOrder::Triple => return BondOrder::Triple,
            BondOrder::Double => max_order = BondOrder::Double,
            _ => {}
        }
    }
    max_order
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_ethane() {
        let mol = parse("CC").expect("ethane");
        let types = assign_dreiding_types(&mol);
        assert_eq!(types.len(), 2);
        assert_eq!(types[0], DREIDINGType::C_3);
        assert_eq!(types[1], DREIDINGType::C_3);
    }

    #[test]
    fn test_benzene() {
        let mol = parse("c1ccccc1").expect("benzene");
        let types = assign_dreiding_types(&mol);
        assert_eq!(types.len(), 6);
        for &t in &types {
            assert_eq!(t, DREIDINGType::C_R);
        }
    }

    #[test]
    fn test_acetylene() {
        let mol = parse("C#C").expect("acetylene");
        let types = assign_dreiding_types(&mol);
        assert_eq!(types.len(), 2);
        assert_eq!(types[0], DREIDINGType::C_1);
        assert_eq!(types[1], DREIDINGType::C_1);
    }

    #[test]
    fn test_acetone() {
        let mol = parse("CC(=O)C").expect("acetone");
        let types = assign_dreiding_types(&mol);
        assert_eq!(types.len(), 4);
        assert_eq!(types[0], DREIDINGType::C_3);
        assert_eq!(types[2], DREIDINGType::O_2);
        assert_eq!(types[3], DREIDINGType::C_3);
    }
}
