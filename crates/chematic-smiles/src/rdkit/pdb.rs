//! `MolToPDBBlock` (RDKit 2026.03.1 `PDBWriter.cpp`, default flavor 0) for
//! a molecule without PDB residue information.

use std::collections::HashMap;
use std::fmt::Write as _;

use super::RdkitSmilesError;
use super::kekulize::kekulize_ranked;
use super::mol::{BondType, Mol};
use super::periodic;
use super::rank::rank_mol_atoms;

/// `GetDefaultAtomNumber`.
fn default_atom_number(anum: u32, elem: &mut HashMap<u32, u32>) -> String {
    let mut ret = *b"  ";
    match elem.get_mut(&anum) {
        None => {
            elem.insert(anum, 1);
            ret[0] = b'1';
        }
        Some(count) => {
            *count += 1;
            let tmp = *count;
            if tmp < 10 {
                ret[0] = b'0' + tmp as u8;
            } else if tmp < 100 {
                ret[0] = b'0' + (tmp / 10) as u8;
                ret[1] = b'0' + (tmp % 10) as u8;
            } else if tmp < 360 {
                ret[0] = b'A' + ((tmp - 100) / 10) as u8;
                ret[1] = b'0' + ((tmp - 100) % 10) as u8;
            } else if tmp < 1036 {
                ret[0] = b'A' + ((tmp - 360) / 26) as u8;
                ret[1] = b'A' + ((tmp - 360) % 26) as u8;
            }
        }
    }
    String::from_utf8_lossy(&ret).into_owned()
}

/// `MolToPDBBlock(mol, confId=-1, flavor=0)`; `coords` stands for the
/// conformer (none: zero coordinates).
pub(crate) fn mol_to_pdb_block(
    mol: &Mol,
    coords: Option<&[[f64; 3]]>,
) -> Result<String, RdkitSmilesError> {
    if let Some(c) = coords
        && c.len() != mol.atoms.len()
    {
        return Err(RdkitSmilesError::Unsupported(format!(
            "expected {} coordinates, got {}",
            mol.atoms.len(),
            c.len()
        )));
    }
    let mut t = mol.clone();
    if !t.bonds.is_empty() {
        let ranks = rank_mol_atoms(&t);
        kekulize_ranked(&mut t, Some(&ranks))?;
    }
    let mut res = String::new();
    let mut elem = HashMap::new();
    for (i, atom) in t.atoms.iter().enumerate() {
        let symb = periodic::symbol(atom.anum).as_bytes();
        let (at1, at2) = match symb.len() {
            0 => (b' ', b'X'),
            1 => (b' ', symb[0]),
            _ => (symb[0], symb[1].to_ascii_uppercase()),
        };
        let (at1, at2) = (at1 as char, at2 as char);
        let atnum = default_atom_number(atom.anum, &mut elem);
        let _ = write!(res, "HETATM{:5} {at1}{at2}{atnum} UNL     1    ", i + 1);
        match coords {
            Some(c) => {
                let p = c[i];
                let _ = write!(res, "{:8.3}{:8.3}{:8.3}", p[0], p[1], p[2]);
            }
            None => res.push_str("   0.000   0.000   0.000"),
        }
        res.push_str("  1.00  0.00          ");
        res.push(at1);
        res.push(at2);
        let charge = atom.charge;
        if charge > 0 && charge < 10 {
            res.push((b'0' + charge as u8) as char);
            res.push('+');
        } else if charge < 0 && charge > -10 {
            res.push((b'0' + (-charge) as u8) as char);
            res.push('-');
        } else {
            res.push_str("  ");
        }
        res.push('\n');
    }
    // GetPDBBondLines with all=true, both=false, mult=true.
    for a in 0..t.atoms.len() {
        let src = a + 1;
        let mut v = Vec::new();
        for &b in &t.atom_bonds[a] {
            let bond = &t.bonds[b];
            let dst = bond.other(a) + 1;
            if dst < src {
                continue;
            }
            let times = match bond.bt {
                BondType::Double => 2,
                BondType::Triple => 3,
                BondType::Quadruple => 4,
                _ => 1,
            };
            for _ in 0..times {
                v.push(dst);
            }
        }
        if v.is_empty() {
            continue;
        }
        v.sort_unstable();
        for (i, d) in v.iter().enumerate() {
            if i & 3 == 0 {
                if i != 0 {
                    res.push('\n');
                }
                let _ = write!(res, "CONECT{src:5}");
            }
            let _ = write!(res, "{d:5}");
        }
        res.push('\n');
    }
    res.push_str("END\n");
    Ok(res)
}
