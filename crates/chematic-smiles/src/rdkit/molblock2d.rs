//! `Chem.MolToMolBlock(mol)` (RDKit 2026.03.1 `MolFileWriter.cpp`,
//! V2000 path) for a molecule carrying the 2D coordinates of
//! [`super::depict`]: `prepareMol` (kekulization), the atom and bond blocks
//! with `Chirality::pickBondsToWedge` / `GetMolFileBondStereoInfo`
//! (`WedgeBonds.cpp`, `Chirality.cpp`) and the `M  CHG`/`RAD`/`ISO` lines.

use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::fmt::Write as _;

use super::RdkitSmilesError;
use super::kekulize::kekulize_ranked;
use super::mol::{BondDir, BondStereo, BondType, ChiralTag, Mol};
use super::periodic;
use super::rank::rank_mol_atoms;

const ZERO_TOLERANCE: f64 = 1.0e-16;
const MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO: usize = 8;
const NO_NBRS: i32 = 100;

fn unsupported(what: &str) -> RdkitSmilesError {
    RdkitSmilesError::Unsupported(format!("MOL block: {what}"))
}

/// `Point3D::directionVector` of two in-plane points (`z = 0`).
fn direction_vector(from: [f64; 2], to: [f64; 2]) -> Option<[f64; 3]> {
    let mut r = [to[0] - from[0], to[1] - from[1], 0.0];
    let l = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
    if l < ZERO_TOLERANCE {
        return None;
    }
    r[0] /= l;
    r[1] /= l;
    r[2] /= l;
    Some(r)
}

/// `Point3D::signedAngleTo`.
fn signed_angle_to(a: [f64; 3], b: [f64; 3]) -> f64 {
    let lsq = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    let mut dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    dot /= lsq.sqrt();
    let mut res = if dot <= -1.0 {
        PI
    } else if dot >= 1.0 {
        0.0
    } else {
        dot.acos()
    };
    if (a[0] * b[1] - a[1] * b[0]) < -ZERO_TOLERANCE {
        res = 2.0 * PI - res;
    }
    res
}

fn is_tetrahedral(tag: ChiralTag) -> bool {
    matches!(tag, ChiralTag::Cw | ChiralTag::Ccw)
}

/// `Chirality::detail::countChiralNbrs` (no bonds are wedged yet).
fn count_chiral_nbrs(mol: &Mol) -> (bool, Vec<i32>) {
    let mut n = vec![NO_NBRS; mol.atoms.len()];
    let mut chi_nbrs = false;
    for a in 0..mol.atoms.len() {
        if !is_tetrahedral(mol.atoms[a].chiral) {
            continue;
        }
        n[a] = 0;
        chi_nbrs = true;
        for nb in mol.nbrs(a) {
            if mol.atoms[nb].anum == 1 {
                n[a] -= 10;
                continue;
            }
            if !is_tetrahedral(mol.atoms[nb].chiral) {
                continue;
            }
            n[a] -= 1;
        }
    }
    (chi_nbrs, n)
}

/// `getDoubleBondPresence`.
fn double_bond_presence(mol: &Mol, a: usize) -> (i32, i32, i32) {
    let (mut has, mut known, mut any) = (0, 0, 0);
    for &b in &mol.atom_bonds[a] {
        let bond = &mol.bonds[b];
        if bond.bt == BondType::Double {
            has += 1;
            if bond.stereo == BondStereo::Any {
                any += 1;
            } else if bond.stereo > BondStereo::Any {
                known += 1;
            }
        }
    }
    (has, known, any)
}

/// `Chirality::detail::pickBondToWedge`.
fn pick_bond_to_wedge(
    mol: &Mol,
    atom: usize,
    n_chiral_nbrs: &[i32],
    wedge: &BTreeMap<usize, usize>,
) -> Option<usize> {
    let ri = mol.ring_info();
    let mut scores: Vec<(i32, usize)> = Vec::new();
    for &bid in &mol.atom_bonds[atom] {
        let bond = &mol.bonds[bid];
        if bond.bt != BondType::Single {
            continue;
        }
        if wedge.contains_key(&bid) {
            continue;
        }
        let o = bond.other(atom);
        let oatom = &mol.atoms[o];
        if oatom.anum == 1 {
            scores.push((-1_000_000, bid));
            continue;
        }
        let mut score = oatom.anum as i32
            + 100 * mol.degree(o) as i32
            + 1000 * i32::from(oatom.chiral != ChiralTag::Unspecified);
        if n_chiral_nbrs[o] < NO_NBRS {
            score -= 100_000 * n_chiral_nbrs[o];
        }
        score += 10_000 * ri.num_atom_rings(o) as i32;
        score += 20_000 * ri.num_bond_rings(bid) as i32;
        let (has, known, any) = double_bond_presence(mol, o);
        score += 11_000 * has;
        score += 12_000 * known;
        score += 23_000 * any;
        scores.push((score, bid));
    }
    scores.iter().min().map(|&(_, b)| b)
}

/// `Chirality::pickBondsToWedge`: bond index -> the chiral atom it wedges.
fn pick_bonds_to_wedge(mol: &Mol) -> BTreeMap<usize, usize> {
    let mut indices: Vec<usize> = (0..mol.atoms.len()).collect();
    let (chi_nbrs, n_chiral) = count_chiral_nbrs(mol);
    if chi_nbrs {
        chematic_perception::rdkit_sssr_order::libstdcxx_sort(&mut indices, |&i1, &i2| {
            n_chiral[i1] < n_chiral[i2]
        });
    }
    let mut wedge = BTreeMap::new();
    for idx in indices {
        if n_chiral[idx] > NO_NBRS {
            continue;
        }
        if !is_tetrahedral(mol.atoms[idx].chiral) {
            break;
        }
        if let Some(b) = pick_bond_to_wedge(mol, idx, &n_chiral, &wedge) {
            wedge.insert(b, idx);
        }
    }
    wedge
}

/// `Chirality::detail::determineBondWedgeState(bond, fromAtomIdx, conf)`;
/// `None` keeps the bond's own direction.
fn determine_bond_wedge_state(mol: &Mol, bid: usize, from: usize, xy: &[[f64; 2]]) -> Option<u8> {
    let bond = &mol.bonds[bid];
    let (atom, bond_atom) = if bond.begin == from {
        (bond.begin, bond.end)
    } else {
        (bond.end, bond.begin)
    };
    let chiral = mol.atoms[atom].chiral;
    let center = xy[atom];
    let ref_vect = direction_vector(center, xy[bond_atom])?;
    let mut idxs: Vec<usize> = vec![bid];
    let mut angles: Vec<f64> = vec![0.0];
    for &nb in &mol.atom_bonds[atom] {
        if nb == bid {
            continue;
        }
        let other = mol.bonds[nb].other(atom);
        let tmp_vect = direction_vector(center, xy[other])?;
        let mut angle = signed_angle_to(ref_vect, tmp_vect);
        if angle < 0.0 {
            angle += 2.0 * PI;
        }
        let mut pos = 0;
        while pos < angles.len() && angle > angles[pos] {
            pos += 1;
        }
        angles.insert(pos, angle);
        idxs.insert(pos, nb);
    }
    let mut n_swaps = usize::from(mol.perturbation_is_odd(atom, &idxs));
    if angles.len() == 3 {
        let angle1 = angles[1];
        let angle2 = angles[2];
        let angle_tol = PI * 1.9 / 180.0;
        if angle2 - angle1 >= (PI - angle_tol) {
            n_swaps += 1;
        }
    }
    let odd = n_swaps % 2 == 1;
    // 1: BEGINWEDGE, 6: BEGINDASH
    Some(match (chiral, odd) {
        (ChiralTag::Ccw, true) | (ChiralTag::Cw, false) => 6,
        _ => 1,
    })
}

/// `Chirality::detail::isBondPotentialStereoBond`.
pub(crate) fn is_bond_potential_stereo_bond(mol: &Mol, b: usize) -> bool {
    let bond = &mol.bonds[b];
    if bond.bt != BondType::Double {
        return false;
    }
    let total_hs = |a: usize| {
        mol.total_num_hs(a) as usize + mol.nbrs(a).filter(|&n| mol.atoms[n].anum == 1).count()
    };
    let beg_deg = mol.total_degree(bond.begin);
    let end_deg = mol.total_degree(bond.end);
    if beg_deg > 1
        && beg_deg < 4
        && end_deg > 1
        && end_deg < 4
        && total_hs(bond.begin) < 2
        && total_hs(bond.end) < 2
    {
        let ri = mol.ring_info();
        !ri.bond_rings
            .iter()
            .any(|r| r.len() < MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO && r.contains(&b))
    } else {
        false
    }
}

/// `Chirality::canBeStereoBond` (legacy perception: `_CIPRank`).
fn can_be_stereo_bond(mol: &Mol, b: usize, cip: &[u32]) -> bool {
    let bond = &mol.bonds[b];
    if !matches!(bond.bt, BondType::Double | BondType::Aromatic) {
        return false;
    }
    for atom in [bond.begin, bond.end] {
        let mut ranks: Vec<i64> = Vec::new();
        for &nb in &mol.atom_bonds[atom] {
            if nb == b {
                continue;
            }
            let nbond = &mol.bonds[nb];
            if nbond.bt == BondType::Single {
                if matches!(nbond.dir, BondDir::EndUpRight | BondDir::EndDownRight) {
                    return false;
                }
                let other = nbond.other(atom);
                let rank = if cip.len() == mol.atoms.len() {
                    i64::from(cip[other])
                } else {
                    -1
                };
                if rank >= 0 {
                    if ranks.contains(&rank) {
                        return false;
                    }
                    ranks.push(rank);
                }
            }
        }
    }
    true
}

/// `Chirality::shouldBeACrossedBond`.
fn should_be_a_crossed_bond(mol: &Mol, b: usize, cip: &[u32]) -> bool {
    let bond = &mol.bonds[b];
    if bond.stereo == BondStereo::Any {
        return true;
    }
    if bond.stereo != BondStereo::None {
        return false;
    }
    if !is_bond_potential_stereo_bond(mol, b) {
        return false;
    }
    let (beg, end) = (bond.begin, bond.end);
    if mol.degree(beg) > 1
        && mol.degree(end) > 1
        && mol.total_valence(beg) - mol.total_degree(beg) as i32 == 1
        && mol.total_valence(end) - mol.total_degree(end) as i32 == 1
        && can_be_stereo_bond(mol, b, cip)
    {
        return true;
    }
    false
}

/// `hasNonDefaultValence`; `bare_dummy` marks a dummy atom written as a
/// bare `*` (RDKit leaves its `noImplicit` flag unset).
fn has_non_default_valence(
    mol: &Mol,
    a: usize,
    bare_dummy: bool,
) -> Result<bool, RdkitSmilesError> {
    let atom = &mol.atoms[a];
    if atom.radicals != 0 {
        return Ok(true);
    }
    const ORGANIC: [u32; 11] = [0, 5, 6, 7, 8, 9, 15, 16, 17, 35, 53];
    if atom.anum == 1 || ORGANIC.contains(&atom.anum) {
        let no_implicit = atom.no_implicit && !bare_dummy;
        if !no_implicit {
            return Ok(false);
        }
        let eff = atom.anum as i32 - atom.charge;
        if eff < 0 || eff > periodic::max_atomic_num() {
            return Err(unsupported("effective atomic number out of range"));
        }
        return Ok(atom.explicit_valence != periodic::default_valence(eff as u32));
    }
    Ok(true)
}

/// `Chem.MolToMolBlock(mol)` (V2000) for `mol` as `MolFromSmiles` leaves
/// it, with conformer coordinates `xy` and `_CIPRank` values `cip`.
/// `bare_dummies[a]` marks dummy atoms read from a bare `*`.
pub(crate) fn mol_block_2d(
    mol: &Mol,
    cip: &[u32],
    xy: &[[f64; 2]],
    bare_dummies: &[bool],
) -> Result<String, RdkitSmilesError> {
    let n_atoms = mol.atoms.len();
    let n_bonds = mol.bonds.len();
    // outputMolToMolBlock picks V3000 for dative bonds, more than 999 atoms
    // or bonds, or coordinates outside the V2000 field widths.
    let v3000 = n_atoms > 999
        || n_bonds > 999
        || mol.bonds.iter().any(|b| b.bt == BondType::Dative)
        || xy.iter().any(|p| {
            p[0] >= 100_000.0 || p[0] <= -10_000.0 || p[1] >= 100_000.0 || p[1] <= -10_000.0
        });
    // prepareMol: remember aromatic bonds, then Kekulize.
    let mut t = mol.clone();
    let mut was_aromatic = vec![false; n_bonds];
    if n_bonds > 0 {
        for (i, b) in t.bonds.iter().enumerate() {
            was_aromatic[i] = b.aromatic;
        }
        let ranks = rank_mol_atoms(&t);
        kekulize_ranked(&mut t, Some(&ranks))?;
    }

    let mut out = String::new();
    out.push('\n');
    out.push_str("     RDKit          2D\n");
    out.push('\n');
    if v3000 {
        out.push_str("  0  0  0  0  0  0  0  0  0  0999 V3000\n");
        out.push_str("M  V30 BEGIN CTAB\n");
        let _ = writeln!(out, "M  V30 COUNTS {n_atoms} {n_bonds} 0 0 0");
        out.push_str("M  V30 BEGIN ATOM\n");
    } else {
        let _ = writeln!(
            out,
            "{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}999 V2000",
            n_atoms, n_bonds, 0, 0, 0, 0, 0, 0, 0, 0
        );
    }
    for a in 0..n_atoms {
        let atom = &t.atoms[a];
        let map = atom.map.unwrap_or(0);
        let mut tot_valence = 0;
        if has_non_default_valence(&t, a, bare_dummies.get(a).copied().unwrap_or(false))? {
            if t.total_degree(a) == 0 {
                tot_valence = 15;
            } else {
                tot_valence = t.total_valence(a) % 15;
            }
        }
        let mut symbol = if atom.anum == 0 {
            "R".to_string()
        } else {
            periodic::symbol(atom.anum).to_string()
        };
        let [x, y] = xy[a];
        if v3000 {
            // GetV3000MolFileAtomLine (precision 6, no parity from SMILES).
            let _ = write!(
                out,
                "M  V30 {} {symbol} {x:.6} {y:.6} {z:.6} {map}",
                a + 1,
                z = 0.0f64
            );
            if atom.charge != 0 {
                let _ = write!(out, " CHG={}", atom.charge);
            }
            if atom.isotope != 0 {
                let _ = write!(out, " MASS={}", atom.isotope);
            }
            if atom.radicals != 0 && t.total_degree(a) != 0 {
                let code = if atom.radicals % 2 == 1 { 2 } else { 3 };
                let _ = write!(out, " RAD={code}");
            }
            if tot_valence != 0 {
                if tot_valence == 15 {
                    out.push_str(" VAL=-1");
                } else {
                    let _ = write!(out, " VAL={tot_valence}");
                }
            }
            out.push('\n');
            continue;
        }
        while symbol.len() < 3 {
            symbol.push(' ');
        }
        let _ = writeln!(
            out,
            "{x:10.4}{y:10.4}{z:10.4} {symbol:>3}{:2}{:3}{:3}{:3}{:3}{tot_valence:3}  0{:3}{:3}{map:3}{:3}{:3}",
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            z = 0.0f64
        );
    }
    if v3000 {
        out.push_str("M  V30 END ATOM\n");
        if n_bonds > 0 {
            out.push_str("M  V30 BEGIN BOND\n");
        }
    }
    let wedge = pick_bonds_to_wedge(&t);
    for (bid, bond) in t.bonds.iter().enumerate() {
        let mut dir_code = 0u8;
        let mut reverse = false;
        if matches!(bond.bt, BondType::Single | BondType::Aromatic) {
            if let Some(&from) = wedge.get(&bid)
                && let Some(code) = determine_bond_wedge_state(&t, bid, from, xy)
            {
                dir_code = code;
                if from != bond.begin {
                    reverse = true;
                }
            }
        } else if bond.bt == BondType::Double && should_be_a_crossed_bond(&t, bid, cip) {
            dir_code = 3;
        }
        if was_aromatic[bid] && dir_code == 3 {
            dir_code = 0;
        }
        let symbol = match bond.bt {
            BondType::Single | BondType::Double if bond.aromatic => 4,
            BondType::Single => 1,
            BondType::Double => 2,
            BondType::Triple => 3,
            BondType::Aromatic => 4,
            BondType::Dative => 9,
            BondType::Quadruple => 0,
        };
        let (b1, b2) = if reverse {
            (bond.end, bond.begin)
        } else {
            (bond.begin, bond.end)
        };
        if v3000 {
            let _ = write!(out, "M  V30 {} {symbol} {} {}", bid + 1, b1 + 1, b2 + 1);
            let cfg = match dir_code {
                1 => 1,
                3 | 4 => 2,
                6 => 3,
                _ => 0,
            };
            if dir_code != 0 {
                let _ = write!(out, " CFG={cfg}");
            }
            out.push('\n');
            continue;
        }
        let _ = writeln!(out, "{:3}{:3}{symbol:3} {dir_code:2}", b1 + 1, b2 + 1);
    }
    if v3000 {
        if n_bonds > 0 {
            out.push_str("M  V30 END BOND\n");
        }
        out.push_str("M  V30 END CTAB\nM  END\n");
        return Ok(out);
    }
    // GetMolFileChargeInfo: full lines of 8 are flushed as they fill.
    let mut chg = String::new();
    let mut rad = String::new();
    let mut iso = String::new();
    let (mut n_chg, mut n_rad, mut n_iso) = (0, 0, 0);
    for a in 0..n_atoms {
        let atom = &t.atoms[a];
        if atom.charge != 0 {
            n_chg += 1;
            let _ = write!(chg, " {:3} {:3}", a + 1, atom.charge);
            if n_chg == 8 {
                let _ = writeln!(out, "M  CHG{n_chg:3}{chg}");
                chg.clear();
                n_chg = 0;
            }
        }
        if atom.radicals != 0 && t.total_degree(a) != 0 {
            n_rad += 1;
            let code = if atom.radicals % 2 == 1 { 2 } else { 3 };
            let _ = write!(rad, " {:3} {:3}", a + 1, code);
            if n_rad == 8 {
                let _ = writeln!(out, "M  RAD{n_rad:3}{rad}");
                rad.clear();
                n_rad = 0;
            }
        }
        if atom.isotope != 0 {
            n_iso += 1;
            let _ = write!(iso, " {:3} {:3}", a + 1, atom.isotope);
            if n_iso == 8 {
                let _ = writeln!(out, "M  ISO{n_iso:3}{iso}");
                iso.clear();
                n_iso = 0;
            }
        }
    }
    if n_chg > 0 {
        let _ = writeln!(out, "M  CHG{n_chg:3}{chg}");
    }
    if n_rad > 0 {
        let _ = writeln!(out, "M  RAD{n_rad:3}{rad}");
    }
    if n_iso > 0 {
        let _ = writeln!(out, "M  ISO{n_iso:3}{iso}");
    }
    out.push_str("M  END\n");
    Ok(out)
}
