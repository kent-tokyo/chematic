//! RDKit's `MolFromMol2Block` (RDKit 2026.03.1
//! `FileParsers/Mol2FileParser.cpp`): Tripos atoms and bonds, Corina-style
//! substructure cleanup and formal-charge guessing (or `UNITY_ATOM_ATTR`
//! charges), chirality from 3D, and the parser's sanitization flow with
//! double-bond stereo detected from the coordinates
//! (`MolOps::detectBondStereochemistry`, `Chirality.cpp`).

use super::RdkitSmilesError;
use super::mol::{Atom, Bond, BondDir, BondStereo, BondType, Mol};
use super::pdb_read::assign_chiral_types_from_3d;
use super::periodic;
use super::sanitize;
use super::stereo;

/// An exception that makes `Chem.MolFromMol2Block` return `None`.
fn fails(what: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Sanitization(what.into())
}

fn unsupported(what: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Unsupported(what.into())
}

/// A string stream read with `std::getline` plus RDKit's `getLine` (a
/// trailing `\r` dropped), with `tellg`/`seekg` and the eof flag.
struct Stream<'a> {
    text: &'a str,
    pos: usize,
    eof: bool,
}

impl<'a> Stream<'a> {
    fn get_line(&mut self) -> &'a str {
        if self.pos >= self.text.len() {
            self.eof = true;
            return "";
        }
        let rest = &self.text[self.pos..];
        let line = match rest.find('\n') {
            Some(i) => {
                self.pos += i + 1;
                &rest[..i]
            }
            None => {
                self.pos = self.text.len();
                self.eof = true;
                rest
            }
        };
        line.strip_suffix('\r').unwrap_or(line)
    }

    fn seek(&mut self, pos: usize) {
        self.pos = pos;
        self.eof = false;
    }
}

/// `boost::char_separator<char>(" \t\n")` tokens.
fn tokens(line: &str) -> Vec<&str> {
    line.split([' ', '\t', '\n'])
        .filter(|t| !t.is_empty())
        .collect()
}

/// `boost::lexical_cast<unsigned int>` (a leading `-` wraps around in
/// boost; such indices and counts only lead to failures, reported here).
fn to_unsigned(s: &str) -> Result<u32, RdkitSmilesError> {
    s.parse::<u32>()
        .map_err(|_| fails(format!("cannot convert {s} to unsigned int")))
}

/// `boost::lexical_cast<double>`.
fn to_double(s: &str) -> Result<f64, RdkitSmilesError> {
    s.parse::<f64>()
        .map_err(|_| fails("Cannot process mol2 coordinates."))
}

/// The parsed molecule with its Tripos atom types.
struct Parsed {
    mol: Mol,
    coords: Vec<[f64; 3]>,
    types: Vec<String>,
}

/// `ParseMol2FileAtomLine`: `None` for a lone pair.
fn parse_atom_line(line: &str) -> Result<Option<(Atom, [f64; 3], String)>, RdkitSmilesError> {
    let toks = tokens(line);
    if toks.is_empty() {
        return Err(fails("no info in mol2 atom line"));
    }
    if toks.len() < 3 {
        return Err(fails("premature end of mol2 atom line"));
    }
    let x = to_double(toks[2])?;
    if toks.len() < 4 {
        return Err(fails("premature end of mol2 atom line"));
    }
    let y = to_double(toks[3])?;
    if toks.len() < 5 {
        return Err(fails("premature end of mol2 atom line"));
    }
    let z = to_double(toks[4])?;
    if toks.len() < 6 {
        return Err(fails("premature end of mol2 atom line"));
    }
    let tat = toks[5];
    let symb = tat.split('.').next().unwrap_or("");
    match symb {
        "LP" => return Ok(None),
        "ANY" | "Du" | "HEV" | "HET" | "HAL" => {
            return Err(unsupported(format!("Tripos query atom type {tat}")));
        }
        _ => {}
    }
    let anum = periodic::atomic_number(symb)
        .ok_or_else(|| fails(format!("Element '{symb}' not found")))?;
    let mut atom = Atom::new(anum);
    atom.no_implicit = true;
    Ok(Some((atom, [x, y, z], tat.to_string())))
}

/// `ParseMol2FileBondLine`: `None` for a bond to a lone pair or of an
/// unsupported type (skipped with a warning).
fn parse_bond_line(
    line: &str,
    idx_corresp: &[i64],
) -> Result<Option<(usize, usize, BondType)>, RdkitSmilesError> {
    let toks = tokens(line);
    if toks.len() < 4 {
        return Err(fails("no info in mol2 bond line"));
    }
    let idx1 = to_unsigned(toks[1])?.wrapping_sub(1) as usize;
    let idx2 = to_unsigned(toks[2])?.wrapping_sub(1) as usize;
    // RDKit only rejects both ends out of range (one end out of range reads
    // past its index table).
    if idx1 >= idx_corresp.len() || idx2 >= idx_corresp.len() {
        return Err(fails("index mismatch"));
    }
    if idx_corresp[idx1] < 0 || idx_corresp[idx2] < 0 {
        return Ok(None);
    }
    let bt = match toks[3] {
        "1" | "am" => BondType::Single,
        "2" => BondType::Double,
        "3" => BondType::Triple,
        "ar" => BondType::Aromatic,
        "du" | "un" => return Err(unsupported("unspecified Tripos bond type")),
        _ => return Ok(None),
    };
    Ok(Some((
        idx_corresp[idx1] as usize,
        idx_corresp[idx2] as usize,
        bt,
    )))
}

/// Plain SSSR ring information (`MolOps::findSSSR`) unless rings are set.
fn ensure_sssr(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    if mol.rings.is_some() {
        return Ok(());
    }
    let bonds: Vec<(usize, usize, bool)> = mol
        .bonds
        .iter()
        .map(|b| (b.begin, b.end, b.bt != BondType::Dative))
        .collect();
    let rings = chematic_perception::rdkit_sssr(mol.atoms.len(), &bonds).ok_or_else(|| {
        unsupported("ring perception falls back to RDKit's approximate ring finder")
    })?;
    mol.set_rings(rings);
    Ok(())
}

fn symbol_is(mol: &Mol, a: usize, s: &str) -> bool {
    periodic::symbol(mol.atoms[a].anum) == s
}

fn set_bond(mol: &mut Mol, a: usize, b: usize, bt: BondType) -> Result<(), RdkitSmilesError> {
    let bi = mol
        .bond_between(a, b)
        .ok_or_else(|| fails("no bond between atoms"))?;
    mol.bonds[bi].bt = bt;
    mol.bonds[bi].aromatic = false;
    Ok(())
}

/// `chkNoHNeighbNOx`: the number of H neighbours of `n`; `to_mod` becomes
/// `n` when it carries a terminal oxygen (an N-oxide).
fn chk_no_h_neighb_nox(mol: &Mol, n: usize, to_mod: &mut Option<usize>) -> u32 {
    let mut n_h = 0;
    for nbr in mol.nbrs(n) {
        if mol.atoms[nbr].anum == 1 {
            n_h += 1;
        } else if mol.atoms[nbr].anum == 8 && mol.degree(nbr) == 1 {
            *to_mod = Some(n);
        }
    }
    n_h
}

/// `cleanUpMol2Substructures`: `false` where RDKit gives up (returns no
/// molecule).
fn clean_up_mol2_substructures(p: &mut Parsed) -> Result<bool, RdkitSmilesError> {
    let n = p.mol.atoms.len();
    let mut is_fixed = vec![false; n];
    for idx in 0..n {
        if is_fixed[idx] {
            continue;
        }
        let tat = p.types[idx].as_str();
        if tat == "N.4" {
            p.mol.atoms[idx].charge = 1;
        } else if tat == "O.co2" {
            if p.mol.degree(idx) != 1 {
                return Ok(false);
            }
            let nbr = p.mol.nbrs(idx).next().expect("degree 1");
            let ntat = p.types[nbr].as_str();
            if ntat == "P.3" {
                set_bond(&mut p.mol, idx, nbr, BondType::Double)?;
                p.mol.atoms[idx].aromatic = false;
                is_fixed[idx] = true;
                let onbrs: Vec<usize> = p.mol.nbrs(nbr).collect();
                for onbr in onbrs {
                    if p.mol.atoms[onbr].anum == 8 && !is_fixed[onbr] && p.types[onbr] == "O.co2" {
                        set_bond(&mut p.mol, nbr, onbr, BondType::Single)?;
                        p.mol.atoms[onbr].charge = -1;
                        p.mol.atoms[onbr].aromatic = false;
                        is_fixed[onbr] = true;
                    }
                }
                p.mol.atoms[nbr].aromatic = false;
                is_fixed[nbr] = true;
            } else if ntat == "C.2" || ntat == "S.o2" {
                if !is_fixed[nbr] {
                    set_bond(&mut p.mol, idx, nbr, BondType::Single)?;
                    p.mol.atoms[idx].charge = -1;
                    p.mol.atoms[idx].aromatic = false;
                    p.mol.atoms[nbr].aromatic = false;
                    is_fixed[idx] = true;
                    is_fixed[nbr] = true;
                } else {
                    set_bond(&mut p.mol, idx, nbr, BondType::Double)?;
                    p.mol.atoms[idx].aromatic = false;
                    is_fixed[idx] = true;
                }
            } else {
                return Ok(false);
            }
        } else if tat == "C.cat" {
            is_fixed[idx] = true;
            let nbrs: Vec<usize> = p.mol.nbrs(idx).collect();
            let n_n = nbrs.iter().filter(|&&x| symbol_is(&p.mol, x, "N")).count();
            if !(2..=3).contains(&n_n) {
                return Ok(false);
            }
            let mut to_mod: Option<usize> = None;
            if n_n == 2 {
                let mut first: Option<usize> = None;
                let mut second: Option<usize> = None;
                for &nb in &nbrs {
                    if symbol_is(&p.mol, nb, "N") {
                        set_bond(&mut p.mol, idx, nb, BondType::Single)?;
                        p.mol.atoms[nb].aromatic = false;
                        is_fixed[nb] = true;
                        if first.is_some() {
                            second = Some(nb);
                        } else {
                            first = Some(nb);
                        }
                    }
                }
                let (n1, n2) = (first.expect("two N"), second.expect("two N"));
                let h1 = chk_no_h_neighb_nox(&p.mol, n1, &mut to_mod);
                let h2 = chk_no_h_neighb_nox(&p.mol, n2, &mut to_mod);
                if to_mod.is_none() {
                    if h1 != h2 {
                        to_mod = Some(if h1 > h2 { n1 } else { n2 });
                    } else {
                        ensure_sssr(&mut p.mol)?;
                        let r1 = p.mol.num_atom_rings(n1);
                        let r2 = p.mol.num_atom_rings(n2);
                        to_mod = Some(if r1 > r2 { n1 } else { n2 });
                    }
                }
            } else {
                let mut lowest_deg = 100u32;
                for &nb in &nbrs {
                    if !is_fixed[nb] {
                        let mut hvy = 0u32;
                        for nn in p.mol.nbrs(nb) {
                            if p.mol.atoms[nn].anum > 1 {
                                hvy += if p.types[nn] == "C.cat" { 2 } else { 1 };
                            }
                        }
                        if hvy < lowest_deg {
                            to_mod = Some(nb);
                            lowest_deg = hvy;
                        }
                        set_bond(&mut p.mol, idx, nb, BondType::Single)?;
                        p.mol.atoms[nb].aromatic = false;
                        is_fixed[nb] = true;
                    } else {
                        set_bond(&mut p.mol, idx, nb, BondType::Single)?;
                    }
                }
            }
            let to_mod = to_mod.ok_or_else(|| fails("C.cat without a nitrogen to charge"))?;
            set_bond(&mut p.mol, idx, to_mod, BondType::Double)?;
            p.mol.atoms[to_mod].charge = 1;
            p.mol.atoms[idx].aromatic = false;
        }
    }
    Ok(true)
}

/// `fixNitroSubstructureAndCharge`.
fn fix_nitro(mol: &mut Mol, a: usize) {
    let mut n_o_dbl = 0;
    let mut to_mod = 0;
    for &b in &mol.atom_bonds[a] {
        let nbr = mol.bonds[b].other(a);
        if mol.atoms[nbr].anum == 8 && mol.bonds[b].bt == BondType::Double {
            n_o_dbl += 1;
            to_mod = nbr;
        }
    }
    if n_o_dbl == 2 {
        let b = mol.bond_between(a, to_mod).expect("bonded");
        mol.bonds[b].bt = BondType::Single;
        mol.atoms[a].charge = 1;
        mol.atoms[to_mod].charge = -1;
    }
}

/// `guessFormalCharges`.
fn guess_formal_charges(p: &mut Parsed) -> Result<(), RdkitSmilesError> {
    for a in 0..p.mol.atoms.len() {
        if p.mol.atoms[a].charge != 0 || symbol_is(&p.mol, a, "C") {
            continue;
        }
        let mut n_arom = 0;
        let mut accum = 0.0f64;
        for &b in &p.mol.atom_bonds[a] {
            accum += p.mol.bonds[b].valence_contrib(a);
            if p.mol.bonds[b].bt == BondType::Aromatic {
                n_arom += 1;
            }
        }
        ensure_sssr(&mut p.mol)?;
        let tatt = p.types[a].as_str();
        if p.mol.atoms[a].aromatic
            && !tatt.contains("ar")
            && p.mol.ring_info().is_atom_in_ring_of_size(a, 5)
        {
            continue;
        }
        if n_arom == 3 && tatt == "N.ar" {
            continue;
        }
        let exp_val = (accum + 0.1).round() as i32;
        let anum = p.mol.atoms[a].anum;
        let valens = periodic::valence_list(anum);
        let n_elec = periodic::n_outer_elecs(anum);
        let mut assign = if n_elec >= 4 {
            exp_val - i32::from(valens[0])
        } else {
            i32::from(valens[0]) - exp_val
        };
        if assign > 0 && n_elec >= 4 {
            for &vi in valens {
                let vi = i32::from(vi);
                assign = exp_val - vi;
                if vi <= exp_val && assign < 2 {
                    break;
                }
            }
        }
        if assign != 0 {
            p.mol.atoms[a].charge = if p.mol.atoms[a].aromatic && assign.abs() > 1 {
                assign.signum()
            } else {
                assign
            };
            if assign == 2 && exp_val == 5 && symbol_is(&p.mol, a, "N") {
                fix_nitro(&mut p.mol, a);
            }
        }
    }
    Ok(())
}

/// `readFormalChargesFromAttr`.
fn read_formal_charges_from_attr(
    s: &mut Stream<'_>,
    mol: &mut Mol,
) -> Result<(), RdkitSmilesError> {
    let mut line = s.get_line();
    if s.eof {
        return Err(fails("premature EOF in readFormalCharges"));
    }
    loop {
        let toks = tokens(line);
        let atom_idx = to_unsigned(toks.first().copied().unwrap_or(""))
            .map_err(|_| fails("Cannot process mol2 UnityAtomAttr."))?;
        let n_attr = to_unsigned(toks.get(1).copied().unwrap_or(""))
            .map_err(|_| fails("Cannot process mol2 UnityAtomAttr."))?;
        for _ in 0..n_attr {
            let attr = s.get_line();
            if s.eof {
                return Err(fails("premature EOF in readFormalCharges"));
            }
            let at = tokens(attr);
            if at.first() == Some(&"AtomExpr") {
                let val = at.get(1).copied().unwrap_or("");
                if !val.contains('=') {
                    let chg: i32 = val
                        .parse()
                        .map_err(|_| fails("Cannot process mol2 formal charge."))?;
                    let a = (atom_idx as usize).wrapping_sub(1);
                    if a >= mol.atoms.len() {
                        return Err(fails("atom index out of range"));
                    }
                    mol.atoms[a].charge = chg;
                }
            }
        }
        if s.eof {
            break;
        }
        line = s.get_line();
        if line.is_empty() || line.starts_with('@') || line.starts_with('#') {
            break;
        }
    }
    Ok(())
}

fn psub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn pdot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn plen_sq(a: [f64; 3]) -> f64 {
    a[0] * a[0] + a[1] * a[1] + a[2] * a[2]
}

/// `Point3D::crossProduct`.
fn pcross(a: [f64; 3], o: [f64; 3]) -> [f64; 3] {
    [
        a[1] * o[2] - a[2] * o[1],
        -a[0] * o[2] + a[2] * o[0],
        a[0] * o[1] - a[1] * o[0],
    ]
}

/// `Point3D::angleTo`.
fn angle_to(a: [f64; 3], o: [f64; 3]) -> f64 {
    let lsq = plen_sq(a) * plen_sq(o);
    let d = pdot(a, o) / lsq.sqrt();
    if d <= -1.0 {
        return std::f64::consts::PI;
    }
    if d >= 1.0 {
        return 0.0;
    }
    d.acos()
}

/// `RDGeom::computeDihedralAngle`.
fn dihedral(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], p4: [f64; 3]) -> f64 {
    let beg_end = psub(p3, p2);
    let beg_nbr = psub(p1, p2);
    let crs1 = pcross(beg_nbr, beg_end);
    let end_nbr = psub(p4, p3);
    let crs2 = pcross(end_nbr, beg_end);
    angle_to(crs1, crs2)
}

/// `isLinearArrangement`.
fn is_linear(v1: [f64; 3], v2: [f64; 3]) -> bool {
    let lsq = plen_sq(v1) * plen_sq(v2);
    if lsq < 1.0e-6 {
        return true;
    }
    let cos178 = -0.999388;
    pdot(v1, v2) < cos178 * lsq.sqrt()
}

/// `Chirality::detail::setStereoForBond(mol, bond, stereo, false)`.
fn set_stereo_for_bond(mol: &mut Mol, b: usize, stereo: BondStereo) {
    let (mut beg, mut end) = (mol.bonds[b].begin, mol.bonds[b].end);
    if beg > end {
        std::mem::swap(&mut beg, &mut end);
    }
    if mol.degree(beg) > 1 && mol.degree(end) > 1 {
        let mut beg_control = mol.atoms.len();
        for nbr in mol.nbrs(beg) {
            if nbr != end {
                beg_control = beg_control.min(nbr);
            }
        }
        let mut end_control = 0;
        for nbr in mol.nbrs(end) {
            if nbr != beg {
                end_control = end_control.max(nbr);
            }
        }
        if beg != mol.bonds[b].begin {
            std::mem::swap(&mut beg_control, &mut end_control);
        }
        mol.bonds[b].stereo_atoms = vec![beg_control, end_control];
        mol.bonds[b].stereo = stereo;
    }
}

/// `setBondDirRelativeToAtom`.
fn set_bond_dir_relative_to_atom(
    mol: &mut Mol,
    b: usize,
    atom: usize,
    dir: BondDir,
    mut reverse: bool,
) -> Result<(), RdkitSmilesError> {
    if dir == BondDir::None {
        return Err(fails("bad bond direction"));
    }
    if mol.bonds[b].begin != atom {
        reverse = !reverse;
    }
    mol.bonds[b].dir = if reverse { dir.flipped() } else { dir };
    Ok(())
}

/// State shared by `updateDoubleBondNeighbors` calls.
struct DirState {
    needs_dir: Vec<bool>,
    single_bond_counts: Vec<u32>,
    single_bond_nbrs: Vec<Vec<usize>>,
}

/// `controllingBondFromAtom`: `(bond, obond)`.
fn controlling_bond_from_atom(
    mol: &Mol,
    st: &DirState,
    dbl: usize,
    atom: usize,
) -> (Option<usize>, Option<usize>) {
    let mut bond: Option<usize> = None;
    let mut obond: Option<usize> = None;
    for &t in &mol.atom_bonds[atom] {
        if t == dbl {
            continue;
        }
        if matches!(mol.bonds[t].bt, BondType::Single | BondType::Aromatic) {
            match bond {
                None => bond = Some(t),
                Some(cur) => {
                    if st.needs_dir[t] {
                        if st.single_bond_counts[t] > st.single_bond_counts[cur] {
                            obond = Some(cur);
                            bond = Some(t);
                        } else {
                            obond = Some(t);
                        }
                    } else {
                        obond = Some(cur);
                        bond = Some(t);
                    }
                }
            }
        }
    }
    (bond, obond)
}

/// `updateDoubleBondNeighbors` with a conformer.
fn update_double_bond_neighbors(
    mol: &mut Mol,
    dbl: usize,
    coords: &[[f64; 3]],
    st: &mut DirState,
) -> Result<(), RdkitSmilesError> {
    if !st.needs_dir[dbl] {
        return Ok(());
    }
    st.needs_dir[dbl] = false;
    let atom1 = mol.bonds[dbl].begin;
    let atom2 = mol.bonds[dbl].end;
    let (Some(mut bond1), mut obond1) = controlling_bond_from_atom(mol, st, dbl, atom1) else {
        return Ok(());
    };
    let (Some(mut bond2), mut obond2) = controlling_bond_from_atom(mol, st, dbl, atom2) else {
        return Ok(());
    };
    let begin_p = coords[atom1];
    let end_p = coords[atom2];
    let mut bond1_p = coords[mol.bonds[bond1].other(atom1)];
    let mut bond2_p = coords[mol.bonds[bond2].other(atom2)];
    let mut linear = false;
    let mut p1 = psub(bond1_p, begin_p);
    let mut p2 = psub(end_p, begin_p);
    if is_linear(p1, p2) {
        match obond1 {
            None => linear = true,
            Some(ob) => {
                obond1 = Some(bond1);
                bond1 = ob;
                bond1_p = coords[mol.bonds[bond1].other(atom1)];
                p1 = psub(bond1_p, begin_p);
                if is_linear(p1, p2) {
                    linear = true;
                }
            }
        }
    }
    if !linear {
        p1 = psub(bond2_p, end_p);
        p2 = psub(begin_p, end_p);
        if is_linear(p1, p2) {
            match obond2 {
                None => linear = true,
                Some(ob) => {
                    obond2 = Some(bond2);
                    bond2 = ob;
                    bond2_p = coords[mol.bonds[bond2].other(atom2)];
                    // RDKit measures this one from the begin atom.
                    p1 = psub(bond2_p, begin_p);
                    if is_linear(p1, p2) {
                        linear = true;
                    }
                }
            }
        }
    }
    if linear {
        set_stereo_for_bond(mol, dbl, BondStereo::Any);
        return Ok(());
    }
    let ang = dihedral(bond1_p, begin_p, end_p, bond2_p);
    let same_torsion_dir = ang >= std::f64::consts::PI / 2.0;
    let mut reverse = same_torsion_dir;

    let mut followup: Vec<usize> = Vec::new();
    for b in [bond1, bond2] {
        if st.needs_dir[b] {
            for &bidx in &st.single_bond_nbrs[b] {
                if st.needs_dir[bidx] {
                    followup.push(bidx);
                }
            }
        }
    }
    if !st.needs_dir[bond1] {
        if st.needs_dir[bond2] {
            if mol.bonds[bond1].begin != atom1 {
                reverse = !reverse;
            }
            let dir = mol.bonds[bond1].dir;
            set_bond_dir_relative_to_atom(mol, bond2, atom2, dir, reverse)?;
        }
    } else if !st.needs_dir[bond2] {
        if mol.bonds[bond2].begin != atom2 {
            reverse = !reverse;
        }
        let dir = mol.bonds[bond2].dir;
        set_bond_dir_relative_to_atom(mol, bond1, atom1, dir, reverse)?;
    } else {
        set_bond_dir_relative_to_atom(mol, bond1, atom1, BondDir::EndDownRight, false)?;
        set_bond_dir_relative_to_atom(mol, bond2, atom2, BondDir::EndDownRight, reverse)?;
    }
    st.needs_dir[bond1] = false;
    st.needs_dir[bond2] = false;
    if let Some(ob) = obond1
        && st.needs_dir[ob]
    {
        let dir = mol.bonds[bond1].dir;
        let rev = mol.bonds[bond1].begin == atom1;
        set_bond_dir_relative_to_atom(mol, ob, atom1, dir, rev)?;
        st.needs_dir[ob] = false;
    }
    if let Some(ob) = obond2
        && st.needs_dir[ob]
    {
        let dir = mol.bonds[bond2].dir;
        let rev = mol.bonds[bond2].begin == atom2;
        set_bond_dir_relative_to_atom(mol, ob, atom2, dir, rev)?;
        st.needs_dir[ob] = false;
    }
    for b in followup {
        update_double_bond_neighbors(mol, b, coords, st)?;
    }
    Ok(())
}

/// `MolOps::detectBondStereochemistry(mol)`
/// (`setDoubleBondNeighborDirections` from the conformer).
fn detect_bond_stereochemistry(mol: &mut Mol, coords: &[[f64; 3]]) -> Result<(), RdkitSmilesError> {
    mol.find_rings()?;
    let nb = mol.bonds.len();
    let mut st = DirState {
        needs_dir: vec![false; nb],
        single_bond_counts: vec![0; nb],
        single_bond_nbrs: vec![Vec::new(); nb],
    };
    let mut dbl_bond_nbrs: Vec<Vec<usize>> = vec![Vec::new(); nb];
    let mut in_play: Vec<usize> = Vec::new();
    for b in 0..nb {
        let bond = &mol.bonds[b];
        if !(bond.bt == BondType::Double
            && bond.stereo != BondStereo::Any
            && mol.degree(bond.begin) > 1
            && mol.degree(bond.end) > 1
            && stereo::should_detect_double_bond_stereo(mol, b))
        {
            continue;
        }
        for atom in [bond.begin, bond.end] {
            for &nbr_bond in &mol.atom_bonds[atom] {
                if matches!(
                    mol.bonds[nbr_bond].bt,
                    BondType::Single | BondType::Aromatic
                ) {
                    st.single_bond_counts[nbr_bond] += 1;
                    st.needs_dir[b] = true;
                    st.needs_dir[nbr_bond] = true;
                    dbl_bond_nbrs[b].push(nbr_bond);
                    if !st.single_bond_nbrs[nbr_bond].contains(&b) {
                        st.single_bond_nbrs[nbr_bond].push(b);
                    }
                }
            }
        }
        in_play.push(b);
    }
    if in_play.is_empty() {
        return Ok(());
    }
    // (count, Bond*) pairs sorted ascending; bond pointers follow creation
    // order, approximated by the bond index.
    let mut ordered: Vec<(u64, usize)> = in_play
        .iter()
        .map(|&b| {
            let mut count: u64 = dbl_bond_nbrs[b].iter().map(|&x| x as u64).sum();
            if mol.num_bond_rings(b) == 0 {
                count *= 10;
            }
            (count, b)
        })
        .collect();
    ordered.sort_unstable();
    for &(_, b) in ordered.iter().rev() {
        update_double_bond_neighbors(mol, b, coords, &mut st)?;
    }
    Ok(())
}

/// `Chem.MolFromMol2Block(text, sanitize, removeHs, cleanupSubstructures)`:
/// the molecule and its positions, `Err(Sanitization)` where RDKit returns
/// no molecule.
pub(crate) fn mol_from_mol2_block(
    text: &str,
    sanitize: bool,
    remove_hs: bool,
    cleanup_substructures: bool,
) -> Result<(Mol, Vec<[f64; 3]>), RdkitSmilesError> {
    let mut s = Stream {
        text,
        pos: 0,
        eof: false,
    };
    let (mut mol_start, mut atom_start, mut bond_start, mut charge_start) = (0, 0, 0, 0);
    while !s.eof {
        let line = s.get_line();
        if s.eof {
            break;
        }
        if line.starts_with('@') {
            let first = tokens(line).first().copied().unwrap_or("");
            match first {
                "@<TRIPOS>MOLECULE" => {
                    if mol_start == 0 {
                        mol_start = s.pos;
                    } else {
                        break;
                    }
                }
                "@<TRIPOS>ATOM" => atom_start = s.pos,
                "@<TRIPOS>BOND" => bond_start = s.pos,
                "@<TRIPOS>UNITY_ATOM_ATTR" => charge_start = s.pos,
                _ => {}
            }
        }
    }
    if mol_start == 0 {
        return Err(fails("No MOLECULE block found in Mol2 data"));
    }
    if atom_start == 0 {
        return Err(fails("No ATOM block found in Mol2 data"));
    }
    s.seek(mol_start);
    let _name = s.get_line();
    let counts = tokens(s.get_line());
    if counts.is_empty() {
        return Err(fails("Empty counts line"));
    }
    let n_atoms = to_unsigned(counts[0])? as usize;
    let n_bonds = match counts.get(1) {
        Some(t) => to_unsigned(t)? as usize,
        None => 0,
    };
    if n_atoms == 0 {
        return Err(fails("molecule has no atoms"));
    }
    // Every atom takes a line: more atoms than bytes hit a premature EOF.
    if n_atoms > text.len() {
        return Err(fails("premature EOF"));
    }

    // ParseMol2AtomBlock
    s.seek(atom_start);
    let mut idx_corresp: Vec<i64> = vec![-1; n_atoms];
    let mut p = Parsed {
        mol: Mol::default(),
        coords: Vec::new(),
        types: Vec::new(),
    };
    for slot in idx_corresp.iter_mut() {
        let line = s.get_line();
        if s.eof {
            return Err(fails("premature EOF"));
        }
        if let Some((atom, pos, tat)) = parse_atom_line(line)? {
            *slot = p.mol.add_atom(atom) as i64;
            p.coords.push(pos);
            p.types.push(tat);
        }
    }
    // ParseMol2BondBlock
    if n_bonds > 0 {
        s.seek(bond_start);
        for _ in 0..n_bonds {
            let line = s.get_line();
            if s.eof {
                return Err(fails("premature EOF"));
            }
            if let Some((a1, a2, bt)) = parse_bond_line(line, &idx_corresp)? {
                if a1 == a2 {
                    return Err(fails("attempt to add self-bond"));
                }
                if p.mol.bond_between(a1, a2).is_some() {
                    return Err(fails("bond already exists"));
                }
                let mut bond = Bond::new(a1, a2, bt);
                if bt == BondType::Aromatic {
                    bond.aromatic = true;
                    p.mol.atoms[a1].aromatic = true;
                    p.mol.atoms[a2].aromatic = true;
                }
                p.mol.add_bond(bond);
            }
        }
    }

    if charge_start == 0 {
        if cleanup_substructures && !clean_up_mol2_substructures(&mut p)? {
            return Err(fails("Mol2 substructure cleanup failed"));
        }
        guess_formal_charges(&mut p)?;
    } else {
        s.seek(charge_start);
        read_formal_charges_from_attr(&mut s, &mut p.mol)?;
    }

    let Parsed {
        mut mol,
        mut coords,
        ..
    } = p;
    assign_chiral_types_from_3d(&mut mol, &coords)?;
    if sanitize {
        sanitize::clean_up(&mut mol)?;
        if remove_hs {
            // sanitizeMol(SANITIZE_CLEANUP)
            sanitize::clean_up(&mut mol)?;
            mol.update_property_cache(false)?;
            detect_bond_stereochemistry(&mut mol, &coords)?;
            let kept = sanitize::remove_hs_with(&mut mol, false, false)?;
            coords = kept.iter().map(|&i| coords[i]).collect();
            sanitize::sanitize_mol_with(&mut mol, false)?;
        } else {
            sanitize::sanitize_mol_with(&mut mol, false)?;
            detect_bond_stereochemistry(&mut mol, &coords)?;
        }
        mol.update_property_cache(false)?;
        stereo::legacy_stereo_perception(&mut mol, true, false);
    }
    Ok((mol, coords))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn block(atoms: &[(&str, [f64; 3])], bonds: &[(usize, usize, &str)]) -> String {
        let mut s = format!(
            "@<TRIPOS>MOLECULE\nx\n{} {} 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n",
            atoms.len(),
            bonds.len()
        );
        for (i, (t, p)) in atoms.iter().enumerate() {
            s += &format!(
                "{} A{} {} {} {} {} 1 LIG 0.0\n",
                i + 1,
                i + 1,
                p[0],
                p[1],
                p[2],
                t
            );
        }
        s += "@<TRIPOS>BOND\n";
        for (i, (a, b, t)) in bonds.iter().enumerate() {
            s += &format!("{} {} {} {}\n", i + 1, a, b, t);
        }
        s
    }

    #[test]
    fn premature_eof_and_bad_symbols_fail() {
        let ok = block(
            &[("C.3", [0.0, 0.0, 0.0]), ("Cl", [1.7, 0.0, 0.0])],
            &[(1, 2, "1")],
        );
        assert!(mol_from_mol2_block(&ok, true, true, true).is_ok());
        let no_newline = ok.trim_end();
        assert!(matches!(
            mol_from_mol2_block(no_newline, true, true, true),
            Err(RdkitSmilesError::Sanitization(_))
        ));
        let bad = ok.replace(" Cl ", " CL ");
        assert!(matches!(
            mol_from_mol2_block(&bad, true, true, true),
            Err(RdkitSmilesError::Sanitization(_))
        ));
    }

    #[test]
    fn carboxylate_gets_charge() {
        // Acetate: C.3-C.2(O.co2)2 with three H.
        let atoms = [
            ("C.3", [0.0, 0.0, 0.0]),
            ("C.2", [1.5, 0.0, 0.0]),
            ("O.co2", [2.1, 1.1, 0.0]),
            ("O.co2", [2.1, -1.1, 0.0]),
            ("H", [-0.4, 1.0, 0.0]),
            ("H", [-0.4, -0.5, 0.9]),
            ("H", [-0.4, -0.5, -0.9]),
        ];
        let bonds = [
            (1, 2, "1"),
            (2, 3, "ar"),
            (2, 4, "ar"),
            (1, 5, "1"),
            (1, 6, "1"),
            (1, 7, "1"),
        ];
        let (m, coords) = mol_from_mol2_block(&block(&atoms, &bonds), true, true, true).unwrap();
        assert_eq!(m.atoms.len(), 4);
        assert_eq!(coords.len(), 4);
        assert_eq!(m.atoms[2].charge, -1);
        assert_eq!(m.atoms[3].charge, 0);
    }
}
