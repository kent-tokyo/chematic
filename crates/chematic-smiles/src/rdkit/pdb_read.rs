//! RDKit's PDB reader (`MolFromPDBBlock`, RDKit 2026.03.1
//! `FileParsers/PDBParser.cpp` and `ProximityBonds.cpp`) on the
//! RDKit-model molecule: ATOM/HETATM/CONECT records, proximity bonding,
//! standard residue bond orders, sanitization and 3D chirality.

use super::RdkitSmilesError;
use super::mol::{Atom, Bond, BondDir, BondType, ChiralTag, Mol};
use super::periodic;
use super::sanitize;

fn parse_error(what: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Unsupported(what.into())
}

/// An error that makes RDKit return no molecule (a file-parse exception,
/// reported like a sanitization failure).
fn rdkit_fails(what: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Sanitization(what.into())
}

/// `AtomPDBResidueInfo` fields the reader uses.
#[derive(Clone, Debug)]
pub(crate) struct ResidueInfo {
    /// The 4-character atom name (columns 13-16).
    pub name: [u8; 4],
    /// The 3-character residue name.
    pub res_name: Vec<u8>,
    pub chain: u8,
    pub insertion_code: u8,
    pub res_no: i32,
    pub hetero: bool,
}

/// A molecule read from PDB records: the RDKit model, one position per
/// atom (the first conformer) and per-atom residue information.
pub(crate) struct PdbMol {
    pub mol: Mol,
    pub coords: Vec<[f64; 3]>,
    pub info: Vec<ResidueInfo>,
    pub is_3d: bool,
}

/// `FileParserUtils::toInt(s, acceptSpaces=true)`.
fn to_int(s: &[u8]) -> Result<i32, ()> {
    for &c in s {
        if c == 0 {
            break;
        }
        if !(c.is_ascii_digit() || c == b' ' || c == b'+' || c == b'-') {
            return Err(());
        }
    }
    let mut i = 0;
    while i < s.len() && s[i] == b' ' {
        i += 1;
    }
    if i == s.len() {
        return Ok(0);
    }
    // std::from_chars: an optional '-' then digits; anything else stops.
    let rest = &s[i..];
    let (neg, digits) = match rest.first() {
        Some(b'-') => (true, &rest[1..]),
        _ => (false, rest),
    };
    let n = digits.iter().take_while(|c| c.is_ascii_digit()).count();
    if n == 0 {
        return Ok(0);
    }
    let text = std::str::from_utf8(&digits[..n]).map_err(|_| ())?;
    // Out-of-range values leave the result 0 (from_chars' error).
    let v: i64 = text.parse().unwrap_or(0);
    let v = if neg { -v } else { v };
    Ok(i32::try_from(v).unwrap_or(0))
}

/// `FileParserUtils::toDouble(s, acceptSpaces=true)`: validation, then
/// `atof`.
fn to_double(s: &[u8]) -> Result<f64, ()> {
    for &c in s {
        if c == 0 {
            break;
        }
        if !(c.is_ascii_digit() || matches!(c, b' ' | b'+' | b'-' | b',' | b'.')) {
            return Err(());
        }
    }
    // atof: leading white space, then the longest [+-]digits[.digits] prefix.
    let mut i = 0;
    while i < s.len() && s[i].is_ascii_whitespace() {
        i += 1;
    }
    let start = i;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') {
        i += 1;
    }
    let int_start = i;
    while i < s.len() && s[i].is_ascii_digit() {
        i += 1;
    }
    let mut has_digits = i > int_start;
    if i < s.len() && s[i] == b'.' {
        let frac_start = i + 1;
        let mut j = frac_start;
        while j < s.len() && s[j].is_ascii_digit() {
            j += 1;
        }
        if j > frac_start || has_digits {
            has_digits = has_digits || j > frac_start;
            i = j;
        }
    }
    if !has_digits {
        return Ok(0.0);
    }
    let text = std::str::from_utf8(&s[start..i]).map_err(|_| ())?;
    let text = text.strip_prefix('+').unwrap_or(text);
    Ok(text.parse::<f64>().unwrap_or(0.0))
}

/// `PDBAtomFromSymbol`: (atomic number, isotope).
fn atom_from_symbol(symb: &str) -> Result<Option<(u32, u32)>, RdkitSmilesError> {
    match symb {
        "D" => return Ok(Some((1, 2))),
        "T" => return Ok(Some((1, 3))),
        _ => {}
    }
    match periodic::atomic_number(symb) {
        Some(0) => Ok(None),
        Some(n) => Ok(Some((n, 0))),
        None => Err(rdkit_fails(format!("Element '{symb}' not found"))),
    }
}

struct Reader<'a> {
    text: &'a [u8],
    flavor: u32,
    pdb: PdbMol,
    /// Serial number -> atom index (`amap`).
    amap: std::collections::HashMap<i32, usize>,
    /// Bond index -> CONECT duplicate bitmap (`bmap`).
    bmap: std::collections::HashMap<usize, u8>,
    /// Whether any atom has a position (RDKit's conformer exists).
    has_conformer: bool,
}

/// Byte `i` of the line starting at `start`, read past the line as
/// RDKit's pointer arithmetic does (0 at the end of the text).
fn at(text: &[u8], start: usize, i: usize) -> u8 {
    text.get(start + i).copied().unwrap_or(0)
}

fn field(text: &[u8], start: usize, off: usize, n: usize) -> Vec<u8> {
    (0..n).map(|k| at(text, start, off + k)).collect()
}

impl Reader<'_> {
    fn atom_line(&mut self, start: usize, len: usize) -> Result<(), RdkitSmilesError> {
        if len < 16 {
            return Ok(());
        }
        let text = self.text;
        let p = |i: usize| at(text, start, i);
        if self.flavor & 1 == 0 {
            if len >= 17 && p(16) != b' ' && p(16) != b'A' && p(16) != b'1' {
                return Ok(());
            }
            if len >= 54 && field(text, start, 30, 24) == b"9999.0009999.0009999.000" {
                return Ok(());
            }
            if p(12) == b' ' && p(13) == b'Q' {
                return Ok(());
            }
            if len >= 20 && field(text, start, 18, 3) == b"DUM" {
                return Ok(());
            }
        }
        let serial_field = field(text, start, 6, 5);
        let serialno = to_int(&serial_field).map_err(|_| {
            rdkit_fails(format!(
                "Non-integer PDB serial number {}",
                String::from_utf8_lossy(&serial_field)
            ))
        })?;

        let mut atom: Option<(u32, u32)> = None;
        // Attempt #1: the element symbol in columns 77-78.
        let mut symb = String::new();
        if len >= 78 {
            if p(76).is_ascii_uppercase() {
                symb.push(char::from(p(76)));
                if p(77).is_ascii_uppercase() {
                    symb.push(char::from(p(77) + 32));
                } else if p(77).is_ascii_lowercase() {
                    symb.push(char::from(p(77)));
                }
            } else if p(76) == b' ' && p(77).is_ascii_uppercase() {
                symb.push(char::from(p(77)));
            }
        } else if len == 77 && p(76).is_ascii_uppercase() {
            symb.push(char::from(p(76)));
        }
        if !symb.is_empty() {
            atom = atom_from_symbol(&symb)?;
        }
        if atom.is_none() {
            // Attempt #2: from the atom name.
            let mut symb = String::new();
            if p(13).is_ascii_uppercase() {
                if p(12) == b' ' {
                    symb.push(char::from(p(13)));
                    if p(14).is_ascii_lowercase() {
                        symb.push(char::from(p(14)));
                    }
                } else if p(12).is_ascii_uppercase() {
                    if p(12) == b'H' && p(0) == b'A' {
                        // No He, Hf, Hg, Ho or Hs in ATOM records.
                        symb.push('H');
                    } else {
                        symb.push(char::from(p(12)));
                        symb.push(char::from(p(13) + 32));
                    }
                } else if p(12).is_ascii_digit() {
                    symb.push(char::from(p(13)));
                }
            }
            if !symb.is_empty() {
                atom = atom_from_symbol(&symb)?;
            }
        }
        let Some((anum, isotope)) = atom else {
            return Err(rdkit_fails(format!(
                "Cannot determine element for PDB atom #{serialno}"
            )));
        };
        let mut a = Atom::new(anum);
        a.isotope = isotope;
        let idx = self.pdb.mol.add_atom(a);
        self.amap.insert(serialno, idx);
        let mut pos = [0.0f64; 3];
        if len >= 38 {
            let coord_err =
                || rdkit_fails(format!("Problem with coordinates for PDB atom #{serialno}"));
            pos[0] = to_double(&field(text, start, 30, 8)).map_err(|_| coord_err())?;
            if len >= 46 {
                pos[1] = to_double(&field(text, start, 38, 8)).map_err(|_| coord_err())?;
            }
            if len >= 54 {
                pos[2] = to_double(&field(text, start, 46, 8)).map_err(|_| coord_err())?;
            }
            if !self.has_conformer {
                self.has_conformer = true;
                self.pdb.is_3d = pos[2] != 0.0;
            } else if pos[2] != 0.0 {
                self.pdb.is_3d = true;
            }
        }
        self.pdb.coords.push(pos);
        if len >= 79 {
            let (c78, c79) = (p(78), p(79));
            let mut charge = 0;
            if (b'1'..=b'9').contains(&c78) {
                if c79 == b'-' {
                    charge = -i32::from(c78 - b'0');
                } else if c79 == b'+' || c79 == b' ' || c79 == 0 {
                    charge = i32::from(c78 - b'0');
                }
            } else if c78 == b'+' {
                if (b'1'..=b'9').contains(&c79) {
                    charge = i32::from(c79 - b'0');
                } else if c79 == b'+' {
                    charge = 2;
                } else if c79 != b'0' {
                    charge = 1;
                }
            } else if c78 == b'-' {
                if (b'1'..=b'9').contains(&c79) {
                    charge = i32::from(c79 - b'0');
                } else if c79 == b'-' {
                    charge = -2;
                } else if c79 != b'0' {
                    charge = -1;
                }
            } else if c78 == b' ' {
                if (b'1'..=b'9').contains(&c79) {
                    charge = i32::from(c79 - b'0');
                } else if c79 == b'+' {
                    charge = 1;
                } else if c79 == b'-' {
                    charge = -1;
                }
            }
            if charge != 0 {
                self.pdb.mol.atoms[idx].charge = charge;
            }
        }
        let name_field = field(text, start, 12, 4);
        let name = [name_field[0], name_field[1], name_field[2], name_field[3]];
        let res_name = if len >= 20 {
            field(text, start, 17, 3)
        } else {
            b"UNL".to_vec()
        };
        let chain = if len >= 22 { p(21) } else { b' ' };
        let insertion_code = if len >= 27 { p(26) } else { b' ' };
        let res_no = if len >= 26 {
            to_int(&field(text, start, 22, 4)).map_err(|_| {
                rdkit_fails(format!(
                    "Problem with residue number for PDB atom #{serialno}"
                ))
            })?
        } else {
            1
        };
        if len >= 60 {
            to_double(&field(text, start, 54, 6)).map_err(|_| {
                rdkit_fails(format!("Problem with occupancy for PDB atom #{serialno}"))
            })?;
        }
        if len >= 66 {
            to_double(&field(text, start, 60, 6)).map_err(|_| {
                rdkit_fails(format!(
                    "Problem with temperature factor for PDB atom #{serialno}"
                ))
            })?;
        }
        self.pdb.info.push(ResidueInfo {
            name,
            res_name,
            chain,
            insertion_code,
            res_no,
            hetero: p(0) == b'H',
        });
        Ok(())
    }

    fn bond_line(&mut self, start: usize, len: usize) -> Result<(), RdkitSmilesError> {
        if len < 16 {
            return Ok(());
        }
        let text = self.text;
        let src_field = field(text, start, 6, 5);
        let conect_err = || {
            rdkit_fails(format!(
                "Problem with CONECT record for PDB atom #{}",
                String::from_utf8_lossy(&src_field)
            ))
        };
        let src = to_int(&src_field).map_err(|_| conect_err())?;
        let Some(&src_idx) = self.amap.get(&src) else {
            return Ok(());
        };
        let len = len.min(41);
        let mut pos = 11;
        while pos + 5 <= len {
            let f = field(text, start, pos, 5);
            pos += 5;
            if f == b"     " {
                break;
            }
            let dst = to_int(&f).map_err(|_| conect_err())?;
            if dst == src {
                continue;
            }
            let Some(&dst_idx) = self.amap.get(&dst) else {
                continue;
            };
            match self.pdb.mol.bond_between(src_idx, dst_idx) {
                Some(b) => {
                    let seen = self.bmap.get(&b).copied().unwrap_or(0);
                    let bond = &mut self.pdb.mol.bonds[b];
                    let (mut new_seen, mut new_bt) = (seen, None);
                    if src < dst {
                        if seen & 0x0f == 0x01 {
                            new_seen = seen | 0x02;
                            if seen & 0x20 == 0 {
                                new_bt = Some(BondType::Double);
                            }
                        } else if seen & 0x0f == 0x03 {
                            new_seen = seen | 0x04;
                            if seen & 0x40 == 0 {
                                new_bt = Some(BondType::Triple);
                            }
                        } else if seen & 0x0f == 0x07 {
                            new_seen = seen | 0x08;
                            if seen & 0x80 == 0 {
                                new_bt = Some(BondType::Quadruple);
                            }
                        }
                    } else if seen & 0xf0 == 0x10 {
                        new_seen = seen | 0x20;
                        if seen & 0x02 == 0 {
                            new_bt = Some(BondType::Double);
                        }
                    } else if seen & 0xf0 == 0x30 {
                        new_seen = seen | 0x40;
                        if seen & 0x04 == 0 {
                            new_bt = Some(BondType::Triple);
                        }
                    } else if seen & 0xf0 == 0x70 {
                        new_seen = seen | 0x80;
                        if seen & 0x08 == 0 {
                            new_bt = Some(BondType::Quadruple);
                        }
                    }
                    if let Some(bt) = new_bt {
                        bond.bt = bt;
                    }
                    if new_seen != seen {
                        self.bmap.insert(b, new_seen);
                    }
                }
                None => {
                    if self.is_blacklisted_pair(src_idx, dst_idx) {
                        return Err(parse_error("zero-order bond (CONECT between residues)"));
                    }
                    let b = self
                        .pdb
                        .mol
                        .add_bond(Bond::new(src_idx, dst_idx, BondType::Single));
                    self.bmap.insert(b, if src < dst { 0x01 } else { 0x10 });
                }
            }
        }
        Ok(())
    }

    fn is_blacklisted_pair(&self, a: usize, b: usize) -> bool {
        is_blacklisted_pair(&self.pdb, a, b)
    }
}

fn same_residue(p: &ResidueInfo, q: &ResidueInfo) -> bool {
    p.res_no == q.res_no
        && p.res_name == q.res_name
        && p.chain == q.chain
        && p.insertion_code == q.insertion_code
}

/// `IsBlacklistedAtom`: metals, noble gases and halogens.
fn is_blacklisted_atom(anum: u32) -> bool {
    !((5..=8).contains(&anum)
        || (14..=16).contains(&anum)
        || (32..=34).contains(&anum)
        || (51..=52).contains(&anum))
}

/// `IsBlacklistedPair`.
fn is_blacklisted_pair(pdb: &PdbMol, a: usize, b: usize) -> bool {
    let (pi, qi) = (&pdb.info[a], &pdb.info[b]);
    if !same_residue(pi, qi) {
        if is_blacklisted_atom(pdb.mol.atoms[a].anum) || is_blacklisted_atom(pdb.mol.atoms[b].anum)
        {
            return true;
        }
        if pi.res_name == b"HOH" || qi.res_name == b"HOH" {
            return true;
        }
    }
    false
}

const EXTDIST: f64 = 0.45;
const MAXDIST: f64 = 5.45;
const MINDIST2: f64 = 0.16;
const MAXDIST2: f64 = 29.7025;
const HASHMASK: i32 = 1023;
const HASHX: i32 = 571;
const HASHY: i32 = 127;
const HASHZ: i32 = 3;

#[derive(Clone, Copy, Default)]
struct ProximityEntry {
    x: f32,
    y: f32,
    z: f32,
    r: f32,
    hash: i32,
    next: i32,
    elem: u32,
}

/// `IsBonded` with `ctdIGNORE_H_H_CONTACTS`.
fn is_bonded(p: &ProximityEntry, q: &ProximityEntry) -> bool {
    if p.elem == 1 && q.elem == 1 {
        return false;
    }
    let dx = f64::from(p.x) - f64::from(q.x);
    let mut dist2 = dx * dx;
    if dist2 > MAXDIST2 {
        return false;
    }
    let dy = f64::from(p.y) - f64::from(q.y);
    dist2 += dy * dy;
    if dist2 > MAXDIST2 {
        return false;
    }
    let dz = f64::from(p.z) - f64::from(q.z);
    dist2 += dz * dz;
    if !(MINDIST2..=MAXDIST2).contains(&dist2) {
        return false;
    }
    let radius = f64::from(p.r) + f64::from(q.r) + EXTDIST;
    dist2 <= radius * radius
}

/// `ConnectTheDots_Large(mol, ctdIGNORE_H_H_CONTACTS)`.
fn connect_the_dots(pdb: &mut PdbMol) {
    let count = pdb.mol.atoms.len();
    let mut table = [-1i32; 1024];
    let mut tmp = vec![ProximityEntry::default(); count];
    for i in 0..count {
        let elem = pdb.mol.atoms[i].anum;
        let p = pdb.coords[i];
        tmp[i] = ProximityEntry {
            x: p[0] as f32,
            y: p[1] as f32,
            z: p[2] as f32,
            r: periodic::RCOV.get(elem as usize).copied().unwrap_or(0.0) as f32,
            hash: 0,
            next: 0,
            elem,
        };
        let hash = HASHX * ((p[0] / MAXDIST) as i32)
            + HASHY * ((p[1] / MAXDIST) as i32)
            + HASHZ * ((p[2] / MAXDIST) as i32);
        for dx in [-HASHX, 0, HASHX] {
            for dy in [-HASHY, 0, HASHY] {
                for dz in [-HASHZ, 0, HASHZ] {
                    let probe = hash + dx + dy + dz;
                    let mut list = table[(probe & HASHMASK) as usize];
                    while list != -1 {
                        let j = list as usize;
                        if tmp[j].hash == probe
                            && is_bonded(&tmp[i], &tmp[j])
                            && pdb.mol.bond_between(i, j).is_none()
                            && !is_blacklisted_pair(pdb, i, j)
                        {
                            pdb.mol.add_bond(Bond::new(i, j, BondType::Single));
                        }
                        list = tmp[j].next;
                    }
                }
            }
        }
        let slot = (hash & HASHMASK) as usize;
        tmp[i].next = table[slot];
        table[slot] = i as i32;
        tmp[i].hash = hash;
    }
    // Cleanup: keep only the shortest bond of a multivalent H.
    for i in 0..count {
        if pdb.mol.atoms[i].anum != 1 || pdb.mol.degree(i) <= 1 {
            continue;
        }
        let p = pdb.coords[i];
        let mut best = 10000.0f32;
        let mut best_idx = count + 1;
        let nbrs: Vec<usize> = pdb.mol.nbrs(i).collect();
        for &nb in &nbrs {
            let q = pdb.coords[nb];
            let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
            let d = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() as f32;
            if d < best && pdb.info[i].res_no == pdb.info[nb].res_no {
                best = d;
                best_idx = nb;
            }
        }
        for &nb in &nbrs {
            let b = pdb.mol.bond_between(i, nb).expect("bonded");
            if nb == best_idx {
                pdb.mol.bonds[b].bt = BondType::Single;
            } else {
                pdb.mol.remove_bond(b);
            }
        }
    }
}

const fn atm(s: &[u8; 4]) -> [u8; 4] {
    *s
}

/// `StandardPDBDoubleBond(rescode, atm1, atm2)`.
fn standard_double_bond(res: &[u8], a1: [u8; 4], a2: [u8; 4]) -> bool {
    let (a1, a2) = if u32::from_be_bytes(a1) > u32::from_be_bytes(a2) {
        (a2, a1)
    } else {
        (a1, a2)
    };
    let is = |x: &[u8; 4], y: &[u8; 4]| a1 == atm(x) && a2 == atm(y);
    let backbone = is(b" C  ", b" O  ");
    match res {
        b"ACE" | b"ALA" | b"CYS" | b"GLY" | b"ILE" | b"LEU" | b"LYS" | b"MET" | b"PRO" | b"SER"
        | b"THR" | b"VAL" => backbone,
        b"ARG" => backbone || is(b" CZ ", b" NH2"),
        b"ASN" | b"ASP" => backbone || is(b" CG ", b" OD1"),
        b"GLN" | b"GLU" => backbone || is(b" CD ", b" OE1"),
        b"HIS" => backbone || is(b" CD2", b" CG ") || is(b" CE1", b" ND1"),
        b"PHE" | b"TYR" => {
            backbone || is(b" CD1", b" CG ") || is(b" CD2", b" CE2") || is(b" CE1", b" CZ ")
        }
        b"TRP" => {
            backbone
                || is(b" CD1", b" CG ")
                || is(b" CD2", b" CE2")
                || is(b" CE3", b" CZ3")
                || is(b" CH2", b" CZ2")
        }
        b"  A" | b" DA" => {
            is(b" C6 ", b" N1 ")
                || is(b" C2 ", b" N3 ")
                || is(b" C4 ", b" C5 ")
                || is(b" C8 ", b" N7 ")
                || is(b" OP1", b" P  ")
        }
        b"  G" | b" DG" => {
            is(b" C6 ", b" O6 ")
                || is(b" C2 ", b" N3 ")
                || is(b" C4 ", b" C5 ")
                || is(b" C8 ", b" N7 ")
                || is(b" OP1", b" P  ")
        }
        b"  C" | b" DC" => {
            is(b" C2 ", b" O2 ")
                || is(b" C4 ", b" N3 ")
                || is(b" C5 ", b" C6 ")
                || is(b" OP1", b" P  ")
        }
        b"  T" | b" DT" | b"  U" | b" DU" => {
            is(b" C2 ", b" O2 ")
                || is(b" C4 ", b" O4 ")
                || is(b" C5 ", b" C6 ")
                || is(b" OP1", b" P  ")
        }
        _ => false,
    }
}

/// `StandardPDBResidueBondOrders`.
fn standard_residue_bond_orders(pdb: &mut PdbMol) {
    for b in 0..pdb.mol.bonds.len() {
        if pdb.mol.bonds[b].bt != BondType::Single {
            continue;
        }
        let (beg, end) = (pdb.mol.bonds[b].begin, pdb.mol.bonds[b].end);
        let (bi, ei) = (&pdb.info[beg], &pdb.info[end]);
        if !same_residue(bi, ei) || bi.hetero || ei.hetero {
            continue;
        }
        if !standard_double_bond(&bi.res_name, bi.name, ei.name) {
            continue;
        }
        let has_double = |a: usize| {
            pdb.mol.atom_bonds[a]
                .iter()
                .any(|&x| pdb.mol.bonds[x].bt == BondType::Double)
        };
        if has_double(beg) || has_double(end) {
            continue;
        }
        pdb.mol.bonds[b].bt = BondType::Double;
    }
}

/// `StandardPDBChiralAtom(resnam, atmnam)`.
fn standard_chiral_atom(res: &[u8], name: [u8; 4]) -> bool {
    match res {
        b"GLY" => false,
        b"ILE" | b"THR" => &name == b" CA " || &name == b" CB ",
        b"ALA" | b"ARG" | b"ASN" | b"ASP" | b"CYS" | b"GLN" | b"GLU" | b"HIS" | b"LEU" | b"LYS"
        | b"MET" | b"PHE" | b"PRO" | b"SER" | b"TRP" | b"TYR" | b"VAL" => &name == b" CA ",
        _ => false,
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        -a[0] * b[2] + a[2] * b[0],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Whether `assignNontetrahedralChiralTypeFrom3D` would tag the atom.
fn nontetrahedral_from_3d(mol: &Mol, coords: &[[f64; 3]], a: usize) -> bool {
    if mol.atoms[a].anum < 15 {
        return false;
    }
    let cen = coords[a];
    let mut v: Vec<[f64; 3]> = Vec::new();
    for nb in mol.nbrs(a) {
        if v.len() == 6 {
            return false;
        }
        let d = sub(coords[nb], cen);
        let l = dot(d, d).sqrt();
        if l < 1.0e-16 {
            // `normalize` throws on a zero-length vector.
            return false;
        }
        v.push([d[0] / l, d[1] / l, d[2] / l]);
    }
    let count = v.len();
    if count < 3 {
        return false;
    }
    let mut pair = [0usize; 6];
    let mut pairs = 0;
    for i in 0..count {
        for j in i + 1..count {
            if dot(v[i], v[j]) < -(1.0 - 0.1) {
                if pair[i] != 0 || pair[j] != 0 {
                    return false;
                }
                pair[i] = j + 1;
                pair[j] = i + 1;
                pairs += 1;
            }
        }
    }
    matches!(
        (pairs, count),
        (1, 3) | (1, 4) | (1, 5) | (2, 4) | (2, 5) | (3, 6)
    )
}

/// `MolOps::assignChiralTypesFrom3D(mol, -1, replaceExistingTags=true)`.
pub(crate) fn assign_chiral_types_from_3d(
    mol: &mut Mol,
    coords: &[[f64; 3]],
) -> Result<(), RdkitSmilesError> {
    const ZERO_VOLUME_TOL: f64 = 0.1;
    for a in 0..mol.atoms.len() {
        mol.atoms[a].chiral = ChiralTag::Unspecified;
        let nz_degree = mol.atom_bonds[a]
            .iter()
            .filter(|&&b| {
                let bond = &mol.bonds[b];
                !(bond.bt == BondType::Dative && bond.begin == a)
            })
            .count();
        let tnz_degree = nz_degree + mol.total_num_hs(a) as usize;
        if nz_degree < 3 || tnz_degree > 6 {
            continue;
        }
        if nontetrahedral_from_3d(mol, coords, a) {
            return Err(parse_error("non-tetrahedral stereo from 3D coordinates"));
        }
        if tnz_degree > 4 {
            continue;
        }
        let anum = mol.atoms[a].anum;
        if anum != 16 && anum != 34 && tnz_degree != 4 {
            continue;
        }
        let p0 = coords[a];
        let nbrs: Vec<[f64; 3]> = mol.atom_bonds[a]
            .iter()
            .filter(|&&b| {
                let bond = &mol.bonds[b];
                !(bond.bt == BondType::Dative && bond.begin == a)
            })
            .map(|&b| coords[mol.bonds[b].other(a)])
            .collect();
        let v1 = sub(nbrs[0], p0);
        let v2 = sub(nbrs[1], p0);
        let v3 = sub(nbrs[2], p0);
        let mut vol = dot(v1, cross(v2, v3));
        if vol < -ZERO_VOLUME_TOL {
            mol.atoms[a].chiral = ChiralTag::Cw;
        } else if vol > ZERO_VOLUME_TOL {
            mol.atoms[a].chiral = ChiralTag::Ccw;
        } else if nbrs.len() == 4 {
            let v4 = sub(nbrs[3], p0);
            vol = -dot(v1, cross(v2, v4));
            if vol < -ZERO_VOLUME_TOL {
                mol.atoms[a].chiral = ChiralTag::Cw;
            } else if vol > ZERO_VOLUME_TOL {
                mol.atoms[a].chiral = ChiralTag::Ccw;
            }
        }
    }
    Ok(())
}

/// `MolFromPDBBlock(text, sanitize, removeHs, flavor, proximityBonding)`:
/// `None` where RDKit returns no molecule (no ATOM/HETATM/COMPND/HEADER
/// record).
pub(crate) fn mol_from_pdb_block(
    text: &str,
    sanitize: bool,
    remove_hs: bool,
    flavor: u32,
    proximity_bonding: bool,
) -> Result<Option<PdbMol>, RdkitSmilesError> {
    let bytes = text.as_bytes();
    let mut r = Reader {
        text: bytes,
        flavor,
        pdb: PdbMol {
            mol: Mol::default(),
            coords: Vec::new(),
            info: Vec::new(),
            is_3d: false,
        },
        amap: std::collections::HashMap::new(),
        bmap: std::collections::HashMap::new(),
        has_conformer: false,
    };
    let mut have_mol = false;
    let mut multi_conformer = false;
    let mut start = 0;
    while start < bytes.len() && bytes[start] != 0 {
        // The line runs to '\r', '\n' or the end (from the second byte).
        let mut next = start + 1;
        let len;
        loop {
            match bytes.get(next) {
                Some(b'\r') => {
                    len = next - start;
                    next += if bytes.get(next + 1) == Some(&b'\n') {
                        2
                    } else {
                        1
                    };
                    break;
                }
                Some(b'\n') => {
                    len = next - start;
                    next += 1;
                    break;
                }
                None | Some(0) => {
                    len = next - start;
                    break;
                }
                _ => next += 1,
            }
        }
        let head = field(bytes, start, 0, 6);
        if head == b"ATOM  " || head == b"HETATM" {
            have_mol = true;
            if !multi_conformer {
                r.atom_line(start, len)?;
            }
        } else if head == b"CONECT" {
            if have_mol && !multi_conformer {
                r.bond_line(start, len)?;
            }
        } else if head == b"COMPND" || head == b"HEADER" {
            have_mol = true;
        } else if head == b"ENDMDL" {
            if !have_mol {
                break;
            }
            multi_conformer = true;
        }
        start = next;
    }
    if !have_mol {
        return Ok(None);
    }
    let mut pdb = r.pdb;
    if proximity_bonding && r.has_conformer {
        connect_the_dots(&mut pdb);
    }
    if proximity_bonding || flavor & 8 != 0 {
        standard_residue_bond_orders(&mut pdb);
    }
    // BasicPDBCleanup: four-valent neutral N -> N+.
    for a in 0..pdb.mol.atoms.len() {
        let ev = pdb.mol.calc_explicit_valence(a, false)?;
        if pdb.mol.atoms[a].anum == 7 && pdb.mol.atoms[a].charge == 0 && ev == 4 {
            pdb.mol.atoms[a].charge = 1;
        }
    }
    if sanitize {
        if remove_hs {
            let kept = sanitize::remove_hs(&mut pdb.mol, false)?;
            pdb.coords = kept.iter().map(|&i| pdb.coords[i]).collect();
            pdb.info = kept.iter().map(|&i| pdb.info[i].clone()).collect();
        } else {
            sanitize::sanitize_keeping_hs(&mut pdb.mol)?;
        }
    } else {
        pdb.mol.update_property_cache(false)?;
    }
    if r.has_conformer && pdb.is_3d {
        assign_chiral_types_from_3d(&mut pdb.mol, &pdb.coords)?;
    }
    // StandardPDBResidueChirality
    for a in 0..pdb.mol.atoms.len() {
        let info = &pdb.info[a];
        if pdb.mol.atoms[a].chiral != ChiralTag::Unspecified
            && !info.hetero
            && !standard_chiral_atom(&info.res_name, info.name)
        {
            pdb.mol.atoms[a].chiral = ChiralTag::Unspecified;
        }
    }
    Ok(Some(pdb))
}

/// The chematic molecule that `parse::from_chematic` turns back into `m`
/// (then sanitized again by the RDKit-compatible writers): every heavy
/// atom carries its total hydrogen count when `m` is sanitized, hydrogen
/// atoms none (so they stay graph atoms), chiral tags are written relative
/// to `m`'s bond order with the hydrogen last.
pub(crate) fn to_chematic(
    m: &Mol,
    sanitized: bool,
) -> Result<chematic_core::Molecule, RdkitSmilesError> {
    use chematic_core::{
        Atom as CAtom, AtomIdx, BondOrder, Chirality, Element, MoleculeBuilder, STEREO_H_SENTINEL,
    };
    let mut b = MoleculeBuilder::new();
    let mut orders: Vec<(AtomIdx, Vec<u32>)> = Vec::new();
    for (i, a) in m.atoms.iter().enumerate() {
        let mut atom = if a.anum == 0 {
            CAtom::wildcard()
        } else {
            CAtom::new(
                Element::from_symbol(periodic::symbol(a.anum))
                    .ok_or_else(|| parse_error(format!("element {}", a.anum)))?,
            )
        };
        atom.isotope = (a.isotope != 0).then_some(a.isotope as u16);
        atom.charge = a.charge as i8;
        atom.aromatic = a.aromatic;
        let n_hs = if a.anum == 1 && a.num_explicit_hs == 0 {
            // A hydrogen graph atom: written without a count, so the
            // writers keep every H atom (the state after `Chem.AddHs`).
            0
        } else if sanitized {
            atom.hydrogen_count = Some(m.total_num_hs(i) as u8);
            m.total_num_hs(i)
        } else if a.no_implicit {
            atom.hydrogen_count = Some(a.num_explicit_hs as u8);
            a.num_explicit_hs
        } else if a.num_explicit_hs == 0 {
            0
        } else {
            return Err(parse_error(
                "explicit hydrogen count on an unsanitized atom",
            ));
        };
        // A bracket atom's radicals follow from its hydrogen count.
        if a.radicals != 0 && atom.hydrogen_count.is_none() {
            return Err(parse_error("radical"));
        }
        if a.chiral != ChiralTag::Unspecified {
            let mut order: Vec<u32> = m.nbrs(i).map(|x| x as u32).collect();
            let is_start = !m.nbrs(i).any(|x| x < i);
            let mut tag = a.chiral;
            if m.degree(i) == 3 && is_start && n_hs == 1 {
                tag = if tag == ChiralTag::Cw {
                    ChiralTag::Ccw
                } else {
                    ChiralTag::Cw
                };
            }
            if n_hs > 0 {
                order.push(STEREO_H_SENTINEL);
            }
            atom.chirality = if tag == ChiralTag::Ccw {
                Chirality::CounterClockwise
            } else {
                Chirality::Clockwise
            };
            orders.push((AtomIdx(i as u32), order));
        }
        b.add_atom(atom);
    }
    let mut aromatic_dirs = Vec::new();
    for (bi, bond) in m.bonds.iter().enumerate() {
        // `/` `\` directions relative to the begin atom (the first atom).
        let dir = match bond.dir {
            BondDir::None => None,
            BondDir::EndUpRight => Some(BondOrder::Up),
            BondDir::EndDownRight => Some(BondOrder::Down),
        };
        let order = match bond.bt {
            BondType::Single if dir.is_some() => dir.expect("set"),
            BondType::Aromatic if dir.is_some() => {
                aromatic_dirs.push((bi, dir.expect("set")));
                BondOrder::Aromatic
            }
            _ if dir.is_some() => return Err(parse_error("direction on a multiple bond")),
            BondType::Single => BondOrder::Single,
            BondType::Double => BondOrder::Double,
            BondType::Triple => BondOrder::Triple,
            BondType::Quadruple => BondOrder::Quadruple,
            BondType::Aromatic => BondOrder::Aromatic,
            BondType::Dative => BondOrder::Dative,
        };
        b.add_bond(AtomIdx(bond.begin as u32), AtomIdx(bond.end as u32), order)
            .map_err(|e| parse_error(e.to_string()))?;
    }
    for (idx, order) in orders {
        b.set_stereo_neighbor_order(idx, order);
    }
    for (bi, dir) in aromatic_dirs {
        b.set_bond_direction(chematic_core::BondIdx(bi as u32), dir);
    }
    Ok(b.build())
}

#[cfg(test)]
mod tests {
    #[test]
    fn reads_rdkit_written_ethanol() {
        // Chem.MolToPDBBlock of an embedded ethanol (heavy atoms).
        let block = "HETATM    1  C1  UNL     1      -0.889  -0.214   0.011  1.00  0.00           C  \n\
HETATM    2  C2  UNL     1       0.457   0.442  -0.252  1.00  0.00           C  \n\
HETATM    3  O1  UNL     1       1.440  -0.506   0.010  1.00  0.00           O  \n\
CONECT    1    2\nCONECT    2    3\nEND\n";
        let res = crate::rdkit_mol_from_pdb_block(block, true, true, 0, true)
            .unwrap()
            .unwrap();
        assert_eq!(res.smiles, "CCO");
        assert_eq!(res.coords[2], [1.440, -0.506, 0.010]);
        assert_eq!(crate::rdkit_canonical_smiles(&res.molecule).unwrap(), "CCO");
    }

    #[test]
    fn conect_duplicates_give_bond_orders() {
        let block = "HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  \n\
HETATM    2  O1  UNL     1       1.200   0.000   0.000  1.00  0.00           O  \n\
CONECT    1    2    2\nCONECT    2    1    1\nEND\n";
        let res = crate::rdkit_mol_from_pdb_block(block, true, true, 0, true)
            .unwrap()
            .unwrap();
        assert_eq!(res.smiles, "C=O");
    }

    #[test]
    fn xyz_block_reads_atoms_without_bonds() {
        let (mol, coords) =
            crate::rdkit_mol_from_xyz_block("2\ncomment\nC 0 0 0\nCL 1.0 2.0 3.5\n").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 0);
        assert_eq!(coords[1], [1.0, 2.0, 3.5]);
        assert!(crate::rdkit_mol_from_xyz_block("1\nc\nC 1.0 2.0\n").is_err());
    }
}
