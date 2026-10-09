//! RDKit's `MolFromXYZBlock` (RDKit 2026.03.1
//! `FileParsers/XYZFileParser.cpp`): atoms and one conformer, no bonds.

use super::RdkitSmilesError;
use super::mol::{Atom, Mol};
use super::periodic;

fn fails(what: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Sanitization(what.into())
}

/// `std::getline` over a string stream plus RDKit's `getLine` (a trailing
/// `\r` dropped), tracking the stream's eof flag.
struct Lines<'a> {
    text: &'a str,
    pos: usize,
    eof: bool,
}

impl<'a> Lines<'a> {
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
}

/// `FileParserUtils::toUnsigned(s, acceptSpaces=true)`.
fn to_unsigned(s: &str) -> Result<u32, ()> {
    for c in s.bytes() {
        if !(c.is_ascii_digit() || c == b' ' || c == b'+') {
            return Err(());
        }
    }
    let t = s.trim_start_matches(' ');
    if t.is_empty() {
        return Ok(0);
    }
    let n = t.bytes().take_while(u8::is_ascii_digit).count();
    Ok(t[..n].parse::<u32>().unwrap_or(0))
}

/// `FileParserUtils::toDouble(s, acceptSpaces=false)`: validation, then
/// `atof` (the longest numeric prefix).
fn to_double(s: &str) -> Result<f64, ()> {
    for c in s.bytes() {
        if !(c.is_ascii_digit() || matches!(c, b'+' | b'-' | b',' | b'.')) {
            return Err(());
        }
    }
    let b = s.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let int_start = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    let mut digits = i > int_start;
    if i < b.len() && b[i] == b'.' {
        let mut j = i + 1;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if j > i + 1 || digits {
            digits = digits || j > i + 1;
            i = j;
        }
    }
    if !digits {
        return Ok(0.0);
    }
    let t = &s[..i];
    Ok(t.strip_prefix('+')
        .unwrap_or(t)
        .parse::<f64>()
        .unwrap_or(0.0))
}

fn is_ws(c: u8) -> bool {
    c == b' ' || c == b'\t'
}

/// `ParseXYZFileAtomLine`.
fn atom_line(line: &str, n: usize) -> Result<(u32, [f64; 3]), RdkitSmilesError> {
    let b = line.as_bytes();
    let mut delims = [0usize; 8];
    let mut prev = 0;
    for (i, d) in delims.iter_mut().enumerate().take(7) {
        let found = if i % 2 == 0 {
            (prev..b.len()).find(|&k| !is_ws(b[k]))
        } else {
            (prev..b.len()).find(|&k| is_ws(b[k]))
        };
        match found {
            Some(k) => *d = k,
            None => return Err(fails(format!("Missing coordinates on line {n}"))),
        }
        prev = *d;
    }
    delims[7] = (0..b.len())
        .rev()
        .find(|&k| !is_ws(b[k]))
        .map_or(0, |k| k + 1);
    let mut pos = [0.0; 3];
    for (k, p) in pos.iter_mut().enumerate() {
        let (s, e) = (delims[2 + 2 * k], delims[3 + 2 * k]);
        let field = line.get(s..e.max(s)).unwrap_or("");
        *p = to_double(field)
            .map_err(|_| fails(format!("Cannot convert '{field}' to double on line {n}")))?;
    }
    let mut symb: Vec<u8> = b[delims[0]..delims[1]].to_vec();
    if symb.len() == 2 && symb[1].is_ascii_uppercase() {
        symb[1] = symb[1].to_ascii_lowercase();
    }
    let symb = String::from_utf8_lossy(&symb).into_owned();
    let anum = periodic::atomic_number(&symb)
        .ok_or_else(|| fails(format!("Element '{symb}' not found")))?;
    Ok((anum, pos))
}

/// `MolFromXYZBlock` as `Chem.MolFromXYZBlock` runs it (an empty block
/// gives an empty molecule).
pub(crate) fn mol_from_xyz_block(text: &str) -> Result<(Mol, Vec<[f64; 3]>), RdkitSmilesError> {
    let mut lines = Lines {
        text,
        pos: 0,
        eof: false,
    };
    let num = lines.get_line();
    let num_atoms = to_unsigned(num).map_err(|_| {
        fails(format!(
            "Unable to recognize the number of atoms: cannot convert '{num}' to unsigned int on line 0"
        ))
    })?;
    let _comment = lines.get_line();
    let mut mol = Mol::default();
    let mut coords = Vec::with_capacity(num_atoms as usize);
    for i in 0..num_atoms as usize {
        if lines.eof {
            return Err(fails("EOF hit while reading atoms"));
        }
        let line = lines.get_line();
        let (anum, pos) = atom_line(line, i + 2)?;
        mol.add_atom(Atom::new(anum));
        coords.push(pos);
    }
    while !lines.eof {
        let extra = lines.get_line();
        if extra.bytes().any(|c| !is_ws(c)) {
            return Err(fails("More lines than expected"));
        }
    }
    Ok((mol, coords))
}
