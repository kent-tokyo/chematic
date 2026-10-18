// Ported from the Avalon Cheminformatics Toolkit, Copyright 2001-2011 Novartis
// Pharma AG, BSD-3-Clause license; see THIRD_PARTY_NOTICES.md.
//! The part of the Avalon toolkit's MOL-file reader (`reaccsio.c`:
//! `MolStr2Mol`, `ReadREACCSAtom`, `ReadREACCSBond`, `ReadProperties`,
//! `ReadV30Atom`, `ReadV30Bond`) that decides the fields the fingerprint
//! reads: atom symbols, charges, radicals, R-atom texts and bond types.

use super::{AvalonAtom, AvalonBond, AvalonMolecule};

/// C `sscanf("%d")` on a fixed-width field: leading blanks skipped, optional
/// sign, digits; `None` when no digits.
fn scan_int(s: &str) -> Option<i32> {
    let t = s.trim_start();
    let bytes = t.as_bytes();
    let mut i = 0;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        i += 1;
    }
    let start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == start {
        return None;
    }
    t[..i].parse().ok()
}

fn field(line: &str, from: usize, len: usize) -> &str {
    let end = (from + len).min(line.len());
    if from >= line.len() {
        ""
    } else {
        &line[from..end]
    }
}

/// `SplitChargeRadical`: MDL's combined charge/radical code.
fn split_charge_radical(code: i32) -> (i32, i32) {
    match code {
        0 => (0, 0),
        4 => (0, 2),
        c => (4 - c, 0),
    }
}

/// Reads a V2000 or V3000 MOL block the way `MolStr2Mol` does. Returns
/// `None` where it fails.
pub fn read_molblock(block: &str) -> Option<AvalonMolecule> {
    let lines: Vec<&str> = block.lines().map(|l| l.trim_end_matches('\r')).collect();
    if lines.len() < 4 {
        return None;
    }
    let counts = lines[3];
    let version = if counts.len() > 34 {
        let v = field(counts, 34, 5);
        if v.contains("V200") { "V2000" } else { v }
    } else {
        ""
    };
    if version == "V3000" {
        return read_v3000(&lines[4..]);
    }
    let n_atoms = scan_int(field(counts, 0, 3)).unwrap_or(0).max(0) as usize;
    let n_bonds = scan_int(field(counts, 3, 3)).unwrap_or(0).max(0) as usize;
    let mut pos = 4;
    let mut atoms = Vec::with_capacity(n_atoms);
    for _ in 0..n_atoms {
        let line = *lines.get(pos)?;
        pos += 1;
        // `sscanf("%10f%10f%10f %s")`: the token starting at column 31.
        let symbol = line
            .get(31..)
            .and_then(|s| s.split_whitespace().next())
            .unwrap_or("?");
        let code = scan_int(field(line, 36, 3)).unwrap_or(0);
        let (charge, radical) = split_charge_radical(code);
        atoms.push(AvalonAtom {
            symbol: symbol.chars().take(3).collect(),
            charge,
            radical,
            atext: String::new(),
        });
    }
    let mut bonds = Vec::with_capacity(n_bonds);
    for _ in 0..n_bonds {
        let line = *lines.get(pos)?;
        pos += 1;
        let a0 = scan_int(&field(line, 0, 3).replace(' ', ""))?;
        let a1 = scan_int(&field(line, 3, 3).replace(' ', ""))?;
        let bt = scan_int(&field(line, 6, 3).replace(' ', "")).unwrap_or(0);
        bonds.push(AvalonBond {
            atoms: [a0, a1],
            bond_type: bt,
        });
    }
    let mut mol = AvalonMolecule { atoms, bonds };
    while pos < lines.len() {
        let line = lines[pos];
        pos += 1;
        if line.starts_with("M  END") {
            break;
        }
        let is_chg = line.starts_with("M  CHG");
        let is_rad = line.starts_with("M  RAD");
        if is_chg || is_rad {
            let n = scan_int(field(line, 6, 3)).unwrap_or(0).max(0) as usize;
            let vals: Vec<i32> = line
                .get(9..)
                .unwrap_or("")
                .split_whitespace()
                .filter_map(|t| t.parse().ok())
                .collect();
            for k in 0..n.min(8) {
                let (Some(&a), Some(&v)) = (vals.get(2 * k), vals.get(2 * k + 1)) else {
                    break;
                };
                if let Some(atom) = mol.atoms.get_mut((a - 1) as usize) {
                    if is_chg {
                        atom.charge = v;
                    } else {
                        atom.radical = v;
                    }
                }
            }
        } else if let Some(rest) = line.strip_prefix("A  ") {
            // atom text line pair; kept only as the label of an "R" atom
            if let Some(a) = scan_int(rest)
                && a >= 1
                && (a as usize) <= mol.atoms.len()
                && mol.atoms[(a - 1) as usize].symbol == "R"
                && let Some(text) = lines.get(pos)
            {
                mol.atoms[(a - 1) as usize].atext = text.chars().take(80).collect();
            }
            pos += 1;
        } else if line.starts_with("G  ") {
            pos += 1;
        }
    }
    Some(mol)
}

fn read_v3000(lines: &[&str]) -> Option<AvalonMolecule> {
    let mut atoms: Vec<AvalonAtom> = Vec::new();
    let mut bonds: Vec<AvalonBond> = Vec::new();
    let mut section = "";
    for line in lines {
        if *line == "M  END" {
            break;
        }
        let Some(rest) = line.strip_prefix("M  V30 ") else {
            continue;
        };
        match rest.trim() {
            "BEGIN ATOM" => section = "atom",
            "BEGIN BOND" => section = "bond",
            "END ATOM" | "END BOND" => section = "",
            _ => {
                let toks: Vec<&str> = rest.split_whitespace().collect();
                if section == "atom" {
                    if toks.len() < 6 {
                        return None;
                    }
                    let mut atom = AvalonAtom {
                        symbol: toks[1].to_string(),
                        charge: 0,
                        radical: 0,
                        atext: String::new(),
                    };
                    for p in toks.iter().skip(6).take(10) {
                        if let Some(v) = p.strip_prefix("CHG=") {
                            atom.charge = scan_int(v).unwrap_or(0);
                        } else if let Some(v) = p.strip_prefix("RAD=") {
                            atom.radical = scan_int(v).unwrap_or(0);
                        }
                    }
                    atoms.push(atom);
                } else if section == "bond" {
                    if toks.len() < 4 {
                        return None;
                    }
                    bonds.push(AvalonBond {
                        atoms: [scan_int(toks[2])?, scan_int(toks[3])?],
                        bond_type: scan_int(toks[1])?,
                    });
                }
            }
        }
    }
    Some(AvalonMolecule { atoms, bonds })
}

#[cfg(test)]
mod boundary_contract_tests {
    use super::*;

    #[test]
    fn fixed_width_integer_and_charge_helpers_follow_c_scanning_rules() {
        assert_eq!(scan_int("  -12tail"), Some(-12));
        assert_eq!(scan_int(" +7 "), Some(7));
        assert_eq!(scan_int(" + "), None);
        assert_eq!(field("abc", 9, 2), "");
        assert_eq!(field("abcdef", 2, 99), "cdef");
        assert_eq!(split_charge_radical(0), (0, 0));
        assert_eq!(split_charge_radical(4), (0, 2));
        assert_eq!(split_charge_radical(1), (3, 0));
    }

    #[test]
    fn v2000_properties_atom_text_and_group_lines_are_bounded() {
        assert!(read_molblock("too\nshort\n").is_none());
        let block = concat!(
            "name\nprogram\ncomment\n",
            "  2  1  0  0  0  0            999 V2000\n",
            "    0.0000    0.0000    0.0000 R   0  0  0  0  0  0\n",
            "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0\n",
            "  1  2  1\n",
            "M  CHG  1   2  -1\n",
            "M  RAD  1   2   2\n",
            "A    1\n",
            "R12\n",
            "G    2\n",
            "ignored group text\n",
            "M  END\n",
        );
        let molecule = read_molblock(block).unwrap();
        assert_eq!(molecule.atoms.len(), 2);
        assert_eq!(molecule.bonds.len(), 1);
        assert_eq!(molecule.atoms[0].atext, "R12");
        assert_eq!(molecule.atoms[1].charge, -1);
        assert_eq!(molecule.atoms[1].radical, 2);
    }

    #[test]
    fn v3000_requires_complete_records_and_preserves_charge_and_radical() {
        let valid = "name\nprogram\ncomment\n  0  0  0  0  0  0            999 V3000\n\
M  V30 BEGIN ATOM\n\
M  V30 1 N 0 0 0 0 CHG=1 RAD=2\n\
M  V30 2 C 1 0 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 END BOND\n\
M  END\n";
        let molecule = read_molblock(valid).unwrap();
        assert_eq!(
            (molecule.atoms[0].charge, molecule.atoms[0].radical),
            (1, 2)
        );
        assert_eq!(molecule.bonds[0].atoms, [1, 2]);

        for body in [
            "M  V30 BEGIN ATOM\nM  V30 1 C 0 0\nM  END\n",
            "M  V30 BEGIN BOND\nM  V30 1 1 1\nM  END\n",
            "M  V30 BEGIN BOND\nM  V30 1 x 1 2\nM  END\n",
        ] {
            let block =
                format!("name\nprogram\ncomment\n  0  0  0  0  0  0            999 V3000\n{body}");
            assert!(read_molblock(&block).is_none(), "{body}");
        }
    }
}
