use super::*;
pub fn parse_gpr(text: &str) -> Result<GprRecord, GprError> {
    parse_gpr_with_limits(text, &GprParseLimits::default())
}
pub fn parse_gpr_with_limits(text: &str, limits: &GprParseLimits) -> Result<GprRecord, GprError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input byte limit exceeded"));
    }
    let mut r = Reader {
        lines: text.lines().enumerate(),
        limits,
    };
    r.exact(&["!Header", "gpr", "100"])?;
    r.exact(&["!Info", "1"])?;
    let n = r.count("!Atoms", limits.max_atoms)?;
    if n == 0 {
        return Err(error(0, "at least one atom required"));
    }
    let atoms = read_atoms(&mut r, n)?;
    let nb = r.count("!Bonds", limits.max_bonds)?;
    let bonds = read_bonds(&mut r, nb)?;
    r.exact(&["!Coord"])?;
    let coordinates = read_coordinates(&mut r, n)?;
    let (line, tag) = r.fields(1)?;
    let partial_charges = match tag[0] {
        "!Charges" => {
            let mut rows = Vec::new();
            for _ in 0..n {
                let (line, f) = r.fields(2)?;
                real(f[1]).map_err(|e| error(line, e.detail))?;
                rows.push(GprPartialCharge {
                    atom_id: integer(f[0], line)?,
                    value: f[1].into(),
                });
            }
            r.exact(&["!End"])?;
            Some(rows)
        }
        "!End" => None,
        _ => {
            return Err(error(
                line,
                "expected !Charges or !End; unsupported section",
            ));
        }
    };
    if let Some((line, _)) = r.next()? {
        return Err(error(line, "trailing record or unsupported section"));
    }
    let record = GprRecord {
        atoms,
        bonds,
        coordinates,
        partial_charges,
    };
    record.validate_with_limits(limits)?;
    Ok(record)
}
fn read_atoms(r: &mut Reader<'_>, n: usize) -> Result<Vec<GprAtom>, GprError> {
    let mut atoms = Vec::new();
    for _ in 0..n {
        let (line, f) = r.fields(2)?;
        let z = integer(f[1], line)?;
        if !(1..=118).contains(&z) {
            return Err(error(line, "atomic number must be 1..118"));
        }
        atoms.push(GprAtom {
            id: integer(f[0], line)?,
            atomic_number: z as u8,
        });
    }
    Ok(atoms)
}
fn read_bonds(r: &mut Reader<'_>, nb: usize) -> Result<Vec<GprBond>, GprError> {
    let mut bonds = Vec::new();
    for _ in 0..nb {
        let (line, f) = r.fields(3)?;
        let order = match f[2] {
            "S" => GprBondOrder::S,
            "D" => GprBondOrder::D,
            "T" => GprBondOrder::T,
            "C" => GprBondOrder::C,
            _ => return Err(error(line, "unsupported bond code; expected S/D/T/C")),
        };
        bonds.push(GprBond {
            atoms: [integer(f[0], line)?, integer(f[1], line)?],
            order,
        });
    }
    Ok(bonds)
}
fn read_coordinates(r: &mut Reader<'_>, n: usize) -> Result<Vec<GprCoordinate>, GprError> {
    let mut coordinates = Vec::new();
    for _ in 0..n {
        let (line, f) = r.fields(4)?;
        let mut position_nm = [String::new(), String::new(), String::new()];
        for i in 0..3 {
            real(f[i + 1]).map_err(|e| error(line, e.detail))?;
            position_nm[i] = f[i + 1].into();
        }
        coordinates.push(GprCoordinate {
            atom_id: integer(f[0], line)?,
            position_nm,
        });
    }
    Ok(coordinates)
}
fn integer(token: &str, line: usize) -> Result<u32, GprError> {
    if token.is_empty() || !token.bytes().all(|b| b.is_ascii_digit()) {
        return Err(error(line, "expected unsigned decimal integer"));
    }
    token
        .parse()
        .map_err(|_| error(line, "integer outside u32 range"))
}
struct Reader<'a> {
    lines: std::iter::Enumerate<std::str::Lines<'a>>,
    limits: &'a GprParseLimits,
}
impl<'a> Reader<'a> {
    fn next(&mut self) -> Result<Option<(usize, Vec<&'a str>)>, GprError> {
        for (i, line) in self.lines.by_ref() {
            if line.len() > self.limits.max_line_len
                || line.chars().any(|c| c.is_control() && c != '\t')
            {
                return Err(error(i + 1, "line limit exceeded or control character"));
            }
            let fields: Vec<_> = line.split_ascii_whitespace().collect();
            if !fields.is_empty() {
                return Ok(Some((i + 1, fields)));
            }
        }
        Ok(None)
    }
    fn fields(&mut self, n: usize) -> Result<(usize, Vec<&'a str>), GprError> {
        let (line, f) = self
            .next()?
            .ok_or_else(|| error(0, "unexpected end of input"))?;
        if f.len() != n {
            return Err(error(line, format!("expected exactly {n} fields")));
        }
        Ok((line, f))
    }
    fn exact(&mut self, want: &[&str]) -> Result<(), GprError> {
        let (line, f) = self.fields(want.len())?;
        if f != want {
            return Err(error(line, format!("expected {}", want.join(" "))));
        }
        Ok(())
    }
    fn count(&mut self, tag: &str, max: usize) -> Result<usize, GprError> {
        let (line, f) = self.fields(2)?;
        if f[0] != tag {
            return Err(error(line, format!("expected {tag}")));
        }
        let n = integer(f[1], line)? as usize;
        if n > max {
            return Err(error(line, "declared count exceeds limit"));
        }
        Ok(n)
    }
}
