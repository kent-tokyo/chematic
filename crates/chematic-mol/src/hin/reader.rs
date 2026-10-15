use super::*;
use std::collections::{BTreeMap, BTreeSet};

type PendingAtom = (HinAtom, Vec<(u32, BondOrder)>, usize);

fn invalid(line: usize, reason: &'static str) -> HinError {
    HinError::InvalidRecord { line, reason }
}
fn positive_id(token: &str, line: usize) -> Result<u32, HinError> {
    token
        .parse::<u32>()
        .ok()
        .filter(|&id| id > 0)
        .ok_or_else(|| invalid(line, "expected a positive integer ID"))
}
fn finite(token: &str, line: usize) -> Result<f64, HinError> {
    token
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| invalid(line, "charge and coordinates must be finite numbers"))
}
fn header(line: &str, number: usize) -> Result<(u32, String), HinError> {
    let rest = line[3..].trim_start();
    let (id, title) = rest
        .split_once(char::is_whitespace)
        .ok_or_else(|| invalid(number, "mol requires an ID and quoted title"))?;
    let title = title
        .trim()
        .strip_prefix('"')
        .and_then(|t| t.strip_suffix('"'))
        .filter(|t| !t.contains('"'))
        .ok_or_else(|| invalid(number, "mol title must be one quoted string"))?;
    Ok((positive_id(id, number)?, title.to_owned()))
}
fn atom(line: &str, number: usize, limits: HinParseLimits) -> Result<PendingAtom, HinError> {
    let fields: Vec<_> = line.split_whitespace().collect();
    if fields.len() < 11 {
        return Err(invalid(
            number,
            "atom requires identity, charge, coordinates and degree",
        ));
    }
    let count = fields[10]
        .parse::<usize>()
        .map_err(|_| invalid(number, "invalid connectivity count"))?;
    if count > limits.max_atoms || count > limits.max_bonds {
        return Err(HinError::ResourceLimit);
    }
    if count.checked_mul(2).and_then(|n| n.checked_add(11)) != Some(fields.len()) {
        return Err(invalid(
            number,
            "connectivity count does not match atom row",
        ));
    }
    let source_id = positive_id(fields[1], number)?;
    let element =
        Element::from_symbol(fields[3]).ok_or_else(|| invalid(number, "unknown element symbol"))?;
    let entry = HinAtom {
        source_id,
        name: fields[2].to_owned(),
        element,
        atom_type: fields[4].to_owned(),
        flags: fields[5].to_owned(),
        partial_charge: finite(fields[6], number)?,
        coordinates: [
            finite(fields[7], number)?,
            finite(fields[8], number)?,
            finite(fields[9], number)?,
        ],
    };
    let mut seen = BTreeSet::new();
    let mut neighbors = Vec::with_capacity(count);
    let (pairs, remainder) = fields[11..].as_chunks::<2>();
    debug_assert!(remainder.is_empty());
    for &[other_token, order_token] in pairs {
        let other = positive_id(other_token, number)?;
        if other == source_id || !seen.insert(other) {
            return Err(invalid(number, "self or duplicate connectivity reference"));
        }
        let order = match order_token {
            "s" => BondOrder::Single,
            "d" => BondOrder::Double,
            "t" => BondOrder::Triple,
            "a" => BondOrder::Aromatic,
            _ => return Err(invalid(number, "unsupported bond token")),
        };
        neighbors.push((other, order));
    }
    Ok((entry, neighbors, number))
}
fn bonds(atoms: &[PendingAtom], limits: HinParseLimits) -> Result<Vec<HinBond>, HinError> {
    let ids: BTreeMap<_, _> = atoms
        .iter()
        .enumerate()
        .map(|(i, a)| (a.0.source_id, i))
        .collect();
    let mut directed = BTreeMap::new();
    for (i, (_, neighbors, line)) in atoms.iter().enumerate() {
        for &(id, order) in neighbors {
            let &j = ids
                .get(&id)
                .ok_or_else(|| invalid(*line, "dangling atom reference"))?;
            directed.insert((i, j), (order, *line));
        }
    }
    let mut result = Vec::new();
    for (&(i, j), &(order, line)) in &directed {
        if directed.get(&(j, i)).map(|&(o, _)| o) != Some(order) {
            return Err(invalid(
                line,
                "connectivity must be reciprocal with matching bond types",
            ));
        }
        if i < j {
            if result.len() >= limits.max_bonds {
                return Err(HinError::ResourceLimit);
            }
            result.push(HinBond {
                atom1: i,
                atom2: j,
                order,
            });
        }
    }
    Ok(result)
}

pub(super) fn parse(source: &str, limits: HinParseLimits) -> Result<HinRecord, HinError> {
    if source.len() > limits.max_bytes {
        return Err(HinError::ResourceLimit);
    }
    let mut identity = None;
    let mut ended = false;
    let mut atoms = Vec::new();
    let mut ids = BTreeSet::new();
    let mut directed_count = 0usize;
    for (offset, row) in source.lines().enumerate() {
        let number = offset + 1;
        let row = row.trim();
        if row.is_empty() || row.starts_with(';') {
            continue;
        }
        if ended {
            return Err(invalid(
                number,
                "content after endmol or multiple molecules",
            ));
        }
        match row.split_whitespace().next().unwrap() {
            "mol" if identity.is_none() => identity = Some(header(row, number)?),
            "atom" if identity.is_some() => {
                if atoms.len() >= limits.max_atoms.min(u32::MAX as usize) {
                    return Err(HinError::ResourceLimit);
                }
                let entry = atom(row, number, limits)?;
                if !ids.insert(entry.0.source_id) {
                    return Err(invalid(number, "duplicate atom ID"));
                }
                directed_count = directed_count.saturating_add(entry.1.len());
                if directed_count > limits.max_bonds.saturating_mul(2) {
                    return Err(HinError::ResourceLimit);
                }
                atoms.push(entry);
            }
            "endmol" if identity.is_some() => {
                let fields: Vec<_> = row.split_whitespace().collect();
                if fields.len() != 2
                    || positive_id(fields[1], number)? != identity.as_ref().unwrap().0
                {
                    return Err(invalid(number, "endmol ID does not match mol"));
                }
                ended = true;
            }
            _ => return Err(invalid(number, "unsupported directive or record order")),
        }
    }
    let (molecule_id, title) = identity.ok_or_else(|| invalid(1, "missing mol header"))?;
    if !ended || atoms.is_empty() {
        return Err(invalid(
            source.lines().count().max(1),
            "missing endmol or empty molecule",
        ));
    }
    let bonds = bonds(&atoms, limits)?;
    Ok(HinRecord {
        source: source.to_owned(),
        molecule_id,
        title,
        atoms: atoms.into_iter().map(|a| a.0).collect(),
        bonds,
    })
}
