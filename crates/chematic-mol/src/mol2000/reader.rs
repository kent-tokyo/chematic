use super::*;

/// Parse a MOL V2000 string, running stereo perception and returning every
/// rejected wedge/hash center as a structured [`StereoDiagnostic`].
///
/// The parser follows the MDL/CTfile fixed-width column layout. `coords[i]`
/// is the `(x, y)` position for atom `i` extracted from the atom block.
/// Z-coordinates are discarded. This is the one parsing core for V2000 MOL
/// text -- [`parse_mol_with_coords`]/[`parse_mol`] are thin wrappers that
/// discard `stereo_diagnostics`, and [`crate::sdf`]'s readers delegate here
/// per record.
///
/// Local tetrahedral parity (`Atom.chirality` + `Molecule::stereo_neighbor_order`)
/// is perceived unconditionally whenever a wedge/hash bond is present --
/// mirroring RDKit's own `assignChiralTypesFromBondDirs`, which runs
/// regardless of a `sanitize`-equivalent flag. It never touches `Atom.cip_code`
/// and never depends on CIP ranking.
pub fn read_mol_with_diagnostics(input: &str) -> Result<MolReadReport, MolParseError> {
    read_mol_internal(input, true)
}

/// Parse a MOL V2000 block without optional stereo/3D diagnostics.
///
/// This is used by the RDKit-compatible `SDMolSupplier` path, whose contract
/// is the molecule graph plus metadata and SD properties. Keeping the parser
/// and fixed-width validation shared with [`read_mol_with_diagnostics`] avoids
/// semantic drift while skipping perception and geometry work.
#[allow(clippy::type_complexity)]
pub(crate) fn parse_mol_fast(input: &str) -> Result<(Molecule, MolMetadata), MolParseError> {
    let report = read_mol_internal(input, false)?;
    Ok((report.mol, report.metadata))
}

type NumberedLines<'a> = std::iter::Enumerate<std::str::Lines<'a>>;

fn next_required<'a>(lines: &mut NumberedLines<'a>) -> Result<(usize, &'a str), MolParseError> {
    lines
        .next()
        .map(|(index, line)| (index + 1, line))
        .ok_or(MolParseError::UnexpectedEnd)
}

fn parse_header(
    lines: &mut NumberedLines<'_>,
) -> Result<(MolMetadata, CoordinateDimension), MolParseError> {
    let name = next_required(lines)?.1.to_string();
    let line2 = next_required(lines)?.1;
    let comment = next_required(lines)?.1.to_string();
    Ok((
        MolMetadata {
            name,
            comment,
            v3000_sgroups: Vec::new(),
            v3000_bond_properties: Vec::new(),
            v3000_atom_properties: Vec::new(),
        },
        parse_dimension_code(line2),
    ))
}

fn parse_counts(lines: &mut NumberedLines<'_>) -> Result<(usize, usize), MolParseError> {
    let (line_number, line) = next_required(lines)?;
    if !line.contains("V2000") {
        return Err(MolParseError::InvalidCountLine {
            line: line_number,
            detail: "missing V2000 version tag".to_string(),
        });
    }
    let make_error = |line, detail| MolParseError::InvalidCountLine { line, detail };
    let atoms = parse_field3(line, 0, line_number, make_error)?;
    let bonds = parse_field3(line, 3, line_number, make_error)?;
    if atoms > MAX_ATOMS || bonds > MAX_BONDS {
        let (kind, count, maximum) = if atoms > MAX_ATOMS {
            ("atom", atoms, MAX_ATOMS)
        } else {
            ("bond", bonds, MAX_BONDS)
        };
        return Err(MolParseError::InvalidCountLine {
            line: line_number,
            detail: format!("{kind} count {count} exceeds maximum allowed {maximum}"),
        });
    }
    Ok((atoms, bonds))
}

fn parse_z_coordinate(line: &str, line_number: usize) -> Result<f64, MolParseError> {
    let Some(raw) = line.get(20..30) else {
        return Ok(0.0);
    };
    let trimmed = raw.trim();
    if trimmed.is_empty() || trimmed == "0.0000" || trimmed == "0" {
        return Ok(0.0);
    }
    let value: f64 = trimmed
        .parse()
        .map_err(|_| MolParseError::InvalidAtomLine {
            line: line_number,
            detail: format!("cannot parse z coordinate from '{trimmed}'"),
        })?;
    value
        .is_finite()
        .then_some(value)
        .ok_or_else(|| MolParseError::InvalidAtomLine {
            line: line_number,
            detail: format!("z coordinate is not finite (NaN/Infinite): '{trimmed}'"),
        })
}

fn parse_atom_identity(
    symbol: &str,
    line_number: usize,
) -> Result<(Atom, Option<RGroupLabel>), MolParseError> {
    let trimmed = symbol.trim();
    if trimmed == "*" {
        return Ok((Atom::wildcard(), None));
    }
    if trimmed == "R" || trimmed == "R#" {
        return Ok((Atom::wildcard(), Some(RGroupLabel::unnumbered())));
    }
    if let Some(number) = trimmed.strip_prefix('R')
        && !number.is_empty()
        && number.bytes().all(|byte| byte.is_ascii_digit())
    {
        let label = number
            .parse::<u16>()
            .ok()
            .and_then(RGroupLabel::numbered)
            .ok_or_else(|| MolParseError::UnknownElement {
                symbol: trimmed.to_string(),
                line: line_number,
            })?;
        return Ok((Atom::wildcard(), Some(label)));
    }
    let element = match symbol {
        "C  " => Some(Element::C),
        "N  " => Some(Element::N),
        "O  " => Some(Element::O),
        "H  " => Some(Element::H),
        "S  " => Some(Element::S),
        "P  " => Some(Element::P),
        "F  " => Some(Element::F),
        "Cl " => Some(Element::CL),
        "Br " => Some(Element::BR),
        _ => Element::from_symbol(trimmed),
    }
    .ok_or_else(|| MolParseError::UnknownElement {
        symbol: trimmed.to_string(),
        line: line_number,
    })?;
    Ok((Atom::new(element), None))
}

struct ParsedAtom {
    atom: Atom,
    r_group: Option<RGroupLabel>,
    coords: (f64, f64),
    z: f64,
    valence: Option<u8>,
    default_valence: bool,
    parity: u8,
}

fn parse_atom_line(
    line: &str,
    line_number: usize,
    atom_index: usize,
    include_diagnostics: bool,
) -> Result<ParsedAtom, MolParseError> {
    let z = parse_z_coordinate(line, line_number)?;
    let coords = if include_diagnostics {
        let parse = |range: std::ops::Range<usize>| {
            line.get(range)
                .and_then(|value| value.trim().parse().ok())
                .unwrap_or(0.0)
        };
        (parse(0..10), parse(10..20))
    } else {
        (0.0, 0.0)
    };
    let symbol = line
        .get(31..34)
        .ok_or_else(|| MolParseError::InvalidAtomLine {
            line: line_number,
            detail: format!("atom line {atom_index} too short for element field"),
        })?;
    let (mut atom, r_group) = parse_atom_identity(symbol, line_number)?;
    atom.charge = line
        .get(36..39)
        .map(|field| {
            if field == "  0" {
                return 0;
            }
            let code = parse_unsigned_ascii(field)
                .and_then(|value| i8::try_from(value).ok())
                .unwrap_or_else(|| field.trim().parse().unwrap_or(0));
            decode_charge(code)
        })
        .unwrap_or(0);
    atom.isotope = line
        .get(34..36)
        .and_then(|field| decode_mass_difference(atom.element, field));
    atom.atom_map = line
        .get(60..63)
        .and_then(parse_unsigned_ascii)
        .and_then(|value| u16::try_from(value).ok())
        .filter(|value| *value != 0);
    let valence_field = line.get(48..51).and_then(parse_unsigned_ascii);
    let valence = valence_field
        .filter(|value| *value != 0 && *value <= 15)
        .map(|value| if value == 15 { 0 } else { value as u8 });
    let parity = line
        .get(39..42)
        .and_then(parse_unsigned_ascii)
        .and_then(|value| u8::try_from(value).ok())
        .unwrap_or(0);
    Ok(ParsedAtom {
        atom,
        r_group,
        coords,
        z,
        valence,
        default_valence: valence_field == Some(0),
        parity,
    })
}

struct AtomBlock {
    builder: MoleculeBuilder,
    coords: Vec<(f64, f64)>,
    raw_z: Vec<f64>,
    valences: Vec<(AtomIdx, u8)>,
    default_valences: Vec<AtomIdx>,
    parities: Vec<(AtomIdx, u8)>,
}

fn parse_atom_block(
    lines: &mut NumberedLines<'_>,
    atom_count: usize,
    bond_count: usize,
    include_diagnostics: bool,
) -> Result<AtomBlock, MolParseError> {
    let mut block = AtomBlock {
        builder: MoleculeBuilder::with_capacity(atom_count, bond_count),
        coords: Vec::with_capacity(if include_diagnostics { atom_count } else { 0 }),
        raw_z: Vec::with_capacity(if include_diagnostics { atom_count } else { 0 }),
        valences: Vec::new(),
        default_valences: Vec::new(),
        parities: Vec::new(),
    };
    for atom_index in 0..atom_count {
        let (line_number, line) = next_required(lines)?;
        let parsed = parse_atom_line(line, line_number, atom_index, include_diagnostics)?;
        if include_diagnostics {
            block.coords.push(parsed.coords);
            block.raw_z.push(parsed.z);
        }
        let index = block.builder.add_atom(parsed.atom);
        if parsed.parity == 1 || parsed.parity == 2 {
            block.parities.push((index, parsed.parity));
        }
        if let Some(valence) = parsed.valence {
            block.valences.push((index, valence));
        } else if parsed.default_valence {
            block.default_valences.push(index);
        }
        if let Some(label) = parsed.r_group {
            block.builder.set_r_group(index, label);
        }
    }
    Ok(block)
}

fn parse_bond_line(
    builder: &mut MoleculeBuilder,
    line: &str,
    line_number: usize,
    bond_index: usize,
) -> Result<(BondIdx, bool), MolParseError> {
    let make_error = |line, detail| MolParseError::InvalidBondLine { line, detail };
    let atom1 = parse_field3(line, 0, line_number, make_error)?;
    let atom2 = parse_field3(line, 3, line_number, make_error)?;
    let bond_type = parse_field3(line, 6, line_number, make_error)?;
    if atom1 == 0 || atom2 == 0 {
        return Err(MolParseError::InvalidBondLine {
            line: line_number,
            detail: format!("bond {bond_index}: atom indices are 1-based; got {atom1}/{atom2}"),
        });
    }
    let stereo = if line.len() >= 12 {
        parse_field3(line, 9, line_number, make_error).unwrap_or(0)
    } else {
        0
    };
    let order = match (bond_type, stereo) {
        (1, 1) => BondOrder::Up,
        (1, 6) => BondOrder::Down,
        (1, _) => BondOrder::Single,
        (2, _) => BondOrder::Double,
        (3, _) => BondOrder::Triple,
        (4, _) => BondOrder::Aromatic,
        (5, _) => BondOrder::QuerySingleOrDouble,
        (6, _) => BondOrder::QuerySingleOrAromatic,
        (7, _) => BondOrder::QueryDoubleOrAromatic,
        (8, _) => BondOrder::QueryAny,
        (9, _) => BondOrder::Dative,
        _ => BondOrder::Single,
    };
    let index = builder
        .add_bond(
            AtomIdx((atom1 - 1) as u32),
            AtomIdx((atom2 - 1) as u32),
            order,
        )
        .map_err(|error| MolParseError::InvalidBondLine {
            line: line_number,
            detail: format!("bond {bond_index}: {error}"),
        })?;
    Ok((index, bond_type == 2 && stereo == 3))
}

fn parse_bond_block(
    lines: &mut NumberedLines<'_>,
    builder: &mut MoleculeBuilder,
    bond_count: usize,
) -> Result<std::collections::HashSet<BondIdx>, MolParseError> {
    let mut unspecified_ez = std::collections::HashSet::new();
    for bond_index in 0..bond_count {
        let (line_number, line) = next_required(lines)?;
        let (index, unspecified) = parse_bond_line(builder, line, line_number, bond_index)?;
        if unspecified {
            unspecified_ez.insert(index);
        }
    }
    Ok(unspecified_ez)
}

#[derive(Default)]
struct PropertyUpdates {
    charges: Vec<(AtomIdx, i8)>,
    isotopes: Vec<(AtomIdx, u16)>,
    r_groups: Vec<(AtomIdx, RGroupLabel)>,
}

fn property_pairs<'a>(
    fields: &'a [&str],
    line_number: usize,
    tag: &str,
) -> Result<impl Iterator<Item = &'a [&'a str]>, MolParseError> {
    let count = fields
        .get(2)
        .and_then(|value| value.parse::<usize>().ok())
        .ok_or_else(|| MolParseError::InvalidPropertyLine {
            line: line_number,
            detail: format!("M {tag} lacks a valid entry count"),
        })?;
    let expected = 3usize.checked_add(count.saturating_mul(2)).ok_or_else(|| {
        MolParseError::InvalidPropertyLine {
            line: line_number,
            detail: format!("M {tag} entry count overflows"),
        }
    })?;
    if fields.len() != expected {
        return Err(MolParseError::InvalidPropertyLine {
            line: line_number,
            detail: format!(
                "M {tag} declares {count} entries but has {} values",
                fields.len().saturating_sub(3)
            ),
        });
    }
    Ok(fields[3..].chunks(2))
}

fn property_atom(
    raw: &str,
    line_number: usize,
    tag: &str,
    atom_count: usize,
) -> Result<AtomIdx, MolParseError> {
    raw.parse::<usize>()
        .ok()
        .filter(|id| *id > 0 && *id <= atom_count)
        .map(|id| AtomIdx((id - 1) as u32))
        .ok_or_else(|| MolParseError::InvalidPropertyLine {
            line: line_number,
            detail: format!("M {tag} atom id '{raw}' is outside 1..={atom_count}"),
        })
}

fn parse_charge_properties(
    fields: &[&str],
    line_number: usize,
    atom_count: usize,
    updates: &mut PropertyUpdates,
) -> Result<(), MolParseError> {
    for pair in property_pairs(fields, line_number, "CHG")? {
        let atom = property_atom(pair[0], line_number, "CHG", atom_count)?;
        let charge = pair[1]
            .parse::<i8>()
            .map_err(|_| MolParseError::InvalidPropertyLine {
                line: line_number,
                detail: format!("M CHG charge '{}' is not an i8", pair[1]),
            })?;
        updates.charges.push((atom, charge));
    }
    Ok(())
}

fn parse_isotope_properties(
    fields: &[&str],
    line_number: usize,
    atom_count: usize,
    updates: &mut PropertyUpdates,
) -> Result<(), MolParseError> {
    for pair in property_pairs(fields, line_number, "ISO")? {
        let atom = property_atom(pair[0], line_number, "ISO", atom_count)?;
        let isotope = pair[1]
            .parse::<u16>()
            .ok()
            .filter(|value| *value != 0)
            .ok_or_else(|| MolParseError::InvalidPropertyLine {
                line: line_number,
                detail: format!("M ISO mass '{}' must be a positive u16", pair[1]),
            })?;
        updates.isotopes.push((atom, isotope));
    }
    Ok(())
}

fn parse_rgroup_properties(
    fields: &[&str],
    line_number: usize,
    atom_count: usize,
    updates: &mut PropertyUpdates,
) -> Result<(), MolParseError> {
    for pair in property_pairs(fields, line_number, "RGP")? {
        let atom = property_atom(pair[0], line_number, "RGP", atom_count)?;
        let label = pair[1]
            .parse::<u16>()
            .ok()
            .and_then(RGroupLabel::numbered)
            .ok_or_else(|| MolParseError::InvalidPropertyLine {
                line: line_number,
                detail: format!("M RGP number '{}' must be in 1..=9999", pair[1]),
            })?;
        updates.r_groups.push((atom, label));
    }
    Ok(())
}

fn parse_property_block(
    lines: &mut NumberedLines<'_>,
    atom_count: usize,
) -> Result<PropertyUpdates, MolParseError> {
    let mut updates = PropertyUpdates::default();
    for (index, line) in lines.by_ref() {
        let line_number = index + 1;
        let fields: Vec<&str> = line.split_whitespace().collect();
        match fields.get(0..2) {
            Some(["M", "CHG"]) => {
                parse_charge_properties(&fields, line_number, atom_count, &mut updates)?
            }
            Some(["M", "ISO"]) => {
                parse_isotope_properties(&fields, line_number, atom_count, &mut updates)?
            }
            Some(["M", "RGP"]) => {
                parse_rgroup_properties(&fields, line_number, atom_count, &mut updates)?
            }
            _ => {}
        }
        if line.trim_start().starts_with("M  END") {
            break;
        }
    }
    Ok(updates)
}

fn ctfile_default_valence(element: Element, bond_sum: u8) -> Option<u8> {
    let allowed: &[u8] = match element {
        Element::NA | Element::K => &[1],
        Element::SE => &[2, 4, 6],
        Element::SN => &[2, 4],
        Element::AS => &[3, 5],
        Element::SI => &[4],
        _ => &[],
    };
    allowed.iter().copied().find(|&valence| valence >= bond_sum)
}

fn apply_atom_updates(
    mol: &mut Molecule,
    valences: Vec<(AtomIdx, u8)>,
    default_valences: Vec<AtomIdx>,
    updates: PropertyUpdates,
) {
    flag_aromatic_bond_atoms(mol);
    for (atom, charge) in updates.charges {
        mol.set_charge(atom, charge);
    }
    for (atom, isotope) in updates.isotopes {
        mol.set_isotope(atom, Some(isotope));
    }
    for (atom, valence) in valences {
        let aromatic = mol
            .neighbors(atom)
            .any(|(_, bond)| mol.bond(bond).order == BondOrder::Aromatic);
        let bonds = chematic_core::valence::bond_order_sum(mol, atom);
        if !aromatic && valence >= bonds {
            mol.set_hydrogen_count(atom, Some(valence - bonds));
        }
    }
    // A CTfile valence code of zero requests the element's default valence.
    // These non-organic-subset elements do not receive implicit H from the
    // SMILES valence model, so materialize the CTfile-inferred count here.
    for atom in default_valences {
        if mol.atom(atom).charge != 0 {
            continue;
        }
        let aromatic = mol
            .neighbors(atom)
            .any(|(_, bond)| mol.bond(bond).order == BondOrder::Aromatic);
        if aromatic {
            continue;
        }
        let bonds = chematic_core::valence::bond_order_sum(mol, atom);
        let Some(valence) = ctfile_default_valence(mol.atom(atom).element, bonds) else {
            continue;
        };
        mol.set_hydrogen_count(atom, Some(valence - bonds));
    }
    for (atom, label) in updates.r_groups {
        mol.set_r_group(atom, label);
    }
}

fn graph_only_report(mol: Molecule, metadata: MolMetadata) -> MolReadReport {
    MolReadReport {
        mol,
        metadata,
        coords: Vec::new(),
        stereo_diagnostics: Vec::new(),
        ez_diagnostics: Vec::new(),
        conformer: None,
        coordinate_dimension: CoordinateDimension::Unknown,
        geometry_rank: GeometryRank::Indeterminate,
        stereo3d_diagnostics: Vec::new(),
        square_planar_diagnostics: Vec::new(),
    }
}

fn geometry_diagnostics(
    mol: &mut Molecule,
    coords: &[(f64, f64)],
    raw_z: &[f64],
    coordinate_dimension: CoordinateDimension,
) -> (
    Option<Coords3D>,
    GeometryRank,
    Vec<Stereo3DDiagnostic>,
    Vec<SquarePlanarPerceptionDiagnostic>,
) {
    let points = coords
        .iter()
        .zip(raw_z)
        .map(|(&(x, y), &z)| Point3::new(x, y, z))
        .collect::<Vec<_>>();
    let geometry_rank = classify_geometry_rank(&points);
    let conformer = match geometry_rank {
        GeometryRank::Coplanar | GeometryRank::ThreeD => Some(Coords3D { points }),
        GeometryRank::FlatZero | GeometryRank::Indeterminate => None,
    };
    let mut diagnostics = match (coordinate_dimension, geometry_rank) {
        (CoordinateDimension::TwoD, GeometryRank::Coplanar | GeometryRank::ThreeD) => {
            vec![Stereo3DDiagnostic::DeclaredTwoDButNonzeroZ {
                observed: geometry_rank,
            }]
        }
        (CoordinateDimension::ThreeD, GeometryRank::FlatZero | GeometryRank::Coplanar) => {
            vec![Stereo3DDiagnostic::DeclaredThreeDButFlat {
                observed: geometry_rank,
            }]
        }
        _ => Vec::new(),
    };
    let square_planar = conformer.as_ref().map_or_else(Vec::new, |coords| {
        let perceived = perceive_square_planar_from_3d(mol, coords);
        let conflicts = wedge_vs_3d_conflicts(mol, coords);
        if geometry_rank == GeometryRank::ThreeD
            && coordinate_dimension != CoordinateDimension::TwoD
        {
            correct_conflicting_wedges_from_3d(mol, &conflicts);
            perceive_tetrahedral_from_3d(mol, coords);
            super::correct_ring_alkene_directions_from_3d(mol, coords);
        }
        diagnostics.extend(conflicts);
        perceived
    });
    (conformer, geometry_rank, diagnostics, square_planar)
}

/// Open Babel's 0D parity numbers non-hydrogen neighbors by atom-block index,
/// with explicit or implicit hydrogen after them. Interpret it only for an
/// all-zero coordinate block:
/// Open Babel labels its 0D output as `2D` in the header, so the header alone
/// cannot decide this. Existing wedge perception keeps precedence.
fn apply_0d_atom_parities(mol: &mut Molecule, parities: &[(AtomIdx, u8)]) {
    for &(center, parity) in parities {
        if mol.atom(center).chirality != Chirality::None
            || has_duplicate_terminal_substituents(mol, center)
        {
            continue;
        }
        let mut order: Vec<u32> = mol.neighbors(center).map(|(atom, _)| atom.0).collect();
        order.sort_by_key(|&index| (mol.atom(AtomIdx(index)).element == Element::H, index));
        // Three-neighbor sulfur can use its lone pair as the fourth
        // tetrahedral direction; its stereo order keeps three entries.
        if order.len() == 3 {
            if chematic_core::implicit_hcount(mol, center) == 1 {
                order.push(STEREO_H_SENTINEL);
            } else if mol.atom(center).element != Element::S {
                continue;
            }
        } else if order.len() != 4 {
            continue;
        }
        mol.set_stereo_neighbor_order(center, order);
        mol.set_chirality(
            center,
            if parity == 1 {
                Chirality::Clockwise
            } else {
                Chirality::CounterClockwise
            },
        );
    }
}

/// Two chemically identical terminal groups make tetrahedral parity
/// meaningless. This local rejection is safe even in bridged systems, where
/// whole-molecule symmetry classifiers can reject valid stereocenters.
fn has_duplicate_terminal_substituents(mol: &Molecule, center: AtomIdx) -> bool {
    let terminal: Vec<(AtomIdx, BondIdx)> = mol
        .neighbors(center)
        .filter(|&(atom, _)| mol.neighbors(atom).count() == 1)
        .collect();
    for (i, &(left, left_bond)) in terminal.iter().enumerate() {
        for &(right, right_bond) in terminal.iter().skip(i + 1) {
            let a = mol.atom(left);
            let b = mol.atom(right);
            let ordinary_order = |bond: BondIdx| match mol.bond(bond).order {
                BondOrder::Up | BondOrder::Down => BondOrder::Single,
                order => order,
            };
            if a.element == b.element
                && a.isotope == b.isotope
                && a.charge == b.charge
                && a.aromatic == b.aromatic
                && a.wildcard == b.wildcard
                && chematic_core::implicit_hcount(mol, left)
                    == chematic_core::implicit_hcount(mol, right)
                && ordinary_order(left_bond) == ordinary_order(right_bond)
            {
                return true;
            }
        }
    }
    false
}

/// Open Babel's 0D MOL extension marks substituent orientation at each end of
/// a double bond. SMILES directions follow bond atom order. Two markers beside
/// the same double bond therefore need opposite orientation conversions. A
/// shared marker in a conjugated system links those conversions together.
fn apply_0d_alkene_directions(mol: &mut Molecule) {
    if !mol
        .bonds()
        .any(|(_, bond)| matches!(bond.order, BondOrder::Up | BondOrder::Down))
    {
        return;
    }
    let aromaticity = chematic_perception::assign_aromaticity(mol);
    let mut links: Vec<Vec<(usize, bool)>> = vec![Vec::new(); mol.bond_count()];
    for (double_idx, double) in mol.bonds() {
        if double.order != BondOrder::Double || aromaticity.is_bond_aromatic(double_idx) {
            continue;
        }
        let markers = |atom| {
            mol.neighbors(atom)
                .filter_map(|(_, idx)| {
                    matches!(mol.bond(idx).order, BondOrder::Up | BondOrder::Down).then_some(idx)
                })
                .collect::<Vec<_>>()
        };
        let left_markers = markers(double.atom1);
        let right_markers = markers(double.atom2);
        for &left in &left_markers {
            for &right in &right_markers {
                let left_reversed = mol.bond(left).atom1 == double.atom1;
                let right_reversed = mol.bond(right).atom1 == double.atom2;
                let opposite = left_reversed ^ right_reversed;
                links[left.0 as usize].push((right.0 as usize, opposite));
                links[right.0 as usize].push((left.0 as usize, opposite));
            }
        }
    }
    let mut graph_flip = vec![None; mol.bond_count()];
    let mut visited = vec![false; mol.bond_count()];
    for root in 0..links.len() {
        if links[root].is_empty() || visited[root] {
            continue;
        }
        let mut stack = vec![root];
        let mut component = Vec::new();
        let mut conflict = false;
        visited[root] = true;
        graph_flip[root] = Some(false);
        while let Some(from) = stack.pop() {
            component.push(from);
            let from_flip = graph_flip[from].expect("linked marker has an orientation");
            for &(to, opposite) in &links[from] {
                let expected = from_flip ^ opposite;
                if !visited[to] {
                    visited[to] = true;
                    graph_flip[to] = Some(expected);
                    stack.push(to);
                } else if graph_flip[to] != Some(expected) {
                    conflict = true;
                }
            }
        }
        if conflict {
            for idx in component {
                graph_flip[idx] = None;
            }
        }
    }
    let flip: Vec<BondIdx> = mol
        .bonds()
        .filter_map(|(idx, bond)| {
            if !matches!(bond.order, BondOrder::Up | BondOrder::Down) {
                return None;
            }
            let has_alkene = |atom| {
                mol.neighbors(atom).any(|(other, adjacent)| {
                    mol.bond(adjacent).order == BondOrder::Double
                        && !aromaticity.is_bond_aromatic(adjacent)
                        && mol
                            .neighbors(other)
                            .any(|(substituent, _)| substituent != atom)
                })
            };
            graph_flip[idx.0 as usize]
                .unwrap_or_else(|| has_alkene(bond.atom1) && !has_alkene(bond.atom2))
                .then_some(idx)
        })
        .collect();
    for idx in flip {
        let reversed = match mol.bond(idx).order {
            BondOrder::Up => BondOrder::Down,
            BondOrder::Down => BondOrder::Up,
            _ => unreachable!("collected directional bonds only"),
        };
        mol.set_bond_order(idx, reversed);
    }
}

fn diagnostic_report(
    mut mol: Molecule,
    metadata: MolMetadata,
    coords: Vec<(f64, f64)>,
    raw_z: Vec<f64>,
    coordinate_dimension: CoordinateDimension,
    unspecified_ez: &std::collections::HashSet<BondIdx>,
    parities: &[(AtomIdx, u8)],
) -> MolReadReport {
    let stereo_diagnostics =
        apply_ctab_local_parity_from_wedges_with_diagnostics(&mut mol, &coords);
    let ez_diagnostics = apply_ez_directions_from_2d_ex(&mut mol, &coords, unspecified_ez);
    let (conformer, geometry_rank, stereo3d_diagnostics, square_planar_diagnostics) =
        geometry_diagnostics(&mut mol, &coords, &raw_z, coordinate_dimension);
    if geometry_rank == GeometryRank::FlatZero
        && coords.iter().all(|&(x, y)| x == 0.0 && y == 0.0)
        && raw_z.iter().all(|&z| z == 0.0)
    {
        apply_0d_atom_parities(&mut mol, parities);
        apply_0d_alkene_directions(&mut mol);
    }
    MolReadReport {
        mol,
        metadata,
        coords,
        stereo_diagnostics,
        ez_diagnostics,
        conformer,
        coordinate_dimension,
        geometry_rank,
        stereo3d_diagnostics,
        square_planar_diagnostics,
    }
}

fn read_mol_internal(
    input: &str,
    include_diagnostics: bool,
) -> Result<MolReadReport, MolParseError> {
    let mut lines = input.lines().enumerate();
    let (metadata, coordinate_dimension) = parse_header(&mut lines)?;
    let (atom_count, bond_count) = parse_counts(&mut lines)?;
    let mut atoms = parse_atom_block(&mut lines, atom_count, bond_count, include_diagnostics)?;
    let unspecified_ez = parse_bond_block(&mut lines, &mut atoms.builder, bond_count)?;
    let properties = parse_property_block(&mut lines, atom_count)?;

    let mut mol = atoms.builder.build();
    apply_atom_updates(&mut mol, atoms.valences, atoms.default_valences, properties);
    if !include_diagnostics {
        return Ok(graph_only_report(mol, metadata));
    }
    Ok(diagnostic_report(
        mol,
        metadata,
        atoms.coords,
        atoms.raw_z,
        coordinate_dimension,
        &unspecified_ez,
        &atoms.parities,
    ))
}

/// Parse a MOL V2000 string into a `(Molecule, MolMetadata, coords)` triple.
///
/// Thin wrapper around [`read_mol_with_diagnostics`] that discards
/// `stereo_diagnostics` -- signature and behavior unchanged from before
/// stereo perception was wired in, except that a wedge/hash bond now
/// populates `Atom.chirality`/`Molecule::stereo_neighbor_order` where it
/// previously left them unset.
#[allow(clippy::type_complexity)]
pub fn parse_mol_with_coords(
    input: &str,
) -> Result<(Molecule, MolMetadata, Vec<(f64, f64)>), MolParseError> {
    read_mol_with_diagnostics(input).map(|r| (r.mol, r.metadata, r.coords))
}

/// Parse a MOL V2000 string into a `(Molecule, MolMetadata)` pair.
///
/// This is a convenience wrapper around [`parse_mol_with_coords`] that discards
/// the 2D coordinate data.
pub fn parse_mol(input: &str) -> Result<(Molecule, MolMetadata), MolParseError> {
    parse_mol_with_coords(input).map(|(mol, meta, _coords)| (mol, meta))
}

/// Parse all molecules from an SDF string, running stereo perception and
/// returning every rejected wedge/hash center per record.
///
/// Stops and returns an error on the first parse failure.
pub fn read_sdf_with_diagnostics(input: &str) -> Result<Vec<MolReadReport>, MolParseError> {
    let mut result = Vec::new();
    let mut remaining = input;
    loop {
        // A blank first line is a legal MOL name line (issue #171) -- do not
        // skip it. A genuinely empty gap between/after `$$$$` delimiters is
        // already handled below via `block.trim().is_empty()`.
        if remaining.is_empty() {
            break;
        }

        // Find the $$$$ delimiter (line-by-line to avoid false matches inside data).
        let mut byte_offset = 0usize;
        let (end_byte, after_delim) = loop {
            let rest = &remaining[byte_offset..];
            match rest.find('\n') {
                Some(nl) => {
                    let line = rest[..nl].trim_end_matches('\r');
                    if line == "$$$$" {
                        break (byte_offset, &remaining[byte_offset + nl + 1..]);
                    }
                    byte_offset += nl + 1;
                }
                None => {
                    if rest.trim_end_matches('\r') == "$$$$" {
                        break (byte_offset, "");
                    }
                    break (remaining.len(), "");
                }
            }
        };

        let block = &remaining[..end_byte];
        remaining = after_delim;
        if block.trim().is_empty() {
            continue;
        }

        result.push(read_mol_with_diagnostics(block)?);
    }
    Ok(result)
}

/// Parse all molecules from an SDF string, returning 2D coordinates.
///
/// Each entry contains the molecule, its metadata, and a `Vec<(x, y)>` of
/// 2D coordinates in atom-insertion order (the same order as `.atoms()`).
///
/// Thin wrapper around [`read_sdf_with_diagnostics`] that discards
/// `stereo_diagnostics`. Stops and returns an error on the first parse failure.
#[allow(clippy::type_complexity)]
pub fn parse_sdf_with_coords(
    input: &str,
) -> Result<Vec<(Molecule, MolMetadata, Vec<(f64, f64)>)>, MolParseError> {
    Ok(read_sdf_with_diagnostics(input)?
        .into_iter()
        .map(|r| (r.mol, r.metadata, r.coords))
        .collect())
}
