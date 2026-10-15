use super::*;

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Write a `Molecule` to MOL V2000 format.
///
/// Coordinates are written as 0.0 because the core `Molecule` type does not
/// store 2D/3D coordinates.  All other atom and bond fields are derived from
/// the molecule graph.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo.** This is a 2D-only
/// writer with no z channel; MOL/CTfile has no other field for a
/// non-tetrahedral stereo tag either (see
/// `docs/rfcs/square_planar_mol_io_rfc.md`). A square-planar-tagged atom is
/// written with no indication anything was dropped. If `mol` may carry
/// `Chirality::SquarePlanar`, use [`write_mol_with_conformer_checked`]
/// instead, which fails closed with a typed error rather than silently
/// discarding the tag.
pub fn write_mol(mol: &Molecule, metadata: &MolMetadata) -> String {
    write_mol_with_coords(mol, metadata, &[])
}

/// Serialize `mol` to a V2000 MOL block, using `coords` for atom positions.
///
/// `coords[i]` is the `(x, y)` position in Ångström for atom index `i`.
/// Atoms beyond `coords.len()` receive `(0.0, 0.0, 0.0)`.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo** -- see
/// [`write_mol`]'s doc comment; the same 2D-only limitation applies here.
/// Use [`write_mol_with_conformer_checked`] if `mol` may carry it.
#[allow(clippy::write_with_newline)]
pub fn write_mol_with_coords(
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
) -> String {
    let mut out = String::new();
    write_mol_with_coords_into(&mut out, mol, metadata, coords);
    out
}

/// Standard CTfile atom symbol for real atoms and the bounded pseudoatom
/// contract supported by this crate.
fn mdl_atom_symbol(mol: &Molecule, idx: AtomIdx, atom: &Atom) -> &'static str {
    if !atom.wildcard {
        return atom.element.symbol();
    }
    match mol.r_group_label(idx).and_then(|label| label.number()) {
        Some(_) => "R#",
        None if mol.r_group_label(idx).is_some() => "R",
        None => "*",
    }
}

fn write_v2000_rgroups(out: &mut String, mol: &Molecule) {
    use std::fmt::Write as _;

    let entries = mol
        .atoms()
        .filter_map(|(idx, _)| {
            mol.r_group_label(idx)
                .and_then(|label| label.number())
                .map(|number| (idx.0 + 1, number))
        })
        .collect::<Vec<_>>();
    // CTfile property records conventionally carry at most eight pairs.
    for chunk in entries.chunks(8) {
        write!(out, "M  RGP{:>3}", chunk.len()).expect("writing to String cannot fail");
        for &(atom_id, number) in chunk {
            write!(out, "{:>4}{:>4}", atom_id, number).expect("writing to String cannot fail");
        }
        out.push('\n');
    }
}

#[inline]
/// Append `x` as `format!("{x:>width$.4}")` would, without the general
/// float formatter (its exact fallback dominated writing laid-out
/// molecules). Values near a rounding tie, very large values and
/// non-finite values take the formatter, so the bytes are always the same.
pub(crate) fn push_fixed4(out: &mut String, x: f64, width: usize) {
    use std::fmt::Write as _;
    let scaled = x.abs() * 10_000.0;
    let frac = scaled - scaled.floor();
    if !x.is_finite() || scaled >= 1e15 || (frac - 0.5).abs() < 1e-6 {
        write!(out, "{x:>width$.4}").expect("writing to String cannot fail");
        return;
    }
    let units = scaled.round() as u64;
    let mut buf = [0u8; 24];
    let mut pos = buf.len();
    let mut push = |b: u8| {
        pos -= 1;
        buf[pos] = b;
    };
    let (mut int, mut dec) = (units / 10_000, units % 10_000);
    for _ in 0..4 {
        push(b'0' + (dec % 10) as u8);
        dec /= 10;
    }
    push(b'.');
    loop {
        push(b'0' + (int % 10) as u8);
        int /= 10;
        if int == 0 {
            break;
        }
    }
    if x.is_sign_negative() {
        push(b'-');
    }
    let digits = &buf[pos..];
    for _ in digits.len()..width {
        out.push(' ');
    }
    out.push_str(std::str::from_utf8(digits).expect("ASCII digits"));
}

pub(super) fn push_right_aligned_u32(out: &mut String, mut value: u32, width: usize) {
    // V2000 counts and bond fields are almost always three columns wide.
    // Append the entire field in one operation, retaining the general path
    // for widths/values outside the fixed-width format.
    if width == 3 && value < 1000 {
        static FIELDS: [[u8; 3]; 1000] = {
            let mut fields = [[b' '; 3]; 1000];
            let mut n = 0;
            while n < fields.len() {
                if n >= 100 {
                    fields[n][0] = b'0' + (n / 100) as u8;
                }
                if n >= 10 {
                    fields[n][1] = b'0' + (n / 10 % 10) as u8;
                }
                fields[n][2] = b'0' + (n % 10) as u8;
                n += 1;
            }
            fields
        };
        // Validate ASCII at compile time, not at every emitted field.
        static TEXT: [&str; 1000] = {
            let mut text = [""; 1000];
            let mut n = 0;
            while n < text.len() {
                text[n] = match std::str::from_utf8(&FIELDS[n]) {
                    Ok(s) => s,
                    Err(_) => panic!("ASCII decimal field"),
                };
                n += 1;
            }
            text
        };
        out.push_str(TEXT[value as usize]);
        return;
    }
    let mut digits = [0_u8; 10];
    let mut start = digits.len();
    loop {
        start -= 1;
        digits[start] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }
    let len = digits.len() - start;
    for _ in len..width {
        out.push(' ');
    }
    for &digit in &digits[start..] {
        out.push(digit as char);
    }
}

#[inline]
pub(super) fn push_right_aligned_i16(out: &mut String, value: i16, width: usize) {
    let negative = value < 0;
    let magnitude = value.unsigned_abs() as u32;
    let digits = if magnitude < 10 {
        1
    } else if magnitude < 100 {
        2
    } else if magnitude < 1_000 {
        3
    } else if magnitude < 10_000 {
        4
    } else {
        5
    };
    for _ in (digits + usize::from(negative))..width {
        out.push(' ');
    }
    if negative {
        out.push('-');
    }
    push_right_aligned_u32(out, magnitude, digits);
}

/// Serialize a MOL V2000 block into a reusable caller-owned buffer.
///
/// This is the allocation-amortized implementation used by streaming
/// writers. The buffer is cleared first; on return it contains exactly the
/// same bytes as [`write_mol_with_coords`].
#[doc(hidden)]
pub fn write_mol_with_coords_into(
    out: &mut String,
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
) {
    write_v2000_reporting(out, mol, metadata, coords, &[]);
}

/// [`write_mol_with_coords`] that also reports the stereo the block does not
/// carry. An empty `coords` lays the molecule out (as [`write_mol`] does).
///
/// ```
/// let mol = chematic_smiles::parse("C[C@H](N)O").unwrap();
/// let (block, loss) =
///     chematic_mol::write_mol_with_stereo_report(&mol, &Default::default(), &[]);
/// assert!(loss.is_empty());
/// assert!(block.contains("V2000"));
/// ```
pub fn write_mol_with_stereo_report(
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
) -> (String, MolStereoLoss) {
    let mut out = String::new();
    let loss = write_v2000_reporting(&mut out, mol, metadata, coords, &[]);
    (out, loss)
}

/// Stereo of a molecule that a written 2D MOL block does not carry. A reader
/// sees each listed centre as unspecified and each listed double bond as
/// "either" (V2000 stereo 3 / V3000 `CFG=2`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MolStereoLoss {
    /// Tetrahedral centres written without a wedge or hash: no wedge on the
    /// coordinates read back as the declared configuration (bridged and cage
    /// layouts, bonds drawn in nearly one direction).
    pub centres: Vec<AtomIdx>,
    /// Declared E/Z double bonds whose geometry the coordinates could not be
    /// given (ring bonds); written "either".
    pub double_bonds: Vec<chematic_core::BondIdx>,
    /// Square-planar centres: a 2D block has no field for them.
    pub non_tetrahedral_centres: Vec<AtomIdx>,
    /// Enhanced stereo groups other than absolute: V2000 has no field for
    /// them (V3000 writes them).
    pub stereo_groups_dropped: bool,
}

impl MolStereoLoss {
    /// Whether the block carries all of the molecule's stereo.
    pub fn is_empty(&self) -> bool {
        self.centres.is_empty()
            && self.double_bonds.is_empty()
            && self.non_tetrahedral_centres.is_empty()
            && !self.stereo_groups_dropped
    }

    pub(crate) fn from_depiction(
        mol: &Molecule,
        depiction: Option<&crate::stereo_depiction::StereoDepiction>,
        format: MolFormat,
    ) -> Self {
        let mut centres: Vec<AtomIdx> = depiction
            .map(|d| d.unexpressed_centres.clone())
            .unwrap_or_default();
        // A tetrahedral tag without a neighbour order has no configuration
        // to draw.
        centres.extend(mol.atoms().filter_map(|(idx, atom)| {
            (atom.chirality.is_tetrahedral() && mol.stereo_neighbor_order(idx).is_none())
                .then_some(idx)
        }));
        centres.sort_by_key(|a| a.0);
        centres.dedup();
        MolStereoLoss {
            centres,
            double_bonds: depiction
                .map(|d| d.unexpressed_double_bonds.clone())
                .unwrap_or_default(),
            non_tetrahedral_centres: mol
                .atoms()
                .filter(|(_, a)| {
                    a.chirality != chematic_core::Chirality::None && !a.chirality.is_tetrahedral()
                })
                .map(|(idx, _)| idx)
                .collect(),
            stereo_groups_dropped: format == MolFormat::V2000
                && mol
                    .stereo_groups()
                    .iter()
                    .any(|g| g.kind != chematic_core::StereoGroupKind::Absolute),
        }
    }
}

impl core::fmt::Display for MolStereoLoss {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let atoms = |v: &[AtomIdx]| v.iter().map(|a| a.0.to_string()).collect::<Vec<_>>();
        let mut parts = Vec::new();
        if !self.centres.is_empty() {
            parts.push(format!(
                "tetrahedral centres without a wedge (atoms {})",
                atoms(&self.centres).join(", ")
            ));
        }
        if !self.double_bonds.is_empty() {
            let bonds: Vec<String> = self.double_bonds.iter().map(|b| b.0.to_string()).collect();
            parts.push(format!(
                "E/Z written as either (bonds {})",
                bonds.join(", ")
            ));
        }
        if !self.non_tetrahedral_centres.is_empty() {
            parts.push(format!(
                "square-planar centres (atoms {})",
                atoms(&self.non_tetrahedral_centres).join(", ")
            ));
        }
        if self.stereo_groups_dropped {
            parts.push("enhanced stereo groups".to_string());
        }
        if parts.is_empty() {
            write!(f, "no stereo lost")
        } else {
            write!(f, "MOL block loses stereo: {}", parts.join("; "))
        }
    }
}

/// [`write_mol_with_stereo_report`] on chematic's own layout (the stereo
/// depiction when the molecule has stereo, the plain layout otherwise), as
/// the Python and WASM `to_mol_block` write it.
///
/// A layout draws every double bond cis or trans, and readers take drawn
/// geometry as E/Z. Double bonds that could carry E/Z but have none
/// declared are therefore written "either" (stereo 3), so the block never
/// gains a configuration the molecule does not have. Caller-supplied
/// coordinates ([`write_mol_with_coords`]) are data and stay unmarked.
///
/// ```
/// let mol = chematic_smiles::parse("CC=CC").unwrap();
/// let (block, loss) = chematic_mol::write_laid_out_mol(&mol, &Default::default());
/// assert!(loss.is_empty());
/// assert!(block.contains("  2  3  2  3\n"));
/// ```
pub fn write_laid_out_mol(mol: &Molecule, metadata: &MolMetadata) -> (String, MolStereoLoss) {
    let mut out = String::new();
    let loss = write_laid_out_mol_into(&mut out, mol, metadata);
    (out, loss)
}

/// [`write_laid_out_mol`] into a reusable buffer (cleared first).
#[doc(hidden)]
pub fn write_laid_out_mol_into(
    out: &mut String,
    mol: &Molecule,
    metadata: &MolMetadata,
) -> MolStereoLoss {
    let coords = crate::stereo_depiction::mol_block_coords(mol);
    let either = if crate::stereo_depiction::needs_stereo_depiction(mol) {
        // The stereo depiction marks its own unspecified double bonds.
        Vec::new()
    } else {
        crate::stereo_depiction::undeclared_stereo_double_bonds(mol)
    };
    write_v2000_reporting(out, mol, metadata, &coords, &either)
}

/// An SDF record written as [`write_laid_out_mol`] writes its MOL block.
pub fn write_laid_out_sdf_record_into(
    out: &mut String,
    mol: &Molecule,
    meta: &MolMetadata,
    props: &std::collections::HashMap<String, String>,
) {
    write_laid_out_mol_into(out, mol, meta);
    append_sd_fields_and_delimiter(out, props);
}

/// The V2000 valence field (vvv) for atom `idx`: 0 (unset) when a reader
/// infers the atom's hydrogens from its bonds, otherwise the total valence
/// (15 for zero), so `[I]`, `[S](=O)=O` or `[SH2]` keep their hydrogen
/// count (they used to read back as `[IH]`, `[SH]`, ...).
fn valence_field(mol: &Molecule, idx: AtomIdx, atom: &Atom) -> u8 {
    let Some(h) = atom.hydrogen_count else {
        return 0;
    };
    if atom.wildcard
        || chematic_core::valence::valence_inferred_hcount(mol, idx) == h
        || mol
            .neighbors(idx)
            .any(|(_, b)| mol.bond(b).order == BondOrder::Aromatic)
    {
        return 0;
    }
    match chematic_core::valence::bond_order_sum(mol, idx).saturating_add(h) {
        0 => 15,
        v if v < 15 => v,
        _ => 0,
    }
}

fn write_atom_row(
    out: &mut String,
    mol: &Molecule,
    idx: AtomIdx,
    atom: &Atom,
    coords: Option<(f64, f64)>,
) {
    let symbol = mdl_atom_symbol(mol, idx, atom);
    let charge = encode_charge(atom.charge);
    let mass = encode_mass_difference(atom.element, atom.isotope).unwrap_or(0);
    let atom_map = atom.atom_map.unwrap_or(0);
    let valence = valence_field(mol, idx, atom);
    if let Some((x, y)) = coords {
        push_fixed4(out, x, 10);
        push_fixed4(out, y, 10);
        out.push_str("    0.0000 ");
    } else {
        out.push_str("    0.0000    0.0000    0.0000 ");
    }
    out.push_str(symbol);
    for _ in symbol.len()..3 {
        out.push(' ');
    }
    if mass == 0 && charge == 0 {
        out.push_str(" 0  0");
    } else {
        push_right_aligned_i16(out, mass, 2);
        push_right_aligned_u32(out, charge as u32, 3);
    }
    out.push_str("  0  0  0");
    push_right_aligned_u32(out, u32::from(valence), 3);
    out.push_str("  0  0  0");
    push_right_aligned_u32(out, atom_map as u32, 3);
    out.push_str("  0\n");
}

fn write_atom_rows(out: &mut String, mol: &Molecule, coords: &[(f64, f64)]) {
    for (idx, atom) in mol.atoms() {
        write_atom_row(out, mol, idx, atom, coords.get(idx.0 as usize).copied());
    }
}

fn v2000_bond_type(order: BondOrder) -> u32 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic | BondOrder::Quadruple => 4,
        BondOrder::QuerySingleOrDouble => 5,
        BondOrder::QuerySingleOrAromatic => 6,
        BondOrder::QueryDoubleOrAromatic => 7,
        BondOrder::QueryAny | BondOrder::Zero => 8,
    }
}

fn v2000_bond_stereo(
    bond_idx: BondIdx,
    bond: &chematic_core::BondEntry,
    depiction: Option<&crate::stereo_depiction::StereoDepiction>,
    either: &[BondIdx],
) -> (u32, u32, u32) {
    let atom1 = bond.atom1.0 + 1;
    let atom2 = bond.atom2.0 + 1;
    match depiction {
        Some(depiction) => match depiction.wedges.get(&bond_idx) {
            Some(wedge) => {
                let (start, end) = if wedge.start == bond.atom1 {
                    (atom1, atom2)
                } else {
                    (atom2, atom1)
                };
                (start, end, if wedge.order == BondOrder::Up { 1 } else { 6 })
            }
            None if depiction.unexpressed_double_bonds.contains(&bond_idx)
                || depiction.unspecified_double_bonds.contains(&bond_idx) =>
            {
                (atom1, atom2, 3)
            }
            None => (atom1, atom2, 0),
        },
        None => {
            let stereo = match bond.order {
                BondOrder::Up => 1,
                BondOrder::Down => 6,
                _ if either.contains(&bond_idx) => 3,
                _ => 0,
            };
            (atom1, atom2, stereo)
        }
    }
}

fn write_bond_rows(
    out: &mut String,
    mol: &Molecule,
    depiction: Option<&crate::stereo_depiction::StereoDepiction>,
    either: &[BondIdx],
) {
    let kekule = kekule_orders_for_writing(mol);
    for (bond_idx, bond) in mol.bonds() {
        let order = kekule
            .as_ref()
            .and_then(|orders| orders.get(&bond_idx).copied())
            .unwrap_or(bond.order);
        let (atom1, atom2, stereo) = v2000_bond_stereo(bond_idx, bond, depiction, either);
        push_right_aligned_u32(out, atom1, 3);
        push_right_aligned_u32(out, atom2, 3);
        push_right_aligned_u32(out, v2000_bond_type(order), 3);
        push_right_aligned_u32(out, stereo, 3);
        out.push('\n');
    }
}

/// `either`: double bonds written with stereo 3 when there is no stereo
/// depiction (see [`write_laid_out_mol`]).
fn write_v2000_reporting(
    out: &mut String,
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
    either: &[chematic_core::BondIdx],
) -> MolStereoLoss {
    // A MOL block is dominated by fixed-width atom and bond rows.  Reserve
    // the common-size output up front so serialization does not repeatedly
    // grow and copy the String for every row (the old empty String was a
    // noticeable cost when writing many small records).
    out.clear();
    out.reserve(128 + mol.atom_count() * 80 + mol.bond_count() * 16);

    // Header lines 1–3
    out.push_str(&metadata.name);
    out.push('\n');
    out.push_str("  chematic\n");
    out.push_str(&metadata.comment);
    out.push('\n');

    // A molecule with stereo is written on coordinates that express it, with
    // one wedge or hash per tetrahedral centre (see `stereo_depiction`).
    let depiction = crate::stereo_depiction::needs_stereo_depiction(mol)
        .then(|| crate::stereo_depiction::stereo_depiction(mol, coords));
    let coords: &[(f64, f64)] = depiction.as_ref().map_or(coords, |d| &d.coords);

    // Counts line (line 4)
    let natoms = mol.atom_count();
    let nbonds = mol.bond_count();
    push_right_aligned_u32(out, natoms as u32, 3);
    push_right_aligned_u32(out, nbonds as u32, 3);
    out.push_str("  0  0  0  0  0  0  0  0999 V2000\n");

    write_atom_rows(out, mol, coords);
    write_bond_rows(out, mol, depiction.as_ref(), either);

    write_v2000_rgroups(out, mol);

    // Terminator
    out.push_str("M  END\n");
    MolStereoLoss::from_depiction(mol, depiction.as_ref(), MolFormat::V2000)
}

/// Serialize `mol` to MOL V2000 format using `conformer`'s real 3D
/// coordinates, stamping the header's line-2 dimensional code as `3D`.
///
/// This is the writer counterpart of [`read_mol_with_diagnostics`]'s new 3D
/// support -- round-tripping this output back through the reader reproduces
/// [`CoordinateDimension::ThreeD`] and never manufactures a fresh
/// [`Stereo3DDiagnostic::WedgeVs3DParityConflict`] on its own output: no
/// wedge/hash stereo field is ever emitted here (always `0`), matching the
/// fact that a real 3D geometry makes a 2D wedge symbol redundant at best,
/// contradictory at worst -- see [`write_mol_with_coords`] (unchanged, still
/// the 2D writer) for the wedge-preserving counterpart. Atoms beyond
/// `conformer.atom_count()` receive `(0.0, 0.0, 0.0)`.
///
/// **Does not validate `Chirality::SquarePlanar` stereo against
/// `conformer`** -- it writes whatever coordinates it is given, trusting the
/// caller (the same "trust the caller" posture this crate's writers have
/// always had for tetrahedral wedge bonds). If `conformer`'s geometry
/// doesn't actually match a declared square-planar tag, this will silently
/// write a self-inconsistent file. Use
/// [`write_mol_with_conformer_checked`] instead to fail closed on that
/// mismatch (or on a missing/flat conformer) rather than trusting it.
pub fn write_mol_with_conformer(
    mol: &Molecule,
    metadata: &MolMetadata,
    conformer: &Coords3D,
) -> String {
    let mut out = String::new();

    out.push_str(&metadata.name);
    out.push('\n');
    // Columns 0..10 "  chematic" match the 2D writer exactly; columns
    // 10..20 are the (blank) date field; columns 20..22 carry the "3D"
    // dimensional code -- the same column `parse_dimension_code` reads,
    // empirically confirmed against RDKit 2026.03.3's own `MolToMolBlock`.
    out.push_str("  chematic          3D\n");
    out.push_str(&metadata.comment);
    out.push('\n');

    let natoms = mol.atom_count();
    let nbonds = mol.bond_count();
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        natoms, nbonds
    ));

    for (idx, atom) in mol.atoms() {
        let sym = mdl_atom_symbol(mol, idx, atom);
        let charge_code = encode_charge(atom.charge);
        let mass_difference = encode_mass_difference(atom.element, atom.isotope).unwrap_or(0);
        let atom_map = atom.atom_map.unwrap_or(0);
        let p = conformer
            .points
            .get(idx.0 as usize)
            .copied()
            .unwrap_or(Point3::zero());
        push_fixed4(&mut out, p.x, 10);
        push_fixed4(&mut out, p.y, 10);
        push_fixed4(&mut out, p.z, 10);
        out.push_str(&format!(
            " {:<3}{:>2}{:>3}  0  0  0  0  0  0  0{:>3}  0\n",
            sym, mass_difference, charge_code, atom_map,
        ));
    }

    let kekule = kekule_orders_for_writing(mol);
    for (bond_idx, bond) in mol.bonds() {
        let a1 = bond.atom1.0 + 1;
        let a2 = bond.atom2.0 + 1;
        let order = kekule
            .as_ref()
            .and_then(|k| k.get(&bond_idx).copied())
            .unwrap_or(bond.order);
        let btype = match order {
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            BondOrder::QuerySingleOrDouble => 5,
            BondOrder::QuerySingleOrAromatic => 6,
            BondOrder::QueryDoubleOrAromatic => 7,
            BondOrder::QueryAny | BondOrder::Zero => 8,
            BondOrder::Quadruple => 4,
        };
        // Stereo field is always 0 -- see doc comment above.
        out.push_str(&format!("{:>3}{:>3}{:>3}{:>3}\n", a1, a2, btype, 0));
    }

    write_v2000_rgroups(&mut out, mol);
    out.push_str("M  END\n");
    out
}

// ---------------------------------------------------------------------------
// Square-planar write-side validation -- see the read-side section above
// and `docs/rfcs/square_planar_mol_io_rfc.md` §8-9 for the design (no
// coordinate fabrication, ever; existing coordinates are validated against
// the declared tag, never trusted silently).
// ---------------------------------------------------------------------------

/// Which MOL/SDF variant a [`MolStereoWriteError`] was raised against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MolFormat {
    V2000,
    V3000,
}

impl core::fmt::Display for MolFormat {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::V2000 => write!(f, "V2000"),
            Self::V3000 => write!(f, "V3000"),
        }
    }
}

/// Why [`validate_square_planar_for_write`] refused to write a molecule
/// carrying `Chirality::SquarePlanar`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnsupportedStereoReason {
    /// No 3D conformer was supplied at all. MDL/CTfile has no symbolic field
    /// for a non-tetrahedral stereo tag (RFC §2) -- geometry is the only
    /// mechanism, and there is none to write here. The writer never
    /// fabricates one (RFC §8).
    NoConformerSupplied,
    /// A conformer was supplied, but this atom has no recorded
    /// `stereo_neighbor_order` to validate/write positions against.
    MissingNeighborOrder,
    /// This atom's `stereo_neighbor_order` contains the implicit-hydrogen
    /// sentinel -- there is no real spatial position for an implicit H to
    /// validate or write against (RFC §11).
    ImplicitHydrogenNeighbor,
    /// A conformer was supplied, but it has no real position for the center
    /// or one of its declared neighbors.
    MissingConformerPosition,
    /// The whole-molecule conformer is flat (`GeometryRank::FlatZero`) or
    /// indeterminate -- indistinguishable from an ordinary 2D depiction with
    /// no 3D data at all (RFC §10); this crate's own reader would never
    /// reperceive a tag from it, so writing one out would silently produce
    /// an unrecoverable file.
    WholeMoleculeConformerFlat,
    /// The conformer's local geometry around this atom does not
    /// unambiguously resolve to *any* square-planar permutation.
    GeometryRejected(SquarePlanarRejectionReason),
    /// The conformer's local geometry resolves to a *different* permutation
    /// than the one declared on the atom.
    GeometryTagMismatch { computed: SquarePlanarPermutation },
}

impl core::fmt::Display for UnsupportedStereoReason {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoConformerSupplied => write!(f, "no 3D conformer supplied"),
            Self::MissingNeighborOrder => write!(f, "atom has no recorded stereo_neighbor_order"),
            Self::ImplicitHydrogenNeighbor => write!(
                f,
                "stereo_neighbor_order contains an implicit-hydrogen slot"
            ),
            Self::MissingConformerPosition => {
                write!(
                    f,
                    "conformer is missing a position for this center or a neighbor"
                )
            }
            Self::WholeMoleculeConformerFlat => {
                write!(f, "whole-molecule conformer is flat/indeterminate")
            }
            Self::GeometryRejected(reason) => write!(f, "geometry rejected: {reason}"),
            Self::GeometryTagMismatch { computed } => write!(
                f,
                "conformer geometry encodes {computed:?}, which does not match the declared tag"
            ),
        }
    }
}

/// A MOL/SDF write was refused because `atom`'s declared
/// `chematic_core::StereoGeometry` could not be represented in `format` --
/// see [`UnsupportedStereoReason`] for why. Never raised for tetrahedral
/// chirality (which round-trips via wedge/hash bonds, already written
/// faithfully by every writer in this file).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolStereoWriteError {
    pub atom: AtomIdx,
    pub geometry: StereoGeometry,
    pub format: MolFormat,
    pub reason: UnsupportedStereoReason,
}

impl core::fmt::Display for MolStereoWriteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "cannot write atom {} ({:?} stereo) to MOL {}: {}",
            self.atom.0, self.geometry, self.format, self.reason
        )
    }
}

impl std::error::Error for MolStereoWriteError {}

/// Validate every `Chirality::SquarePlanar` atom in `mol` against `conformer`
/// before writing `format`. `Ok(())` when `mol` has no square-planar atom at
/// all (the common case -- this is a no-op for every existing caller), or
/// when every square-planar atom's declared tag matches what its
/// `conformer` positions actually encode. `pub`, not `pub(crate)`: callers
/// building their own write pipeline can use this as a pre-flight check
/// even without going through this crate's `*_checked` writer wrappers.
///
/// `conformer: None` means "the 2D-only writer path was used" -- always an
/// error (`NoConformerSupplied`) when `mol` has any square-planar atom, since
/// no 2D-only writer in this crate has a z channel to encode one (RFC §2).
pub fn validate_square_planar_for_write(
    mol: &Molecule,
    conformer: Option<&Coords3D>,
    format: MolFormat,
) -> Result<(), MolStereoWriteError> {
    // Whole-molecule flatness check once, up front: if ANY square-planar
    // atom exists and the conformer is flat/indeterminate, every one of them
    // fails the same way (RFC §10) -- report the first by atom order for a
    // deterministic, non-HashMap-dependent error.
    let whole_mol_flat_or_missing = match conformer {
        None => true,
        Some(conf) => matches!(
            classify_geometry_rank(&conf.points),
            GeometryRank::FlatZero | GeometryRank::Indeterminate
        ),
    };

    for (idx, atom) in mol.atoms() {
        let Chirality::SquarePlanar(declared) = atom.chirality else {
            continue;
        };
        let geometry = StereoGeometry::SquarePlanar;

        let Some(conf) = conformer else {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::NoConformerSupplied,
            });
        };
        if whole_mol_flat_or_missing {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::WholeMoleculeConformerFlat,
            });
        }

        let Some(order) = mol.stereo_neighbor_order(idx) else {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::MissingNeighborOrder,
            });
        };
        if order.len() != 4 {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::MissingNeighborOrder,
            });
        }
        if order.contains(&STEREO_H_SENTINEL) {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::ImplicitHydrogenNeighbor,
            });
        }

        let Some(center) = conf.points.get(idx.0 as usize).copied() else {
            return Err(MolStereoWriteError {
                atom: idx,
                geometry,
                format,
                reason: UnsupportedStereoReason::MissingConformerPosition,
            });
        };
        let mut pts = [Point3::zero(); 4];
        for (slot, &n) in order.iter().enumerate() {
            match conf.points.get(n as usize).copied() {
                Some(p) => pts[slot] = p,
                None => {
                    return Err(MolStereoWriteError {
                        atom: idx,
                        geometry,
                        format,
                        reason: UnsupportedStereoReason::MissingConformerPosition,
                    });
                }
            }
        }

        match classify_square_planar_geometry(center, pts) {
            Ok(computed) if computed == declared => {}
            Ok(computed) => {
                return Err(MolStereoWriteError {
                    atom: idx,
                    geometry,
                    format,
                    reason: UnsupportedStereoReason::GeometryTagMismatch { computed },
                });
            }
            Err(reason) => {
                return Err(MolStereoWriteError {
                    atom: idx,
                    geometry,
                    format,
                    reason: UnsupportedStereoReason::GeometryRejected(reason),
                });
            }
        }
    }

    Ok(())
}

/// [`write_mol_with_conformer`], but fails closed with a typed
/// [`MolStereoWriteError`] instead of silently writing coordinates that
/// don't actually match a molecule's declared square-planar stereo (or
/// don't exist at all). See [`validate_square_planar_for_write`].
///
/// For a molecule with no `Chirality::SquarePlanar` atom, this is exactly
/// [`write_mol_with_conformer`] wrapped in `Ok`.
pub fn write_mol_with_conformer_checked(
    mol: &Molecule,
    metadata: &MolMetadata,
    conformer: &Coords3D,
) -> Result<String, MolStereoWriteError> {
    validate_square_planar_for_write(mol, Some(conformer), MolFormat::V2000)?;
    Ok(write_mol_with_conformer(mol, metadata, conformer))
}

// ---------------------------------------------------------------------------
// SDF writer
// ---------------------------------------------------------------------------

/// Serialise one or more molecules to SDF format.
///
/// `records` — slice of `(molecule, metadata, coords)` tuples.
/// `coords` is optional; pass an empty slice to write zero coordinates.
/// Each molecule block is terminated with `$$$$`.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo** -- see
/// [`write_mol`]'s doc comment; this delegates to [`write_mol_with_coords`]
/// per record. There is no `_checked` SDF-multi-record entry point today;
/// call [`write_sdf_record_with_conformer_checked`] once per record and
/// concatenate if `records` may carry it.
#[allow(clippy::type_complexity)]
pub fn write_sdf(records: &[(&Molecule, &MolMetadata, &[(f64, f64)])]) -> String {
    let mut out = String::new();
    for (mol, meta, coords) in records {
        out.push_str(&write_mol_with_coords(mol, meta, coords));
        out.push_str("$$$$\n");
    }
    out
}

/// Serialise one or more molecules to SDF format, appending per-atom partial
/// charges as an SD property `<PARTIAL_CHARGES>`.
///
/// `records` — slice of `(molecule, metadata, 2D-coords, charges)` tuples.
/// `charges[i]` is the partial charge for atom `i` (heavy atoms only).
/// Pass an empty charges slice to omit the property block.
///
/// Example SD block appended after `M  END`:
/// ```text
/// > <PARTIAL_CHARGES>
/// -0.2359 0.1076 -0.4500 0.1806
///
/// $$$$
/// ```
///
/// **Does not preserve `Chirality::SquarePlanar` stereo** -- see
/// [`write_mol`]'s doc comment; this delegates to [`write_mol_with_coords`]
/// per record.
#[allow(clippy::type_complexity)]
pub fn write_sdf_with_charges(
    records: &[(&Molecule, &MolMetadata, &[(f64, f64)], &[f64])],
) -> String {
    let mut out = String::new();
    for (mol, meta, coords, charges) in records {
        out.push_str(&write_mol_with_coords(mol, meta, coords));
        if !charges.is_empty() {
            out.push_str("> <PARTIAL_CHARGES>\n");
            let vals: Vec<String> = charges.iter().map(|q| format!("{q:.4}")).collect();
            out.push_str(&vals.join(" "));
            out.push_str("\n\n");
        }
        out.push_str("$$$$\n");
    }
    out
}

/// Serialise a single molecule to one SDF record with arbitrary SD data fields.
///
/// Keys starting with `_` are treated as internal/computed properties and are
/// omitted from the SD block (e.g. `_Name` is written into the MOL header, not
/// as an SD field).  The record is terminated with `$$$$`.
///
/// With supplied 2D coordinates, unspecified stereogenic double bonds are
/// marked "either" so the drawing does not invent an E/Z assignment.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo** -- see
/// [`write_mol`]'s doc comment.
/// Use [`write_sdf_record_with_conformer_checked`] if `mol` may carry it.
pub fn write_sdf_record(
    mol: &Molecule,
    meta: &MolMetadata,
    coords: &[(f64, f64)],
    props: &std::collections::HashMap<String, String>,
) -> String {
    let mut out = String::new();
    write_sdf_record_into(&mut out, mol, meta, coords, props);
    out
}

/// Serialize one SDF record into a reusable caller-owned buffer.
///
/// The buffer is cleared first. Reusing it across records amortizes the MOL
/// block allocation while preserving byte-for-byte output compatibility with
/// [`write_sdf_record`].
#[doc(hidden)]
pub fn write_sdf_record_into(
    out: &mut String,
    mol: &Molecule,
    meta: &MolMetadata,
    coords: &[(f64, f64)],
    props: &std::collections::HashMap<String, String>,
) {
    // Preserve the molecule's unspecified double-bond stereo when coordinates
    // are supplied. Otherwise an SDF consumer can infer a definite E/Z solely
    // from the drawing, including azo bonds whose explicit "either" marker
    // was read without any directional markers in the molecule graph.
    let either = if !coords.is_empty() && !crate::stereo_depiction::needs_stereo_depiction(mol) {
        crate::stereo_depiction::undeclared_stereo_double_bonds(mol)
    } else {
        Vec::new()
    };
    write_v2000_reporting(out, mol, meta, coords, &either);
    append_sd_fields_and_delimiter(out, props);
}

/// Append visible SD fields and the record delimiter without allocating a
/// temporary formatted `String` for every property.  Reserving the aggregate
/// payload once also avoids repeated growth when records carry many fields.
fn append_sd_fields_and_delimiter(
    out: &mut String,
    props: &std::collections::HashMap<String, String>,
) {
    use std::fmt::Write as _;

    let additional = props
        .iter()
        .filter(|(key, _)| !key.starts_with('_'))
        .map(|(key, value)| key.len().saturating_add(value.len()).saturating_add(7))
        .sum::<usize>()
        .saturating_add(5);
    out.reserve(additional);
    for (key, value) in props {
        if !key.starts_with('_') {
            write!(out, "> <{key}>\n{value}\n\n").expect("writing to String cannot fail");
        }
    }
    out.push_str("$$$$\n");
}

/// Like [`write_sdf_record`] but emits a MOL V3000 (Extended Ctab) block
/// instead of V2000 — required for molecules with more than 999 atoms or
/// bonds, which don't fit V2000's fixed-width count fields.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo** -- see
/// [`write_mol`]'s doc comment; this delegates to [`crate::mol3000::write_mol_v3000`],
/// itself 2D-only.
pub fn write_sdf_record_v3000(
    mol: &Molecule,
    meta: &MolMetadata,
    coords: &[(f64, f64)],
    props: &std::collections::HashMap<String, String>,
) -> String {
    let mut out = crate::mol3000::write_mol_v3000(mol, meta, coords);
    append_sd_fields_and_delimiter(&mut out, props);
    out
}

/// Like [`write_sdf_record`] but writes `conformer`'s real 3D coordinates
/// (V2000, via [`write_mol_with_conformer`]) instead of a 2D `(x, y)` slice
/// -- the 3D counterpart for a single SDF record. To write several
/// conformers of the same molecule as separate, repeated records (readable
/// back as one [`crate::sdf::ConformerEnsemble`] by
/// [`crate::sdf::read_sdf_conformer_ensembles`]), call this once per
/// conformer and concatenate the results.
///
/// **Does not validate `Chirality::SquarePlanar` stereo against
/// `conformer`** -- see [`write_mol_with_conformer`]'s doc comment for why.
/// Use [`write_sdf_record_with_conformer_checked`] instead to fail closed on
/// a geometry/tag mismatch (or a missing/flat conformer) rather than
/// trusting it.
pub fn write_sdf_record_with_conformer(
    mol: &Molecule,
    meta: &MolMetadata,
    conformer: &Coords3D,
    props: &std::collections::HashMap<String, String>,
) -> String {
    let mut out = write_mol_with_conformer(mol, meta, conformer);
    append_sd_fields_and_delimiter(&mut out, props);
    out
}

/// [`write_sdf_record_with_conformer`], but fails closed with a typed
/// [`MolStereoWriteError`] instead of silently writing coordinates that
/// don't actually match a molecule's declared square-planar stereo (or
/// don't exist at all) -- the SDF-record counterpart of
/// [`write_mol_with_conformer_checked`]. See
/// [`validate_square_planar_for_write`].
///
/// For a molecule with no `Chirality::SquarePlanar` atom, this is exactly
/// [`write_sdf_record_with_conformer`] wrapped in `Ok`. To write a
/// multi-record SDF (e.g. several distinct molecules, or several conformers
/// of the same one), call this once per record and concatenate the
/// results -- there is no separate multi-record entry point, matching
/// [`write_sdf`]'s own per-record-call convention.
pub fn write_sdf_record_with_conformer_checked(
    mol: &Molecule,
    meta: &MolMetadata,
    conformer: &Coords3D,
    props: &std::collections::HashMap<String, String>,
) -> Result<String, MolStereoWriteError> {
    validate_square_planar_for_write(mol, Some(conformer), MolFormat::V2000)?;
    Ok(write_sdf_record_with_conformer(mol, meta, conformer, props))
}
