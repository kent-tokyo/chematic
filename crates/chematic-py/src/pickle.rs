//! Versioned Python-pickle state for [`crate::Mol`].
//!
//! The payload is deliberately independent of Rust struct layout. It stores
//! graph state and semantic sidecars, but not derived caches, so a restored
//! molecule cannot inherit stale computed values.

use crate::{Mol, PropertyStore, PropertyType};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
#[cfg(test)]
use std::collections::HashMap;
use std::sync::Arc;

use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Chirality, CipCode, Element, MoleculeBuilder, RGroupLabel,
    STEREO_H_SENTINEL, SquarePlanarPermutation, StereoGroup, StereoGroupKind,
};

const MAGIC: &[u8; 8] = b"CHMPKL\0\0";
const VERSION: u16 = 2;
const LEGACY_VERSION: u16 = 1;
const MAX_PAYLOAD: usize = 64 * 1024 * 1024;
const MAX_ATOMS: usize = 1_000_000;
const MAX_BONDS: usize = 4_000_000;
const MAX_GROUPS: usize = 1_000_000;
const MAX_PROPERTIES: usize = 1_000_000;
const MAX_STRING: usize = 16 * 1024 * 1024;

struct Writer(Vec<u8>);

impl Writer {
    fn u8(&mut self, value: u8) {
        self.0.push(value);
    }
    fn u16(&mut self, value: u16) {
        self.0.extend_from_slice(&value.to_le_bytes());
    }
    fn u32(&mut self, value: u32) {
        self.0.extend_from_slice(&value.to_le_bytes());
    }
    fn bytes(&mut self, value: &[u8]) {
        self.u32(value.len() as u32);
        self.0.extend_from_slice(value);
    }
    fn option_u8(&mut self, value: Option<u8>) {
        self.u8(value.is_some() as u8);
        if let Some(value) = value {
            self.u8(value);
        }
    }
    fn option_u16(&mut self, value: Option<u16>) {
        self.u8(value.is_some() as u8);
        if let Some(value) = value {
            self.u16(value);
        }
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> Reader<'a> {
    fn new(bytes: &'a [u8]) -> Result<Self, String> {
        if bytes.len() > MAX_PAYLOAD {
            return Err(format!("pickle payload exceeds {MAX_PAYLOAD} bytes"));
        }
        Ok(Self { bytes, offset: 0 })
    }

    fn take(&mut self, len: usize) -> Result<&'a [u8], String> {
        let end = self
            .offset
            .checked_add(len)
            .ok_or_else(|| "pickle payload length overflow".to_string())?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or_else(|| "truncated chematic molecule pickle".to_string())?;
        self.offset = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, String> {
        Ok(self.take(1)?[0])
    }
    fn bool(&mut self, field: &str) -> Result<bool, String> {
        match self.u8()? {
            0 => Ok(false),
            1 => Ok(true),
            value => Err(format!("invalid {field} flag {value}")),
        }
    }
    fn u16(&mut self) -> Result<u16, String> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> Result<u32, String> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn bounded_count(&mut self, field: &str, max: usize) -> Result<usize, String> {
        let value = self.u32()? as usize;
        if value > max {
            return Err(format!("{field} count {value} exceeds limit {max}"));
        }
        Ok(value)
    }
    fn string(&mut self, field: &str) -> Result<String, String> {
        let len = self.bounded_count(field, MAX_STRING)?;
        let bytes = self.take(len)?;
        std::str::from_utf8(bytes)
            .map(str::to_owned)
            .map_err(|_| format!("{field} is not valid UTF-8"))
    }
    fn option_u8(&mut self, field: &str) -> Result<Option<u8>, String> {
        if self.bool(field)? {
            Ok(Some(self.u8()?))
        } else {
            Ok(None)
        }
    }
    fn option_u16(&mut self, field: &str) -> Result<Option<u16>, String> {
        if self.bool(field)? {
            Ok(Some(self.u16()?))
        } else {
            Ok(None)
        }
    }
    fn finish(self) -> Result<(), String> {
        if self.offset == self.bytes.len() {
            Ok(())
        } else {
            Err("trailing bytes in chematic molecule pickle".into())
        }
    }
}

fn chirality_to_parts(value: Chirality, is_allene: bool) -> (u8, u8) {
    if is_allene {
        return match value {
            Chirality::CounterClockwise => (6, 1),
            Chirality::Clockwise => (6, 2),
            _ => (0, 0),
        };
    }
    match value {
        Chirality::None => (0, 0),
        Chirality::CounterClockwise => (1, 0),
        Chirality::Clockwise => (2, 0),
        Chirality::SquarePlanar(SquarePlanarPermutation::SP1) => (3, 1),
        Chirality::SquarePlanar(SquarePlanarPermutation::SP2) => (3, 2),
        Chirality::SquarePlanar(SquarePlanarPermutation::SP3) => (3, 3),
        Chirality::SquarePlanarUnnumbered => (3, 0),
        Chirality::TrigonalBipyramidal(value) => (4, value),
        Chirality::Octahedral(value) => (5, value),
    }
}

fn chirality_from_parts(kind: u8, value: u8) -> Result<Chirality, String> {
    match (kind, value) {
        (0, 0) => Ok(Chirality::None),
        (1, 0) => Ok(Chirality::CounterClockwise),
        (2, 0) => Ok(Chirality::Clockwise),
        (6, 1) => Ok(Chirality::Allene1),
        (6, 2) => Ok(Chirality::Allene2),
        (3, 0) => Ok(Chirality::SquarePlanarUnnumbered),
        (3, 1) => Ok(Chirality::SquarePlanar(SquarePlanarPermutation::SP1)),
        (3, 2) => Ok(Chirality::SquarePlanar(SquarePlanarPermutation::SP2)),
        (3, 3) => Ok(Chirality::SquarePlanar(SquarePlanarPermutation::SP3)),
        (4, 0..=20) => Ok(Chirality::TrigonalBipyramidal(value)),
        (5, 0..=30) => Ok(Chirality::Octahedral(value)),
        _ => Err(format!("invalid chirality kind {kind} value {value}")),
    }
}

fn cip_to_u8(value: Option<CipCode>) -> u8 {
    match value {
        None => 0,
        Some(CipCode::R) => 1,
        Some(CipCode::S) => 2,
        Some(CipCode::E) => 3,
        Some(CipCode::Z) => 4,
        Some(CipCode::LowerR) => 5,
        Some(CipCode::LowerS) => 6,
    }
}

fn property_type_to_u8(value: PropertyType) -> u8 {
    match value {
        PropertyType::String => 0,
        PropertyType::Int => 1,
        PropertyType::Unsigned => 2,
        PropertyType::Double => 3,
        PropertyType::Bool => 4,
    }
}

fn property_type_from_u8(value: u8) -> Result<PropertyType, String> {
    match value {
        0 => Ok(PropertyType::String),
        1 => Ok(PropertyType::Int),
        2 => Ok(PropertyType::Unsigned),
        3 => Ok(PropertyType::Double),
        4 => Ok(PropertyType::Bool),
        _ => Err(format!("invalid property type {value}")),
    }
}

fn cip_from_u8(value: u8) -> Result<Option<CipCode>, String> {
    match value {
        0 => Ok(None),
        1 => Ok(Some(CipCode::R)),
        2 => Ok(Some(CipCode::S)),
        3 => Ok(Some(CipCode::E)),
        4 => Ok(Some(CipCode::Z)),
        5 => Ok(Some(CipCode::LowerR)),
        6 => Ok(Some(CipCode::LowerS)),
        _ => Err(format!("invalid CIP code {value}")),
    }
}

fn bond_to_u8(value: BondOrder) -> u8 {
    match value {
        BondOrder::Single => 0,
        BondOrder::Double => 1,
        BondOrder::Triple => 2,
        BondOrder::Quadruple => 3,
        BondOrder::Aromatic => 4,
        BondOrder::Up => 5,
        BondOrder::Down => 6,
        BondOrder::Zero => 7,
        BondOrder::Dative => 8,
        BondOrder::QueryAny => 9,
        BondOrder::QuerySingleOrDouble => 10,
        BondOrder::QuerySingleOrAromatic => 11,
        BondOrder::QueryDoubleOrAromatic => 12,
    }
}

fn bond_from_u8(value: u8) -> Result<BondOrder, String> {
    match value {
        0 => Ok(BondOrder::Single),
        1 => Ok(BondOrder::Double),
        2 => Ok(BondOrder::Triple),
        3 => Ok(BondOrder::Quadruple),
        4 => Ok(BondOrder::Aromatic),
        5 => Ok(BondOrder::Up),
        6 => Ok(BondOrder::Down),
        7 => Ok(BondOrder::Zero),
        8 => Ok(BondOrder::Dative),
        9 => Ok(BondOrder::QueryAny),
        10 => Ok(BondOrder::QuerySingleOrDouble),
        11 => Ok(BondOrder::QuerySingleOrAromatic),
        12 => Ok(BondOrder::QueryDoubleOrAromatic),
        _ => Err(format!("invalid bond order {value}")),
    }
}

fn checked_atom(value: u32, atom_count: usize, field: &str) -> Result<AtomIdx, String> {
    if value as usize >= atom_count {
        Err(format!("{field} atom index {value} is out of range"))
    } else {
        Ok(AtomIdx(value))
    }
}

fn checked_bond(value: u32, bond_count: usize, field: &str) -> Result<BondIdx, String> {
    if value as usize >= bond_count {
        Err(format!("{field} bond index {value} is out of range"))
    } else {
        Ok(BondIdx(value))
    }
}

fn checked_payload_size(mol: &Mol) -> Result<usize, String> {
    let inner = mol.inner.as_ref();
    if inner.atom_count() > MAX_ATOMS {
        return Err(format!("atom count exceeds pickle limit {MAX_ATOMS}"));
    }
    if inner.bond_count() > MAX_BONDS {
        return Err(format!("bond count exceeds pickle limit {MAX_BONDS}"));
    }
    if inner.stereo_groups().len() > MAX_GROUPS {
        return Err(format!(
            "stereo group count exceeds pickle limit {MAX_GROUPS}"
        ));
    }
    if mol.props.len() > MAX_PROPERTIES {
        return Err(format!(
            "property count exceeds pickle limit {MAX_PROPERTIES}"
        ));
    }

    let mut size = MAGIC.len() + 2 + 4 * 4;
    let mut add = |amount: usize| -> Result<(), String> {
        size = size
            .checked_add(amount)
            .ok_or_else(|| "pickle payload size overflow".to_string())?;
        if size > MAX_PAYLOAD {
            return Err(format!("pickle payload exceeds {MAX_PAYLOAD} bytes"));
        }
        Ok(())
    };
    for (idx, atom) in inner.atoms() {
        let option_bytes = usize::from(atom.isotope.is_some()) * 2
            + usize::from(atom.hydrogen_count.is_some())
            + usize::from(atom.atom_map.is_some()) * 2;
        add(18 + option_bytes)?;
        add(inner
            .stereo_neighbor_order(idx)
            .map_or(0, |order| order.len().saturating_mul(4)))?;
    }
    for (idx, _) in inner.bonds() {
        add(16)?;
        if inner.is_smiles_ring_closure(idx) {
            add(9)?;
        }
    }
    for (group_index, group) in inner.stereo_groups().iter().enumerate() {
        add(13usize.saturating_add(group.atom_indices.len().saturating_mul(4)))?;
        add(inner
            .stereo_group_bond_indices(group_index)
            .map_or(0, |bonds| bonds.len().saturating_mul(4)))?;
    }
    for (key, value) in &mol.props {
        if key.len() > MAX_STRING || value.len() > MAX_STRING {
            return Err(format!(
                "property {key:?} exceeds pickle string limit {MAX_STRING}"
            ));
        }
        add(9 + key.len() + value.len())?;
    }
    Ok(size)
}

pub(crate) fn encode(mol: &Mol) -> Result<Vec<u8>, String> {
    encode_version(mol, VERSION)
}

fn encode_version(mol: &Mol, version: u16) -> Result<Vec<u8>, String> {
    if !matches!(version, LEGACY_VERSION | VERSION) {
        return Err(format!("unsupported pickle version {version}"));
    }
    let inner = mol.inner.as_ref();
    let mut writer = Writer(Vec::with_capacity(checked_payload_size(mol)?));
    writer.0.extend_from_slice(MAGIC);
    writer.u16(version);
    writer.u32(inner.atom_count() as u32);
    writer.u32(inner.bond_count() as u32);
    writer.u32(inner.stereo_groups().len() as u32);
    writer.u32(mol.props.len() as u32);

    for (idx, atom) in inner.atoms() {
        writer.u8(atom.element.atomic_number());
        writer.option_u16(atom.isotope);
        writer.u8(atom.charge as u8);
        writer.option_u8(atom.hydrogen_count);
        writer.u8(atom.aromatic as u8);
        let (kind, value) = chirality_to_parts(atom.chirality, inner.is_allene_stereo(idx));
        writer.u8(kind);
        writer.u8(value);
        writer.u8(atom.wildcard as u8);
        writer.option_u16(atom.atom_map);
        writer.u8(cip_to_u8(atom.cip_code));
        writer.u16(inner.atom_tag(idx).map_or(0, |tag| tag.get()));
        let r_group = match inner.r_group_label(idx) {
            None => 0,
            Some(label) => label.number().map_or(1, |number| number + 1),
        };
        writer.u16(r_group);
        let order = inner.stereo_neighbor_order(idx).unwrap_or(&[]);
        writer.u32(order.len() as u32);
        for &neighbor in order {
            writer.u32(neighbor);
        }
    }

    for (idx, bond) in inner.bonds() {
        writer.u32(bond.atom1.0);
        writer.u32(bond.atom2.0);
        writer.u8(bond_to_u8(bond.order));
        writer.u8(inner.bond_direction(idx).map_or(255, bond_to_u8));
        writer.u32(
            inner
                .bond_direction_anchor(idx)
                .map_or(u32::MAX, |atom| atom.0),
        );
        if let Some((label, occurrence, begins)) = inner.smiles_ring_closure_metadata(idx) {
            writer.u8(1);
            writer.u32(label);
            writer.u32(occurrence);
            writer.u8(match begins {
                None => 0,
                Some(false) => 1,
                Some(true) => 2,
            });
        } else {
            writer.u8(0);
        }
    }

    for (group_index, group) in inner.stereo_groups().iter().enumerate() {
        let (kind, number) = match group.kind {
            StereoGroupKind::Absolute => (0, 0),
            StereoGroupKind::Or(number) => (1, number),
            StereoGroupKind::And(number) => (2, number),
        };
        writer.u8(kind);
        writer.u32(number);
        writer.u32(group.atom_indices.len() as u32);
        for atom in &group.atom_indices {
            writer.u32(atom.0);
        }
        let bonds = inner.stereo_group_bond_indices(group_index).unwrap_or(&[]);
        writer.u32(bonds.len() as u32);
        for bond in bonds {
            writer.u32(bond.0);
        }
    }

    let mut properties: Vec<_> = mol.props.iter().collect();
    properties.sort_unstable_by(|left, right| left.0.cmp(right.0));
    for (key, value) in properties {
        if version >= 2 {
            writer.u8(property_type_to_u8(mol.props.property_type(key)));
        }
        writer.bytes(key.as_bytes());
        writer.bytes(value.as_bytes());
    }
    debug_assert!(writer.0.len() <= MAX_PAYLOAD);
    Ok(writer.0)
}

struct PickleCounts {
    version: u16,
    atoms: usize,
    bonds: usize,
    groups: usize,
    properties: usize,
}

fn decode_header(reader: &mut Reader<'_>) -> Result<PickleCounts, String> {
    if reader.take(MAGIC.len())? != MAGIC {
        return Err("invalid chematic molecule pickle header".into());
    }
    let version = reader.u16()?;
    if !matches!(version, LEGACY_VERSION | VERSION) {
        return Err(format!(
            "unsupported chematic molecule pickle version {version} (expected {LEGACY_VERSION} or {VERSION})"
        ));
    }
    Ok(PickleCounts {
        version,
        atoms: reader.bounded_count("atom", MAX_ATOMS)?,
        bonds: reader.bounded_count("bond", MAX_BONDS)?,
        groups: reader.bounded_count("stereo group", MAX_GROUPS)?,
        properties: reader.bounded_count("property", MAX_PROPERTIES)?,
    })
}

fn decode_stereo_order(reader: &mut Reader<'_>, atom_count: usize) -> Result<Vec<u32>, String> {
    let neighbor_count = reader.bounded_count("stereo neighbor", atom_count + 1)?;
    let mut order = Vec::with_capacity(neighbor_count);
    for _ in 0..neighbor_count {
        let neighbor = reader.u32()?;
        if neighbor != STEREO_H_SENTINEL && neighbor as usize >= atom_count {
            return Err(format!(
                "stereo neighbor atom index {neighbor} is out of range"
            ));
        }
        order.push(neighbor);
    }
    Ok(order)
}

fn apply_r_group(builder: &mut MoleculeBuilder, atom: AtomIdx, encoded: u16) -> Result<(), String> {
    match encoded {
        0 => Ok(()),
        1 => {
            builder.set_r_group(atom, RGroupLabel::unnumbered());
            Ok(())
        }
        2..=10_000 => {
            let label = RGroupLabel::numbered(encoded - 1)
                .ok_or_else(|| format!("invalid R-group label {encoded}"))?;
            builder.set_r_group(atom, label);
            Ok(())
        }
        _ => Err(format!("invalid R-group label encoding {encoded}")),
    }
}

fn decode_atoms(
    reader: &mut Reader<'_>,
    builder: &mut MoleculeBuilder,
    atom_count: usize,
) -> Result<(), String> {
    for atom_index in 0..atom_count {
        let atomic_number = reader.u8()?;
        let element = Element::from_atomic_number(atomic_number)
            .ok_or_else(|| format!("invalid atomic number {atomic_number}"))?;
        let isotope = reader.option_u16("isotope")?;
        let charge = reader.u8()? as i8;
        let hydrogen_count = reader.option_u8("hydrogen count")?;
        let aromatic = reader.bool("aromatic")?;
        let chirality = chirality_from_parts(reader.u8()?, reader.u8()?)?;
        let wildcard = reader.bool("wildcard")?;
        let atom_map = reader.option_u16("atom map")?;
        let cip_code = cip_from_u8(reader.u8()?)?;
        let tag = reader.u16()?;
        let r_group = reader.u16()?;
        let stereo_order = decode_stereo_order(reader, atom_count)?;
        let idx = builder.add_atom(Atom {
            element,
            isotope,
            charge,
            hydrogen_count,
            aromatic,
            chirality,
            wildcard,
            atom_map,
            cip_code,
        });
        debug_assert_eq!(idx.0 as usize, atom_index);
        if tag != 0 {
            builder.set_tag(idx, Some(tag));
        }
        apply_r_group(builder, idx, r_group)?;
        if !stereo_order.is_empty() {
            builder.set_stereo_neighbor_order(idx, stereo_order);
        }
    }
    Ok(())
}

type RingClosure = Option<(u32, u32, Option<bool>)>;

fn decode_ring_closure(reader: &mut Reader<'_>) -> Result<RingClosure, String> {
    if !reader.bool("ring closure")? {
        return Ok(None);
    }
    let label = reader.u32()?;
    let occurrence = reader.u32()?;
    let begins = match reader.u8()? {
        0 => None,
        1 => Some(false),
        2 => Some(true),
        value => return Err(format!("invalid ring-closure endpoint flag {value}")),
    };
    Ok(Some((label, occurrence, begins)))
}

fn decode_bonds(
    reader: &mut Reader<'_>,
    builder: &mut MoleculeBuilder,
    atom_count: usize,
    bond_count: usize,
) -> Result<(), String> {
    for _ in 0..bond_count {
        let atom1 = checked_atom(reader.u32()?, atom_count, "bond endpoint")?;
        let atom2 = checked_atom(reader.u32()?, atom_count, "bond endpoint")?;
        let order = bond_from_u8(reader.u8()?)?;
        let direction = match reader.u8()? {
            255 => None,
            value => Some(bond_from_u8(value)?),
        };
        if direction.is_some_and(|value| !matches!(value, BondOrder::Up | BondOrder::Down)) {
            return Err("bond direction must be Up or Down".into());
        }
        let anchor = match reader.u32()? {
            u32::MAX => None,
            value => Some(checked_atom(value, atom_count, "bond direction anchor")?),
        };
        let ring = decode_ring_closure(reader)?;
        let bond = builder
            .add_bond(atom1, atom2, order)
            .map_err(|error| format!("invalid bond in pickle: {error}"))?;
        if let Some(direction) = direction {
            builder.set_bond_direction(bond, direction);
        }
        if let Some(anchor) = anchor {
            if anchor != atom1 && anchor != atom2 {
                return Err("bond direction anchor is not a bond endpoint".into());
            }
            builder.set_bond_direction_anchor(bond, anchor);
        }
        if let Some((label, occurrence, begins)) = ring {
            builder.set_smiles_ring_closure(bond, label, occurrence);
            if let Some(begins) = begins {
                builder.set_smiles_ring_closure_begins_at_open(bond, begins);
            }
        }
    }
    Ok(())
}

fn decode_group_kind(reader: &mut Reader<'_>) -> Result<StereoGroupKind, String> {
    let kind = reader.u8()?;
    let number = reader.u32()?;
    match (kind, number) {
        (0, 0) => Ok(StereoGroupKind::Absolute),
        (1, number) if number != 0 => Ok(StereoGroupKind::Or(number)),
        (2, number) if number != 0 => Ok(StereoGroupKind::And(number)),
        _ => Err(format!("invalid stereo group kind {kind} number {number}")),
    }
}

fn decode_stereo_groups(
    reader: &mut Reader<'_>,
    builder: &mut MoleculeBuilder,
    counts: &PickleCounts,
) -> Result<(), String> {
    for _ in 0..counts.groups {
        let kind = decode_group_kind(reader)?;
        let atom_members = reader.bounded_count("stereo group atom", counts.atoms)?;
        let mut atoms = Vec::with_capacity(atom_members);
        for _ in 0..atom_members {
            atoms.push(checked_atom(reader.u32()?, counts.atoms, "stereo group")?);
        }
        let bond_members = reader.bounded_count("stereo group bond", counts.bonds)?;
        let mut bonds = Vec::with_capacity(bond_members);
        for _ in 0..bond_members {
            bonds.push(checked_bond(reader.u32()?, counts.bonds, "stereo group")?);
        }
        builder.add_stereo_group_with_bonds(StereoGroup::new(kind, atoms), bonds);
    }
    Ok(())
}

fn decode_properties(
    reader: &mut Reader<'_>,
    count: usize,
    version: u16,
) -> Result<PropertyStore, String> {
    let mut props = PropertyStore::default();
    for _ in 0..count {
        let property_type = if version >= 2 {
            property_type_from_u8(reader.u8()?)?
        } else {
            PropertyType::String
        };
        let key = reader.string("property key")?;
        let value = reader.string("property value")?;
        if props.contains_key(&key) {
            return Err(format!("duplicate property key {key:?}"));
        }
        props.insert_typed(key, value, property_type);
    }
    Ok(props)
}

pub(crate) fn decode(payload: &[u8]) -> Result<Mol, String> {
    let mut reader = Reader::new(payload)?;
    let counts = decode_header(&mut reader)?;
    let mut builder = MoleculeBuilder::with_capacity(counts.atoms, counts.bonds);
    decode_atoms(&mut reader, &mut builder, counts.atoms)?;
    decode_bonds(&mut reader, &mut builder, counts.atoms, counts.bonds)?;
    decode_stereo_groups(&mut reader, &mut builder, &counts)?;
    let props = decode_properties(&mut reader, counts.properties, counts.version)?;
    reader.finish()?;
    Ok(Mol {
        inner: Arc::new(builder.build()),
        props,
    })
}

#[pyfunction(name = "_restore_mol_pickle")]
fn restore_mol_pickle(payload: &Bound<'_, PyBytes>) -> PyResult<Mol> {
    decode(payload.as_bytes()).map_err(PyValueError::new_err)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(restore_mol_pickle, module)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::{Atom, Element};

    fn rich_mol() -> Mol {
        let mut builder = MoleculeBuilder::new();
        let mut carbon = Atom::new(Element::C);
        carbon.isotope = Some(13);
        carbon.chirality = Chirality::Clockwise;
        carbon.atom_map = Some(17);
        carbon.cip_code = Some(CipCode::R);
        let carbon = builder.add_atom(carbon);
        let r_group = builder.add_atom(Atom::wildcard());
        let bond = builder.add_bond(carbon, r_group, BondOrder::Up).unwrap();
        builder.set_tag(carbon, Some(42));
        builder.set_r_group(r_group, RGroupLabel::numbered(7).unwrap());
        builder.set_stereo_neighbor_order(carbon, vec![r_group.0, STEREO_H_SENTINEL]);
        builder.set_bond_direction(bond, BondOrder::Down);
        builder.set_bond_direction_anchor(bond, carbon);
        builder.set_smiles_ring_closure(bond, 12, 3);
        builder.set_smiles_ring_closure_begins_at_open(bond, true);
        builder.add_stereo_group_with_bonds(
            StereoGroup::new(StereoGroupKind::Or(2), vec![carbon]),
            vec![bond],
        );
        Mol {
            inner: Arc::new(builder.build()),
            props: HashMap::from([
                ("name".into(), "molecule".into()),
                ("unicode".into(), "β-lactam".into()),
            ])
            .into(),
        }
    }

    #[test]
    fn round_trip_preserves_graph_sidecars_and_properties() {
        let original = rich_mol();
        let payload = encode(&original).unwrap();
        let restored = decode(&payload).unwrap();
        assert_eq!(
            original.inner.atoms().collect::<Vec<_>>(),
            restored.inner.atoms().collect::<Vec<_>>()
        );
        assert_eq!(
            original.inner.bonds().collect::<Vec<_>>(),
            restored.inner.bonds().collect::<Vec<_>>()
        );
        assert_eq!(
            original.inner.atom_tag(AtomIdx(0)),
            restored.inner.atom_tag(AtomIdx(0))
        );
        assert_eq!(
            original.inner.r_group_label(AtomIdx(1)),
            restored.inner.r_group_label(AtomIdx(1))
        );
        assert_eq!(
            original.inner.stereo_neighbor_order(AtomIdx(0)),
            restored.inner.stereo_neighbor_order(AtomIdx(0))
        );
        assert_eq!(
            original.inner.stereo_groups(),
            restored.inner.stereo_groups()
        );
        assert_eq!(
            original.inner.stereo_group_bond_indices(0),
            restored.inner.stereo_group_bond_indices(0)
        );
        assert_eq!(
            original.inner.bond_direction(BondIdx(0)),
            restored.inner.bond_direction(BondIdx(0))
        );
        assert_eq!(
            original.inner.bond_direction_anchor(BondIdx(0)),
            restored.inner.bond_direction_anchor(BondIdx(0))
        );
        assert_eq!(
            original.inner.smiles_ring_closure_metadata(BondIdx(0)),
            restored.inner.smiles_ring_closure_metadata(BondIdx(0))
        );
        assert_eq!(original.props, restored.props);
    }

    #[test]
    fn encoding_is_deterministic_and_decoder_is_strict() {
        let original = rich_mol();
        assert_eq!(encode(&original).unwrap(), encode(&original).unwrap());
        assert!(decode(b"bad").is_err());
        let mut trailing = encode(&original).unwrap();
        trailing.push(0);
        assert!(decode(&trailing).is_err());
        let mut wrong_version = encode(&original).unwrap();
        wrong_version[MAGIC.len()] = 3;
        assert!(decode(&wrong_version).is_err());

        let legacy = encode_version(&original, LEGACY_VERSION).unwrap();
        let legacy_restored = decode(&legacy).unwrap();
        assert_eq!(legacy_restored.props, original.props);
    }

    #[test]
    fn round_trip_preserves_property_types() {
        let mut original = Mol::bare(MoleculeBuilder::new().build());
        original
            .props
            .insert_typed("double".into(), "42".into(), PropertyType::Double);
        original
            .props
            .insert_typed("bool".into(), "1".into(), PropertyType::Bool);

        let restored = decode(&encode(&original).unwrap()).unwrap();
        assert_eq!(restored.props.property_type("double"), PropertyType::Double);
        assert_eq!(restored.props.property_type("bool"), PropertyType::Bool);
    }
}
