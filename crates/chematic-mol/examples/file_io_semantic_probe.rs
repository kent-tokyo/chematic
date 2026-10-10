//! Produce an engine-neutral semantic snapshot or a same-format rewrite.
//!
//! This helper is intentionally an example binary, not public API. The
//! Open-Babel comparison gate uses it to judge both engines through the same
//! chematic reader after conversion.

use std::env;
use std::error::Error;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use chematic_core::Molecule;
use chematic_mol::{
    Mol2Record, MolMetadata, parse_cdxml_all, parse_cml, parse_mol_v3000_with_coords,
    parse_mol2_record, write_cdxml, write_cml, write_mol_v3000, write_mol2_record,
};
use serde::Serialize;
use serde_json::{Value, json};

#[derive(Serialize)]
struct Snapshot {
    schema_version: u8,
    format: String,
    records: Vec<RecordSnapshot>,
}

#[derive(Serialize)]
struct RecordSnapshot {
    graph_key: String,
    native_key: String,
    atom_count: usize,
    bond_count: usize,
    total_formal_charge: i32,
    isotopes: Vec<(usize, u16)>,
    pairwise_distances: Vec<f64>,
    metadata: Value,
}

enum ParsedRecord {
    V3000 {
        molecule: Molecule,
        metadata: MolMetadata,
        coords: Vec<(f64, f64)>,
    },
    Mol2(Mol2Record),
    Cml {
        molecule: Molecule,
        coords: Vec<(f64, f64)>,
    },
    Cdxml {
        molecule: Molecule,
        coords: Vec<(f64, f64)>,
    },
}

impl ParsedRecord {
    fn molecule(&self) -> &Molecule {
        match self {
            Self::V3000 { molecule, .. }
            | Self::Cml { molecule, .. }
            | Self::Cdxml { molecule, .. } => molecule,
            Self::Mol2(record) => &record.molecule,
        }
    }

    fn coords3(&self) -> Vec<(f64, f64, f64)> {
        match self {
            Self::V3000 { coords, .. } | Self::Cml { coords, .. } | Self::Cdxml { coords, .. } => {
                coords.iter().map(|&(x, y)| (x, y, 0.0)).collect()
            }
            Self::Mol2(record) => record.coords.clone(),
        }
    }

    fn metadata(&self) -> Value {
        match self {
            Self::V3000 { metadata, .. } => json!({
                "name": metadata.name,
                "comment": metadata.comment,
                "sgroups": metadata.v3000_sgroups,
                "bond_properties": metadata.v3000_bond_properties,
                "atom_properties": metadata.v3000_atom_properties,
            }),
            Self::Mol2(record) => json!({
                "name": record.name,
                "molecule_type": record.molecule_type,
                "charge_type": record.charge_type,
                "atoms": record.atoms.iter().map(|atom| json!({
                    "name": atom.atom_name,
                    "type": atom.atom_type,
                    "subst_id": atom.subst_id,
                    "subst_name": atom.subst_name,
                    "partial_charge": atom.partial_charge,
                    "status_bits": atom.status_bits,
                })).collect::<Vec<_>>(),
                "bonds": record.bonds.iter().map(|bond| json!({
                    "type": bond.bond_type,
                    "status_bits": bond.status_bits,
                })).collect::<Vec<_>>(),
                "duplicate_bonds": record.duplicate_bonds.iter().map(|bond| json!({
                    "atom1_id": bond.atom1_id,
                    "atom2_id": bond.atom2_id,
                    "type": bond.bond_type,
                    "status_bits": bond.status_bits,
                })).collect::<Vec<_>>(),
                "unity_atom_attributes": record.unity_atom_attributes.iter().map(|group| json!({
                    "atom_id": group.atom_id,
                    "attributes": group.attributes,
                })).collect::<Vec<_>>(),
                "opaque_sections": record.opaque_sections.iter().map(|section| json!({
                    "name": section.name,
                    "lines": section.lines,
                })).collect::<Vec<_>>(),
            }),
            Self::Cml { .. } | Self::Cdxml { .. } => json!({}),
        }
    }

    fn rewrite(&self) -> String {
        match self {
            Self::V3000 {
                molecule,
                metadata,
                coords,
            } => write_mol_v3000(molecule, metadata, coords),
            Self::Mol2(record) => write_mol2_record(record),
            Self::Cml { molecule, coords } => write_cml(molecule, Some(coords)),
            Self::Cdxml { molecule, coords } => write_cdxml(molecule, coords),
        }
    }
}

fn argument(name: &str) -> Option<String> {
    let mut args = env::args().skip(1);
    while let Some(value) = args.next() {
        if value == name {
            return args.next();
        }
    }
    None
}

fn benchmark(
    format: &str,
    text: &str,
    operation: &str,
    repeats: usize,
) -> Result<(), Box<dyn Error>> {
    if repeats == 0 {
        return Err("--repeats must be positive".into());
    }
    let mut bytes = 0usize;
    let mut records = 0usize;
    let parsed = (operation == "write")
        .then(|| parse_records(format, text))
        .transpose()?;
    let started = Instant::now();
    for _ in 0..repeats {
        match operation {
            "parse" => {
                let current = parse_records(format, text)?;
                records += current.len();
                bytes = bytes.wrapping_add(
                    current
                        .iter()
                        .map(|record| record.molecule().atom_count())
                        .sum::<usize>(),
                );
            }
            "write" => {
                let current = parsed.as_ref().expect("write input is parsed above");
                records += current.len();
                for record in current {
                    bytes = bytes.wrapping_add(record.rewrite().len());
                }
            }
            "roundtrip" => {
                let current = parse_records(format, text)?;
                records += current.len();
                for record in &current {
                    bytes = bytes.wrapping_add(record.rewrite().len());
                }
            }
            _ => return Err(format!("unsupported benchmark operation {operation}").into()),
        }
    }
    println!(
        "{}",
        json!({
            "engine": "chematic",
            "format": format,
            "operation": operation,
            "repeats": repeats,
            "records": records,
            "output_units": bytes,
            "elapsed_ns": started.elapsed().as_nanos(),
        })
    );
    Ok(())
}

fn parse_records(format: &str, text: &str) -> Result<Vec<ParsedRecord>, Box<dyn Error>> {
    match format {
        "v3000" => {
            let (molecule, metadata, coords) = parse_mol_v3000_with_coords(text)?;
            Ok(vec![ParsedRecord::V3000 {
                molecule,
                metadata,
                coords,
            }])
        }
        "mol2" => Ok(vec![ParsedRecord::Mol2(parse_mol2_record(text)?)]),
        "cml" => {
            let (molecule, coords) = parse_cml(text)?;
            Ok(vec![ParsedRecord::Cml { molecule, coords }])
        }
        "cdxml" => Ok(parse_cdxml_all(text)?
            .into_iter()
            .map(|(molecule, coords)| ParsedRecord::Cdxml { molecule, coords })
            .collect()),
        _ => Err(format!("unsupported --format {format}").into()),
    }
}

fn round_float(value: f64) -> f64 {
    (value * 1_000_000.0).round() / 1_000_000.0
}

fn snapshot(record: &ParsedRecord) -> RecordSnapshot {
    let molecule = record.molecule();
    let (native_key, canonical_order) = chematic_smiles::canonical_smiles_with_atom_order(molecule);
    let graph_key =
        chematic_smiles::rdkit_canonical_smiles(molecule).unwrap_or_else(|_| native_key.clone());
    let coords = record.coords3();
    let ordered = canonical_order
        .iter()
        .filter_map(|idx| coords.get(idx.0 as usize).copied())
        .collect::<Vec<_>>();
    let mut pairwise_distances = Vec::with_capacity(ordered.len().saturating_pow(2) / 2);
    for left in 0..ordered.len() {
        for right in left + 1..ordered.len() {
            let (ax, ay, az) = ordered[left];
            let (bx, by, bz) = ordered[right];
            pairwise_distances.push(round_float(
                ((ax - bx).powi(2) + (ay - by).powi(2) + (az - bz).powi(2)).sqrt(),
            ));
        }
    }
    RecordSnapshot {
        graph_key,
        native_key,
        atom_count: molecule.atom_count(),
        bond_count: molecule.bond_count(),
        total_formal_charge: molecule
            .atoms()
            .map(|(_, atom)| i32::from(atom.charge))
            .sum(),
        isotopes: canonical_order
            .iter()
            .enumerate()
            .filter_map(|(canonical_idx, idx)| {
                molecule
                    .atom(*idx)
                    .isotope
                    .map(|isotope| (canonical_idx, isotope))
            })
            .collect(),
        pairwise_distances,
        metadata: record.metadata(),
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let format = argument("--format").ok_or("missing --format")?;
    let path = PathBuf::from(argument("--path").ok_or("missing --path")?);
    let text = fs::read_to_string(&path)?;
    if let Some(operation) = argument("--benchmark-operation") {
        let repeats = argument("--repeats")
            .unwrap_or_else(|| "1000".into())
            .parse::<usize>()?;
        return benchmark(&format, &text, &operation, repeats);
    }
    let records = parse_records(&format, &text)?;
    if let Some(output) = argument("--rewrite-output") {
        let rewritten = records
            .iter()
            .map(ParsedRecord::rewrite)
            .collect::<Vec<_>>()
            .join("\n");
        fs::write(output, rewritten)?;
        return Ok(());
    }
    let result = Snapshot {
        schema_version: 1,
        format,
        records: records.iter().map(snapshot).collect(),
    };
    println!("{}", serde_json::to_string_pretty(&result)?);
    Ok(())
}
