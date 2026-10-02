//! Output-only operation matrix from the locked crates.io `chematic` dependency.
//! This lane maps the published Python operation matrix onto the published
//! Rust crate. It captures outputs only; it is not a timing comparison.

use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

use chematic::{chem, core, depict, fp, inchi, mol as mol_io, perception, smarts, smiles, threed};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

const SPECS: &[(&str, usize)] = &[
    ("parse_smiles", 5000),
    ("canonical_smiles(prepared)", 5000),
    ("parse+canonical_smiles", 5000),
    ("mol_block_write", 2000),
    ("mol_block_read", 2000),
    ("inchi", 1000),
    ("inchikey", 1000),
    ("mw", 5000),
    ("exact_mass", 5000),
    ("logp", 5000),
    ("mr", 5000),
    ("tpsa", 5000),
    ("rdkit_tpsa", 5000),
    ("hbd", 5000),
    ("hba", 5000),
    ("rotatable_bonds", 5000),
    ("fsp3", 5000),
    ("ring_count", 5000),
    ("aromatic_ring_count", 5000),
    ("qed", 5000),
    ("kappa2", 5000),
    ("chi1v", 5000),
    ("bertz_ct", 5000),
    ("formula", 5000),
    ("num_stereocenters", 5000),
    ("lipinski_bundle(mw,logp,hbd,hba)", 5000),
    ("parse+logp", 5000),
    ("parse+tpsa", 5000),
    ("parse+ring_count", 5000),
    ("parse+aromatic_ring_count", 5000),
    ("morgan_r2_2048(native ecfp4)", 5000),
    ("morgan_r2_2048(rdkit-compatible)", 5000),
    ("maccs", 5000),
    ("atom_pair(rdkit-compatible)", 5000),
    ("torsion(rdkit-compatible)", 5000),
    ("rdkit_fp(rdkit-compatible)", 2000),
    ("pattern_fp(rdkit-compatible)", 2000),
    ("tanimoto_1xN(per target, compatible Morgan)", 50),
    ("has_substructure c1ccccc1", 2000),
    ("has_substructure [OH]", 2000),
    ("has_substructure C(=O)N", 2000),
    ("has_substructure [#7;R]", 2000),
    ("has_substructure c1ccc2ccccc2c1", 2000),
    ("has_substructure [CX3](=O)[OX2H1]", 2000),
    ("has_substructure [NX3;H2,H1;!$(NC=O)]", 2000),
    ("has_substructure *~*~*~*~*~*", 2000),
    ("find_matches [#6]~[#7]", 2000),
    ("murcko_scaffold", 5000),
    ("add_hydrogens", 5000),
    ("remove_hydrogens", 5000),
    ("largest_fragment", 2000),
    ("neutralize", 2000),
    ("standardize vs Cleanup", 300),
    ("standardize vs Cleanup+Uncharge+Canonicalize", 300),
    ("canonical_tautomer", 300),
    ("brics_fragments", 500),
    ("cip_labels", 1000),
    ("sssr_rings", 5000),
    ("svg_depiction", 300),
    ("2d_layout", 500),
    ("embed_3d", 100),
    ("embed+minimize_mmff94", 10),
    ("mcs_pair", 50),
];

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn fp_hex(bits: &fp::BitVec2048, bytes: usize) -> String {
    let mut hex = String::with_capacity(bytes * 2);
    for byte_index in 0..bytes {
        let mut byte = 0u8;
        for bit in 0..8 {
            if bits.get(byte_index * 8 + bit) {
                byte |= 1 << bit;
            }
        }
        hex.push_str(&format!("{byte:02x}"));
    }
    hex
}

fn product(mol: &core::Molecule) -> Value {
    json!({"mol_smiles": smiles::canonical_smiles(mol)})
}

fn cip_labels(mol: &core::Molecule) -> Result<Value, String> {
    let assignments = chem::assign_cip(mol).assignments;
    let mut ez_bonds = chem::assign_ez_bonds_with_mode(mol, chem::CipMode::LegacyFast);
    let rows: Vec<_> = assignments
        .iter()
        .map(|(atom, code)| {
            let descriptor = match code {
                core::CipCode::R => "R",
                core::CipCode::S => "S",
                core::CipCode::E => "E",
                core::CipCode::Z => "Z",
                core::CipCode::LowerR => "r",
                core::CipCode::LowerS => "s",
            };
            if matches!(code, core::CipCode::E | core::CipCode::Z) {
                let Some(pos) = ez_bonds
                    .iter()
                    .position(|&(bond_idx, c)| c == *code && mol.bond(bond_idx).atom1 == *atom)
                else {
                    return Err(format!("missing E/Z bond for atom {}", atom.0));
                };
                let (bond_idx, _) = ez_bonds.remove(pos);
                let bond = mol.bond(bond_idx);
                Ok(json!({"atom_idx": atom.0, "bond_idx": bond_idx.0,
                    "bond_atoms": [bond.atom1.0, bond.atom2.0], "descriptor": descriptor}))
            } else {
                Ok(json!({"atom_idx": atom.0, "descriptor": descriptor}))
            }
        })
        .collect::<Result<_, _>>()?;
    Ok(json!(rows))
}

fn depict_data(mol: &core::Molecule) -> Value {
    let data = depict::compute_depict_data(mol);
    let atoms: Vec<_> = data
        .atoms
        .iter()
        .map(|a| {
            json!({
                "idx": a.idx.0, "element": a.element.symbol(), "x": a.pos.x, "y": a.pos.y,
                "label": a.label, "color": a.color, "charge": a.charge,
            })
        })
        .collect();
    let bonds: Vec<_> = data
        .bonds
        .iter()
        .map(|b| {
            json!({
                "idx": b.idx.0, "atom1": b.atom1.0, "atom2": b.atom2.0,
                "kind": format!("{:?}", b.kind),
            })
        })
        .collect();
    json!({"atoms": atoms, "bonds": bonds})
}

fn coords_json(coords: &threed::Coords3D) -> Value {
    json!(
        coords
            .points
            .iter()
            .map(|p| [p.x, p.y, p.z])
            .collect::<Vec<_>>()
    )
}

// The Python binding reconstructs a concrete Molecule from an MCS query.
// Mirror that public-binding adaptation, not a new chemistry algorithm.
fn mcs_molecule(query: &smarts::QueryMolecule) -> Option<core::Molecule> {
    use core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
    use smarts::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery};
    if query.atoms.is_empty() {
        return None;
    }
    fn atomic_num(q: &AtomQuery) -> Option<u8> {
        match q {
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)) => Some(*n),
            AtomQuery::And(left, right) => atomic_num(left).or_else(|| atomic_num(right)),
            _ => None,
        }
    }
    let mut aromatic = vec![false; query.atoms.len()];
    for (i, neighbors) in query.adj.iter().enumerate() {
        for &(bond_idx, neighbor_idx) in neighbors {
            if matches!(
                query.bonds[bond_idx].query,
                BondQuery::Primitive(BondPrimitive::Aromatic)
            ) {
                aromatic[i] = true;
                aromatic[neighbor_idx] = true;
            }
        }
    }
    let mut builder = MoleculeBuilder::new();
    for (i, query_atom) in query.atoms.iter().enumerate() {
        let element = atomic_num(&query_atom.query)
            .and_then(Element::from_atomic_number)
            .unwrap_or(Element::C);
        let mut atom = Atom::new(element);
        atom.aromatic = aromatic[i];
        builder.add_atom(atom);
    }
    for (i, neighbors) in query.adj.iter().enumerate() {
        for &(bond_idx, neighbor_idx) in neighbors {
            if i >= neighbor_idx {
                continue;
            }
            let order = match &query.bonds[bond_idx].query {
                BondQuery::Primitive(BondPrimitive::Double) => BondOrder::Double,
                BondQuery::Primitive(BondPrimitive::Triple) => BondOrder::Triple,
                BondQuery::Primitive(BondPrimitive::Aromatic) => BondOrder::Aromatic,
                _ => BondOrder::Single,
            };
            let _ = builder.add_bond(AtomIdx(i as u32), AtomIdx(neighbor_idx as u32), order);
        }
    }
    Some(builder.build())
}

fn evaluate(
    name: &str,
    mol: &core::Molecule,
    index: usize,
    molecules: &[core::Molecule],
    compatible_fps: &[fp::BitVec2048],
) -> Result<Value, String> {
    let value = match name {
        "parse_smiles" => json!({"mol_smiles": smiles::canonical_smiles(mol)}),
        "canonical_smiles(prepared)" | "parse+canonical_smiles" => {
            json!(smiles::canonical_smiles(mol))
        }
        "mol_block_write" => json!(mol_io::write_mol(mol, &mol_io::MolMetadata::default())),
        "mol_block_read" => {
            let block = mol_io::write_mol(mol, &mol_io::MolMetadata::default());
            product(&mol_io::parse_mol(&block).map_err(|e| e.to_string())?.0)
        }
        "inchi" => json!(inchi::inchi(mol)),
        "inchikey" => json!(inchi::inchi_key(&inchi::inchi(mol))),
        "mw" => json!(chem::rdkit_molecular_weight(mol)),
        "exact_mass" => json!(chem::exact_mass(mol)),
        "logp" | "parse+logp" => json!(chem::logp_crippen(mol)),
        "mr" => json!(chem::molar_refractivity(mol)),
        "tpsa" | "parse+tpsa" => json!(chem::tpsa(mol)),
        "rdkit_tpsa" => json!(chem::rdkit_tpsa(mol)),
        "hbd" => json!(chem::hbd_count(mol)),
        "hba" => json!(chem::rdkit_hba_count(mol)),
        "rotatable_bonds" => json!(chem::rotatable_bond_count(mol)),
        "fsp3" => json!(chem::fsp3(mol)),
        "ring_count" | "parse+ring_count" => json!(chem::ring_count(mol)),
        "aromatic_ring_count" | "parse+aromatic_ring_count" => {
            json!(chem::aromatic_ring_count(mol))
        }
        "qed" => json!(chem::qed(mol)),
        "kappa2" => json!(chem::kappa2(mol)),
        "chi1v" => json!(chem::chi1v(mol)),
        "bertz_ct" => json!(chem::bertz_ct(mol)),
        "formula" => json!(chem::calc_mol_formula(mol)),
        "num_stereocenters" => json!(chem::num_stereocenters(mol)),
        "lipinski_bundle(mw,logp,hbd,hba)" => json!([
            chem::rdkit_molecular_weight(mol),
            chem::logp_crippen(mol),
            chem::hbd_count(mol) as f64,
            chem::rdkit_hba_count(mol) as f64,
        ]),
        "morgan_r2_2048(native ecfp4)" => json!({"bytes_hex": fp_hex(&fp::ecfp4(mol), 256)}),
        "morgan_r2_2048(rdkit-compatible)" => json!({"bytes_hex": fp_hex(
            &fp::rdkit_morgan_ecfp4_bitvec(mol).map_err(|e| e.to_string())?, 256)}),
        "maccs" => json!({"bytes_hex": fp_hex(&fp::maccs(mol), 21)}),
        "atom_pair(rdkit-compatible)" => {
            json!({"bytes_hex": fp_hex(&fp::rdkit_atom_pair_fp(mol), 256)})
        }
        "torsion(rdkit-compatible)" => {
            json!({"bytes_hex": fp_hex(&fp::rdkit_torsion_fp(mol), 256)})
        }
        "rdkit_fp(rdkit-compatible)" => json!({"bytes_hex": fp_hex(&fp::rdkit_rdk_fp(mol), 256)}),
        "pattern_fp(rdkit-compatible)" => {
            json!({"bytes_hex": fp_hex(&fp::rdkit_pattern_fp(mol), 256)})
        }
        "tanimoto_1xN(per target, compatible Morgan)" => {
            json!(fp::tanimoto_slice(&compatible_fps[index], compatible_fps))
        }
        "find_matches [#6]~[#7]" => {
            let query = smarts::parse_smarts("[#6]~[#7]").map_err(|e| e.to_string())?;
            json!(smarts::find_match_atom_sets_perceived(
                &query,
                mol,
                &smarts::MatchConfig::default()
            ))
        }
        "murcko_scaffold" => product(&chem::murcko_scaffold(mol)),
        "add_hydrogens" => product(&chem::add_hydrogens(mol)),
        "remove_hydrogens" => product(&chem::remove_hydrogens(mol)),
        "largest_fragment" => product(&chem::largest_fragment(mol)),
        "neutralize" => product(&chem::neutralize_charges(mol)),
        "standardize vs Cleanup" | "standardize vs Cleanup+Uncharge+Canonicalize" => product(
            &chem::standardize(mol, &chem::StandardizeOptions::default()),
        ),
        "canonical_tautomer" => product(&chem::canonical_tautomer(mol)),
        "brics_fragments" => json!(
            chem::brics_fragments(mol)
                .iter()
                .map(product)
                .collect::<Vec<_>>()
        ),
        "cip_labels" => cip_labels(mol)?,
        "sssr_rings" => json!(perception::ring_membership(mol)),
        "svg_depiction" => json!(depict::depict_svg(mol)),
        "2d_layout" => depict_data(mol),
        "embed_3d" => coords_json(&threed::generate_coords_etkdg(mol)),
        "embed+minimize_mmff94" => {
            let coords = threed::generate_coords_etkdg(mol);
            coords_json(&threed::minimize_mmff94(mol, coords))
        }
        "mcs_pair" => {
            let pair = [&molecules[2 * index], &molecules[2 * index + 1]];
            let query = smarts::find_mcs_with_config(
                &pair,
                &smarts::McsConfig {
                    timeout_ms: Some(2000),
                    ..smarts::McsConfig::default()
                },
            );
            mcs_molecule(&query).map_or(Value::Null, |m| product(&m))
        }
        name if name.starts_with("has_substructure ") => {
            let pattern = name.strip_prefix("has_substructure ").unwrap_or("");
            let query = smarts::parse_smarts(pattern).map_err(|e| e.to_string())?;
            let config = smarts::MatchConfig {
                max_matches: Some(1),
                uniquify: false,
                ..smarts::MatchConfig::default()
            };
            json!(smarts::has_match_perceived(&query, mol, &config))
        }
        _ => return Err(format!("unmapped operation: {name}")),
    };
    Ok(value)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 5 || !matches!(args[4].as_str(), "1.0.29" | "1.0.30") {
        return Err("usage: operations CORPUS ROWS_JSONL SUMMARY_JSON 1.0.29|1.0.30".into());
    }
    let corpus_bytes = fs::read(&args[1])?;
    let smiles: Vec<_> = String::from_utf8(corpus_bytes.clone())?
        .lines()
        .take(5000)
        .map(|line| line.split_whitespace().next().unwrap_or("").to_string())
        .collect();
    if smiles.len() != 5000 || smiles.iter().any(String::is_empty) {
        return Err("corpus must contain at least 5,000 indexed SMILES".into());
    }
    let molecules: Vec<_> = smiles
        .iter()
        .map(|s| smiles::parse(s))
        .collect::<Result<_, _>>()?;
    let compatible_fps: Vec<_> = molecules
        .iter()
        .map(fp::rdkit_morgan_ecfp4_bitvec)
        .collect::<Result<_, _>>()?;
    let mut raw = Vec::new();
    for &(name, count) in SPECS {
        let rows: Vec<_> = molecules
            .iter()
            .take(count)
            .enumerate()
            .map(|(index, mol)| match evaluate(name, mol, index, &molecules, &compatible_fps) {
                Ok(value) => json!({"input_index": index, "value": value}),
                Err(message) => json!({"input_index": index, "error": {"type": "Error", "message": message}}),
            })
            .collect();
        raw.extend_from_slice(
            serde_json::to_string(&json!({"op": name, "rows": rows}))?.as_bytes(),
        );
        raw.push(b'\n');
    }
    let rows_path = Path::new(&args[2]);
    let mut writer = BufWriter::new(fs::File::create(rows_path)?);
    writer.write_all(&raw)?;
    writer.flush()?;
    let summary = json!({
        "schema": "published-rust-python-63op-output-slice/v1",
        "scope": "output only; not a timing or RDKit parity claim",
        "crate": format!("chematic {} from crates.io", args[4]),
        "corpus_sha256": digest(&corpus_bytes),
        "input_count": smiles.len(),
        "adapted_operations": SPECS.len(),
        "rows_sha256": digest(&raw),
    });
    fs::write(&args[3], serde_json::to_vec_pretty(&summary)?)?;
    println!(
        "published Rust operation output: {} operations",
        SPECS.len()
    );
    Ok(())
}
