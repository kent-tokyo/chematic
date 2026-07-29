//! IO-3 acceptance-gate dump: parses every fixture named in
//! `scripts/gen_rdkit_mrv_oracle.py`'s manifest with chematic's own
//! `parse_mrv`, dumps each fixture's extracted atoms/bonds/coordinates plus
//! chematic's own isomeric SMILES output as JSONL for
//! `scripts/mrv_io_parity.py`, which re-canonicalizes both chematic's and
//! RDKit's output through RDKit itself (never a direct chematic-vs-RDKit
//! canonicalizer string comparison -- ring-closure digits, atom traversal
//! order, and branch representation legitimately differ between
//! canonicalizers for the same graph).
//!
//! Also writes each parsed molecule back out via `write_mrv` into
//! `<out_dir>/chematic_written/<id>.mrv` (chematic-write -> RDKit-read leg),
//! and performs a chematic-only parse -> write -> parse round trip,
//! reporting per-fixture identity (chematic-write -> chematic-read leg,
//! entirely independent of RDKit).
//!
//! Usage:
//! ```text
//! cargo run -p chematic-mol --release --example mrv_dump -- \
//!     <manifest.json> <fixtures_dir> <out.jsonl> <written_out_dir>
//! ```

use chematic_mol::{MrvWriteOptions, parse_mrv, write_mrv};
use serde_json::{Value, json};
use std::fs;
use std::io::Write;

fn atom_json(mol: &chematic_core::Molecule) -> Vec<Value> {
    mol.atoms()
        .map(|(_, a)| {
            json!({
                "symbol": a.element.symbol(),
                "charge": a.charge,
                "isotope": a.isotope,
                "atom_map": a.atom_map,
                "aromatic": a.aromatic,
            })
        })
        .collect()
}

fn bond_json(mol: &chematic_core::Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(_, b)| {
            json!({
                "begin": b.atom1.0,
                "end": b.atom2.0,
                "order": format!("{:?}", b.order),
            })
        })
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let manifest_path = args
        .get(1)
        .expect("usage: mrv_dump <manifest.json> <fixtures_dir> <out.jsonl> <written_out_dir>");
    let fixtures_dir = args.get(2).expect("usage: ...");
    let out_path = args.get(3).expect("usage: ...");
    let written_dir = args.get(4).expect("usage: ...");

    fs::create_dir_all(written_dir).unwrap_or_else(|e| panic!("create {written_dir}: {e}"));

    let manifest_text =
        fs::read_to_string(manifest_path).unwrap_or_else(|e| panic!("read manifest: {e}"));
    let manifest: Value =
        serde_json::from_str(&manifest_text).unwrap_or_else(|e| panic!("parse manifest: {e}"));

    let mut out = fs::File::create(out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));
    let fixtures = manifest["fixtures"].as_array().expect("fixtures array");

    let mut total = 0usize;
    let mut errors = 0usize;

    for fixture in fixtures {
        let id = fixture["id"].as_str().unwrap();
        let category = fixture["category"].as_str().unwrap();
        let file_name = fixture["file"].as_str().unwrap();
        let path = format!("{fixtures_dir}/{file_name}");
        let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path}: {e}"));

        total += 1;
        let row = match parse_mrv(&text) {
            Ok(rec) => {
                // chematic's own isomeric SMILES output -- this is what the
                // comparator feeds back into RDKit (Chem.MolFromSmiles) as
                // the "B" side of the same-tool recanonicalization gate.
                // Re-perceive aromaticity first: RDKit's MolToMrvBlock
                // defaults to kekulize=True, so fixture bonds are plain
                // Kekule single/double (order != "A"), and chematic's
                // reader correctly leaves those atoms non-aromatic -- a
                // bare `chematic_smiles::write` here would emit a
                // structurally-correct but needlessly Kekule SMILES.
                let aromatic_view = chematic_perception::apply_aromaticity(&rec.mol);
                let chematic_isomeric_smiles = chematic_smiles::write(&aromatic_view);

                // Write chematic's own MRV output for the chematic-write ->
                // RDKit-read leg (kekulize=false: preserves the aromatic
                // bond order token as parsed, avoiding a needless
                // representation change on top of the one already being
                // measured).
                let write_opts = MrvWriteOptions {
                    kekulize: false,
                    ..Default::default()
                };
                let mut write_error = None;
                match write_mrv(&rec, &write_opts) {
                    Ok(block) => {
                        let written_path = format!("{written_dir}/{id}.mrv");
                        fs::write(&written_path, &block)
                            .unwrap_or_else(|e| panic!("write {written_path}: {e}"));
                    }
                    Err(e) => write_error = Some(e.to_string()),
                }

                // chematic-only round trip: parse -> write -> parse again,
                // compare atom/bond arrays for identity. Independent of
                // RDKit and of the kekulize option above (round trip uses
                // the same non-kekulized write so structure is preserved
                // exactly, not just re-derivable via re-aromatization).
                let round_trip_ok = match write_mrv(&rec, &write_opts) {
                    Ok(block) => match parse_mrv(&block) {
                        Ok(rec2) => {
                            atom_json(&rec.mol) == atom_json(&rec2.mol)
                                && bond_json(&rec.mol) == bond_json(&rec2.mol)
                        }
                        Err(_) => false,
                    },
                    Err(_) => false,
                };

                json!({
                    "id": id,
                    "category": category,
                    "status": "success",
                    "atom_count": rec.mol.atom_count(),
                    "bond_count": rec.mol.bond_count(),
                    "atoms": atom_json(&rec.mol),
                    "bonds": bond_json(&rec.mol),
                    "coordinates_2d": rec.coordinates_2d,
                    "coordinates_3d": rec.coordinates_3d,
                    "chematic_isomeric_smiles": chematic_isomeric_smiles,
                    "write_error": write_error,
                    "round_trip_ok": round_trip_ok,
                })
            }
            Err(e) => {
                errors += 1;
                json!({
                    "id": id,
                    "category": category,
                    "status": "error",
                    "error": e.to_string(),
                })
            }
        };
        writeln!(out, "{row}").unwrap();
    }

    eprintln!("total={total} errors={errors} out={out_path}");
}
