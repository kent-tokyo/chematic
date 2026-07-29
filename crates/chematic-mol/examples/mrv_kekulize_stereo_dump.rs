//! IO-3 dedicated kekulize/stereo option dump -- independent of the main
//! 206-fixture oracle pool (`mrv_dump.rs`, which always uses the default
//! `kekulize=false, include_stereo=true` options). Verifies the
//! *non-default* option values specifically:
//!
//! - kekulize=True/False: does the written bond-order token shape match
//!   the option (order="A" absent/present), and does RDKit read both
//!   variants back to the same canonical structure (kekulize is a
//!   representation choice, not a structural change)?
//! - include_stereo=False: given a REAL RDKit-generated MRV fixture that
//!   already encodes tetrahedral/E-Z stereo via a native wedge/dash bond
//!   (reused from the main oracle pool's tetrahedral_stereo_*/ez_stereo_*
//!   fixtures -- not synthesized from SMILES chirality, which tests an
//!   unsupported, different direction: chematic has no converter from
//!   Atom.chirality to a wedge bond on write, only the reverse read path;
//!   see mrv_io_parity.py's `tetrahedral_or_ez_stereo_lost_...` finding),
//!   does turning the option off correctly drop the stereo assignment
//!   when RDKit re-reads the output (documented, expected loss)?
//!   `include_stereo=True`'s round trip is already covered by the main
//!   pool's own phase2 (chematic-write -> RDKit-read) check.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-mol --release --example mrv_kekulize_stereo_dump -- \
//!     <fixtures_dir> <out.json>
//! ```

use chematic_mol::{MrvWriteOptions, parse_mrv};
use serde_json::json;
use std::fs;

const AROMATIC_CASES: &[(&str, &str)] = &[("benzene", "c1ccccc1"), ("pyridine", "c1ccncc1")];
const STEREO_FIXTURE_IDS: &[&str] = &[
    "tetrahedral_stereo_0",
    "tetrahedral_stereo_1",
    "ez_stereo_0",
    "ez_stereo_1",
];

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let fixtures_dir = args
        .get(1)
        .expect("usage: mrv_kekulize_stereo_dump <fixtures_dir> <out.json>");
    let out_path = args.get(2).expect("usage: ...");

    let mut out = Vec::new();

    for (id, smi) in AROMATIC_CASES {
        let mol = chematic_smiles::parse(smi).unwrap_or_else(|e| panic!("parse {smi}: {e}"));
        let aromatic_mol = chematic_perception::apply_aromaticity(&mol);
        let record = chematic_mol::MoleculeRecord::new(aromatic_mol);

        let kekulized = chematic_mol::write_mrv(
            &record,
            &MrvWriteOptions {
                kekulize: true,
                ..Default::default()
            },
        )
        .unwrap();
        let non_kekulized = chematic_mol::write_mrv(
            &record,
            &MrvWriteOptions {
                kekulize: false,
                ..Default::default()
            },
        )
        .unwrap();

        out.push(json!({
            "id": id,
            "kind": "kekulize",
            "known_smiles": smi,
            "kekulize_true_mrv": kekulized,
            "kekulize_false_mrv": non_kekulized,
            "kekulize_true_has_aromatic_token": kekulized.contains("order=\"A\""),
            "kekulize_false_has_aromatic_token": non_kekulized.contains("order=\"A\""),
        }));
    }

    for id in STEREO_FIXTURE_IDS {
        let path = format!("{fixtures_dir}/{id}.mrv");
        let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path}: {e}"));
        let rec = parse_mrv(&text).unwrap_or_else(|e| panic!("parse {id}: {e}"));

        let with_stereo = chematic_mol::write_mrv(
            &rec,
            &MrvWriteOptions {
                include_stereo: true,
                kekulize: false,
                ..Default::default()
            },
        )
        .unwrap();
        let without_stereo = chematic_mol::write_mrv(
            &rec,
            &MrvWriteOptions {
                include_stereo: false,
                kekulize: false,
                ..Default::default()
            },
        )
        .unwrap();

        out.push(json!({
            "id": id,
            "kind": "stereo",
            "original_mrv": text,
            "with_stereo_mrv": with_stereo,
            "without_stereo_mrv": without_stereo,
        }));
    }

    fs::write(out_path, serde_json::to_string_pretty(&out).unwrap())
        .unwrap_or_else(|e| panic!("write {out_path}: {e}"));
    eprintln!("wrote {} cases to {out_path}", out.len());
}
