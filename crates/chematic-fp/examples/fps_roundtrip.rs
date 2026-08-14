//! Minimal FPS write + read demo: compute ECFP4 for a couple of molecules,
//! write them to an in-memory `.fps` buffer, then read the buffer back.
//!
//! Run with: `cargo run -p chematic-fp --example fps_roundtrip`

use chematic_fp::{FpsHeader, FpsReader, FpsWriter, ecfp4};
use chematic_smiles::parse;
use std::io::Cursor;

fn main() {
    let molecules = [("benzene", "c1ccccc1"), ("ethane", "CC"), ("methane", "C")];

    let header = FpsHeader::for_chematic(2048, "ECFP4").with_source("fps_roundtrip.rs example");

    let mut buf = Vec::new();
    {
        let mut writer = FpsWriter::new(&mut buf, &header).expect("write header");
        for (id, smiles) in molecules {
            let mol = parse(smiles).expect("valid SMILES");
            let fp = ecfp4(&mol);
            writer.write_record_2048(id, &fp).expect("write record");
        }
    }

    print!("{}", String::from_utf8_lossy(&buf));

    let reader = FpsReader::new(Cursor::new(&buf)).expect("parse header");
    println!("\n# num_bits={}", reader.header().num_bits);
    for record in reader {
        let record = record.expect("valid record");
        println!(
            "read back: {} (popcount={})",
            record.id,
            record.fingerprint.popcount()
        );
    }
}
