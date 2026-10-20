//! Bounded stdin adapter for the CJSON → MOL2 development comparison.
//! This example exercises the shared checked converter, not a separate writer.
use chematic_mol::{ConversionLimits, convert_with_report_json};
use std::io::{self, Read};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let limits = ConversionLimits::default();
    let mut source = String::new();
    io::stdin()
        .take(limits.mol2.max_input_bytes as u64 + 1)
        .read_to_string(&mut source)?;
    println!(
        "{}",
        convert_with_report_json(&source, "cjson", "mol2", "allow_loss", &limits)?
    );
    Ok(())
}
