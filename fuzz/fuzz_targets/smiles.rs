#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(mol) = chematic_smiles::parse(input) {
        let canonical = chematic_smiles::canonical_smiles(&mol);
        // Empty component input such as "." can produce an empty molecule;
        // its canonical representation is intentionally not reparsed here.
        if !canonical.is_empty() {
            let _ = chematic_smiles::parse(&canonical);
        }
    }
});
