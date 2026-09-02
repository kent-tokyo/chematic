#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    let limits = chematic_mol::MolJsonParseLimits {
        max_input_bytes: 1 << 20,
        max_json_depth: 64,
        max_array_items: 1024,
        max_string_bytes: 100_000,
        max_atoms: 10_000,
        max_bonds: 20_000,
    };
    let _ = chematic_mol::parse_moljson_with_limits(input, &limits);
});
