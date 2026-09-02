#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    let limits = chematic_3d::XyzParseLimits {
        max_input_bytes: 1 << 20,
        max_atoms: 10_000,
        max_line_bytes: 4096,
    };
    let _ = chematic_3d::parse_xyz_with_limits(input, limits);
});
