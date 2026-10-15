//! Decimal formatter oracle for padded V2000 and unpadded V3000 fields.
use super::*;

#[test]
fn decimal_table_and_fallback_match_integer_formatter() {
    let mut values: Vec<u32> = (0..=1024).collect();
    for power in [
        10u32, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000,
    ] {
        values.extend([power - 1, power, power + 1]);
    }
    values.extend([u32::MAX - 1, u32::MAX]);
    let mut seed = 0x72646974u32;
    for _ in 0..10000 {
        seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
        values.push(seed);
    }
    for width in 0..=12 {
        for &value in &values {
            let mut actual = String::from("prefix:");
            push_right_aligned_u32(&mut actual, value, width);
            assert_eq!(actual, format!("prefix:{value:>width$}"));
        }
    }
}
