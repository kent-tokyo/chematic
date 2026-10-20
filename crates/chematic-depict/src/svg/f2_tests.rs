use super::F2;

#[test]
fn f2_matches_the_float_formatter() {
    let mut values = vec![
        0.0,
        -0.0,
        0.005,
        -0.005,
        0.015,
        0.025,
        1.005,
        2.675,
        123.455,
        -0.004,
        -0.0049,
        0.0049,
        99.995,
        1e14,
        3.0e15,
        f64::NAN,
        f64::INFINITY,
        -1.5,
    ];
    // Exact ties (multiples of 1/8) and their neighbours.
    for m in -4000..4000 {
        let v = m as f64 / 8.0;
        values.push(v);
        values.push(f64::from_bits(v.to_bits() + 1));
        values.push(f64::from_bits(v.to_bits().saturating_sub(1)));
        values.push(m as f64 * 0.005);
        values.push(m as f64 * 0.125 + 0.005);
    }
    let mut seed = 0x2545_f491_4f6c_dd1du64;
    for _ in 0..200_000 {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        let v = (seed % 2_000_000_000) as f64 / 1_000_000.0 - 1000.0;
        values.push(v);
        values.push(v / 997.0);
    }
    for v in values {
        assert_eq!(format!("{:.2}", F2(v)), format!("{v:.2}"), "{v:e}");
    }
}
