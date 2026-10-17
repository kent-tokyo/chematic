use super::push_fixed4;
#[inline]
/// Append `x` as `format!("{x:>width$.4}")` would, without the general
/// float formatter (its exact fallback dominated writing laid-out
/// molecules). Values near a rounding tie, very large values and
/// non-finite values take the formatter, so the bytes are always the same.
fn frozen_push_fixed4(out: &mut String, x: f64, width: usize) {
    use std::fmt::Write as _;
    let scaled = x.abs() * 10_000.0;
    let frac = scaled - scaled.floor();
    if !x.is_finite() || scaled >= 1e15 || (frac - 0.5).abs() < 1e-6 {
        write!(out, "{x:>width$.4}").expect("writing to String cannot fail");
        return;
    }
    let units = scaled.round() as u64;
    let mut buf = [0u8; 24];
    let mut pos = buf.len();
    let mut push = |b: u8| {
        pos -= 1;
        buf[pos] = b;
    };
    let (mut int, mut dec) = (units / 10_000, units % 10_000);
    for _ in 0..4 {
        push(b'0' + (dec % 10) as u8);
        dec /= 10;
    }
    push(b'.');
    loop {
        push(b'0' + (int % 10) as u8);
        int /= 10;
        if int == 0 {
            break;
        }
    }
    if x.is_sign_negative() {
        push(b'-');
    }
    let digits = &buf[pos..];
    for _ in digits.len()..width {
        out.push(' ');
    }
    out.push_str(std::str::from_utf8(digits).expect("ASCII digits"));
}

fn check(x: f64) {
    for width in [0, 1, 6, 10, 24, 30] {
        let mut actual = String::from("prefix:");
        let mut frozen = actual.clone();
        push_fixed4(&mut actual, x, width);
        frozen_push_fixed4(&mut frozen, x, width);
        assert_eq!(actual, frozen, "bits={:016x}, width={width}", x.to_bits());
        assert_eq!(
            actual,
            format!("prefix:{x:>width$.4}"),
            "bits={:016x}",
            x.to_bits()
        );
    }
}
#[test]
fn decimal_pairs_match_frozen_and_standard_formatter() {
    for n in 0..10_000 {
        for x in [
            n as f64 / 10_000.0,
            -(n as f64 / 10_000.0),
            n as f64 + 0.0001,
        ] {
            check(x);
        }
    }
    let mut bits = 0x5678_1234_abcd_9876u64;
    for _ in 0..100_000 {
        bits ^= bits << 13;
        bits ^= bits >> 7;
        bits ^= bits << 17;
        check(f64::from_bits(bits));
    }
}
#[test]
fn rounding_edges_signs_and_padding_match() {
    for x in [
        0.0,
        -0.0,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        -f64::NAN,
        f64::MIN_POSITIVE,
        f64::MAX,
        -f64::MAX,
        f64::from_bits(1),
    ] {
        check(x);
    }
    for n in [
        0u64,
        1,
        9,
        10,
        99,
        100,
        999,
        1000,
        9999,
        10000,
        99999,
        1_000_000,
        99_999_999_999,
    ] {
        for fraction in [0.0, 0.00005, 0.00995, 0.09995, 0.99995] {
            let x = n as f64 + fraction;
            for value in [
                x,
                x.next_up(),
                x.next_down(),
                -x,
                (-x).next_up(),
                (-x).next_down(),
            ] {
                check(value);
            }
        }
    }
}
