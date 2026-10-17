// Each pair is ASCII, and the table occupies 200 bytes.
const DECIMAL_PAIRS: &[u8; 200] = b"00010203040506070809101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899";

/// Append `x` exactly as `format!("{x:>width$.4}")` would. Keep the general
/// formatter for rounding ties, very large coordinates and nonfinite values.
#[inline]
pub(crate) fn push_fixed4(out: &mut String, x: f64, width: usize) {
    use std::fmt::Write as _;
    let scaled = x.abs() * 10_000.0;
    let frac = scaled - scaled.floor();
    if !x.is_finite() || scaled >= 1e15 || (frac - 0.5).abs() < 1e-6 {
        write!(out, "{x:>width$.4}").expect("writing to String cannot fail");
        return;
    }
    let units = scaled.round() as u64;
    let (mut integer, decimal) = (units / 10_000, units % 10_000);
    let mut buf = [b' '; 24];
    let mut pos = buf.len() - 4;
    let high = (decimal / 100) as usize * 2;
    let low = (decimal % 100) as usize * 2;
    buf[pos..pos + 2].copy_from_slice(&DECIMAL_PAIRS[high..high + 2]);
    buf[pos + 2..].copy_from_slice(&DECIMAL_PAIRS[low..low + 2]);
    pos -= 1;
    buf[pos] = b'.';
    while integer >= 100 {
        let pair = (integer % 100) as usize * 2;
        integer /= 100;
        pos -= 2;
        buf[pos..pos + 2].copy_from_slice(&DECIMAL_PAIRS[pair..pair + 2]);
    }
    if integer < 10 {
        pos -= 1;
        buf[pos] = b'0' + integer as u8;
    } else {
        pos -= 2;
        let pair = integer as usize * 2;
        buf[pos..pos + 2].copy_from_slice(&DECIMAL_PAIRS[pair..pair + 2]);
    }
    if x.is_sign_negative() {
        pos -= 1;
        buf[pos] = b'-';
    }
    if width <= buf.len() {
        pos = pos.min(buf.len() - width);
    } else {
        for _ in (buf.len() - pos)..width {
            out.push(' ');
        }
    }
    out.push_str(std::str::from_utf8(&buf[pos..]).expect("ASCII digits and padding"));
}

#[cfg(test)]
#[path = "fixed4_tests.rs"]
mod tests;
