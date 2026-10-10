//! Transcendental functions from the pure-Rust `libm` crate instead of the
//! platform C library (#739).
//!
//! `f64::sin` and friends call the host libm, whose last-bit results differ
//! between glibc, Apple and MSVC; a force-field minimisation amplifies one
//! ULP into a different trajectory, so the same wheel gave a stereo-clean
//! conformer on macOS and a `FinalStereoViolation` on Linux for the same
//! seed. `libm` is what `wasm32` already uses, so native and WASM builds now
//! compute the same values.

#[allow(dead_code)]
pub(crate) trait DetMath: Sized {
    fn dsin(self) -> Self;
    fn dcos(self) -> Self;
    fn dsin_cos(self) -> (Self, Self);
    fn dacos(self) -> Self;
    fn dasin(self) -> Self;
    fn datan(self) -> Self;
    fn datan2(self, x: Self) -> Self;
    fn dexp(self) -> Self;
    fn dln(self) -> Self;
    fn dpowf(self, e: Self) -> Self;
    fn dcbrt(self) -> Self;
    fn dhypot(self, other: Self) -> Self;
}

impl DetMath for f64 {
    fn dsin(self) -> f64 {
        libm::sin(self)
    }
    fn dcos(self) -> f64 {
        libm::cos(self)
    }
    fn dsin_cos(self) -> (f64, f64) {
        libm::sincos(self)
    }
    fn dacos(self) -> f64 {
        libm::acos(self)
    }
    fn dasin(self) -> f64 {
        libm::asin(self)
    }
    fn datan(self) -> f64 {
        libm::atan(self)
    }
    fn datan2(self, x: f64) -> f64 {
        libm::atan2(self, x)
    }
    fn dexp(self) -> f64 {
        libm::exp(self)
    }
    fn dln(self) -> f64 {
        libm::log(self)
    }
    fn dpowf(self, e: f64) -> f64 {
        libm::pow(self, e)
    }
    fn dcbrt(self) -> f64 {
        libm::cbrt(self)
    }
    fn dhypot(self, other: f64) -> f64 {
        libm::hypot(self, other)
    }
}

impl DetMath for f32 {
    fn dsin(self) -> f32 {
        libm::sinf(self)
    }
    fn dcos(self) -> f32 {
        libm::cosf(self)
    }
    fn dsin_cos(self) -> (f32, f32) {
        libm::sincosf(self)
    }
    fn dacos(self) -> f32 {
        libm::acosf(self)
    }
    fn dasin(self) -> f32 {
        libm::asinf(self)
    }
    fn datan(self) -> f32 {
        libm::atanf(self)
    }
    fn datan2(self, x: f32) -> f32 {
        libm::atan2f(self, x)
    }
    fn dexp(self) -> f32 {
        libm::expf(self)
    }
    fn dln(self) -> f32 {
        libm::logf(self)
    }
    fn dpowf(self, e: f32) -> f32 {
        libm::powf(self, e)
    }
    fn dcbrt(self) -> f32 {
        libm::cbrtf(self)
    }
    fn dhypot(self, other: f32) -> f32 {
        libm::hypotf(self, other)
    }
}

#[cfg(test)]
mod numerical_contract_tests {
    use super::DetMath;

    macro_rules! contracts {
        ($name:ident, $ty:ty, $pi:expr, $tol:expr) => {
            #[test]
            fn $name() {
                let close = |actual: $ty, expected: $ty| {
                    assert!(
                        (actual - expected).abs() <= $tol * (1.0 + expected.abs()),
                        "{actual} != {expected}"
                    );
                };
                let pi: $ty = $pi;
                close((pi / 6.0).dsin(), 0.5);
                close((pi / 3.0).dcos(), 0.5);
                let (s, c) = (pi / 4.0).dsin_cos();
                close(s, c);
                close(s * s + c * c, 1.0);
                close((0.5 as $ty).dasin(), pi / 6.0);
                close((0.5 as $ty).dacos(), pi / 3.0);
                close((1.0 as $ty).datan(), pi / 4.0);
                for (y, x, expected) in [
                    (1.0, 1.0, pi / 4.0),
                    (1.0, -1.0, 3.0 * pi / 4.0),
                    (-1.0, -1.0, -3.0 * pi / 4.0),
                    (-1.0, 1.0, -pi / 4.0),
                ] {
                    close((y as $ty).datan2(x), expected);
                }
                for x in [-3.0 as $ty, -0.5, 0.0, 0.5, 3.0] {
                    close(x.dexp().dln(), x);
                }
                close((9.0 as $ty).dpowf(0.5), 3.0);
                close((-2.0 as $ty).dpowf(3.0), -8.0);
                close((-27.0 as $ty).dcbrt(), -3.0);
                close((3.0 as $ty).dhypot(4.0), 5.0);
                close((3e30 as $ty).dhypot(4e30) / (5e30 as $ty), 1.0);
                close((3e-30 as $ty).dhypot(4e-30) / (5e-30 as $ty), 1.0);
                assert!((2.0 as $ty).dasin().is_nan());
                assert!((2.0 as $ty).dacos().is_nan());
                assert!((-1.0 as $ty).dln().is_nan());
                assert!((-2.0 as $ty).dpowf(0.5).is_nan());
                assert_eq!((0.0 as $ty).dln(), <$ty>::NEG_INFINITY);
                assert_eq!(<$ty>::NEG_INFINITY.dexp(), 0.0);
                assert_eq!((-0.0 as $ty).dsin().to_bits(), (-0.0 as $ty).to_bits());
                assert_eq!((-0.0 as $ty).dcbrt().to_bits(), (-0.0 as $ty).to_bits());
                assert_eq!((-0.0 as $ty).datan2(1.0).to_bits(), (-0.0 as $ty).to_bits());
            }
        };
    }

    contracts!(
        f32_math_identities_and_ieee_domains,
        f32,
        std::f32::consts::PI,
        2e-6
    );
    contracts!(
        f64_math_identities_and_ieee_domains,
        f64,
        std::f64::consts::PI,
        2e-13
    );
}
