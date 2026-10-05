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
