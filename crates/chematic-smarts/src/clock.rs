//! Monotonic clock, portable across native and `wasm32-unknown-unknown`.
//!
//! `std::time::Instant::now()` panics with "time not implemented on this
//! platform" under real `wasm32-unknown-unknown` (no host time source) --
//! see <https://github.com/kent-tokyo/chematic/issues/221> (same underlying
//! trap as chematic-3d's <https://github.com/kent-tokyo/chematic/issues/219>,
//! independently fixed here since chematic-3d depends on this crate, not the
//! other way around). `web_time::Instant` is API-compatible and backed by
//! `Performance.now()` there, while transparently re-exporting
//! `std::time::Instant` everywhere else.

#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
pub(crate) use web_time::Instant;

#[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
pub(crate) use std::time::Instant;
