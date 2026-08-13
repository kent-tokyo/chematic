//! Baseline performance measurements. Correctness comes first in `v0.1`
//! (see `docs/crystal_scope.md`'s "no spatial-partitioning optimization"
//! non-goal) -- these benchmarks exist to record where the naive O(n^2 *
//! search-box) neighbor search currently stands, not to defend a
//! performance target. Recorded numbers (CPU, rustc version, commit SHA,
//! site count, cell shape, cutoff, neighbor count, method) are kept in
//! this crate's README under "Benchmarks".

use chematic_core::Element;
use chematic_crystal::{FractionalCoord, Lattice, PeriodicSite, PeriodicStructure, SiteSpecies};
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

/// `n` sites on a simple grid inside a cubic cell sized so the average
/// nearest-neighbor spacing stays ~3 Angstrom regardless of `n` (density
/// held roughly constant across site counts, so cutoff-relative behavior
/// is comparable).
fn grid_structure_cubic(n_per_axis: usize) -> PeriodicStructure {
    let spacing = 3.0;
    let a = spacing * n_per_axis as f64;
    let lattice = Lattice::cubic(a).unwrap();
    let mut sites = Vec::new();
    for i in 0..n_per_axis {
        for j in 0..n_per_axis {
            for k in 0..n_per_axis {
                let f = FractionalCoord::new([
                    i as f64 / n_per_axis as f64,
                    j as f64 / n_per_axis as f64,
                    k as f64 / n_per_axis as f64,
                ]);
                sites
                    .push(PeriodicSite::new(vec![SiteSpecies::full(Element::C)], f, None).unwrap());
            }
        }
    }
    PeriodicStructure::new(lattice, sites).unwrap()
}

/// Same site count/spacing as [`grid_structure_cubic`], but a skewed
/// triclinic cell instead of cubic.
fn grid_structure_triclinic(n_per_axis: usize) -> PeriodicStructure {
    let spacing = 3.0;
    let a = spacing * n_per_axis as f64;
    let lattice = Lattice::from_parameters(a, a * 1.05, a * 0.95, 75.0, 100.0, 60.0).unwrap();
    let mut sites = Vec::new();
    for i in 0..n_per_axis {
        for j in 0..n_per_axis {
            for k in 0..n_per_axis {
                let f = FractionalCoord::new([
                    i as f64 / n_per_axis as f64,
                    j as f64 / n_per_axis as f64,
                    k as f64 / n_per_axis as f64,
                ]);
                sites
                    .push(PeriodicSite::new(vec![SiteSpecies::full(Element::C)], f, None).unwrap());
            }
        }
    }
    PeriodicStructure::new(lattice, sites).unwrap()
}

fn bench_frac_to_cart(c: &mut Criterion) {
    let lattice = Lattice::from_parameters(5.0, 6.0, 7.0, 75.0, 100.0, 60.0).unwrap();
    let points: Vec<FractionalCoord> = (0..1000)
        .map(|i| {
            let t = i as f64 / 1000.0;
            FractionalCoord::new([t, (t * 1.3) % 1.0, (t * 2.7) % 1.0])
        })
        .collect();
    c.bench_function("frac_to_cart_1000", |b| {
        b.iter(|| {
            for p in &points {
                let _ = black_box(lattice.frac_to_cart(black_box(*p)));
            }
        })
    });
}

fn bench_cart_to_frac(c: &mut Criterion) {
    let lattice = Lattice::from_parameters(5.0, 6.0, 7.0, 75.0, 100.0, 60.0).unwrap();
    let points: Vec<_> = (0..1000)
        .map(|i| {
            let t = i as f64;
            chematic_crystal::CartesianCoord::new([t % 5.0, (t * 1.3) % 6.0, (t * 2.7) % 7.0])
        })
        .collect();
    c.bench_function("cart_to_frac_1000", |b| {
        b.iter(|| {
            for p in &points {
                let _ = black_box(lattice.cart_to_frac(black_box(*p)));
            }
        })
    });
}

fn bench_neighbor_search(c: &mut Criterion) {
    // n_per_axis^3 sites: 5^3=125 (~100), 10^3=1000. The naive O(n^2 *
    // search-box) baseline (see this file's module doc -- spatial
    // partitioning is an explicit non-goal for v0.1) makes a ~10000-site
    // structure (22^3=10648) take on the order of a minute per *single*
    // call; criterion's minimum sample_size(10) would make that an
    // impractically slow `cargo bench` run for what v0.1 already documents
    // as "not optimized yet". That size is instead measured once via plain
    // wall-clock timing and recorded directly in README.md's Benchmarks
    // section rather than run through criterion here.
    for (label, n_per_axis) in [("100", 5usize), ("1000", 10)] {
        let cubic = grid_structure_cubic(n_per_axis);
        let cutoff = 6.5; // a few shells at 3 Angstrom grid spacing
        c.bench_function(&format!("neighbor_search_cubic_{label}"), |b| {
            b.iter(|| black_box(cubic.neighbors_within(black_box(cutoff)).unwrap()))
        });

        let triclinic = grid_structure_triclinic(n_per_axis);
        c.bench_function(&format!("neighbor_search_triclinic_{label}"), |b| {
            b.iter(|| black_box(triclinic.neighbors_within(black_box(cutoff)).unwrap()))
        });
    }
}

fn bench_supercell(c: &mut Criterion) {
    let base = grid_structure_cubic(5); // 125 sites
    c.bench_function("supercell_generation_2x2x2_from_125_sites", |b| {
        b.iter(|| black_box(base.make_supercell(black_box([2, 2, 2])).unwrap()))
    });
}

criterion_group!(
    benches,
    bench_frac_to_cart,
    bench_cart_to_frac,
    bench_neighbor_search,
    bench_supercell
);
criterion_main!(benches);
