//! Ranked-screening EF boundary regressions motivated by RDKit #7981.
use chematic_chem::{ScreeningError, enrichment_factors};

#[test]
fn issue_7981_top_ninety_percent_includes_last_rank() {
    // The issue's descending (score, label) order is [true,true,false,false,false].
    assert_eq!(
        enrichment_factors(&[true, true, false, false, false], &[0.9]).unwrap(),
        [1.0]
    );
}

#[test]
fn full_library_and_singleton_return_one() {
    for labels in [
        vec![true],
        vec![false, true],
        vec![true, true, true],
        vec![true, false, true, false],
    ] {
        assert_eq!(enrichment_factors(&labels, &[1.0]).unwrap(), [1.0]);
    }
}

#[test]
fn fractions_keep_order_duplicates_and_shared_cutoffs() {
    let actual = enrichment_factors(
        &[true, false, true, false, false],
        &[1.0, 0.2, 0.01, 0.2, 0.5, 0.9],
    )
    .unwrap();
    let expected = [1.0, 2.5, 2.5, 2.5, 5.0 / 3.0, 1.0];
    for (a, b) in actual.iter().zip(expected) {
        assert!((a - b).abs() < 1e-14);
    }
}

#[test]
fn selected_count_not_nominal_fraction_sets_denominator() {
    // f=0.21 selects 2 of 5 records, both active: EF = 1 / (2/5) = 2.5.
    assert_eq!(
        enrichment_factors(&[true, true, false, false, false], &[0.21]).unwrap(),
        [2.5]
    );
}

#[test]
fn no_top_hits_can_have_zero_enrichment() {
    assert_eq!(
        enrichment_factors(&[false, false, true, true], &[0.5, 1.0]).unwrap(),
        [0.0, 1.0]
    );
}

#[test]
fn invalid_inputs_are_typed_errors() {
    assert_eq!(
        enrichment_factors(&[], &[1.]),
        Err(ScreeningError::EmptyRanking)
    );
    assert_eq!(
        enrichment_factors(&[true], &[]),
        Err(ScreeningError::EmptyFractions)
    );
    assert_eq!(
        enrichment_factors(&[false, false], &[1.]),
        Err(ScreeningError::NoActives)
    );
    for value in [
        0.,
        -0.0,
        -0.1,
        1.0000001,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        let error = enrichment_factors(&[true], &[1., value]).unwrap_err();
        assert!(matches!(
            error,
            ScreeningError::InvalidFraction { index: 1, .. }
        ));
        assert!(error.to_string().contains("index 1"));
    }
}

#[test]
fn tiny_positive_fraction_selects_one_record() {
    assert_eq!(
        enrichment_factors(&[true, false], &[f64::from_bits(1)]).unwrap(),
        [2.0]
    );
}

#[test]
fn exhaustive_binary_rankings_match_integer_count_oracle() {
    // Exact quarter fractions allow an integer-only rank/count oracle.
    let quarters = [4usize, 1, 3, 2, 1, 4];
    let fractions: Vec<_> = quarters.iter().map(|&q| q as f64 / 4.).collect();
    for n in 1..=8 {
        for mask in 1..1usize << n {
            let labels: Vec<_> = (0..n).map(|i| mask & (1 << i) != 0).collect();
            let total = mask.count_ones() as usize;
            let actual = enrichment_factors(&labels, &fractions).unwrap();
            assert_eq!(actual.len(), fractions.len());
            for (&q, &value) in quarters.iter().zip(&actual) {
                let k = (n * q).div_ceil(4);
                let hits = (mask & ((1 << k) - 1)).count_ones() as usize;
                let expected = (n * hits) as f64 / (k * total) as f64;
                assert!(
                    (value - expected).abs() < 1e-14,
                    "n={n}, mask={mask}, q={q}"
                );
            }
        }
    }
}
