use super::RdkitEmbedError;

/// Preserve the single seed-zero result, then skip seed one: minstd_rand
/// aliases zero and one. Positive base seeds keep the RDKit schedule.
pub(super) fn conformer_seed(base: i32, index: usize) -> Result<i32, RdkitEmbedError> {
    let multiplier = index
        .checked_add(1)
        .and_then(|n| i32::try_from(n).ok())
        .ok_or_else(|| {
            RdkitEmbedError::Unsupported("numConfs exceeds the reproducible seed range".into())
        })?;
    if base == 0 && index > 0 {
        return Ok(multiplier);
    }
    base.checked_mul(multiplier).ok_or_else(|| {
        RdkitEmbedError::Unsupported("conformer seed overflows a signed 32-bit integer".into())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_and_positive_seed_boundaries_are_checked() {
        assert_eq!(
            (0..5)
                .map(|i| conformer_seed(0, i).unwrap())
                .collect::<Vec<_>>(),
            [0, 2, 3, 4, 5]
        );
        assert_eq!(conformer_seed(42, 2).unwrap(), 126);
        assert_eq!(conformer_seed(0, i32::MAX as usize - 1).unwrap(), i32::MAX);
        assert!(conformer_seed(0, i32::MAX as usize).is_err());
        assert!(conformer_seed(0, usize::MAX).is_err());
        assert!(conformer_seed(i32::MAX, 1).is_err());
    }
}
