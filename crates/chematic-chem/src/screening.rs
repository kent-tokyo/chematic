//! Evaluation of an already ranked virtual-screening library.
//!
//! Labels are supplied in best-to-worst order. Sorting scores and resolving
//! ties belong to the caller, so an activity label never breaks a score tie
//! implicitly. This module does not reproduce RDKit's boundary bugs (#7981).

/// Invalid input to [`enrichment_factors`].
#[derive(Debug, Clone, PartialEq)]
pub enum ScreeningError {
    EmptyRanking,
    EmptyFractions,
    /// Fractions must be finite and in `(0, 1]`.
    InvalidFraction {
        index: usize,
        value: f64,
    },
    /// With no actives, the background active rate and EF denominator are zero.
    NoActives,
}

impl core::fmt::Display for ScreeningError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyRanking => write!(f, "screening ranking is empty"),
            Self::EmptyFractions => write!(f, "screening fraction list is empty"),
            Self::InvalidFraction { index, value } => write!(
                f,
                "fraction at index {index} must be finite and in (0, 1], got {value}"
            ),
            Self::NoActives => write!(
                f,
                "enrichment factor is undefined for a ranking with no actives"
            ),
        }
    }
}

impl std::error::Error for ScreeningError {}

/// Compute enrichment factors for an already ranked library, best first.
///
/// For `N` labels, `A` total actives and fraction `f`, selects the first
/// `k = ceil(N * f)` entries and returns `(actives_in_top_k / k) / (A / N)`.
/// The selected count, rather than the nominal fraction, supplies the
/// denominator. A fraction of 1 therefore returns 1, even for a single entry.
/// Fractions may be unsorted or repeated; output has exactly the same order
/// and length as `fractions`.
///
/// Returns an error for an empty ranking, empty fractions, fractions outside
/// finite `(0, 1]`, or no actives (EF is undefined). Labels are binary, not
/// weights. Scores and ties are not inspected or reordered. Work is
/// `O(N + F log F)`, with `O(F)` auxiliary memory for `F` requested fractions.
///
/// ```
/// use chematic_chem::enrichment_factors;
/// let labels = [true, true, false, false, false];
/// assert_eq!(enrichment_factors(&labels, &[0.9, 1.0]).unwrap(), [1.0, 1.0]);
/// ```
pub fn enrichment_factors(
    ranked_actives: &[bool],
    fractions: &[f64],
) -> Result<Vec<f64>, ScreeningError> {
    if ranked_actives.is_empty() {
        return Err(ScreeningError::EmptyRanking);
    }
    if fractions.is_empty() {
        return Err(ScreeningError::EmptyFractions);
    }
    let n = ranked_actives.len();
    let mut cutoffs = Vec::with_capacity(fractions.len());
    for (index, &value) in fractions.iter().enumerate() {
        if !value.is_finite() || value <= 0.0 || value > 1.0 {
            return Err(ScreeningError::InvalidFraction { index, value });
        }
        let cutoff = ((n as f64 * value).ceil() as usize).clamp(1, n);
        cutoffs.push((cutoff, index));
    }
    let total_actives = ranked_actives.iter().filter(|&&active| active).count();
    if total_actives == 0 {
        return Err(ScreeningError::NoActives);
    }
    cutoffs.sort_unstable();
    let background_rate = total_actives as f64 / n as f64;
    let mut result = vec![0.0; fractions.len()];
    let mut scanned = 0;
    let mut hits = 0;
    for (cutoff, index) in cutoffs {
        hits += ranked_actives[scanned..cutoff]
            .iter()
            .filter(|&&active| active)
            .count();
        scanned = cutoff;
        result[index] = (hits as f64 / cutoff as f64) / background_rate;
    }
    Ok(result)
}
