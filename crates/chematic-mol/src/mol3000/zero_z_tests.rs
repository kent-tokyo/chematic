//! Differential oracle; see the v3000-zero-z fixture recipe and manifest.
use super::*;

fn reference_coords(tokens: &[&str], line: usize) -> Result<(f64, f64, f64), MolParseError> {
    let parse = |index: usize, axis: &str| {
        tokens[index]
            .parse::<f64>()
            .map_err(|_| MolParseError::InvalidAtomLine {
                line,
                detail: format!("cannot parse {axis} coordinate from '{}'", tokens[index]),
            })
            .and_then(|value| {
                value
                    .is_finite()
                    .then_some(value)
                    .ok_or_else(|| MolParseError::InvalidAtomLine {
                        line,
                        detail: format!(
                            "{axis} coordinate is not finite (NaN/Infinite): '{}'",
                            tokens[index]
                        ),
                    })
            })
    };
    Ok((parse(2, "x")?, parse(3, "y")?, parse(4, "z")?))
}

fn compare(x: &str, y: &str, z: &str, line: usize) {
    let tokens = ["1", "C", x, y, z];
    let expected = reference_coords(&tokens, line);
    let actual = parse_v3000_coords(&tokens, line);
    match (expected, actual) {
        (Ok((x, y, z)), Ok((a, b, c))) => assert_eq!(
            [x.to_bits(), y.to_bits(), z.to_bits()],
            [a.to_bits(), b.to_bits(), c.to_bits()]
        ),
        (Err(left), Err(right)) => assert_eq!(format!("{left:?}"), format!("{right:?}")),
        (left, right) => panic!("coordinate result differs: {left:?} vs {right:?}"),
    }
}

fn recipe() -> serde_json::Value {
    serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/benchmark_corpora/v3000-zero-z-speed-v2.recipe.json"
    )))
    .unwrap()
}

#[test]
fn zero_spellings_and_failure_precedence_match_original() {
    let recipe = recipe();
    let values: Vec<_> = recipe["lexemes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    for line in recipe["lines"].as_array().unwrap() {
        // Bound the portable test location to the platform's usize range.
        let line = line.as_u64().unwrap().min(usize::MAX as u64) as usize;
        for &x in &values {
            for &y in &values {
                for &z in &values {
                    compare(x, y, z, line);
                }
            }
        }
    }
    assert_eq!(values.len(), 25);
}

#[test]
fn deterministic_f64_fallback_and_signed_zero_match_original() {
    let recipe = recipe();
    let random = &recipe["random_f64"];
    let mut bits = u64::from_str_radix(random["seed_hex"].as_str().unwrap(), 16).unwrap();
    for _ in 0..random["count"].as_u64().unwrap() {
        bits ^= bits << 13;
        bits ^= bits >> 7;
        bits ^= bits << 17;
        let value = f64::from_bits(bits).to_string();
        compare("1", "-1", &value, 42);
    }
}
