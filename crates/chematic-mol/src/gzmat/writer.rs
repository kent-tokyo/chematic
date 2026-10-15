use super::{ZMatrixError, ZMatrixRow, error, validate_values, validate_zmatrix};

/// Explicit electronic state, single-line header fields and decimal precision.
/// The route defaults to a template marker (#); supply calculation keywords
/// before submitting a job. No method, charge or multiplicity is inferred.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GzmatWriteOptions {
    pub charge: i32,
    pub multiplicity: u32,
    /// 1..=17 significant decimal digits, in scientific notation.
    /// 17 preserves the computed f64 internal values on parsing.
    pub significant_digits: u8,
    pub title: String,
    pub route: String,
}
impl GzmatWriteOptions {
    pub fn new(charge: i32, multiplicity: u32) -> Self {
        Self {
            charge,
            multiplicity,
            significant_digits: 17,
            title: "chematic internal coordinates".into(),
            route: "#".into(),
        }
    }
    fn validate(&self) -> Result<(), ZMatrixError> {
        if !(1..=17).contains(&self.significant_digits) {
            return Err(error(None, "significant_digits must be 1..=17"));
        }
        if self.multiplicity == 0 {
            return Err(error(None, "multiplicity must be positive"));
        }
        for value in [&self.title, &self.route] {
            if value.trim().is_empty() || value.len() > 1024 || value.chars().any(char::is_control)
            {
                return Err(error(
                    None,
                    "header fields must be nonempty single lines of at most 1024 bytes",
                ));
            }
        }
        if !self.route.starts_with('#') {
            return Err(error(None, "Gaussian route must start with #"));
        }
        Ok(())
    }
}

fn formatted_values(i: usize, row: &ZMatrixRow, digits: u8) -> Result<Vec<String>, ZMatrixError> {
    let values: Vec<String> = row
        .values
        .iter()
        .map(|v| {
            // Normalize negative zero without changing a nonzero signed torsion.
            let v = if *v == 0.0 { 0.0 } else { *v };
            format!("{:.*e}", usize::from(digits) - 1, v)
        })
        .collect();
    let rounded: Vec<f64> = values
        .iter()
        .map(|s| {
            s.parse::<f64>()
                .map_err(|_| error(Some(i), "formatted value is invalid"))
        })
        .collect::<Result<_, _>>()?;
    // A low-precision request must not silently turn a valid angle into 180.
    validate_values(i, &rounded)?;
    Ok(values)
}

fn rounded_values(rows: &[ZMatrixRow], digits: u8) -> Result<Vec<Vec<String>>, ZMatrixError> {
    let formatted = rows
        .iter()
        .enumerate()
        .map(|(i, row)| formatted_values(i, row, digits))
        .collect::<Result<Vec<_>, _>>()?;
    let rounded = rows
        .iter()
        .zip(&formatted)
        .enumerate()
        .map(|(i, (row, values))| {
            let values = values
                .iter()
                .map(|v| {
                    v.parse::<f64>()
                        .map_err(|_| error(Some(i), "formatted value is invalid"))
                })
                .collect::<Result<_, _>>()?;
            Ok(ZMatrixRow {
                element: row.element,
                references: row.references.clone(),
                values,
            })
        })
        .collect::<Result<Vec<_>, ZMatrixError>>()?;
    // Also check that rounding did not collapse an earlier reference plane.
    super::zmatrix_to_cartesian(&rounded)?;
    Ok(formatted)
}

/// Emit a Gaussian-style Z-matrix with independent r/a/d variables per atom.
/// Serializes references one-based; does not emit bonds, isotopes or periodicity.
/// Rounding that destroys a defined coordinate is rejected.
pub fn write_gzmat(
    rows: &[ZMatrixRow],
    options: &GzmatWriteOptions,
) -> Result<String, ZMatrixError> {
    options.validate()?;
    validate_zmatrix(rows)?;
    let formatted = rounded_values(rows, options.significant_digits)?;
    let mut output = format!(
        "{}\n\n{}\n\n{} {}\n",
        options.route, options.title, options.charge, options.multiplicity
    );
    let mut variables = String::new();
    for (i, row) in rows.iter().enumerate() {
        let values = &formatted[i];
        output.push_str(row.element.symbol());
        for (k, (&reference, value)) in row.references.iter().zip(values).enumerate() {
            let name = format!("{}{}", ['r', 'a', 'd'][k], i + 1);
            output.push_str(&format!(" {} {}", reference + 1, name));
            variables.push_str(&format!("{name}= {value}\n"));
        }
        output.push('\n');
    }
    if !variables.is_empty() {
        output.push_str("Variables:\n");
        output.push_str(&variables);
    }
    output.push('\n');
    Ok(output)
}
