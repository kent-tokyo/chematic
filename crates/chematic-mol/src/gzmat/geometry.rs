use super::{ZMATRIX_MIN_SINE, ZMatrixError, error};

type V = [f64; 3];
fn sub(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn scale(a: V, s: f64) -> V {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(v: V) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn direction(i: usize, a: V, b: V) -> Result<(V, f64), ZMatrixError> {
    let delta = sub(a, b);
    let length = norm(delta);
    if !length.is_finite() || length <= 0.0 {
        return Err(error(
            Some(i),
            "coincident reference or overflowing distance",
        ));
    }
    // Divide componentwise: computing 1/length first can overflow for tiny values.
    Ok((
        [delta[0] / length, delta[1] / length, delta[2] / length],
        length,
    ))
}

pub(super) fn internal_values(i: usize, p: &[V], r: &[usize]) -> Result<Vec<f64>, ZMatrixError> {
    if r.is_empty() {
        return Ok(Vec::new());
    }
    let (u, distance) = direction(i, p[i], p[r[0]])?;
    let mut values = vec![distance];
    if r.len() == 1 {
        return Ok(values);
    }
    let (axis, _) = direction(i, p[r[1]], p[r[0]])?;
    let sine = norm(cross(u, axis));
    if sine <= ZMATRIX_MIN_SINE {
        return Err(error(Some(i), "angle is linear or near-linear"));
    }
    values.push(sine.atan2(dot(u, axis)).to_degrees());
    if r.len() == 2 {
        return Ok(values);
    }
    let (w, _) = direction(i, p[r[2]], p[r[1]])?;
    if norm(cross(axis, w)) <= ZMATRIX_MIN_SINE {
        return Err(error(Some(i), "reference plane is linear or near-linear"));
    }
    let v0 = sub(u, scale(axis, dot(u, axis)));
    let v2 = sub(w, scale(axis, dot(w, axis)));
    values.push(dot(cross(axis, v0), v2).atan2(dot(v0, v2)).to_degrees());
    Ok(values)
}

/// Reconstruct a canonical Cartesian frame (angstroms), preserving handedness.
/// Absolute translation and rotation are not part of an internal-coordinate set.
pub fn zmatrix_to_cartesian(rows: &[super::ZMatrixRow]) -> Result<Vec<V>, ZMatrixError> {
    super::validate_zmatrix(rows)?;
    let mut positions = Vec::with_capacity(rows.len());
    for (i, row) in rows.iter().enumerate() {
        let position = place_row(i, row, &positions)?;
        if position.iter().any(|v| !v.is_finite()) {
            return Err(error(Some(i), "reconstructed coordinate overflows"));
        }
        positions.push(position);
    }
    Ok(positions)
}
fn place_row(i: usize, row: &super::ZMatrixRow, p: &[V]) -> Result<V, ZMatrixError> {
    if i == 0 {
        return Ok([0.0; 3]);
    }
    let a = p[row.references[0]];
    let distance = row.values[0];
    if i == 1 {
        return Ok([a[0], a[1], a[2] + distance]);
    }
    let b = p[row.references[1]];
    let (axis, _) = direction(i, b, a)?;
    let perpendicular = if i == 2 {
        let basis = if axis[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        sub(basis, scale(axis, dot(basis, axis)))
    } else {
        let (w, _) = direction(i, p[row.references[2]], b)?;
        sub(w, scale(axis, dot(w, axis)))
    };
    let length = norm(perpendicular);
    if length <= ZMATRIX_MIN_SINE {
        return Err(error(Some(i), "reference plane is linear or near-linear"));
    }
    let e2 = scale(perpendicular, 1.0 / length);
    let e3 = cross(axis, e2);
    let angle = row.values[1].to_radians();
    let torsion = row.values.get(2).copied().unwrap_or(0.0).to_radians();
    let mut result = [0.0; 3];
    for k in 0..3 {
        let unit =
            angle.cos() * axis[k] + angle.sin() * (torsion.cos() * e2[k] - torsion.sin() * e3[k]);
        result[k] = a[k] + distance * unit;
    }
    Ok(result)
}
