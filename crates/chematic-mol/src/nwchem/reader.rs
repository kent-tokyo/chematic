use super::*;

pub fn parse_nwchem_geometry(text: &str) -> Result<NwchemGeometry, NwchemError> {
    parse_nwchem_geometry_with_limits(text, &NwchemLimits::default())
}
pub fn parse_nwchem_geometry_with_limits(
    text: &str,
    limits: &NwchemLimits,
) -> Result<NwchemGeometry, NwchemError> {
    if text.len() > limits.max_input_bytes || !text.is_ascii() {
        return Err(error(0, "input exceeds byte limit or is not ASCII"));
    }
    let mut geometry: Option<NwchemGeometry> = None;
    let mut comments = Vec::new();
    let mut ended = false;
    for (i, raw) in text.lines().enumerate() {
        let line = i + 1;
        if raw.len() > limits.max_line_len
            || raw.bytes().any(|b| b.is_ascii_control() && b != b'\t')
        {
            return Err(error(line, "line length/control character limit"));
        }
        let (code, comment) = raw
            .split_once('#')
            .map_or((raw, None), |(a, b)| (a, Some(b)));
        if let Some(c) = comment {
            comments.push(c.to_owned());
        }
        let code = code.trim();
        if code.is_empty() {
            continue;
        }
        if code.bytes().any(|b| b";\\\"'".contains(&b)) {
            return Err(error(
                line,
                "semicolon, continuation and quoted/escaped syntax are unsupported",
            ));
        }
        if ended {
            return Err(error(
                line,
                "content after geometry END; full job decks are unsupported",
            ));
        }
        let words: Vec<_> = code.split_whitespace().collect();
        if let Some(g) = &mut geometry {
            if words[0].eq_ignore_ascii_case("end") {
                if words.len() != 1 {
                    return Err(error(line, "unexpected END arguments"));
                }
                ended = true;
            } else {
                if g.centers.len() >= limits.max_centers {
                    return Err(error(line, "center count limit"));
                }
                g.centers.push(center(&words, line)?);
            }
        } else {
            geometry = Some(header(&words, line)?);
        }
    }
    let mut g = geometry.ok_or_else(|| error(0, "missing GEOMETRY block"))?;
    if !ended {
        return Err(error(0, "missing geometry END"));
    }
    g.comments = comments;
    g.validate_with_limits(limits)?;
    Ok(g)
}
fn header(words: &[&str], line: usize) -> Result<NwchemGeometry, NwchemError> {
    if !words[0].eq_ignore_ascii_case("geometry") {
        return Err(error(
            line,
            "expected a GEOMETRY block, not a full job deck",
        ));
    }
    let mut g = NwchemGeometry {
        name: None,
        units: None,
        angstrom_to_au: None,
        center: None,
        autosym: None,
        symmetry_tolerance: None,
        autoz: None,
        print: None,
        bqbq: false,
        nucleus: None,
        centers: Vec::new(),
        comments: Vec::new(),
    };
    let mut pos = 1;
    if let Some(name) = words.get(pos)
        && !OPTIONS.contains(&name.to_ascii_lowercase().as_str())
    {
        if UNSUPPORTED.contains(&name.to_ascii_lowercase().as_str()) {
            return Err(error(line, "unsupported geometry option"));
        }
        g.name = Some((*name).to_owned());
        pos += 1;
    }
    while pos < words.len() {
        option(&mut g, words, &mut pos, line)?;
    }
    Ok(g)
}
fn next<'a>(words: &[&'a str], pos: &mut usize, line: usize) -> Result<&'a str, NwchemError> {
    let word = words
        .get(*pos)
        .ok_or_else(|| error(line, "missing option value"))?;
    *pos += 1;
    Ok(word)
}
fn set<T>(target: &mut Option<T>, value: T, line: usize) -> Result<(), NwchemError> {
    if target.replace(value).is_some() {
        return Err(error(line, "duplicate or conflicting option"));
    }
    Ok(())
}
fn option(
    g: &mut NwchemGeometry,
    w: &[&str],
    p: &mut usize,
    line: usize,
) -> Result<(), NwchemError> {
    let keyword = next(w, p, line)?.to_ascii_lowercase();
    match keyword.as_str() {
        "units" => set(&mut g.units, next(w, p, line)?.to_ascii_lowercase(), line)?,
        "ang2au" | "angstrom_to_au" => set(
            &mut g.angstrom_to_au,
            number(next(w, p, line)?, line)?,
            line,
        )?,
        "center" | "nocenter" => set(&mut g.center, keyword == "center", line)?,
        "autoz" | "noautoz" => set(&mut g.autoz, keyword == "autoz", line)?,
        "autosym" | "noautosym" => {
            set(&mut g.autosym, keyword == "autosym", line)?;
            if keyword == "autosym"
                && let Some(v) = w.get(*p)
                && let Ok(v) = number(v, line)
            {
                g.symmetry_tolerance = Some(v);
                *p += 1;
            }
        }
        "print" => {
            let mode = if w.get(*p).is_some_and(|s| s.eq_ignore_ascii_case("xyz")) {
                *p += 1;
                NwchemPrint::Xyz
            } else {
                NwchemPrint::Normal
            };
            set(&mut g.print, mode, line)?;
        }
        "noprint" => set(&mut g.print, NwchemPrint::None, line)?,
        "bqbq" if !g.bqbq => g.bqbq = true,
        "nuc" | "nucl" | "nucleus" => set(&mut g.nucleus, nucleus(next(w, p, line)?, line)?, line)?,
        _ => return Err(error(line, "unsupported or duplicate geometry option")),
    }
    Ok(())
}
fn center(w: &[&str], line: usize) -> Result<NwchemCenter, NwchemError> {
    if w.len() < 4 {
        return Err(error(line, "expected tag and three Cartesian coordinates"));
    }
    let mut c = NwchemCenter {
        tag: w[0].to_owned(),
        position: vector(&w[1..4], line)?,
        velocity_au: None,
        nuclear_charge: None,
        mass_amu: None,
        nucleus: None,
    };
    let mut p = 4;
    if let Some(word) = w.get(p)
        && !["charge", "mass", "nuc", "nucl", "nucleus"]
            .contains(&word.to_ascii_lowercase().as_str())
    {
        c.velocity_au = Some(vector(
            w.get(p..p + 3)
                .ok_or_else(|| error(line, "incomplete velocity"))?,
            line,
        )?);
        p += 3;
    }
    while p < w.len() {
        let key = next(w, &mut p, line)?.to_ascii_lowercase();
        let value = next(w, &mut p, line)?;
        match key.as_str() {
            "charge" => set(&mut c.nuclear_charge, number(value, line)?, line)?,
            "mass" => set(&mut c.mass_amu, number(value, line)?, line)?,
            "nuc" | "nucl" | "nucleus" => set(&mut c.nucleus, nucleus(value, line)?, line)?,
            _ => return Err(error(line, "unsupported center attribute")),
        }
    }
    Ok(c)
}
fn vector(w: &[&str], line: usize) -> Result<[f64; 3], NwchemError> {
    Ok([
        number(w[0], line)?,
        number(w[1], line)?,
        number(w[2], line)?,
    ])
}
fn number(s: &str, line: usize) -> Result<f64, NwchemError> {
    let v: f64 = s
        .replace(['D', 'd'], "e")
        .parse()
        .map_err(|_| error(line, "invalid real number"))?;
    if v == 0.0
        && s.split(['e', 'E', 'd', 'D'])
            .next()
            .unwrap_or("")
            .bytes()
            .any(|b| matches!(b, b'1'..=b'9'))
    {
        return Err(error(line, "nonzero real number underflows f64"));
    }
    if !v.is_finite() {
        return Err(error(line, "nonfinite real number"));
    }
    Ok(v)
}
fn nucleus(s: &str, line: usize) -> Result<NwchemNucleus, NwchemError> {
    match s.to_ascii_lowercase().as_str() {
        "point" | "pt" => Ok(NwchemNucleus::Point),
        "finite" | "fi" => Ok(NwchemNucleus::Finite),
        _ => Err(error(line, "unknown nuclear model")),
    }
}
