//! Read a single archive restart deck without interpreting its calculation results.
use super::{Element, MopacError, MopacOptimizationFlag, MopacParseLimits, error, reader};

/// One restart row. `references == None` means Cartesian values in Å.
/// With references present, values are MOPAC distance/angle/dihedral fields
/// (Å/degrees/degrees), including unused fields on the first three atoms.
/// Reference numbers are retained one-based; zero means an unused reference.
#[derive(Debug, Clone, PartialEq)]
pub struct MopacArchiveAtom {
    pub element: Element,
    pub label: Option<String>,
    pub values: [f64; 3],
    pub optimization_flags: [MopacOptimizationFlag; 3],
    pub references: Option<[usize; 3]>,
    pub partial_charge: Option<f64>,
}

/// A single MOPAC archive. Summary and ancillary restart sections remain raw;
/// reading this type does not establish convergence, topology or bond orders.
#[derive(Debug, Clone, PartialEq)]
pub struct MopacArchive {
    pub keywords: String,
    pub title: String,
    pub comment: String,
    pub atoms: Vec<MopacArchiveAtom>,
    /// Exact bytes before the FINAL GEOMETRY OBTAINED marker line.
    pub summary_text: String,
    /// Exact restart input after the marker line, including ancillary sections.
    pub restart_text: String,
    /// Exact bytes after the atom block (e.g. SYMMETRY constraints).
    pub trailing_text: String,
    pub source_text: String,
}

pub fn parse_mopac_archive(text: &str) -> Result<MopacArchive, MopacError> {
    parse_mopac_archive_with_limits(text, &MopacParseLimits::default())
}

/// Accept exactly one final-geometry marker and a three-line restart header.
/// Flagged Cartesian and numeric internal-coordinate rows are supported,
/// optionally with a CHARGE column announced on the marker line. Periodic Tv,
/// dummy/sparkle atoms, continuation headers and atomic units are refused.
pub fn parse_mopac_archive_with_limits(
    text: &str,
    limits: &MopacParseLimits,
) -> Result<MopacArchive, MopacError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "archive input byte limit exceeded"));
    }
    let mut offset = 0;
    let mut marker = None;
    for (i, raw) in text.split_inclusive('\n').enumerate() {
        let line = raw.trim_end_matches(['\r', '\n']);
        if line.len() > limits.max_line_len {
            return Err(error(i + 1, "archive line length limit exceeded"));
        }
        if line.chars().any(|c| c.is_control() && c != '\t') {
            return Err(error(i + 1, "unsupported control character"));
        }
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.starts_with(&["FINAL", "GEOMETRY", "OBTAINED"]) {
            if marker.is_some() {
                return Err(error(i + 1, "multiple final geometries are ambiguous"));
            }
            let charge = match fields.as_slice() {
                ["FINAL", "GEOMETRY", "OBTAINED"] => false,
                ["FINAL", "GEOMETRY", "OBTAINED", "CHARGE"] => true,
                _ => return Err(error(i + 1, "unsupported final geometry marker")),
            };
            marker = Some((offset, offset + raw.len(), i + 1, charge));
        }
        offset += raw.len();
    }
    let (start, end, line, charge) =
        marker.ok_or_else(|| error(0, "missing final geometry marker"))?;
    parse_restart(text, start, end, line, charge, limits)
}

fn parse_restart(
    text: &str,
    start: usize,
    end: usize,
    line: usize,
    charge: bool,
    limits: &MopacParseLimits,
) -> Result<MopacArchive, MopacError> {
    let restart = &text[end..];
    let mut lines = restart.split_inclusive('\n');
    let mut header = Vec::new();
    let mut consumed = 0;
    for i in 1..=3 {
        let raw = lines
            .next()
            .ok_or_else(|| error(line + i, "three restart header lines required"))?;
        consumed += raw.len();
        header.push(raw.trim_end_matches(['\r', '\n']).to_string());
    }
    for token in header[0].split_whitespace() {
        let key = token.split('=').next().unwrap_or("").to_ascii_uppercase();
        if matches!(key.as_str(), "A0" | "+" | "++" | "&") {
            return Err(error(
                line + 1,
                "atomic units and keyword continuations unsupported",
            ));
        }
    }
    let mut atoms = Vec::new();
    for (i, raw) in lines.enumerate() {
        if raw.trim().is_empty() {
            break;
        }
        if atoms.len() >= limits.max_atoms.min(MopacParseLimits::default().max_atoms) {
            return Err(error(line + i + 4, "archive atom limit exceeded"));
        }
        atoms.push(parse_atom(raw, line + i + 4, charge)?);
        consumed += raw.len();
    }
    if atoms.is_empty() {
        return Err(error(line + 4, "archive requires at least one atom"));
    }
    for (i, atom) in atoms.iter().enumerate() {
        if let Some(refs) = atom.references {
            let mut used = Vec::new();
            for reference in refs.into_iter().filter(|r| *r != 0) {
                if reference > atoms.len() || reference == i + 1 || used.contains(&reference) {
                    return Err(error(line + i + 4, "invalid internal coordinate reference"));
                }
                used.push(reference);
            }
        }
    }
    Ok(MopacArchive {
        keywords: header.remove(0),
        title: header.remove(0),
        comment: header.remove(0),
        atoms,
        summary_text: text[..start].into(),
        restart_text: restart.into(),
        trailing_text: restart[consumed..].into(),
        source_text: text.into(),
    })
}

fn parse_atom(raw: &str, line: usize, charge: bool) -> Result<MopacArchiveAtom, MopacError> {
    let (symbol, label, body) = reader::identity(raw, line)?;
    let element = reader::element(symbol)
        .filter(|e| !matches!(e.atomic_number(), 85 | 87 | 99))
        .ok_or_else(|| {
            error(
                line,
                "unsupported archive element, dummy/sparkle or translation vector",
            )
        })?;
    if let Some(label) = &label {
        super::validate_label(label).map_err(|e| error(line, e.detail))?;
    }
    let fields: Vec<_> = body.split_whitespace().collect();
    let n = fields.len().saturating_sub(usize::from(charge));
    if !matches!(n, 6 | 9) {
        return Err(error(
            line,
            "expected three flagged values, optional three references, and announced charge",
        ));
    }
    let row = reader::parse_row(&fields[..6].join(" "), line)?;
    let references = if n == 9 {
        let parse = |s: &str| {
            s.parse::<usize>()
                .map_err(|_| error(line, "invalid reference number"))
        };
        Some([parse(fields[6])?, parse(fields[7])?, parse(fields[8])?])
    } else {
        None
    };
    Ok(MopacArchiveAtom {
        element,
        label,
        values: row.position_angstrom,
        optimization_flags: row
            .optimization_flags
            .expect("six fields always contain flags"),
        references,
        partial_charge: if charge {
            Some(reader::number(fields[n], line)?)
        } else {
            None
        },
    })
}
