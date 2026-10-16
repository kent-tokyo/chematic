use super::*;
use std::io::{self, BufRead, Read};

fn reference_fill_next_block<R: BufRead>(
    reader: &mut SdfFileReader<R>,
) -> Result<bool, MolParseError> {
    reader.block.clear();
    loop {
        let line_start = reader.block.len();
        let bytes = reader
            .reader
            .read_until(b'\n', &mut reader.block)
            .map_err(|error| MolParseError::Io(error.to_string()))?;
        if bytes == 0 {
            reader.done = true;
            return Ok(!reader.block.iter().all(u8::is_ascii_whitespace));
        }
        let line = &reader.block[line_start..];
        if line.len() > reader.limits.max_line_bytes {
            return Err(MolParseError::ResourceLimit {
                resource: "line bytes",
                actual: line.len(),
                limit: reader.limits.max_line_bytes,
            });
        }
        reader.bytes_read = reader.bytes_read.saturating_add(line.len());
        if reader.bytes_read > reader.limits.max_input_bytes {
            return Err(MolParseError::ResourceLimit {
                resource: "input bytes",
                actual: reader.bytes_read,
                limit: reader.limits.max_input_bytes,
            });
        }
        let trimmed = line.strip_suffix(b"\n").unwrap_or(line);
        let trimmed = trimmed.strip_suffix(b"\r").unwrap_or(trimmed);
        if trimmed == b"$$$$" {
            reader.block.truncate(line_start);
            return Ok(true);
        }
    }
}

// The old iterator and independently frozen line reader. The limit validator
// is unchanged; the shared molecule/SD-field parser only accepts a slice now.
fn reference_next<R: BufRead>(
    reader: &mut SdfFileReader<R>,
) -> Option<Result<SdfRecord, MolParseError>> {
    if reader.done {
        return None;
    }
    loop {
        match reference_fill_next_block(reader) {
            Ok(false) => return None,
            Ok(true) if reader.block.iter().all(u8::is_ascii_whitespace) => continue,
            Ok(true) => break,
            Err(error) => {
                reader.done = true;
                return Some(Err(error));
            }
        }
    }
    if let Err(error) = reader.validate_record_limits() {
        reader.done = true;
        return Some(Err(error));
    }
    Some(SdfFileReader::<R>::parse_block(
        &reader.block,
        reader.diagnostics,
    ))
}

#[derive(Clone)]
struct Chunks<'a> {
    input: &'a [u8],
    position: usize,
    chunk: usize,
    fault: Option<(usize, ErrorKind)>,
    interrupted: bool,
}

impl Read for Chunks<'_> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        let available = self.fill_buf()?;
        let n = available.len().min(out.len());
        out[..n].copy_from_slice(&available[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for Chunks<'_> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        let mut end = self
            .input
            .len()
            .min(self.position.saturating_add(self.chunk));
        if let Some((position, kind)) = self.fault {
            if self.position == position && !self.interrupted {
                self.interrupted = kind == ErrorKind::Interrupted;
                return Err(io::Error::new(kind, "injected SDF read fault"));
            }
            if self.position < position {
                end = end.min(position);
            }
        }
        Ok(&self.input[self.position..end])
    }

    fn consume(&mut self, amount: usize) {
        self.position += amount;
        assert!(self.position <= self.input.len());
    }
}

fn assert_record_equal(a: &SdfRecord, b: &SdfRecord) {
    assert_eq!(a.properties, b.properties);
    assert_eq!(format!("{:?}", a.meta), format!("{:?}", b.meta));
    assert_eq!(a.mol.atom_count(), b.mol.atom_count());
    assert_eq!(a.mol.bond_count(), b.mol.bond_count());
    for ((i, a_atom), (j, b_atom)) in a.mol.atoms().zip(b.mol.atoms()) {
        assert_eq!(i, j);
        assert_eq!(format!("{a_atom:?}"), format!("{b_atom:?}"));
        assert_eq!(
            a.mol.neighbors(i).collect::<Vec<_>>(),
            b.mol.neighbors(j).collect::<Vec<_>>()
        );
        assert_eq!(
            a.mol.stereo_neighbor_order(i),
            b.mol.stereo_neighbor_order(j)
        );
        assert_eq!(a.mol.r_group_label(i), b.mol.r_group_label(j));
        assert_eq!(a.mol.atom_tag(i), b.mol.atom_tag(j));
    }
    assert_eq!(
        format!("{:?}", a.mol.stereo_groups()),
        format!("{:?}", b.mol.stereo_groups())
    );
    for ((i, a_bond), (j, b_bond)) in a.mol.bonds().zip(b.mol.bonds()) {
        assert_eq!(i, j);
        assert_eq!(format!("{a_bond:?}"), format!("{b_bond:?}"));
        assert_eq!(a.mol.bond_direction(i), b.mol.bond_direction(j));
        assert_eq!(
            a.mol.bond_direction_anchor(i),
            b.mol.bond_direction_anchor(j)
        );
    }
    let bits = |points: &[(f64, f64)]| {
        points
            .iter()
            .map(|(x, y)| (x.to_bits(), y.to_bits()))
            .collect::<Vec<_>>()
    };
    assert_eq!(bits(&a.coords), bits(&b.coords));
    assert_eq!(format!("{:?}", a.conformer), format!("{:?}", b.conformer));
    assert_eq!(a.stereo_diagnostics, b.stereo_diagnostics);
    assert_eq!(a.ez_diagnostics, b.ez_diagnostics);
    assert_eq!(a.coordinate_dimension, b.coordinate_dimension);
    assert_eq!(a.geometry_rank, b.geometry_rank);
    assert_eq!(a.stereo3d_diagnostics, b.stereo3d_diagnostics);
}

fn compare(
    input: &[u8],
    chunk: usize,
    limits: SdfParseLimits,
    diagnostics: bool,
    fault: Option<(usize, ErrorKind)>,
) -> usize {
    let source = Chunks {
        input,
        position: 0,
        chunk,
        fault,
        interrupted: false,
    };
    let mut candidate =
        SdfFileReader::with_limits_and_diagnostics(source.clone(), limits, diagnostics);
    let mut reference = SdfFileReader::with_limits_and_diagnostics(source, limits, diagnostics);
    reference.block = Vec::with_capacity(2048);
    let mut records = 0;
    loop {
        let a = candidate.next();
        let b = reference_next(&mut reference);
        assert_eq!(
            candidate.reader.position, reference.reader.position,
            "source position"
        );
        assert_eq!(
            candidate.bytes_read, reference.bytes_read,
            "byte accounting"
        );
        assert_eq!(
            candidate.records_read, reference.records_read,
            "record accounting"
        );
        assert_eq!(candidate.done, reference.done, "termination state");
        match (a, b) {
            (None, None) => break,
            (Some(Ok(a)), Some(Ok(b))) => {
                assert_record_equal(&a, &b);
                records += 1;
            }
            (Some(Err(a)), Some(Err(b))) => assert_eq!(a, b),
            _ => panic!("different iterator outcome"),
        }
    }
    assert!(candidate.next().is_none());
    assert!(reference_next(&mut reference).is_none());
    records
}

fn small_inputs() -> Vec<Vec<u8>> {
    let v3000 = include_str!("../../../../benchmarks/fixtures/ethanol.v3000");
    let (mol, meta) = crate::mol3000::parse_mol_v3000(v3000).unwrap();
    let v2000 = crate::mol2000::write_mol(&mol, &meta);
    let full = format!("{v2000}> <DATA>\ninner $$$$ text\n\n$$$$\n{v3000}$$$$\n");
    vec![
        full.as_bytes().to_vec(),
        full.replace('\n', "\r\n").into_bytes(),
        format!("$$$$\n \t\n$$$$\n{v2000}$$$$\n{v3000}").into_bytes(),
        format!("bad\n$$$$\n{v2000}$$$$\n").into_bytes(),
        b"bad\xfftitle\n$$$$\n\n$$$$\r\n$$$$\r\r\n".to_vec(),
        b"$$$$\n$$$$\r\n$$$$".to_vec(),
        Vec::new(),
    ]
}

#[test]
fn boundaries_limits_errors_and_source_position_match_line_reader() {
    let mut cases = 0;
    for input in small_inputs() {
        let sizes = [
            0,
            1,
            4,
            5,
            64,
            input.len().saturating_sub(1),
            input.len(),
            input.len() + 1,
            usize::MAX,
        ];
        let mut limits = vec![SdfParseLimits::default()];
        for cap in sizes {
            for resource in 0..4 {
                let mut limit = SdfParseLimits::default();
                match resource {
                    0 => limit.max_input_bytes = cap,
                    1 => limit.max_line_bytes = cap,
                    2 => limit.max_record_bytes = cap,
                    _ => limit.max_records = cap,
                }
                limits.push(limit);
            }
            limits.push(SdfParseLimits {
                max_input_bytes: cap,
                max_line_bytes: cap,
                max_record_bytes: cap,
                max_records: cap,
            });
        }
        for limit in limits {
            for chunk in [1, 2, 4, 5, 7, 16, 63, 128, usize::MAX] {
                for diagnostics in [false, true] {
                    compare(&input, chunk, limit, diagnostics, None);
                    cases += 1;
                }
            }
        }
    }
    eprintln!("buffer differential boundary cases: {cases}");
}

#[test]
fn io_failures_and_interruption_match_line_reader() {
    let input = &small_inputs()[0];
    let mut cases = 0;
    for position in 0..=input.len() {
        for chunk in [1, 17, usize::MAX] {
            for kind in [ErrorKind::Other, ErrorKind::Interrupted] {
                compare(
                    input,
                    chunk,
                    SdfParseLimits::default(),
                    true,
                    Some((position, kind)),
                );
                cases += 1;
            }
        }
    }
    eprintln!("buffer differential I/O cases: {cases}");
}

#[test]
fn complete_records_borrow_without_filling_owned_block() {
    let input = &small_inputs()[0];
    let mut reader = SdfFileReader::new(io::Cursor::new(input));
    assert!(reader.next().unwrap().is_ok());
    assert!(reader.block.is_empty());
    assert_eq!(reader.block.capacity(), 0);
    assert!(reader.next().unwrap().is_ok());
    assert!(reader.block.is_empty());
    assert_eq!(reader.block.capacity(), 0);
    assert!(reader.next().is_none());
}

#[test]
fn buffer_boundary_and_resource_limit_cases_use_original_fallback() {
    let input = b"record\n$$$$\nremaining\n";
    let defaults = SdfParseLimits::default();
    assert_eq!(record_span(input, defaults, 0, 0), Some((7, 12)));
    assert_eq!(record_span(b"record\n$$$$", defaults, 0, 0), None);
    for limits in [
        SdfParseLimits {
            max_line_bytes: 6,
            ..defaults
        },
        SdfParseLimits {
            max_input_bytes: 11,
            ..defaults
        },
        SdfParseLimits {
            max_record_bytes: 6,
            ..defaults
        },
        SdfParseLimits {
            max_records: 0,
            ..defaults
        },
    ] {
        assert_eq!(record_span(input, limits, 0, 0), None);
    }
    assert_eq!(
        record_span(
            b"$$$$\n",
            SdfParseLimits {
                max_record_bytes: 0,
                max_records: 0,
                ..defaults
            },
            0,
            0
        ),
        Some((0, 5))
    );
}

#[test]
fn pinned_sdf_batches_match_reference_at_buffer_boundaries() {
    use flate2::read::GzDecoder;
    let corpora: [&[u8]; 2] = [
        include_bytes!(
            "../../../../validation/benchmark_corpora/openbabel-sdf-speed-1000-v1.sdf.gz"
        ),
        include_bytes!(
            "../../../../validation/benchmark_corpora/openbabel-sdf-speed-slice5000-1000-v1.sdf.gz"
        ),
    ];
    let mut comparisons = 0;
    for corpus in corpora {
        let mut input = Vec::new();
        GzDecoder::new(corpus).read_to_end(&mut input).unwrap();
        for chunk in [1, 17, 8192, 65536, usize::MAX] {
            for diagnostics in [false, true] {
                comparisons += compare(&input, chunk, SdfParseLimits::default(), diagnostics, None);
            }
        }
    }
    assert_eq!(comparisons, 20_000);
    eprintln!("buffer differential pinned record comparisons: {comparisons}");
}

#[test]
fn actual_bufreader_refills_match_line_reader() {
    use flate2::read::GzDecoder;
    use std::io::{BufReader, Cursor};
    let mut input = Vec::new();
    GzDecoder::new(
        include_bytes!(
            "../../../../validation/benchmark_corpora/openbabel-sdf-speed-1000-v1.sdf.gz"
        )
        .as_slice(),
    )
    .read_to_end(&mut input)
    .unwrap();
    let mut comparisons = 0;
    let mut borrowed = 0;
    for capacity in [1, 17, 256, 8192, 65536] {
        for diagnostics in [false, true] {
            let source = || BufReader::with_capacity(capacity, Cursor::new(&input));
            let mut candidate = SdfFileReader::with_limits_and_diagnostics(
                source(),
                SdfParseLimits::default(),
                diagnostics,
            );
            let mut reference = SdfFileReader::with_limits_and_diagnostics(
                source(),
                SdfParseLimits::default(),
                diagnostics,
            );
            reference.block = Vec::with_capacity(2048);
            loop {
                borrowed += usize::from(
                    record_span(
                        candidate.reader.fill_buf().unwrap(),
                        candidate.limits,
                        candidate.bytes_read,
                        candidate.records_read,
                    )
                    .is_some(),
                );
                let a = candidate.next();
                let b = reference_next(&mut reference);
                let position = |reader: &BufReader<Cursor<&Vec<u8>>>| {
                    reader.get_ref().position() - reader.buffer().len() as u64
                };
                assert_eq!(position(&candidate.reader), position(&reference.reader));
                assert_eq!(candidate.bytes_read, reference.bytes_read);
                assert_eq!(candidate.records_read, reference.records_read);
                assert_eq!(candidate.done, reference.done);
                match (a, b) {
                    (None, None) => break,
                    (Some(Ok(a)), Some(Ok(b))) => {
                        assert_record_equal(&a, &b);
                        comparisons += 1;
                    }
                    _ => panic!("different BufReader outcome"),
                }
            }
        }
    }
    assert_eq!(comparisons, 10_000);
    assert!(borrowed > 0 && borrowed < comparisons);
    eprintln!(
        "actual BufReader record comparisons: {comparisons}; complete buffered spans: {borrowed}"
    );
}
