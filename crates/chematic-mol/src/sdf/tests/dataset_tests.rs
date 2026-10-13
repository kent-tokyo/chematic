use super::*;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_FILE: AtomicU64 = AtomicU64::new(0);

struct TempSdf(PathBuf);

impl TempSdf {
    fn new(contents: &[u8]) -> Self {
        let sequence = NEXT_FILE.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "chematic-sdf-dataset-{}-{sequence}.sdf",
            std::process::id()
        ));
        std::fs::write(&path, contents).unwrap();
        Self(path)
    }
}

impl Drop for TempSdf {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

fn record(title: &str, symbol: &str) -> String {
    format!(
        "{title}\n  chematic\n\n  1  0  0  0  0  0  0  0  0  0  0 V2000\n    0.0000    0.0000    0.0000 {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n>  <ID>\n{title}-id\n\n"
    )
}

#[test]
fn indexes_metadata_and_reads_records_in_any_order() {
    let first = record("carbon", "C");
    let second = record("oxygen", "O");
    let text = format!("{first}$$$$\n{second}$$$$\n");
    let file = TempSdf::new(text.as_bytes());
    let dataset = SdfDataset::open(&file.0).unwrap();

    assert_eq!(dataset.len(), 2);
    assert_eq!(dataset.metadata(0).unwrap().title, "carbon");
    assert_eq!(dataset.metadata(0).unwrap().byte_offset, 0);
    assert_eq!(dataset.metadata(0).unwrap().byte_len, first.len());
    assert_eq!(dataset.metadata(1).unwrap().title, "oxygen");
    assert_eq!(
        dataset.metadata(1).unwrap().byte_offset,
        (first.len() + "$$$$\n".len()) as u64
    );

    let oxygen = dataset.get(1).unwrap().unwrap();
    assert_eq!(oxygen.meta.name, "oxygen");
    assert_eq!(oxygen.properties.get("ID").unwrap(), "oxygen-id");
    assert_eq!(
        oxygen.mol.atoms().next().unwrap().1.element.atomic_number(),
        8
    );
    let carbon = dataset.get(0).unwrap().unwrap();
    assert_eq!(
        carbon.mol.atoms().next().unwrap().1.element.atomic_number(),
        6
    );
    assert!(dataset.get(2).unwrap().is_none());
}

#[test]
fn keeps_malformed_record_index_and_allows_later_access() {
    let valid = record("valid", "N");
    let text = format!("broken\nnot a mol\n$$$$\n{valid}$$$$\n");
    let file = TempSdf::new(text.as_bytes());
    let dataset = SdfDataset::open(&file.0).unwrap();

    assert_eq!(dataset.len(), 2);
    assert!(dataset.get(0).is_err());
    assert_eq!(dataset.get(1).unwrap().unwrap().meta.name, "valid");
}

#[test]
fn supports_crlf_unterminated_last_record_and_skips_empty_records() {
    let record = record("last", "F").replace('\n', "\r\n");
    let text = format!("$$$$\r\n\r\n$$$$\r\n{record}");
    let file = TempSdf::new(text.as_bytes());
    let dataset = SdfDataset::open(&file.0).unwrap();

    assert_eq!(dataset.len(), 1);
    assert_eq!(dataset.metadata(0).unwrap().title, "last");
    assert_eq!(
        dataset
            .get(0)
            .unwrap()
            .unwrap()
            .mol
            .atoms()
            .next()
            .unwrap()
            .1
            .element
            .atomic_number(),
        9
    );
}

#[test]
fn range_reads_are_bounded_and_preserve_parse_failures() {
    let first = record("one", "C");
    let third = record("three", "O");
    let text = format!("{first}$$$$\nbroken\n$$$$\n{third}$$$$\n");
    let file = TempSdf::new(text.as_bytes());
    let dataset = SdfDataset::open(&file.0).unwrap();

    let range = dataset.read_range(1, 99).unwrap();
    assert_eq!(range.len(), 2);
    assert!(range[0].is_err());
    assert_eq!(range[1].as_ref().unwrap().meta.name, "three");
    assert!(dataset.read_range(9, 12).unwrap().is_empty());
}

#[test]
fn indexing_enforces_record_and_count_limits() {
    let text = format!("{}$$$$\n{}$$$$\n", record("a", "C"), record("b", "O"));
    let file = TempSdf::new(text.as_bytes());
    let limits = SdfParseLimits {
        max_records: 1,
        ..SdfParseLimits::default()
    };
    assert!(matches!(
        SdfDataset::open_with_limits(&file.0, limits),
        Err(MolParseError::ResourceLimit {
            resource: "records",
            ..
        })
    ));
}

#[test]
fn empty_dataset_and_raw_record_access_are_explicit() {
    let file = TempSdf::new(b"$$$$\n\r\n$$$$\r\n");
    let dataset = SdfDataset::open(&file.0).unwrap();

    assert!(dataset.is_empty());
    assert_eq!(dataset.len(), 0);
    assert!(dataset.metadata(0).is_none());
    assert!(dataset.record_text(0).unwrap().is_none());
    assert!(dataset.read_range(0, 0).unwrap().is_empty());
}

#[test]
fn opening_and_indexing_report_io_encoding_and_size_failures() {
    let missing = std::env::temp_dir().join(format!(
        "chematic-sdf-dataset-missing-{}-{}.sdf",
        std::process::id(),
        NEXT_FILE.fetch_add(1, Ordering::Relaxed)
    ));
    assert!(matches!(
        SdfDataset::open(missing),
        Err(MolParseError::Io(_))
    ));

    let invalid_utf8 = TempSdf::new(b"bad\xfftitle\n$$$$\n");
    let dataset = SdfDataset::open(&invalid_utf8.0).unwrap();
    assert!(matches!(dataset.record_text(0), Err(MolParseError::Io(_))));

    let text = record("limited", "C");
    let file = TempSdf::new(text.as_bytes());
    for (limits, resource) in [
        (
            SdfParseLimits {
                max_input_bytes: text.len() - 1,
                ..SdfParseLimits::default()
            },
            "input bytes",
        ),
        (
            SdfParseLimits {
                max_record_bytes: text.len() - 1,
                ..SdfParseLimits::default()
            },
            "record bytes",
        ),
        (
            SdfParseLimits {
                max_line_bytes: 2,
                ..SdfParseLimits::default()
            },
            "line bytes",
        ),
    ] {
        assert!(matches!(
            SdfDataset::open_with_limits(&file.0, limits),
            Err(MolParseError::ResourceLimit {
                resource: actual,
                ..
            }) if actual == resource
        ));
    }
}
