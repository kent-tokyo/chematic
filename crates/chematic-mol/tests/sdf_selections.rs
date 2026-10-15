//! Original MIT OR Apache-2.0 fixtures; provenance and oracle in manifest.json.
use chematic_mol::{SdfDataset, SdfSelectionError, SdfSelectionLimits};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
const INPUT: &str =
    include_str!("../../../validation/cosmolkit_comparison/sdf_selections/minimal.sdf");
static NEXT: AtomicU64 = AtomicU64::new(0);
struct Temp(PathBuf);
impl Temp {
    fn new(text: &str) -> Self {
        let path = std::env::temp_dir().join(format!(
            "chematic-selected-sdf-{}-{}.sdf",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, text).unwrap();
        Self(path)
    }
    fn dataset(&self) -> SdfDataset {
        SdfDataset::open(&self.0).unwrap()
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}
#[test]
fn ordered_duplicate_records_and_errors_keep_their_positions() {
    let f = Temp::new(INPUT);
    let d = f.dataset();
    let b = d.read_indices(&[4, 0, 1, 4, 2]).unwrap();
    assert_eq!(b.indices, [4, 0, 1, 4, 2]);
    assert_eq!(b.valid_mask(), [true, true, false, true, true]);
    assert_eq!(b.records[0].as_ref().unwrap().meta.name, "fluorine");
    assert_eq!(b.records[4].as_ref().unwrap().properties["ID"], "oxygen-id");
    assert_eq!(
        b.records[4].as_ref().unwrap().properties["note"],
        "line1\nline2"
    );
    let mut records = b.records;
    records[0]
        .as_mut()
        .unwrap()
        .properties
        .insert("ID".into(), "changed".into());
    assert_eq!(records[3].as_ref().unwrap().properties["ID"], "fluorine-id");
    assert_eq!(d.get(4).unwrap().unwrap().properties["ID"], "fluorine-id");
}
#[test]
fn empty_selection_and_all_bad_selection() {
    let f = Temp::new(INPUT);
    let d = f.dataset();
    assert!(d.read_indices(&[]).unwrap().is_empty());
    let b = d.read_indices(&[1, 1]).unwrap();
    assert_eq!(b.valid_mask(), [false, false]);
    assert_eq!(b.len(), 2);
}
#[test]
fn index_validation_happens_before_reading_any_record() {
    let f = Temp::new(INPUT);
    let d = f.dataset();
    std::fs::write(&f.0, "").unwrap();
    assert!(matches!(
        d.read_indices(&[0, 5]),
        Err(SdfSelectionError::Index { index: 5, len: 5 })
    ));
}
#[test]
fn record_and_cumulative_byte_caps_include_duplicate_reads() {
    let f = Temp::new(INPUT);
    let d = f.dataset();
    let bytes = d.metadata(0).unwrap().byte_len;
    let mut limits = SdfSelectionLimits {
        max_records: 2,
        max_total_record_bytes: 2 * bytes,
    };
    assert!(d.read_indices_with_limits(&[0, 0], &limits).is_ok());
    limits.max_total_record_bytes -= 1;
    assert!(matches!(
        d.read_indices_with_limits(&[0, 0], &limits),
        Err(SdfSelectionError::ResourceLimit {
            resource: "record bytes",
            ..
        })
    ));
    limits.max_total_record_bytes = usize::MAX;
    limits.max_records = 1;
    assert!(matches!(
        d.read_indices_with_limits(&[0, 0], &limits),
        Err(SdfSelectionError::ResourceLimit {
            resource: "records",
            ..
        })
    ));
}
#[test]
fn byte_validation_precedes_reads_and_io_errors_abort() {
    let f = Temp::new(INPUT);
    let d = f.dataset();
    std::fs::write(&f.0, "").unwrap();
    let limits = SdfSelectionLimits {
        max_records: 2,
        max_total_record_bytes: 0,
    };
    assert!(matches!(
        d.read_indices_with_limits(&[0], &limits),
        Err(SdfSelectionError::ResourceLimit { .. })
    ));
    assert!(matches!(
        d.read_indices(&[0]),
        Err(SdfSelectionError::Read(chematic_mol::MolParseError::Io(_)))
    ));
}
#[test]
fn selected_reads_match_scalar_reads_after_crlf_and_unterminated_input() {
    for text in [
        INPUT.to_string(),
        INPUT.replace('\n', "\r\n"),
        INPUT.trim_end_matches("$$$$\n").to_string(),
    ] {
        let f = Temp::new(&text);
        let d = f.dataset();
        assert_eq!(d.len(), 5);
        let b = d.read_indices(&[3, 4, 2, 0]).unwrap();
        for (index, record) in b.indices.iter().zip(b.records) {
            let scalar = d.get(*index).unwrap().unwrap();
            let selected = record.unwrap();
            assert_eq!(selected.meta.name, scalar.meta.name);
            assert_eq!(selected.properties, scalar.properties);
        }
    }
}
#[test]
fn concurrent_selections_do_not_share_seek_positions() {
    let f = Temp::new(INPUT);
    let d = std::sync::Arc::new(f.dataset());
    let threads: Vec<_> = (0..4)
        .map(|_| {
            let d = d.clone();
            std::thread::spawn(move || {
                for _ in 0..20 {
                    let b = d.read_indices(&[4, 0, 3, 2]).unwrap();
                    assert_eq!(b.records[0].as_ref().unwrap().meta.name, "fluorine");
                    assert_eq!(b.records[3].as_ref().unwrap().meta.name, "oxygen");
                }
            })
        })
        .collect();
    for t in threads {
        t.join().unwrap();
    }
}

#[test]
fn all_sdf_record_readers_dispatch_mixed_ctab_versions() {
    let graph: Vec<_> = chematic_mol::SdfReader::new(INPUT).collect();
    let records: Vec<_> = chematic_mol::SdfRecordReader::new(INPUT).collect();
    let stream: Vec<_> = chematic_mol::SdfFileReader::new(std::io::Cursor::new(INPUT)).collect();
    let fast: Vec<_> = chematic_mol::SdfFileReader::fast(std::io::Cursor::new(INPUT)).collect();
    for valid in [
        graph.iter().map(Result::is_ok).collect::<Vec<_>>(),
        records.iter().map(Result::is_ok).collect(),
        stream.iter().map(Result::is_ok).collect(),
        fast.iter().map(Result::is_ok).collect(),
    ] {
        assert_eq!(valid, [true, false, true, true, true]);
    }
    assert_eq!(records[4].as_ref().unwrap().properties["ID"], "fluorine-id");
}

#[test]
fn exact_m_end_title_or_comment_is_not_a_ctab_terminator() {
    for title in ["carbon", "fluorine"] {
        let text = INPUT.replacen(
            &format!("{title}\n  chematic synthetic\n\n"),
            "M  END\n  chematic synthetic\nM  END\n",
            1,
        );
        let f = Temp::new(&text);
        let d = f.dataset();
        let index = if title == "carbon" { 0 } else { 4 };
        let b = d.read_indices(&[index]).unwrap();
        let record = b.records[0].as_ref().unwrap();
        assert_eq!(record.meta.name, "M  END");
        assert_eq!(record.properties["ID"], format!("{title}-id"));
    }
}

#[test]
fn v2000_permissive_terminator_still_keeps_data_fields() {
    let text = INPUT.replacen("M  END\n", "  M  END trailing\n", 1);
    let f = Temp::new(&text);
    let d = f.dataset();
    let b = d.read_indices(&[0]).unwrap();
    assert_eq!(b.records[0].as_ref().unwrap().properties["ID"], "carbon-id");
}
