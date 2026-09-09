{
  "schema_version": 1,
  "target_version": "1.0.9",
  "format": "mmcif",
  "compression": "gzip",
  "fixture": {
    "path": "benchmarks/fixtures/minimal.mmcif",
    "decompressed_bytes": 698,
    "compressed_bytes": 265,
    "sha256": "d583d8047758d0e804c98729796d76eea3fa87f3b5b8ed85962717e784e98998"
  },
  "repeats": 20,
  "expected_records": 20,
  "rows": {
    "chematic": {
      "format": "mmcif",
      "compression": "gzip",
      "execution_mode": "materialized_one_shot",
      "path": "/var/folders/_h/f8_99c4s0vn6f7lb4p3tq6km0000gn/T/chematic-gzip-openbabel-dh7q9n29/fixture.mmcif.gz",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 5300,
      "limits": {
        "max_input_bytes": 1073741824,
        "max_record_bytes": 16777216,
        "max_line_bytes": 16777216,
        "max_records": 100000,
        "max_atoms": 1000000
      },
      "seconds": 0.007526,
      "records_per_second": 2657.54,
      "bytes_per_second": 704248.75
    },
    "openbabel": {
      "engine": "openbabel",
      "format": "mmcif",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 13960,
      "comparison_boundary": "Open Babel CLI over the identical decompressed payload; gzip wrapper excluded"
    }
  },
  "comparison_boundary": {
    "chematic": "Rust gzip-decompressing file-backed MMCIF reader",
    "openbabel": "Open Babel CLI over the identical decompressed payload; process startup included"
  },
  "tool_versions": {
    "openbabel": "Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27"
  }
}
