{
  "schema_version": 1,
  "target_version": "1.0.9",
  "format": "cdxml",
  "compression": "gzip",
  "fixture": {
    "path": "benchmarks/fixtures/ethanol.cdxml",
    "decompressed_bytes": 298,
    "compressed_bytes": 175,
    "sha256": "788130936f9a9cf46bb2c64d54e05bb9da39a707e04199f8ec4ae97416436857"
  },
  "repeats": 20,
  "expected_records": 20,
  "rows": {
    "chematic": {
      "format": "cdxml",
      "compression": "gzip",
      "execution_mode": "materialized_one_shot",
      "path": "/var/folders/_h/f8_99c4s0vn6f7lb4p3tq6km0000gn/T/chematic-gzip-openbabel-7yx3jty8/fixture.cdxml.gz",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 3500,
      "limits": {
        "max_input_bytes": 1073741824,
        "max_record_bytes": 16777216,
        "max_line_bytes": 16777216,
        "max_records": 100000,
        "max_atoms": 1000000
      },
      "seconds": 0.007929,
      "records_per_second": 2522.43,
      "bytes_per_second": 441424.54
    },
    "openbabel": {
      "engine": "openbabel",
      "format": "cdxml",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 5960,
      "comparison_boundary": "Open Babel CLI over the identical decompressed payload; gzip wrapper excluded"
    }
  },
  "comparison_boundary": {
    "chematic": "Rust gzip-decompressing file-backed CDXML reader",
    "openbabel": "Open Babel CLI over the identical decompressed payload; process startup included"
  },
  "tool_versions": {
    "openbabel": "Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27"
  }
}
