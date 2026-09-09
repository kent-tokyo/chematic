{
  "schema_version": 1,
  "target_version": "1.0.9",
  "format": "cml",
  "compression": "gzip",
  "fixture": {
    "path": "benchmarks/fixtures/ethanol.cml",
    "decompressed_bytes": 383,
    "compressed_bytes": 204,
    "sha256": "0d89c3e1250890147debd2ae097b27f8128f73b8d44bb60cbe3279242053e06c"
  },
  "repeats": 20,
  "expected_records": 20,
  "rows": {
    "chematic": {
      "format": "cml",
      "compression": "gzip",
      "execution_mode": "materialized_one_shot",
      "path": "/var/folders/_h/f8_99c4s0vn6f7lb4p3tq6km0000gn/T/chematic-gzip-openbabel-by5ietip/fixture.cml.gz",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 4080,
      "limits": {
        "max_input_bytes": 1073741824,
        "max_record_bytes": 16777216,
        "max_line_bytes": 16777216,
        "max_records": 100000,
        "max_atoms": 1000000
      },
      "seconds": 0.006532,
      "records_per_second": 3061.81,
      "bytes_per_second": 624609.24
    },
    "openbabel": {
      "engine": "openbabel",
      "format": "cml",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 7660,
      "comparison_boundary": "Open Babel CLI over the identical decompressed payload; gzip wrapper excluded"
    }
  },
  "comparison_boundary": {
    "chematic": "Rust gzip-decompressing file-backed CML reader",
    "openbabel": "Open Babel CLI over the identical decompressed payload; process startup included"
  },
  "tool_versions": {
    "openbabel": "Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27"
  }
}
