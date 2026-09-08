{
  "schema_version": 1,
  "target_version": "1.0.9",
  "format": "sdf",
  "compression": "gzip",
  "fixture": {
    "path": "benchmarks/fixtures/streaming.sdf",
    "decompressed_bytes": 633,
    "compressed_bytes": 134,
    "sha256": "627e7814694f7d00a9212712ec697f2174d30327b459d983d2a10f717cae0465"
  },
  "repeats": 20,
  "expected_records": 40,
  "rows": {
    "chematic": {
      "format": "sdf",
      "compression": "gzip",
      "execution_mode": "file_backed_bufread",
      "path": "/var/folders/_h/f8_99c4s0vn6f7lb4p3tq6km0000gn/T/chematic-gzip-openbabel-cju2_t1_/fixture.sdf.gz",
      "repeats": 20,
      "records": 40,
      "failures": 0,
      "input_bytes": 2680,
      "limits": {
        "max_input_bytes": 1073741824,
        "max_record_bytes": 16777216,
        "max_line_bytes": 16777216,
        "max_records": 100000,
        "max_atoms": 1000000
      },
      "seconds": 0.004613,
      "records_per_second": 8670.52,
      "bytes_per_second": 580924.9
    },
    "openbabel": {
      "engine": "openbabel",
      "format": "sdf",
      "repeats": 20,
      "records": 40,
      "failures": 0,
      "input_bytes": 12660,
      "comparison_boundary": "Open Babel CLI over the identical decompressed payload; gzip wrapper excluded"
    }
  },
  "comparison_boundary": {
    "chematic": "Rust gzip-decompressing file-backed SDF reader",
    "openbabel": "Open Babel CLI over the identical decompressed payload; process startup included"
  },
  "tool_versions": {
    "openbabel": "Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27"
  }
}
