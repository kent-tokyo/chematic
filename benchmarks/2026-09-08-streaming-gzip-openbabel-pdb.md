{
  "schema_version": 1,
  "target_version": "1.0.9",
  "format": "pdb",
  "compression": "gzip",
  "fixture": {
    "path": "benchmarks/fixtures/minimal.pdb",
    "decompressed_bytes": 286,
    "compressed_bytes": 123,
    "sha256": "69d1e5a2d6cbdb39377325cbbf077face011ef0893718e28b42d0563901ad15b"
  },
  "repeats": 20,
  "expected_records": 20,
  "rows": {
    "chematic": {
      "format": "pdb",
      "compression": "gzip",
      "execution_mode": "materialized_one_shot",
      "path": "/var/folders/_h/f8_99c4s0vn6f7lb4p3tq6km0000gn/T/chematic-gzip-openbabel-xqhx4sii/fixture.pdb.gz",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 2460,
      "limits": {
        "max_input_bytes": 1073741824,
        "max_record_bytes": 16777216,
        "max_line_bytes": 16777216,
        "max_records": 100000,
        "max_atoms": 1000000
      },
      "seconds": 0.003798,
      "records_per_second": 5265.24,
      "bytes_per_second": 647624.06
    },
    "openbabel": {
      "engine": "openbabel",
      "format": "pdb",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 5720,
      "comparison_boundary": "Open Babel CLI over the identical decompressed payload; gzip wrapper excluded"
    }
  },
  "comparison_boundary": {
    "chematic": "Rust gzip-decompressing file-backed PDB reader",
    "openbabel": "Open Babel CLI over the identical decompressed payload; process startup included"
  },
  "tool_versions": {
    "openbabel": "Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27"
  }
}
