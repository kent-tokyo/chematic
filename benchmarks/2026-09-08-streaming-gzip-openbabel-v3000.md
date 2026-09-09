{
  "schema_version": 1,
  "target_version": "1.0.9",
  "format": "v3000",
  "compression": "gzip",
  "fixture": {
    "path": "benchmarks/fixtures/ethanol.v3000",
    "decompressed_bytes": 298,
    "compressed_bytes": 150,
    "sha256": "64ac11b3de1c186deae555d54d6a47e7f3c72607ef0a90b0d93f467496bbd3ad"
  },
  "repeats": 20,
  "expected_records": 20,
  "rows": {
    "chematic": {
      "format": "v3000",
      "compression": "gzip",
      "execution_mode": "materialized_one_shot",
      "path": "/var/folders/_h/f8_99c4s0vn6f7lb4p3tq6km0000gn/T/chematic-gzip-openbabel-vlju5g_w/fixture.v3000.gz",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 3000,
      "limits": {
        "max_input_bytes": 1073741824,
        "max_record_bytes": 16777216,
        "max_line_bytes": 16777216,
        "max_records": 100000,
        "max_atoms": 1000000
      },
      "seconds": 0.005008,
      "records_per_second": 3993.41,
      "bytes_per_second": 599011.63
    },
    "openbabel": {
      "engine": "openbabel",
      "format": "mol",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 5960,
      "comparison_boundary": "Open Babel CLI over the identical decompressed payload; gzip wrapper excluded"
    }
  },
  "comparison_boundary": {
    "chematic": "Rust gzip-decompressing file-backed V3000 reader",
    "openbabel": "Open Babel CLI over the identical decompressed payload; process startup included"
  },
  "tool_versions": {
    "openbabel": "Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27"
  }
}
