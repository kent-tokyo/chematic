{
  "schema_version": 1,
  "target_version": "1.0.9",
  "format": "mol2",
  "compression": "gzip",
  "fixture": {
    "path": "benchmarks/fixtures/ethanol.mol2",
    "decompressed_bytes": 361,
    "compressed_bytes": 161,
    "sha256": "80144aad97332e2f58eb42bb2da1472587a9dfc2466f9d2fe081eb9834c2224d"
  },
  "repeats": 20,
  "expected_records": 20,
  "rows": {
    "chematic": {
      "format": "mol2",
      "compression": "gzip",
      "execution_mode": "materialized_one_shot",
      "path": "/var/folders/_h/f8_99c4s0vn6f7lb4p3tq6km0000gn/T/chematic-gzip-openbabel-0xxccr1c/fixture.mol2.gz",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 3220,
      "limits": {
        "max_input_bytes": 1073741824,
        "max_record_bytes": 16777216,
        "max_line_bytes": 16777216,
        "max_records": 100000,
        "max_atoms": 1000000
      },
      "seconds": 0.004509,
      "records_per_second": 4435.45,
      "bytes_per_second": 714107.5
    },
    "openbabel": {
      "engine": "openbabel",
      "format": "mol2",
      "repeats": 20,
      "records": 20,
      "failures": 0,
      "input_bytes": 7220,
      "comparison_boundary": "Open Babel CLI over the identical decompressed payload; gzip wrapper excluded"
    }
  },
  "comparison_boundary": {
    "chematic": "Rust gzip-decompressing file-backed MOL2 reader",
    "openbabel": "Open Babel CLI over the identical decompressed payload; process startup included"
  },
  "tool_versions": {
    "openbabel": "Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27"
  }
}
