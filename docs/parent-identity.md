# Parent identity reports

Parent operations reduce one source of molecular variation at a time. The
fixed composition order is `fragment → charge → isotope → stereo → tautomer`.

Python callers can inspect the composed result and every intermediate stage:

```python
report = mol.super_parent_report()
assert report["status"] == "Completed"
for stage in report["stages"]:
    print(stage["name"], stage["smiles"])
```

WASM callers use `super_parent_report_json(mol, max_transforms, max_tautomers,
timeout_ms)`. The JSON shape is intentionally identical: `smiles`, `status`,
and five ordered `stages` entries. Validate captured reports with
`validation/validate_parent_report.py` and the checked-in
`validation/parent_identity_schema.json` contract.

Budget-limited tautomer stages retain their explicit status. Callers must not
interpret a non-`Completed` result as a silently canonical identity.
