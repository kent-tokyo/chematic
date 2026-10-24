# Extended XYZ writer integrity packet

Nine original positive frames and 19 refusal inputs. Python/WASM receive the
structured data; the comparator receives independently authored text. ASE 3.26.0
and a separate raw-geometry observer assess the result, with fixture assertions
for source bits, cell convention, atom-data shapes and metadata.

Run `scripts/validation/check_extxyz_output_integrity.py` with a fresh Python
extension and Node-target WASM binary. `--baseline` records earlier bindings
without enforcing candidate gates. Dated gzip files preserve complete outputs,
errors, commands, hashes and reader/version identities. Pure observer controls
are in `scripts/tests/test_extxyz_output_integrity.py`.

This measures an ASE-consumer workflow. Open Babel is invoked through plain XYZ;
its separate EXYZ syntax is not treated as ASE ExtXYZ. The nine frames do not
establish chemical graph/stereo preservation, performance or overall superiority.
See [contract and evidence](../../docs/reference/extxyz-output-integrity.md).
