# Node/WASM contract smoke: v1.0.10 workspace

The checked-in Node/WASM package passed the shared contract, semantic
expansion, and pipeline error-envelope tests:

```sh
node crates/chematic-wasm/tests/cross_binding_contract.test.mjs
node crates/chematic-wasm/tests/semantic_expansion_contract.test.mjs
node crates/chematic-wasm/tests/pipeline_v2.test.mjs
```

The runtime was Node `v24.5.0`. The contract assertions passed, including
malformed-record continuation, input limits, descriptor/fingerprint fixtures,
semantic expansion, and typed pipeline failures.

An explicit version probe reported `chematic_version() == 1.0.9`, while the
workspace is `1.0.10`. Therefore this is contract evidence against the
checked-in generated artifact only; it is not current v1.0.10 browser or Node
artifact evidence. Regenerating a current Web/Node package requires the missing
`wasm-bindgen` CLI and remains an environment-dependent open item.
