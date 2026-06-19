---
name: Bug report
about: Something isn't working correctly
labels: bug
---

## Description

<!-- A clear description of the bug. -->

## Reproduction

```python
# Minimal example that reproduces the bug
import chematic
mol = chematic.from_smiles("...")
```

or in Rust:

```rust
use chematic_smiles::parse;
let mol = parse("...").unwrap();
```

## Expected behavior

<!-- What you expected to happen. -->

## Actual behavior

<!-- What actually happened. Include error messages and stack traces. -->

## Environment

- chematic version: <!-- `pip show chematic` or `cargo metadata | grep chematic` -->
- Python version: <!-- if using Python bindings -->
- OS: <!-- macOS / Linux / Windows -->
