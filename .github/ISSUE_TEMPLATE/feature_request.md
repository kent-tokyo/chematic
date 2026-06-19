---
name: Feature request
about: Suggest a new feature or improvement
labels: enhancement
---

## Summary

<!-- One sentence: what feature do you want? -->

## Motivation

<!-- Why is this useful? What use case does it solve? -->

## Proposed API

```python
# Python
import chematic
mol = chematic.from_smiles("CCO")
result = mol.new_feature()
```

or in Rust:

```rust
use chematic_chem::new_feature;
let result = new_feature(&mol);
```

## RDKit / OpenBabel equivalent

<!-- If this exists in RDKit or another toolkit, link to it. -->

## Additional context

<!-- Anything else relevant — papers, references, edge cases. -->
