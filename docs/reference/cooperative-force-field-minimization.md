# Cooperative force-field minimization

Prepared MMFF94, MMFF94s and UFF fields can stop a minimization when a time
budget expires:

```python
import chematic

mol = chematic.from_smiles("CCO").add_hydrogens()
field = mol.mmff_force_field(coords)  # one XYZ coordinate per explicit atom
result = field.minimize_with_timeout_(0.1, max_iters=200)
print(result["status"], result["iterations"], result["energy"])
coordinates = field.positions
```

`timeout_seconds` must be finite and positive. `status` is `converged`,
`iteration_limit`, `line_search_failed`, or `timeout`. `force_field_status`
is the native status (0, 1, or -1) for a completed call and `None` on timeout.
`iterations` counts completed iterations, and `elapsed_seconds` records the
elapsed time within the method. On timeout, coordinates remain at the last
accepted point; `energy` is the energy at that point if it has already been
evaluated, otherwise `None`. Finished calls retain the existing minimizer's
energy semantics. The coordinates are not checked for stereochemical validity
or convergence after a timeout.

The clock starts after the field has been constructed. Stopping is cooperative:
one energy or gradient evaluation, allocation, scheduling, and final cleanup
can exceed the requested budget. This is not a hard deadline and does not
limit molecule parsing, embedding, or force-field construction. The separate
timed conformer-generation API has its own process-level time limit; this
method does not change its behavior.

Rust callers can use `RdkitMmffField::optimize_interruptible` or
`RdkitUffField::optimize_interruptible` with a cancellation predicate, or
`rdkit_bfgs::bfgs_minimize_interruptible` directly. Their interrupted outcomes
report the last accepted coordinates' completed-iteration count and optional
energy. A predicate that never cancels follows the existing minimization path.

A local selected-development-set check is recorded in the
[cooperative minimization benchmark](https://github.com/kent-tokyo/chematic/blob/b7bd7e57db8bc597bd92ba8e713a7e454712bfd9/benchmarks/2026-11-27-cooperative-minimization.md).
