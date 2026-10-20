# Shared script helpers

This directory contains modules imported by maintenance scripts. They are not
standalone commands. `rdkit_standardization_classification.py` supplies the
same residual taxonomy to the source-built and Python cleanup comparisons.
`same_process_report.py` builds the common provenance envelope and writes
format-contract reports; each caller retains its own chemistry checks.
