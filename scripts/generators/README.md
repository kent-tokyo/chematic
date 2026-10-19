# Boundary fixture generators

`generate_*_boundary_fixtures.py` writes the committed RDKit reference
fixtures under `validation/`. Use the exact RDKit version and input noted
in the corresponding validation record. Run a generator from the repository
root, then review the fixture diff before committing it. Moving these
generators did not regenerate or alter any fixture.
