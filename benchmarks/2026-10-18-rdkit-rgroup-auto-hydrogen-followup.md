# RDKit labelled R-group: implicit H at an automatic site

Date: 2026-10-18. Candidate source builds on `2ab1e555` with this branch's changes. The RDKit comparator is `2026.09.1`. The prior [8-case Python source result](2026-10-18-rdkit-rgroup-hydrogen-followup.json) is retained unchanged.

For core `c1ccncc1[*:1]` and the ordered molecules pyridine, methylpyridine, and chloropyridine, RDKit assigns an automatic `R2` column. Its three R2 groups are H, methyl, and chloro. Previously CheMatic omitted R2 on the first row, so its core graph did not align with the other two rows. The updated Rust implementation fills an implicit H at an automatically labelled site only when another row establishes that column. A lone unsubstituted row still has no R2 column.

The Rust test `automatically_labelled_site_fills_hydrogen_in_established_column` checks the three row graphs and labels against the RDKit outputs; the full `chematic-chem` unit suite passed 981 tests with 1 ignored. The Python source gate has been extended from 8 to 9 cases, but **9/9 Python-extension agreement is not yet verified locally**. Source wheel builds with Rust 1.97 and 1.88 both completed on this host, but macOS refused to load either extension with a misaligned Mach-O LINKEDIT string pool. The local test therefore cannot establish wheel behavior; the 9-case Python gate should be run on a host that can load the built extension.

This is a bounded R-group compatibility correction. It does not establish general R-group accuracy superiority over RDKit.
