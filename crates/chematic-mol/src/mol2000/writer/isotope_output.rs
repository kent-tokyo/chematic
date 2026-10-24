//! Absolute isotope labels cannot be represented by a zero mass difference.
use chematic_core::Molecule;
use std::fmt::Write as _;

pub(super) fn write_isotopes(out: &mut String, mol: &Molecule) {
    let entries: Vec<_> = mol
        .atoms()
        .filter_map(|(idx, atom)| {
            atom.isotope
                .filter(|&mass| mass != 0)
                .map(|mass| (idx.0 + 1, mass))
        })
        .collect();
    // M ISO supersedes atom-line differences. Include every labelled atom,
    // even one whose legacy difference is nonzero, in at most eight pairs.
    for chunk in entries.chunks(8) {
        write!(out, "M  ISO{:>3}", chunk.len()).expect("writing to String cannot fail");
        for &(idx, mass) in chunk {
            write!(out, "{idx:>4}{mass:>4}").expect("writing to String cannot fail");
        }
        out.push('\n');
    }
}
