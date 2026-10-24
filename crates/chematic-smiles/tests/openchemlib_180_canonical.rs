//! OpenChemLib #180 cycles among these spellings while its IDCode stays fixed.
//! A canonical SMILES writer must converge on one spelling and reach a fixed
//! point after reparsing.

use chematic_smiles::{canonical_smiles, parse};

#[test]
fn methylquinoline_spellings_converge_to_one_canonical_fixed_point() {
    let spellings = ["Cc(cc1)cc2c1nccc2", "Cc1ccc2ncccc2c1", "Cc1cc2cccnc2cc1"];
    let canonical = canonical_smiles(&parse(spellings[0]).unwrap());
    for spelling in spellings {
        let output = canonical_smiles(&parse(spelling).unwrap());
        assert_eq!(output, canonical, "{spelling}");
        assert_eq!(canonical_smiles(&parse(&output).unwrap()), output);
    }
}
