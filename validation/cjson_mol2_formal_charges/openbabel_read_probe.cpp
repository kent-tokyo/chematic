// Original MIT diagnostic: indexed MOL2 fields, no charge calculation.
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
int main() {
    const std::string text{std::istreambuf_iterator<char>(std::cin), {}};
    OpenBabel::OBConversion conversion;
    if (!conversion.SetInFormat("mol2")) return 2;
    OpenBabel::OBMol mol;
    if (!conversion.ReadString(&mol, text)) return 3;
    // Read stored source values, preventing GetPartialCharge's lazy calculation.
    mol.SetPartialChargesPerceived(true);
    std::cout << mol.NumAtoms() << ' ' << mol.NumBonds() << '\n' << std::setprecision(17);
    for (unsigned i=1; i<=mol.NumAtoms(); ++i) {
        const auto* a=mol.GetAtom(i);
        std::cout << a->GetAtomicNum() << ' ' << a->GetFormalCharge() << ' '
                  << mol.GetAtom(i)->GetPartialCharge() << ' ' << a->GetX() << ' '
                  << a->GetY() << ' ' << a->GetZ() << '\n';
    }
    for (unsigned i=0; i<mol.NumBonds(); ++i) {
        const auto* b=mol.GetBond(i);
        std::cout << b->GetBeginAtomIdx()-1 << ' ' << b->GetEndAtomIdx()-1 << ' '
                  << b->GetBondOrder() << '\n';
    }
}
