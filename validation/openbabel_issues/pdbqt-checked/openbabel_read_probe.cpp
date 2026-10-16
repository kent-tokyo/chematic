// Original MIT OR Apache-2.0 diagnostic; uses the public Open Babel API.
// Read one rigid PDBQT block from stdin without automatic bonding. Emits
// element, coordinates and stored partial charge in source atom order.
#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
int main() {
    const std::string text{std::istreambuf_iterator<char>(std::cin), {}};
    OpenBabel::OBConversion conversion;
    if (!conversion.SetInFormat("pdbqt")) return 2;
    conversion.AddOption("b", OpenBabel::OBConversion::INOPTIONS);
    OpenBabel::OBMol mol;
    if (!conversion.ReadString(&mol, text)) return 3;
    std::cout << mol.NumAtoms() << '\n' << std::setprecision(17);
    for (unsigned i = 1; i <= mol.NumAtoms(); ++i) {
        auto* a = mol.GetAtom(i);
        std::cout << a->GetAtomicNum() << ' ' << a->GetX() << ' ' << a->GetY()
                  << ' ' << a->GetZ() << ' ' << a->GetPartialCharge() << '\n';
    }
}
