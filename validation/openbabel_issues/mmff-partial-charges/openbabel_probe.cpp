#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/forcefield.h>
#include <iomanip>
#include <iostream>
#include <memory>

int main() {
    using namespace OpenBabel;
    std::cout << std::setprecision(17) << "[";
    bool first = true;
    for (const char* kind : {"GAFF", "MMFF94"}) {
        for (double scale : {0.0, 1.0, 2.0}) {
            OBMol mol;
            OBConversion conv;
            if (!conv.SetInFormat("smi") || !conv.ReadString(&mol, "CO")) return 2;
            mol.AddHydrogens();
            if (mol.NumAtoms() != 6) return 3;
            const double xyz[6][3] = {{0,0,0},{1.4,0,0},{-.5,.9,0},{-.5,-.45,.8},{-.5,-.45,-.8},{1.9,.8,0}};
            const double charges[6] = {-.3,.1,.05,.05,.05,.05};
            for (unsigned i = 0; i < 6; ++i) {
                mol.GetAtom(i+1)->SetVector(xyz[i][0], xyz[i][1], xyz[i][2]);
                mol.GetAtom(i+1)->SetPartialCharge(charges[i]*scale);
            }
            mol.SetPartialChargesPerceived();
            auto* prototype = OBForceField::FindForceField(kind);
            if (!prototype) return 4;
            std::unique_ptr<OBForceField> field(prototype->MakeNewInstance());
            if (!field->Setup(mol)) return 5;
            if (!first) std::cout << ",";
            first = false;
            std::cout << "{\"field\":\"" << kind << "\",\"supplied_scale\":" << scale
                      << ",\"electrostatic_energy\":" << field->E_Electrostatic(false)
                      << ",\"total_energy\":" << field->Energy(false) << "}";
        }
    }
    std::cout << "]\n";
}
