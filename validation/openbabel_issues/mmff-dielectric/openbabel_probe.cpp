// Original two-ion reproduction of Open Babel #2811; no issue attachment needed.
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <iomanip>
#include <iostream>
#include <memory>

int main() {
  OpenBabel::OBConversion conversion;
  conversion.SetInFormat("smi");
  OpenBabel::OBMol molecule;
  if (!conversion.ReadString(&molecule, "[Na+].[Cl-]")) return 1;
  molecule.GetAtom(1)->SetVector(0, 0, 0);
  molecule.GetAtom(2)->SetVector(4, 0, 0);
  auto* field = OpenBabel::OBForceField::FindForceField("MMFF94");
  if (!field) return 2;
  std::cout << std::setprecision(17) << '[';
  bool first = true;
  auto row = [&](const char* mode, double epsilon, OpenBabel::OBForceField* f) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << "{\"mode\":\"" << mode << "\",\"epsilon\":" << epsilon
      << ",\"reported_epsilon\":" << f->GetDielectricConstant()
      << ",\"energy\":" << f->Energy()
      << ",\"electrostatic\":" << f->E_Electrostatic() << '}';
  };
  field->SetDielectricConstant(1);
  if (!field->Setup(molecule)) return 3;
  for (double e : {1., 4., 78., 1.}) {
    field->SetDielectricConstant(e);
    row("live", e, field);
  }
  for (double e : {1., 4., 78., 1.}) {
    field->SetDielectricConstant(e);
    if (!field->Setup(molecule)) return 4;
    row("setup_again", e, field);
  }
  for (double e : {1., 4., 78., 1.}) {
    std::unique_ptr<OpenBabel::OBForceField> fresh(field->MakeNewInstance());
    fresh->SetDielectricConstant(e);
    if (!fresh->Setup(molecule)) return 5;
    row("fresh", e, fresh.get());
  }
  std::cout << "]\n";
}
