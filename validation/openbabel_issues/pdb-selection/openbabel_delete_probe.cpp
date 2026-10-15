// Original MIT OR Apache-2.0 validation harness; Open Babel is an external oracle.
// This file is never linked into chematic. Uses only public Open Babel APIs.
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/residue.h>
#include <openbabel/obconversion.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static std::string snapshot(OpenBabel::OBMol& mol) {
  std::ostringstream out;
  out << "[";
  for (unsigned i = 1; i <= mol.NumAtoms(); ++i) {
    auto* atom = mol.GetAtom(i);
    auto* res = atom->GetResidue();
    if (i > 1) out << ",";
    out << "{\"id\":" << atom->GetId();
    if (res) {
      out << ",\"atom_name\":" << std::quoted(res->GetAtomID(atom))
          << ",\"residue_name\":" << std::quoted(res->GetName())
          << ",\"chain\":" << std::quoted(std::string(1, res->GetChain()))
          << ",\"residue_number\":" << std::quoted(res->GetNumString());
    }
    out << "}";
  }
  out << "]";
  return out.str();
}

int main(int argc, char** argv) {
  if (argc != 4) return 2;
  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv;
  if (!conv.SetInAndOutFormats("pdb", "pdb") || !conv.ReadFile(&mol, argv[1])) return 3;
  auto before = snapshot(mol);
  std::vector<unsigned> indices;
  std::istringstream input(argv[2]);
  std::string token;
  while (std::getline(input, token, ',')) {
    unsigned index = std::stoul(token);
    if (index >= mol.NumAtoms()) return 4;
    indices.push_back(index);
  }
  std::sort(indices.rbegin(), indices.rend());
  if (std::adjacent_find(indices.begin(), indices.end()) != indices.end()) return 5;
  for (auto index : indices) if (!mol.DeleteAtom(mol.GetAtom(index+1))) return 6;
  auto after = snapshot(mol);
  if (!conv.WriteFile(&mol, argv[3])) return 7;
  std::cout << "{\"version\":" << std::quoted(OpenBabel::OBReleaseVersion())
            << ",\"before\":" << before << ",\"after\":" << after
            << ",\"after_write\":" << snapshot(mol) << "}\n";
}
