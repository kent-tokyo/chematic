// Original four-cycle title reuse control; no upstream source copied.
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <iostream>
#include <string>

int main() {
    OpenBabel::OBConversion reader;
    reader.SetInFormat("smi");
    OpenBabel::OBMol reused;
    std::cout << "{\"rows\":[";
    for (int cycle = 0; cycle < 4; ++cycle) {
        if (cycle) std::cout << ",";
        const std::string title = "original title " + std::to_string(cycle);
        if (!reader.ReadString(&reused, "CCO " + title + "\n")) return 2;
        reused.GetFormula(); // Warm derived data before clearing.
        const bool cleared = reused.Clear();
        const auto empty_atoms = reused.NumAtoms();
        const std::string empty_title = reused.GetTitle();
        if (!reader.ReadString(&reused, "CN\n")) return 3;
        OpenBabel::OBMol fresh;
        if (!reader.ReadString(&fresh, "CN\n")) return 4;
        std::cout << "{\"cycle\":" << cycle << ",\"clear_success\":" << cleared
                  << ",\"after_clear_atoms\":" << empty_atoms
                  << ",\"after_clear_title\":\"" << empty_title
                  << "\",\"reused_title\":\"" << reused.GetTitle()
                  << "\",\"fresh_title\":\"" << fresh.GetTitle()
                  << "\",\"reused_formula\":\"" << reused.GetFormula()
                  << "\",\"fresh_formula\":\"" << fresh.GetFormula() << "\"}";
    }
    std::cout << "]}\n";
}
