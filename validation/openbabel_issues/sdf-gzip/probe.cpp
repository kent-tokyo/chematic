// Original stream-replacement control for Open Babel issue #3030.
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <fstream>
#include <iostream>
#include <sstream>

static void report(const char *label, OpenBabel::OBConversion &conversion,
                   bool ok, const OpenBabel::OBMol &mol) {
  std::cout << label << "\tok=" << ok << "\tatoms=" << mol.NumAtoms()
            << "\tgzip=" << conversion.GetInGzipped() << "\n";
}

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  OpenBabel::OBConversion fresh;
  fresh.SetInFormat("smi");
  OpenBabel::OBMol control;
  bool ok = fresh.ReadString(&control, "CCO\n");
  report("fresh_plain_string", fresh, ok, control);

  OpenBabel::OBConversion shared;
  shared.SetInFormat("smi");
  std::ifstream compressed(argv[1], std::ios::binary);
  if (!compressed) return 3;
  OpenBabel::OBMol first;
  ok = shared.Read(&first, &compressed);
  report("gzip_input", shared, ok, first);
  OpenBabel::OBMol string_mol;
  ok = shared.ReadString(&string_mol, "CCO\n");
  report("replacement_plain_string", shared, ok, string_mol);
  std::istringstream plain("CN\n");
  OpenBabel::OBMol stream_mol;
  ok = shared.Read(&stream_mol, &plain);
  report("replacement_plain_stream", shared, ok, stream_mol);

  OpenBabel::OBConversion another;
  another.SetInFormat("smi");
  std::istringstream another_plain("CN\n");
  OpenBabel::OBMol another_mol;
  ok = another.Read(&another_mol, &another_plain);
  report("fresh_plain_stream", another, ok, another_mol);
  return 0;
}
