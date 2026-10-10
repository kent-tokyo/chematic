#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

std::string value_after(int argc, char** argv, const std::string& name,
                        const std::string& fallback = "") {
  for (int index = 1; index + 1 < argc; ++index) {
    if (argv[index] == name) return argv[index + 1];
  }
  return fallback;
}

std::string read_file(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot open " + path);
  std::ostringstream contents;
  contents << input.rdbuf();
  return contents.str();
}

std::string openbabel_format(const std::string& format) {
  return format == "v3000" ? "mol" : format;
}

void configure(OpenBabel::OBConversion& conversion, const std::string& format) {
  const auto identifier = openbabel_format(format);
  if (!conversion.SetInAndOutFormats(identifier.c_str(), identifier.c_str())) {
    throw std::runtime_error("Open Babel does not support format " + format);
  }
  if (format == "v3000") {
    conversion.AddOption("3", OpenBabel::OBConversion::OUTOPTIONS);
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto format = value_after(argc, argv, "--format");
    const auto path = value_after(argc, argv, "--path");
    const auto operation = value_after(argc, argv, "--benchmark-operation");
    const auto repeats = std::stoull(value_after(argc, argv, "--repeats", "1000"));
    if (format.empty() || path.empty() || operation.empty() || repeats == 0) {
      throw std::runtime_error(
          "required: --format, --path, --benchmark-operation, positive --repeats");
    }
    const auto input = read_file(path);
    OpenBabel::OBConversion conversion;
    configure(conversion, format);
    OpenBabel::OBMol prepared;
    if (operation == "write" && !conversion.ReadString(&prepared, input)) {
      throw std::runtime_error("Open Babel failed to parse write fixture");
    }

    std::uint64_t output_units = 0;
    std::uint64_t records = 0;
    const auto started = std::chrono::steady_clock::now();
    for (std::uint64_t iteration = 0; iteration < repeats; ++iteration) {
      if (operation == "parse") {
        OpenBabel::OBMol molecule;
        if (!conversion.ReadString(&molecule, input)) {
          throw std::runtime_error("Open Babel parse failed");
        }
        output_units += molecule.NumAtoms();
        ++records;
      } else if (operation == "write") {
        const auto output = conversion.WriteString(&prepared, true);
        if (output.empty()) throw std::runtime_error("Open Babel write failed");
        output_units += output.size();
        ++records;
      } else if (operation == "roundtrip") {
        OpenBabel::OBMol molecule;
        if (!conversion.ReadString(&molecule, input)) {
          throw std::runtime_error("Open Babel round-trip parse failed");
        }
        const auto output = conversion.WriteString(&molecule, true);
        if (output.empty()) throw std::runtime_error("Open Babel round-trip write failed");
        output_units += output.size();
        ++records;
      } else {
        throw std::runtime_error("unsupported operation " + operation);
      }
    }
    const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
                             std::chrono::steady_clock::now() - started)
                             .count();
    std::cout << "{\"engine\":\"openbabel\",\"format\":\"" << format
              << "\",\"operation\":\"" << operation << "\",\"repeats\":"
              << repeats << ",\"records\":" << records << ",\"output_units\":"
              << output_units << ",\"elapsed_ns\":" << elapsed << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
