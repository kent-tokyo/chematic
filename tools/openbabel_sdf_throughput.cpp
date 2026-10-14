#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Molecules = std::vector<std::unique_ptr<OpenBabel::OBMol>>;

// Match Rust's borrowed Cursor: constructing an input stream does not copy
// the entire SDF buffer on every parse. Only the get area is exposed.
class InputBuffer : public std::streambuf {
 public:
  explicit InputBuffer(const std::string& input) {
    auto first = const_cast<char*>(input.data());
    setg(first, first, first + input.size());
  }
};

Molecules read_records(const std::string& input, std::size_t expected) {
  InputBuffer buffer(input);
  std::istream stream(&buffer);
  OpenBabel::OBConversion conversion(&stream, nullptr);
  if (!conversion.SetInFormat("sdf")) throw std::runtime_error("SDF reader unavailable");
  Molecules molecules;
  molecules.reserve(expected);
  while (true) {
    auto molecule = std::make_unique<OpenBabel::OBMol>();
    if (!conversion.Read(molecule.get())) break;
    molecules.push_back(std::move(molecule));
  }
  if (molecules.size() != expected) throw std::runtime_error("SDF record count differs");
  return molecules;
}

std::size_t write_records(const Molecules& molecules, OpenBabel::OBConversion& writer) {
  std::size_t bytes = 0;
  for (const auto& molecule : molecules) {
    const auto text = writer.WriteString(molecule.get(), false);
    if (text.empty()) throw std::runtime_error("SDF write failed");
    bytes += text.size();
  }
  return bytes;
}

int main(int argc, char** argv) {
  try {
    if (argc != 5) throw std::runtime_error("usage: PATH MODE REPEATS EXPECTED_RECORDS");
    std::ifstream file(argv[1], std::ios::binary | std::ios::ate);
    if (!file) throw std::runtime_error("input unavailable");
    const auto size = file.tellg();
    if (size < 0) throw std::runtime_error("input length unavailable");
    std::string input(static_cast<std::size_t>(size), '\0');
    file.seekg(0);
    if (!file.read(input.data(), input.size())) throw std::runtime_error("input read failed");
    const std::string mode(argv[2]);
    const auto repeats = std::stoull(argv[3]);
    const auto expected = std::stoull(argv[4]);
    if (repeats == 0 || expected == 0) throw std::runtime_error("positive counts required");
    auto prepared = read_records(input, expected);
    if (mode == "smiles") {
      OpenBabel::OBConversion observer;
      if (!observer.SetOutFormat("can")) throw std::runtime_error("canonical SMILES unavailable");
      for (const auto& molecule : prepared) std::cout << observer.WriteString(molecule.get(), false);
      return 0;
    }
    OpenBabel::OBConversion writer;
    if (!writer.SetOutFormat("sdf")) throw std::runtime_error("SDF writer unavailable");
    write_records(prepared, writer);
    if (mode == "dump") {
      for (const auto& molecule : prepared) std::cout << writer.WriteString(molecule.get(), false);
      return 0;
    }
    std::size_t units = 0;
    const auto started = std::chrono::steady_clock::now();
    for (std::uint64_t iteration = 0; iteration < repeats; ++iteration) {
      if (mode == "parse") {
        const auto molecules = read_records(input, expected);
        for (const auto& molecule : molecules) units += molecule->NumAtoms();
      } else if (mode == "write") {
        units += write_records(prepared, writer);
      } else if (mode == "roundtrip") {
        units += write_records(read_records(input, expected), writer);
      } else {
        throw std::runtime_error("unknown mode");
      }
    }
    const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - started).count();
    std::cout << "{\"engine\":\"openbabel\",\"mode\":\"" << mode
              << "\",\"records\":" << expected * repeats << ",\"repeats\":" << repeats
              << ",\"output_units\":" << units << ",\"elapsed_ns\":" << elapsed << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
