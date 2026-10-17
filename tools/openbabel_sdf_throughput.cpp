#include "openbabel_steady_timing.hpp"
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
      if (chematic_benchmark::warmed_dump()) {
        chematic_benchmark::measure(chematic_benchmark::dump_repeats(), [&]() { return write_records(prepared, writer); });
      }
      for (const auto& molecule : prepared) std::cout << writer.WriteString(molecule.get(), false);
      return 0;
    }
    const auto timing = chematic_benchmark::measure(repeats, [&]() -> std::uint64_t {
      if (mode == "parse") {
        const auto molecules = read_records(input, expected);
        std::uint64_t units = 0;
        for (const auto& molecule : molecules) units += molecule->NumAtoms();
        return units;
      }
      if (mode == "write") return write_records(prepared, writer);
      if (mode == "roundtrip") return write_records(read_records(input, expected), writer);
      throw std::runtime_error("unknown mode");
    });
    std::cout << "{\"engine\":\"openbabel\",\"mode\":\"" << mode
              << "\",\"records\":" << expected * timing.repeats << ",\"repeats\":" << timing.repeats
              << ",\"output_units\":" << timing.units << ",\"elapsed_ns\":" << timing.elapsed_ns << ",\"warmup_repeats\":" << timing.warmup_repeats
              << ",\"warmup_elapsed_ns\":" << timing.warmup_elapsed_ns
              << ",\"minimum_warmup_ns\":" << timing.minimum_warmup_ns
              << ",\"minimum_timed_ns\":" << timing.minimum_timed_ns << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
