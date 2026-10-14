#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

using Molecules = std::vector<std::unique_ptr<OpenBabel::OBMol>>;

class InputBuffer : public std::streambuf {
 public:
  explicit InputBuffer(std::string_view input) {
    auto first = const_cast<char*>(input.data());
    setg(first, first, first + input.size());
  }
};

std::vector<std::string_view> split_blocks(const std::string& input, std::size_t expected) {
  std::vector<std::string_view> blocks;
  std::size_t start = 0, position = 0;
  while (position < input.size()) {
    auto end = input.find('\n', position);
    if (end == std::string::npos) end = input.size();
    auto line = std::string_view(input).substr(position, end - position);
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    const auto next = end < input.size() ? end + 1 : end;
    if (line == "$$$$") {
      blocks.emplace_back(input.data() + start, position - start);
      start = next;
    }
    position = next;
  }
  if (start != input.size() || blocks.size() != expected) {
    throw std::runtime_error("V3000 framed record count differs");
  }
  return blocks;
}

Molecules read(const std::vector<std::string_view>& blocks) {
  OpenBabel::OBConversion reader;
  if (!reader.SetInFormat("mol")) throw std::runtime_error("MOL reader unavailable");
  Molecules result;
  result.reserve(blocks.size());
  for (auto block : blocks) {
    InputBuffer buffer(block);
    std::istream stream(&buffer);
    auto mol = std::make_unique<OpenBabel::OBMol>();
    if (!reader.Read(mol.get(), &stream)) throw std::runtime_error("MOL parse failed");
    result.push_back(std::move(mol));
  }
  return result;
}

std::size_t write(const Molecules& molecules, OpenBabel::OBConversion& writer) {
  std::size_t bytes = 0;
  for (const auto& mol : molecules) {
    const auto output = writer.WriteString(mol.get(), false);
    if (output.empty()) throw std::runtime_error("V3000 write failed");
    bytes += output.size();
  }
  return bytes;
}

int main(int argc, char** argv) {
  try {
    if (argc != 5) throw std::runtime_error("usage: PATH MODE REPEATS EXPECTED_RECORDS");
    std::ifstream file(argv[1], std::ios::binary | std::ios::ate);
    if (!file) throw std::runtime_error("input unavailable");
    const auto size = file.tellg();
    if (size < 0) throw std::runtime_error("input size unavailable");
    std::string input(static_cast<std::size_t>(size), '\0');
    file.seekg(0);
    if (!file.read(input.data(), input.size())) throw std::runtime_error("input read failed");
    const auto repeats = std::stoull(argv[3]), expected = std::stoull(argv[4]);
    if (!repeats || !expected) throw std::runtime_error("positive counts required");
    const auto blocks = split_blocks(input, expected);
    auto prepared = read(blocks);
    const std::string mode(argv[2]);
    if (mode == "smiles") {
      OpenBabel::OBConversion observer;
      if (!observer.SetOutFormat("can")) throw std::runtime_error("SMILES unavailable");
      for (const auto& mol : prepared) std::cout << observer.WriteString(mol.get(), false);
      return 0;
    }
    OpenBabel::OBConversion writer;
    if (!writer.SetOutFormat("mol")) throw std::runtime_error("MOL writer unavailable");
    writer.AddOption("3", OpenBabel::OBConversion::OUTOPTIONS);
    write(prepared, writer);
    if (mode == "dump") {
      for (const auto& mol : prepared) std::cout << writer.WriteString(mol.get(), false) << "$$$$\n";
      return 0;
    }
    std::size_t units = 0;
    const auto started = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repeats; ++i) {
      if (mode == "parse") {
        const auto molecules = read(blocks);
        for (const auto& mol : molecules) units += mol->NumAtoms();
      } else if (mode == "write") {
        units += write(prepared, writer);
      } else if (mode == "roundtrip") {
        units += write(read(blocks), writer);
      } else {
        throw std::runtime_error("unknown operation");
      }
    }
    const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - started).count();
    std::cout << "{\"engine\":\"openbabel\",\"mode\":\"" << mode << "\",\"records\":"
              << expected * repeats << ",\"repeats\":" << repeats << ",\"output_units\":"
              << units << ",\"elapsed_ns\":" << elapsed << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
