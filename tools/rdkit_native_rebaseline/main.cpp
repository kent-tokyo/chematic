// Independent RDKit C++ oracle for the exposed 10k rebaseline corpus.
// This is output evidence, not a speed benchmark or a Python-wrapper substitute.
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string escape(const std::string &value) {
  std::string result;
  result.reserve(value.size() + 2);
  result.push_back('"');
  constexpr char hex[] = "0123456789abcdef";
  for (unsigned char ch : value) {
    switch (ch) {
      case '"': result += "\\\""; break;
      case '\\': result += "\\\\"; break;
      case '\n': result += "\\n"; break;
      case '\r': result += "\\r"; break;
      case '\t': result += "\\t"; break;
      default:
        if (ch < 0x20) {
          result += "\\u00";
          result += hex[ch >> 4];
          result += hex[ch & 0x0f];
        } else {
          result.push_back(static_cast<char>(ch));
        }
    }
  }
  result.push_back('"');
  return result;
}

std::vector<std::string> read_lines(const std::string &path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open " + path);
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty()) lines.push_back(line);
  }
  return lines;
}

template <typename T>
void emit_indices(std::ostream &out, const T &indices) {
  out << '[';
  bool first = true;
  for (auto index : indices) {
    if (!first) out << ',';
    first = false;
    out << index;
  }
  out << ']';
}

void emit_cip(std::ostream &out, RDKit::ROMol &mol) {
  RDKit::CIPLabeler::assignCIPLabels(mol);
  out << "\"cip_atoms\":{ ";
  bool first = true;
  for (const auto atom : mol.atoms()) {
    if (!atom->hasProp("_CIPCode")) continue;
    if (!first) out << ',';
    first = false;
    out << escape(std::to_string(atom->getIdx())) << ':'
        << escape(atom->getProp<std::string>("_CIPCode"));
  }
  out << "},\"cip_bonds\":{ ";
  first = true;
  for (const auto bond : mol.bonds()) {
    if (!bond->hasProp("_CIPCode")) continue;
    if (!first) out << ',';
    first = false;
    auto begin = bond->getBeginAtomIdx();
    auto end = bond->getEndAtomIdx();
    if (begin > end) std::swap(begin, end);
    out << escape(std::to_string(begin) + "-" + std::to_string(end)) << ':'
        << escape(bond->getProp<std::string>("_CIPCode"));
  }
  out << '}';
}

void emit_morgan(std::ostream &out, const RDKit::ROMol &mol,
                 const RDKit::FingerprintGenerator<std::uint32_t> &generator) {
  std::unique_ptr<ExplicitBitVect> bits(generator.getFingerprint(mol));
  if (!bits || bits->getNumBits() != 2048) {
    throw std::runtime_error("unexpected Morgan fingerprint size");
  }
  out << "\"morgan_on_bits\":[";
  bool first = true;
  for (unsigned int bit = 0; bit < 2048; ++bit) {
    if (!bits->getBit(bit)) continue;
    if (!first) out << ',';
    first = false;
    out << bit;
  }
  out << ']';
}

void emit_smarts(std::ostream &out, const RDKit::ROMol &mol,
                 const std::vector<std::unique_ptr<RDKit::RWMol>> &queries) {
  out << "\"smarts\":[";
  RDKit::SubstructMatchParameters params;
  params.uniquify = true;
  params.maxMatches = 1000;  // Python GetSubstructMatches default.
  for (std::size_t query_index = 0; query_index < queries.size(); ++query_index) {
    if (query_index) out << ',';
    std::set<std::vector<int>> sets;
    for (const auto &match : RDKit::SubstructMatch(mol, *queries[query_index], params)) {
      std::vector<int> atoms;
      atoms.reserve(match.size());
      for (const auto &pair : match) atoms.push_back(pair.second);
      std::sort(atoms.begin(), atoms.end());
      sets.insert(std::move(atoms));
    }
    out << '[';
    bool first = true;
    for (const auto &atoms : sets) {
      if (!first) out << ',';
      first = false;
      emit_indices(out, atoms);
    }
    out << ']';
  }
  out << ']';
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "usage: rdkit_native_rebaseline CORPUS QUERIES_TXT OUTPUT_JSONL\n";
    return 2;
  }
  try {
    const auto smiles = read_lines(argv[1]);
    const auto query_strings = read_lines(argv[2]);
    if (smiles.size() != 10000 || query_strings.size() != 31) {
      throw std::runtime_error("expected 10000 molecules and 31 SMARTS queries");
    }
    std::vector<std::unique_ptr<RDKit::RWMol>> queries;
    queries.reserve(query_strings.size());
    for (const auto &query : query_strings) {
      auto parsed = std::unique_ptr<RDKit::RWMol>(RDKit::SmartsToMol(query));
      if (!parsed) throw std::runtime_error("SMARTS parse failed: " + query);
      queries.push_back(std::move(parsed));
    }
    std::unique_ptr<RDKit::FingerprintGenerator<std::uint32_t>> generator(
        RDKit::MorganFingerprint::getMorganGenerator<std::uint32_t>(2));
    std::ofstream out(argv[3]);
    if (!out) throw std::runtime_error("cannot open output file");
    for (std::size_t index = 0; index < smiles.size(); ++index) {
      out << "{\"input_index\":" << index << ",\"smiles\":" << escape(smiles[index]);
      try {
        auto mol = std::unique_ptr<RDKit::RWMol>(RDKit::SmilesToMol(smiles[index]));
        if (!mol) {
          out << ",\"status\":\"parse_failure\"}";
        } else {
          out << ",\"status\":\"ok\",\"canonical\":"
              << escape(RDKit::MolToSmiles(*mol));
          out << ',';
          emit_cip(out, *mol);
          out << ',';
          emit_morgan(out, *mol, *generator);
          out << ',';
          emit_smarts(out, *mol, queries);
          out << '}';
        }
      } catch (const std::exception &exc) {
        // A failed operation must never turn an incomplete JSON row into evidence.
        throw std::runtime_error("row " + std::to_string(index) + ": " + exc.what());
      }
      out << '\n';
    }
  } catch (const std::exception &exc) {
    std::cerr << "native rebaseline failed: " << exc.what() << '\n';
    return 1;
  }
  return 0;
}
