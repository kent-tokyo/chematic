cmake_minimum_required(VERSION 3.18)
project(rdkit_native_rebaseline LANGUAGES CXX)

find_package(RDKit REQUIRED CONFIG)

add_executable(rdkit_native_rebaseline main.cpp)
target_compile_features(rdkit_native_rebaseline PRIVATE cxx_std_17)
target_link_libraries(rdkit_native_rebaseline PRIVATE
  RDKit::SmilesParse
  RDKit::GraphMol
  RDKit::SubstructMatch
  RDKit::Fingerprints
  RDKit::CIPLabeler
)
