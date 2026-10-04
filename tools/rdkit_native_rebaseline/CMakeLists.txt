cmake_minimum_required(VERSION 3.18)
project(rdkit_native_rebaseline LANGUAGES CXX)

# RDKit's installed rdkit_base target can reference Eigen even with the
# descriptors3D feature disabled; load the imported target first.
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(RDKit REQUIRED CONFIG)

add_executable(rdkit_native_rebaseline main.cpp)
target_compile_features(rdkit_native_rebaseline PRIVATE cxx_std_20)
target_link_libraries(rdkit_native_rebaseline PRIVATE
  RDKit::SmilesParse
  RDKit::GraphMol
  RDKit::SubstructMatch
  RDKit::Fingerprints
  RDKit::CIPLabeler
)
