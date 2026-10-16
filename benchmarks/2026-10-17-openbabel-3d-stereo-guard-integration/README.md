# Main integration validation

These checks belong to combined source head 33ad98f2, containing main 4a97a599 (PR #960 writer and PR #963 MOL2 preservation) plus the measured tetrahedral guard. All 1,356 native tests and 4,256 untimed profile snapshots/rewrites pass. The measured guard function is unchanged. These are correctness observations only, not timing samples. The separate frozen performance packet names measured commit a093ce86. Compiler-cache cleaning and a sparse-fixture setup failure are explicit in the retained logs.
