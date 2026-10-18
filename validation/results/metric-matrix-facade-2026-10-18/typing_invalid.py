from chematic.rdkit_compat import rdMetricMatrixCalc as M
M.GetEuclideanDistMat(1)
M.GetEuclideanDistMat([["feature"]])
M.GetTanimotoSimMat([[1]])
wrong: list[float] = M.GetTanimotoDistMat([])
