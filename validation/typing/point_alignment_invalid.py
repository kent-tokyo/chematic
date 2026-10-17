from chematic.rdkit_compat import rdAlignment

rows = [[0.0, 0.0, 0.0]]
rdAlignment.GetAlignmentTransform("bad", rows)
rdAlignment.GetAlignmentTransform(rows, rows, weights="bad")
rdAlignment.GetAlignmentTransform(rows, rows, reflect="yes")
rdAlignment.GetAlignmentTransform(rows, rows, maxIterations=2.5)
