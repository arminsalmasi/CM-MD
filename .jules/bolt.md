## 2024-05-18 - Vectorization and Newton's Third Law in MATLAB
**Learning:** MATLAB loops are exceptionally slow compared to vectorized operations. Furthermore, in pairwise computations (like MD force calculations), calculating interactions twice ((N^2)$) instead of using Newton's Third Law ((N^2/2)$) is a common anti-pattern in academic codebases.
**Action:** Always look for nested spatial loops in MATLAB and replace them with array operations. For any pairwise interaction loop, apply symmetry (Newton's 3rd law) to halve the work.
