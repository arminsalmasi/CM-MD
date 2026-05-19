## 2024-05-24 - Avoid sqrt() and explicit powers in MD inner loops
**Learning:** The codebase contains manual Lennard-Jones force calculations that compute `sqrt()` for distance and explicit powers (e.g., `^13` and `^7`) in the inner N^2 loops. This is a severe anti-pattern for Molecular Dynamics and also caused a calculation bug in scaling.
**Action:** Always refactor MD distance and potential/force calculations to use squared distance (`r^2`) and multiplicative inverses instead of square roots and fractional/large exponentiation.
