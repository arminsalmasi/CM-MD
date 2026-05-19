## 2024-05-19 - [Missing Input Validation in Fortran MD Monitor]
**Vulnerability:** The `fortran/md_monitor.f90` file read variables like time, timestep, temperature, volume, and number of atoms directly from standard input without any validation.
**Learning:** This could lead to a variety of crashes, including division by zero when calculating `nTstps = floor(t/dt)` if `dt <= 0`, or memory corruption/allocation errors if negative numbers of atoms were provided.
**Prevention:** Always validate user-provided numerical inputs, ensuring physical parameters have sensible non-negative or strictly positive constraints before using them in mathematical operations or array sizing.
