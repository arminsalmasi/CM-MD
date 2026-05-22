## 2024-05-24 - File I/O Error Leakage
**Vulnerability:** The Fortran code reading potential data (`fortran/md_monitor.f90`) used standard `OPEN` and `READ` without checking for file I/O errors. Unhandled file exceptions can cause the application to crash, exposing system paths, OS configurations, or full stack traces to end users (Information Leakage).
**Learning:** Legacy scientific Fortran applications often lack basic fail-safe handling in I/O operations, trusting that the specific file paths and structures exist.
**Prevention:** In Fortran, always use the `IOSTAT=` specifier with `OPEN` and `READ`/`WRITE` statements, and check for `IOSTAT /= 0` to fail securely and handle the issue gracefully, ensuring that system internals are not leaked on failure.
