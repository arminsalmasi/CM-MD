#!/bin/bash
cd "$(dirname "$0")/.."

# Compile the modules
gfortran -c md_global.f90 md_utilities.f90

# Compile and run the tests
gfortran -o tests/test_velVerlet tests/test_velVerlet.f90 md_global.o md_utilities.o
echo "Running test_velVerlet..."
./tests/test_velVerlet

# Cleanup
rm -f *.mod md_global.o md_utilities.o tests/test_velVerlet
