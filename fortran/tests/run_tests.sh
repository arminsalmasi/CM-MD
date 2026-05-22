#!/bin/bash
set -e

cd $(dirname "$0")/..

# Ensure mod files are clean
rm -f *.mod

echo "Compiling tests..."
gfortran -o tests/test_eam md_global.f90 md_utilities.f90 tests/test_calcEAMPot.f90

echo "Running test_calcEAMPot..."
./tests/test_eam
