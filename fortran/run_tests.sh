#!/bin/bash
set -e

echo "Compiling tests..."
cd $(dirname $0)
gfortran -o test_do_genRand.out md_global.f90 md_utilities.f90 test_do_genRand.f90

echo "Running tests..."
./test_do_genRand.out

echo "Cleaning up..."
rm -f test_do_genRand.out *.mod
