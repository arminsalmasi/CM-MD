#!/bin/bash

# Exit on error
set -e

echo "Compiling test_do_alloc..."
gfortran -c md_global.f90
gfortran -c md_utilities.f90
gfortran -o test_do_alloc md_global.o md_utilities.o test_do_alloc.f90

echo "Running test_do_alloc..."
./test_do_alloc
