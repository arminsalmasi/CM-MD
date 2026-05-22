#!/bin/bash
cd $(dirname $0)/..
gfortran -I. -o tests/test_md_utilities md_global.f90 md_utilities.f90 tests/test_md_utilities.f90
./tests/test_md_utilities
