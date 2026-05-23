#!/bin/bash
set -e
gfortran -o test_cm.out md_global.f90 md_utilities.f90 test_center_of_mass.f90
./test_cm.out
