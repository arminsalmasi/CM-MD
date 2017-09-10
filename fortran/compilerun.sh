#!/bin/bash
gfortran -o md.out md_global.f90 md_monitor.f90 md_utilities.f90 main.f90 
rm *.mod

