#!/bin/bash

gfortran md_datastructure.f90 md_monitor.f90  md_utilities.f90 md_md.f90 main.f90 -o out.x

rm *.mod

