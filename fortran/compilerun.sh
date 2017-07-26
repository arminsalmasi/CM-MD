#!/bin/bash

gfortran md_datastructure.f95 md_monitor.f95 md_utilities.f95 main.f95 -o out.x

rm *.mod

