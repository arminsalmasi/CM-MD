#!/bin/bash
gfortran -o test_md.out md_global.f90 md_utilities.f90 test_md_utilities.f90
if [ $? -eq 0 ]; then
    echo "Compilation successful. Running tests..."
    ./test_md.out
    if [ $? -eq 0 ]; then
        echo "Tests finished successfully."
    else
        echo "Tests failed."
        exit 1
    fi
else
    echo "Compilation failed."
    exit 1
fi
rm *.mod
