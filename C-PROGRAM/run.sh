#!/bin/bash

# base directory
base_dir="$( cd "$( dirname "$0" )" && pwd )"

# output directory
output_dir="$base_dir/OUTPUTS"

# input directory
input_dir="$base_dir/INPUTS"

# log file
log_file="$output_dir/simulator_log.txt"

# go into the inputs directory and generate input file
cd INPUTS/
./bubbleSizeDistribution.py
./column.py
./convergence.py
./gas.py
./grid.py
./interfacial.py
./liquid.py
./separator.py
./solids.py
cd ..


# clean up output directory
cd OUTPUTS/
rm *.dat *.txt
cd ..

# clean up figure directory
cd FIGURES/
rm *.pdf
cd ..

# compile source files
cd SOURCES/
make clean
make all

# run executable
./simulator $input_dir $output_dir
cd ..

# create figures
./createGraphs.py