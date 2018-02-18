#!/bin/bash

OUTPUT_ROOT=`cat output_root.txt|xargs`

echo "Removing old input..."
rm -r /tmp/final_interp_input
echo "Copying input..."
cp -RH "${OUTPUT_ROOT}/uncertainty_output"  /tmp/final_interp_input
echo "Removing old output..."
rm -r  "${OUTPUT_ROOT}/final_output"
echo "Running..."
R --slave --no-save < final_interp.R

