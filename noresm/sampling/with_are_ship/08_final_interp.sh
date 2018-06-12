#!/bin/bash

track="$1"

OUTPUT_ROOT="`cat output_root.txt|xargs`/$track/with_are_ship"

echo "Removing old input..."
rm -r /tmp/final_interp_input
echo "Copying input..."
cp -RH "${OUTPUT_ROOT}/uncertainty_output"  /tmp/final_interp_input
echo "Removing old output..."
rm -r  "${OUTPUT_ROOT}/final_output"
echo "Running..."
R --slave --no-save track="$track" < final_interp.R

