#!/bin/bash
echo "Removing old input..."
rm -r /tmp/final_interp_input
echo "Copying input..."
cp -RH /Data/Scratch/science/bradshaw-tracks/interpolations/SOCATv5_full/uncertainty_output  /tmp/final_interp_input
echo "Removing old output..."
rm -r  /Data/Scratch/science/bradshaw-tracks/interpolations/SOCATv5_full/final_output
echo "Running..."
R --slave --no-save < final_interp.R

