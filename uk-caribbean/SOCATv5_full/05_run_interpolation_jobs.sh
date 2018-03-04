#!/bin/bash

indir=$1
outdir=$2
threads=$3

if [ -z $indir ] || [ -z $outdir ] || [ -z $threads ]
then
	echo "Usage: 05_run_interpolation_jobs.sh <indir> <outdir> <threads>"
	exit
fi

./make_interpolation_jobs.tcl $indir $outdir
make -j $threads
./count_success_cells.tcl $outdir
