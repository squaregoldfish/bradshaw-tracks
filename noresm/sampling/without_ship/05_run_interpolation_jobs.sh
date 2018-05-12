#!/bin/bash

track=$1
indir=$2
outdir=$3
threads=$4

if [ -z $track ] || [ -z $indir ] || [ -z $outdir ] || [ -z $threads ]
then
	echo "Usage: 05_run_interpolation_jobs.sh <track> <indir> <outdir> <threads>"
	exit
fi

./make_interpolation_jobs.tcl $track $indir $outdir
make -j $threads
./count_success_cells.tcl $track $outdir
