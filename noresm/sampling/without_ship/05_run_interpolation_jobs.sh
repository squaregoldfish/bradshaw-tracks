#!/bin/bash

threads=$1
track=$2
indir=$3
outdir=$4

if [ -z $track ] || [ -z $indir ] || [ -z $outdir ] || [ -z $threads ]
then
	echo "Usage: 05_run_interpolation_jobs.sh <threads> <track> <indir> <outdir>"
	exit
fi

./make_interpolation_jobs.tcl $track $indir $outdir
make -j $threads
./count_success_cells.tcl $track $outdir
