#!/bin/bash

# Regrid NorESM files onto a 0.01°x0.01° grid
# for spatial autocorrelation analysis

indir=$1
outfile=$2

rm -rf /tmp/regrid
mkdir "/tmp/regrid/"

if [ -z $indir ] || [ -z $outfile ]
then
	echo "Usage: 2.5_grid.sh <indir> <outfile>"
	exit
fi

for f in `find $indir`
do
    filename=`basename $f`
    echo $filename
	cdo remapbil,r144x72 $f "/tmp/regrid/${filename}"
done

ncrcat /tmp/regrid/*nc "${outfile}"

rm -rf /tmp/regrid
