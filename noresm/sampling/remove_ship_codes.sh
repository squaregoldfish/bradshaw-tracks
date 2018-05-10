#!/bin/bash

INFILE="/Data/Scratch/science/bradshaw-tracks/noresm/sampling/SOCATv5_indices.tsv"
OUTDIR="/Data/Scratch/science/bradshaw-tracks/noresm/sampling"

ship_track_dir=$1


grep -v -f "./ship_tracks/$ship_track_dir/removal_codes.csv" $INFILE > "$OUTDIR/$ship_track_dir/SOCATv5_indices.tsv"
