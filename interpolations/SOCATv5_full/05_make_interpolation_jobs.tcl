#!/usr/bin/tclsh

set outputDirRoot "/Data/Scratch/science/bradshaw-tracks/interpolations/SOCATv5_full/interpolation_outputs"

set indir [lindex $argv 0]
set outdir [lindex $argv 1]

if {[string length $indir] == 0 || [string length $outdir] == 0} {
	puts "Missing directories"
	exit
}

set jobList ""
set jobText ""

for {set lon 1} {$lon <= 144} {incr lon} {
    for {set lat 1} {$lat <= 72} {incr lat} {
        append jobList "${lon}_${lat} "

        append jobText "${lon}_${lat} :\n"
        append jobText "\tR --no-save --slave lon=\"${lon}\" lat=\"${lat}\" indir=\"${outputDirRoot}/${indir}\" outdir=\"${outputDirRoot}/${outdir}\" < do_interpolation.R\n\n"
    }
}

set outChan [open "Makefile" w]
puts $outChan "all : $jobList\n"
puts $outChan $jobText
close $outChan
