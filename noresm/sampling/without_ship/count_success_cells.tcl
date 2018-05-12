#!/usr/bin/tclsh

set track [lindex $argv 0]
set outdir [lindex $argv 1]

set chan [open "output_root.txt"]
set outputRoot [string trim [read $chan]]
close $chan
set dir "${outputRoot}/${track}/without_ship/interpolation_outputs/${outdir}"

set cellCount 0
set cellTotal 10368.0

set successCount 0

for {set lon 1} {$lon <= 144} {incr lon} {
    for {set lat 1} {$lat <= 72} {incr lat} {
        incr cellCount
        puts -nonewline "\r [format "%.2f" [expr $cellCount / $cellTotal * 100]]%    "
        flush stdout

        set logfile "${dir}/${lon}_${lat}.log"

        if {[file exists $logfile]} {

            set chan [open $logfile r]
            set data [read $chan]
            close $chan

            if {[string first "SUCCESS" $data] > -1} {
                incr successCount
            }
        }
    }
}    

puts "\n$successCount"
