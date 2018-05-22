#!/usr/bin/tclsh

# Extracts the ship code, day (as day index from 1985-01-01)
# and lon/lat cell index to allow easy sampling of the NorESM model

set IN_FILE "/Data/science/SOCAT/v5/SOCATv5.tsv"
set OUT_FILE "/Data/Scratch/science/bradshaw-tracks/noresm/sampling/SOCATv5_indices.tsv"

variable START_YEAR 1985
variable END_YEAR 2016

variable MONTH_STARTS [list -1 0 31 59 90 120 151 181 212 243 273 304 334]
variable LEAP_MONTH_STARTS [list -1 0 31 60 91 121 152 182 213 244 274 305 335]

variable DAY_STARTS [list]
for {set i 0} {$i <= 365} {incr i} {
    lappend DAY_STARTS $i
}

variable LEAP_DAY_STARTS [list]

for {set i 1} {$i <= 365} {incr i} {
    lappend LEAP_DAY_STARTS [expr $i + (((366 / 365) / 365) * ($i - 1))]
}

proc calcJDate {year month day} {
    set result 0

    if {[isLeapYear $year]} {
        set result [lindex $::LEAP_MONTH_STARTS $month]
    } else {
        set result [lindex $::MONTH_STARTS $month]
    }

    set result [expr $result + $day]

    return $result
}

proc getDateIndex {year month day} {
    set jdate [calcJDate $year $month $day]

    set index ""
    set dayIndex 0

    if {$year >= $::START_YEAR && $year <= $::END_YEAR} {
        if {[isLeapYear $year]} {
            set dayIndex [getDayIndex $jdate $::DAY_STARTS]
        } else {
            set dayIndex [getDayIndex $jdate $::LEAP_DAY_STARTS]
        }
    
        set index [expr (($year - $::START_YEAR) * 365) + $dayIndex]
    }

    return $index
}

proc isLeapYear {year} {
    set leapYear 0
    if {[expr $year % 4] == 0} {
        set leapYear 1
        if {[expr $year % 100] == 0} {
            if {[expr $year % 400] != 0} {
                set leapYear 0
            }
        }
    }

    return $leapYear
}

proc getDayIndex {jdate days} {
    return [lindex $days [expr $jdate - 1]]
}


proc getLatCell {lat} {
    set boundary [expr floor($lat / 2.5)]
    set cell [expr $boundary + 36]
    if {$lat > 0} {
        set cell [expr $cell + 1]
    }
    return [expr int($cell)]
}

proc getLonCell {lon} {
    set result [expr (floor($lon / 2.5) + 1)]
    if {$result == 145} {
        set result 144
    }
    return [expr int($result)]
}


proc writeCell {outChan shipCode dateIndex lonCell latCell} {
    if {[string length $shipCode] > 0} {
        puts $outChan "$shipCode\t$dateIndex\t$lonCell\t$latCell"
    }
}

set outChan [open $OUT_FILE w]
puts $outChan "Ship Code\tDay Index\tLon Index\tLat Index"


set inChan [open $IN_FILE]

set currentShipCode ""
set currentDateIndex -1
set currentLonCell -1
set currentLatCell -1
set cellIndices [list]

set lineCount 0

while {[gets $inChan line] > 0} {

    incr lineCount
    if {[expr $lineCount % 1000] == 0} {
        puts -nonewline "\033\[2K\r$lineCount"
        flush stdout
    }

    set fields [split $line "\t"]
    
    # Check that we're in the right time range
    set year [lindex $fields 4]

    if {$year >= $START_YEAR && $year <= $END_YEAR} {
        set changed 0

        set expoCode [lindex $fields 0]
        set shipCode [string range $expoCode 0 [expr [string length $expoCode] - 9]]

        if {[string compare $currentShipCode $shipCode] != 0} {
            set changed 1
        }

        scan [lindex $fields 5] %d month
        scan [lindex $fields 6] %d day


        set dateIndex [getDateIndex $year $month $day]
        if {$dateIndex != $currentDateIndex} {
            set changed 1
        }

        set lonCell [getLonCell [lindex $fields 10]]
        if {$lonCell != $currentLonCell} {
            set changed 1
        }

        set latCell [getLatCell [lindex $fields 11]]
        if {$latCell != $currentLatCell} {
            set changed 1
        }

        if {$changed} {
            writeCell $outChan $shipCode $dateIndex $lonCell $latCell

            set currentShipCode $shipCode
            set currentDateIndex $dateIndex
            set currentLonCell $lonCell
            set currentLatCell $latCell

        }
    }
}

# Write the last cell data
writeCell $outChan $shipCode $dateIndex $lonCell $latCell

close $inChan
close $outChan
puts ""

