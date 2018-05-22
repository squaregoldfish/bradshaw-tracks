#!/usr/bin/tclsh

# This must be run after the without_ship run

set START_YEAR 1985
set END_YEAR 2016

set track [lindex $argv 0]
set areFile [lindex $argv 1]

set IN_FILE "/Data/Scratch/science/bradshaw-tracks/noresm/sampling/${track}/SOCATv5_indices.tsv"
set OUT_FILE "/Data/Scratch/science/bradshaw-tracks/noresm/sampling/${track}/SOCATv5_indices_with_are_ship.tsv"



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

proc writeCell {outChan dateIndex lonCell latCell} {
    puts $outChan "AREO\t$dateIndex\t$lonCell\t$latCell"
}

# Duplicate the existing file, which doesn't have anything from our ship track in it
file copy -force $IN_FILE $OUT_FILE

# ...and open it for writing
set outChan [open $OUT_FILE a]

# Now load Are's track file
set areChan [open $areFile]
set areData [read $areChan]
close $areChan

# Skip the header
set lines [lrange [split $areData "\n"] 1 end]


set currentDateIndex -1
set currentLonCell -1
set currentLatCell -1

foreach line $lines {

    set fields [split $line "\t"]
    if {[llength $fields] > 1} {
        set changed 0

        set year [expr int([lindex $fields 1])]
        set month [expr int([lindex $fields 2])]
        set day [expr int([lindex $fields 3])]

        set lon [lindex $fields 6]
        if {$lon < 0} {
            set lon [expr 360 - abs($lon)]
        }
        set lat [lindex $fields 7]

        set dateIndex [getDateIndex $year $month $day]

        set dateIndex [getDateIndex $year $month $day]
        if {$dateIndex != $currentDateIndex} {
            set changed 1
        }

        set lonCell [getLonCell $lon]
        if {$lonCell != $currentLonCell} {
            set changed 1
        }

        set latCell [getLatCell $lat]
        if {$latCell != $currentLatCell} {
            set changed 1
        }

        if {$changed} {
            writeCell $outChan $dateIndex $lonCell $latCell

            set currentDateIndex $dateIndex
            set currentLonCell $lonCell
            set currentLatCell $latCell
        }
    }

}

# Write the last cell data
writeCell $outChan $dateIndex $lonCell $latCell

close $outChan
