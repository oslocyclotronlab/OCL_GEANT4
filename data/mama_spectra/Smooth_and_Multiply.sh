#!/usr/bin/expect

set SleepTime 0.5

set multFactor [lindex $argv 2];
set filename [lindex $argv 1];
set prefix [lindex $argv 0];
set outfile "$prefix$filename"

set timeout 2

spawn mama
sleep 1.5
expect "mama>"
sleep $SleepTime
send "re\r"
expect "Destination spectrum <1>:"
send "1\r"
sleep $SleepTime
expect "Filename          <TEST>:"
send "$filename\r"
#expect "mama>"
#send "ds2\r"

# Applying the smoothing
expect "mama>"
# sleep $SleepTime
send "sm\r"
expect "Destination spectrum <2>:"
# sleep $SleepTime
send "2\r"
expect "Source spectrum      <1>:"
# sleep $SleepTime
send "1\r"
expect "Dimension of spectrum"
# sleep $SleepTime
send "\r"
expect "Write FWHMx (ch) around ch "
# sleep $SleepTime
send "20\r"
expect "Write FWHMx (ch) around ch "
# sleep $SleepTime
send "30\r"
expect "Smoothing-matrix OK? (y/n) <y>:"
sleep $SleepTime
send "y\r"
expect "mama>"

# multiplication by "arb." factor
# sleep $SleepTime
# send "ar\r"
# expect "Write your expression:"
# sleep $SleepTime
# send "2=2*$multFactor \r"
# expect "mama>"

# write file
sleep $SleepTime
send "wr\r"
expect "Spectrum to write            <2>:"
# sleep $SleepTime
send "2\r"
expect "Please, choose your type     <1>:"
# sleep $SleepTime
send "1\r"
expect "Cal. coeff. a0 (keV) on x-axis     <        0.0>:"
# sleep $SleepTime
send "\r"
expect "Cal. coeff. a1 (keV/ch) on x-axis "
# sleep $SleepTime
send "\r"
expect "Cal. coeff. a2 (keV/ch2) on x-axis < 0.0000E+00>:"
# sleep $SleepTime
send "\r"
expect "Length of output-spectrum"
# sleep $SleepTime
send "\r"
expect "Filename                  <SPEC>:"
# sleep $SleepTime
send "$outfile\r"
expect "mama>"
sleep $SleepTime
send "st\r"
expect "Are you sure you want to exit? (y/n)"
send "y\r"