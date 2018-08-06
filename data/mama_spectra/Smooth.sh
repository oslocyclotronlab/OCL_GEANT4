#!/usr/bin/expect



set filename [lindex $argv 1];
set prefix [lindex $argv 0];
set outfile "$prefix$filename"

spawn mama
sleep 1.5
expect "mama>"
send "re\r"
expect "Destination spectrum <1>:"
send "1\r"
expect "Filename          <TEST>:"
send "$filename\r"
#expect "mama>"
#send "ds2\r"

# Applying the smoothing
expect "mama>"
send "sm\r"
expect "Destination spectrum <2>:"
send "2\r"
expect "Source spectrum      <1>:"
send "1\r"
expect "Dimension of spectrum <4200>:"
send "4200\r"
expect "Write FWHMx (ch) around ch x=  420 <   1.0>:"
send "20\r"
expect "Write FWHMx (ch) around ch x= 3780 <   3.0>:"
send "30\r"
expect "Smoothing-matrix OK? (y/n) <y>:"
send "y\r"
expect "mama>"
send "wr\r"
expect "Spectrum to write            <2>:"
send "2\r"
expect "Please, choose your type     <1>:"
send "1\r"
expect "Cal. coeff. a0 (keV) on x-axis     <        0.0>:"
send "0.0\r"
expect "Cal. coeff. a1 (keV/ch) on x-axis  <      5.000>:"
send "5.000\r"
expect "Cal. coeff. a2 (keV/ch2) on x-axis < 0.0000E+00>:"
send "0.0000E+00\r"
expect "Length of output-spectrum <4200>:"
send "4200\r"
expect "Filename                  <SPEC>:"
send "$outfile\r"
expect "mama>"
send "st\r"
expect "Are you sure you want to exit? (y/n)"
send "y\r"