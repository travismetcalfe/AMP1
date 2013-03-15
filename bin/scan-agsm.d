#!/bin/csh
#  script to call scan-agsm.d
if($#argv == 0 || $1 == "-help") then
	echo "Usage: scan-agsm.d <grand summary file> [ishort]"
	echo "If ishort is set and equal to 1 or 2, print only first mode"
	echo "for each degree"
	echo "If ishort is set and equal to -1, print detailed results"
	exit(1)
endif
if($#argv > 1) then
	set ib = $2
else
	set ib = 0
endif

(echo $1; echo $ib) | $eprgdir/auxprg/scan-agsm.d.x
