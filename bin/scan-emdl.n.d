#!/bin/csh
#  script to call  scan-emdl.n.d.x

if($#argv == 0) then
	$eprgdir/auxprg/scan-emdl.n.d.x
else if($#argv == 1) then
	(echo $1; echo "0") | $eprgdir/auxprg/scan-emdl.n.d.x
else if($#argv == 2) then
	(echo $1; echo $2) | $eprgdir/auxprg/scan-emdl.n.d.x
endif
