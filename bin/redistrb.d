#!/bin/csh
#  script to call  $prgdir/adiajobs/redistrb.d.x
if($#argv == 0) then
	$eprgdir/auxprg/redistrb.d.x
else if($1 == "-help") then
	echo "Usage: redistrb.d [control file]"
	exit(1)
else
	(get-input $1) | $eprgdir/auxprg/redistrb.d.x
endif
