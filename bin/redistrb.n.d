#!/bin/csh
#  script to call  $prgdir/adiajobs/redistrb.n.d.x
if($#argv == 0) then
	$eprgdir/auxprg/redistrb.n.d.x
else
	get-input $1 | $eprgdir/auxprg/redistrb.n.d.x
endif
