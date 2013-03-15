#!/bin/csh
#  script to call  sel-emdl.n.d.x
#  usage: sel-emdl.n.d <input file> [<output file> <model number>]
#  if no model number is given, the first model is used
#  if no output model is given, it is set to <input file>.0

if($#argv == 0) then
	$eprgdir/auxprg/sel-emdl.n.d.x
else if($1 == "-help") then
	echo "Usage: sel-emdl.n.d <input file> [<output file> <model number>]"
	exit(1)
else
	if($#argv == 1) then
		set out = $1.0
		set num = 1
	else if($#argv == 2) then
		set out = $2
		set num = 1
	else 
		set out = $2
		set num = $3
	endif
	( echo $1; echo $out; echo $num ) | \
	$eprgdir/auxprg/sel-emdl.n.d.x
endif
