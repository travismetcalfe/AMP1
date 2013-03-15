#!/bin/csh
#  script to call  $prgdir/adiajobs/set-obs.d.x
#  usage (with arguments) set-obs.d <case> <input file> [<output file>]
#  If output file is not given, assumes that <input file> contains 
#  trailer to be applied to agsm on input and obs on output.
#  with no arguments: prompts

if($#argv == 1 && $1 == "-help") then
  echo "usage (with arguments) "
  echo "set-obs.d <case> <input file> [<output file>] [precision]"
  echo "If output file is not given, assumes that <input file> contains "
  echo "trailer to be applied to agsm on input and obs on output."
  echo "with no arguments: prompts"
  echo "case:"
  echo "1: grand summary, variational frequency."
  echo "2: short summary."
  echo "4: grand summary, from eigenfrequency in cs(20)."
  echo "   Note that this allows setting Cowling approximation frequency"
  echo "5: grand summary, from Richardson extrapolation frequency"
  echo "6: grand summary, from (possibly corrected) eigenfrequency in cs(21)"
  echo "If icasein gt 10, set according to icasein-10, including mode energy"
  echo "precision:"
  echo "1: f8.2"
  echo "2: f10.4 (default)"
  echo "3: f12.6"
  exit(1)
endif

if($#argv > 3) then
	set prec = $4
else
	set prec = "2"
endif

if($#argv < 2) then
	$eprgdir/auxprg/set-obs.d.x
else if($#argv == 2) then
	(echo $1; echo agsm.$2; echo obs.$2; echo $prec) | \
	$eprgdir/auxprg/set-obs.d.x
else
	(echo $1; echo $2; echo $3; echo $prec) | \
	$eprgdir/auxprg/set-obs.d.x
endif
