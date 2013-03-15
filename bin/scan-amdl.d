#  script to call  $prgdir/adiajobs/scan-amdl.d.x
if($#argv == 0 || $1 == "-help") then
	echo "Usage: scan-amdl.d <input file>"
	exit(1)
else
	echo $1 | $eprgdir/auxprg/scan-amdl.d.x
endif
