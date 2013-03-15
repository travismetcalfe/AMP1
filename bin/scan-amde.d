#  script to call  $prgdir/adiajobs/scan-amde.d.x
if($#argv == 0 || $1 == "-help") then
	echo "Usage: scan-amde.d <file> [case]"
	echo "       case = 1: full set of eigenfunctions"
	echo "       case = 2: restricted set of eigenfunctions"
	echo "       default case: 2"
	exit(1)
endif 

if($#argv == 2) then
	set case = $2
else
	set case = 2
endif
(echo $1; echo $case) | $eprgdir/auxprg/scan-amde.d.x
