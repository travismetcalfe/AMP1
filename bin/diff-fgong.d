#  script to call  $eprgdir/auxprg/diff-fgong.d.x
if($#argv == 0) then
	$eprgdir/auxprg/diff-fgong.d.x
	exit(0)
else if ($#argv < 3) then
	echo "Usage: diff-fgong.d <first file> <second file> <output file> \"
	echo "                    [case]" 
	echo "Here case is an optional case number (default 1):"
	echo "1: differences at fixed r/R "
	echo "2: differences at fixed q"
	echo "3: differences at fixed r/r(last point)"
	echo "4: differences at fixed p"
	echo "5: differences at fixed (R-r)/R"
	exit(1)
endif

set in1 = $1
set in2 = $2
set out = $3

if($#argv < 4) then
	set case = 1
else
	set case = $4
endif

(echo $in1; echo $in2; echo "1 1"; echo $out; echo $case; echo "0"; echo "0") | \
	$eprgdir/auxprg/diff-fgong.d.x
