#  script to call  $eprgdir/auxprg/diff-gong.d.x
#  No arguments: prompts for input.
#  Alternatively may be called by
#  diff-gong.n.d <first file> <second file> <output file>
#  Here first model on each file is used, and differences are at fixed r/R

if($#argv == 0) then
	$eprgdir/auxprg/diff-gong.n.d.x
	exit(0)
else if($#argv < 3) then
	echo "Usage: diff-gong.n.d <first file> <second file> <output file> \"
	echo "                     [case [(r/R)_f]] [no. in first file] " \
	     "                     [no. in second file]"
	echo "Here case is an optional case number (default 1):"
	echo "1: differences at fixed r/R "
	echo "2: differences at fixed q"
	echo "3: differences at fixed r/r(last point)"
	echo "4: differences at fixed r"
	echo "5: differences at fixed p"
	echo "6: differences at fixed m/M(last point)"
	echo "7: differences at fixed m/M(fixed r/R)"
	echo \
	"If case = 7, the following argument must be the surface value of r/R"
	exit(1)
endif

set in1 = $1
set in2 = $2
set out = $3

if($#argv < 4) then
	set case = 1
else
	set case = $4
	if($case == "7") then 
	  set rfix = $5
	  shift
        endif
endif

if($#argv < 5) then
	set num1 = 1
else
	set num1 = $5
endif

if($#argv < 6) then
	set num2 = 1
else
	set num2 = $6
endif


if($case == 7) then
  (echo $in1; echo $in2; echo $num1 $num2; echo "2 2"; echo $out; echo $case; \
	 echo $rfix; echo "0"; echo "0" ) |  $eprgdir/auxprg/diff-gong.n.d.x
else
  (echo $in1; echo $in2; echo $num1 $num2; echo "2 2"; echo $out; echo $case; \
	 echo "0"; echo "0" ) |  $eprgdir/auxprg/diff-gong.n.d.x
endif
