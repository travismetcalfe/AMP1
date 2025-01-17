#!/bin/csh
#  driver for adipls.d.x Usage:
#  run-adipls  <input> <output>
#  uses the file <input> as input (assumed to be in current directory,
#  if complete path is not given)
#  produces output file on <output>, 
#  or on scratch/<user>/adipls/<output>
#  if complete path is not given

# root scratch directory

set scratch = scratch
if(-d scratch) then
else
  mkdir scratch
  echo "Making scratch"
endif


if($#argv == 0) then
  echo "Usage:"
  echo "run-adipls  <input> [output]"
  echo "If output file is not specified, output goes to standard output"
  exit(1)
endif

set user = `whoami`

#  test for creating scratch directory

if(-d $scratch/$user) then
else
	mkdir $scratch/$user
	echo Making $scratch/$user
endif

if(-d $scratch/$user/adipls) then
else
	mkdir $scratch/$user/adipls
	echo Making $scratch/$user/adipls
endif

set prog = adipls.c.d.x

set ain = $1
if($#argv > 1) then
  set aout = $2
endif

set dirin=`test-dir $ain`
if ($dirin == "n") then
	set in = $cwd/$ain
else
	set in = $in
endif

if($?aout) then
  set dirout=`test-dir $aout`
  if ($dirout == "n") then
	set out = $scratch/$user/adipls/$aout
  else
	set out = $aout
  endif
endif

if($?out) then
  get-input $in | $eprgdir/adi/$prog > $out
else
  get-input $in | $eprgdir/adi/$prog 
endif

