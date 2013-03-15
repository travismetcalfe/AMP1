#!/bin/csh
#  driver for for evolution programme, varying mass
#  usage: 
#  run-evol.mass <new mass> <trial mass> <100 Z> <case>
#  [<no. of timesteps>]

if(-d ttt) then
else
	echo "Making ttt"
	mkdir ttt
endif

if(-d log) then
else
	echo "Making log"
	mkdir log
endif


if($#argv == 0) then
  echo "Usage:"
  echo "run-evol.mass <new mass> <trial mass> <100 Z> <case>"
  echo "[<no. of timesteps>]"
  exit(1)
endif


set sdir = ./ttt


if($#argv <= 4) then
	set nt = 200
else
	set nt = $5
endif

set iqfit = "5"

echo "s/#mass/"$1"/" > tmp$$
echo "s/#tmass/"$2"/" >> tmp$$
echo "s/#z/"$3"/" >> tmp$$
echo "s/#case/"$4"/" >> tmp$$
echo "s/#iqfit/"$iqfit"/" >> tmp$$
echo "s/#nt/"$nt"/" >> tmp$$
#echo "s?#dir?"$dir"?" >> tmp$$
#echo "s?#root?"$root"?" >> tmp$$
echo "s?#evolfile?evol-file"$$".log?" >> tmp$$
sed -f tmp$$ evol.mass.d.tm.rin > tmp1$$
run-evol -dcliv0511z tmp1$$ $sdir/ttt.evol.mass.d.$1.Z$3.$4.s run-evol$$.d.log

if(-f ttt.rhs) cp -p ttt.rhs ttt.rhs.$4
if(-f ttt.crhs) cp -p ttt.crhs ttt.crhs.$4
if(-f ttt.emdl) cp -p ttt.emdl ttt.emdl.$4

#rm tmp$$ tmp1$$

cat run-evol$$.d.log evol-file$$.log > log/log.d.$1.Z$3.$4.s

rm run-evol$$.log evol-file$$.log

set-bsum csum.d.$1.Z$3.$4.s

sel-emdl.n.d emdl.d.$1.Z$3.$4.s

