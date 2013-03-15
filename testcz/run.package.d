#!/bin/csh
#  driver for combined evolution and pulsation package
#  usage: 
#  run-evol.mass <trial mass> <100 Z> <case>

if($#argv == 0) then
  echo "Usage:"
  echo "run.package.d <trial mass> <100 Z> <case>"
  exit(1)
endif

if(-d emdl) then
else
  mkdir emdl
  echo "Making directory emdl"
endif

if(-d osc) then
else
  mkdir osc
  echo "Making directory osc"
endif

if(-d amdl) then
else
  mkdir amdl
  echo "Making directory amdl"
endif

if(-d gong) then
else
  mkdir gong
  echo "Making directory gong"
endif

if(-d ttt) then
else
  mkdir ttt
  echo "Making directory ttt"
endif

if(-d log) then
else
  mkdir log
  echo "Making directory log"
endif

set sdir = ./ttt

echo "s/#tmass/"$1"/" > tmp_sed$$
echo "s/#z/"$2"/" >> tmp_sed$$
echo "s?#evolfile?evol-file"$$".log?" >> tmp_sed$$
sed -f tmp_sed$$ evol.d.cp.rin | get-input > tmp_evol$$
get-input adipls.d.cp.in > tmp_adi$$
get-input redistrb.cp.in > tmp_redistrb$$

echo "s/#tmass/"$1"/" > tmp_sedp$$
echo "s/#case/"$3"/" >> tmp_sedp$$
echo "s?#in_evol?"tmp_evol$$"?" >> tmp_sedp$$
echo "s?#in_adi?"tmp_adi$$"?" >> tmp_sedp$$
echo "s?#in_rdist?"tmp_redistrb$$"?" >> tmp_sedp$$

sed -f tmp_sedp$$ package.d.rin > tmp_pack$$

run-evol -dcap9z tmp_pack$$ $sdir/ttt.evol.mass.$2.$3.s run-evol$$.log


cat run-evol$$.log evol-file$$.log.* > log/log.Z$2.$3.$$.s

rm -f run-evol$$.log evol-file$$.log.*
rm -f tmp_sed$$ tmp_evol$$ tmp_adi$$ tmp_redistrb$$ tmp_sedp$$ tmp_pack$$
