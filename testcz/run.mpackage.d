#!/bin/csh
#  driver for combined evolution and pulsation package
#  usage: 
#  run.package <case>

if($#argv == 0) then
  echo "Usage:"
  echo "run.mpackage <case>"
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

echo "s?#evolfile?evol-file"$$".log?" > tmp_sed$$
sed -f tmp_sed$$ evol.d.mcp.rin | get-input > tmp_evol$$
get-input adipls.cp.in > tmp_adi$$
get-input redistrb.cp.in > tmp_redistrb$$

echo "s/#case/"$1"/" > tmp_sedp$$
echo "s?#in_evol?"tmp_evol$$"?" >> tmp_sedp$$
echo "s?#in_adi?"tmp_adi$$"?" >> tmp_sedp$$
echo "s?#in_rdist?"tmp_redistrb$$"?" >> tmp_sedp$$

sed -f tmp_sedp$$ mpackage.d.rin > tmp_pack$$

run-evol -dcamliv0511z tmp_pack$$ $sdir/ttt.evol.mass.$1.s run-evol$$.log

#rm tmp$$ tmp1$$

cat run-evol$$.log evol-file$$.log.* > log/log.$1.$$.s

rm -f run-evol$$.log evol-file$$.log.*
#rm -f tmp_sed$$ tmp_evol$$ tmp_adi$$ tmp_redistrb$$ tmp_sedp$$ tmp_pack$$
