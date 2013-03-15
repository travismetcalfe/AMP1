#  script to call  $prgdir/auxprg/form-gong.n.d.x
if($#argv == 0) then
  $eprgdir/auxprg/form-gong.n.d.x
  exit(0)
else if($1 == "-help") then
  echo "Usage: form-gong.n.d <trailer> [case] [model no.]"
  echo "   case 0: old (pre-1995) format"
  echo "   case 1: including He3 and CNO abundances"
  echo "   case 2: including He3 and CNO abundances and Z (default)"
  echo "   case 3: including He3 and CNO abundances and Z,"
  echo "           as well as Gamma_1 derivatives"
  echo "Note: It is assumed that the header file is in head.<trailer>"
  exit(1)
endif

if($#argv < 2) then
    set case = "2"
else
    set case = $2
endif

if($#argv < 3) then
    set number = "1"
    set fgong = fgong.$1
    set head = head.$1
else
    set number = $3
    set fgong = fgong.$1.$number
    set head = head.$1.$number
endif

echo $case

(echo gong.$1; echo $fgong; echo $number; echo $head; echo $case; \
  echo dgamma1.$1) |  \
  $eprgdir/auxprg/form-gong.n.d.x
