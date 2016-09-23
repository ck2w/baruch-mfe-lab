#!/usr/bin/ksh

lag=$1
spx=$2

./regress.sh SP15COND $lag $spx
./regress.sh SP15CONS $lag $spx
./regress.sh SP15FINL $lag $spx
./regress.sh SP15HLTH $lag $spx
./regress.sh SP15INDU $lag $spx
./regress.sh SP15INFT $lag $spx
./regress.sh SP15BANK $lag $spx
./regress.sh SP15RETL $lag $spx
./regress.sh SP15AUCO $lag $spx
./regress.sh SP15MEDA $lag $spx
./regress.sh SP15INSU $lag $spx
./regress.sh SP15TRAN $lag $spx
./regress.sh SP15ENRS $lag $spx
./regress.sh SP15REAL $lag $spx
./regress.sh SP15SFTW $lag $spx
./regress.sh SP15TECH $lag $spx
./regress.sh SP15PHRM $lag $spx
./regress.sh SP15DIVF $lag $spx
./regress.sh SP15CPGS $lag $spx
./regress.sh SP15UTIL $lag $spx
./regress.sh SP15FDSR $lag $spx
./regress.sh SP15HCES $lag $spx
./regress.sh SP15CODU $lag $spx
./regress.sh SP15HOTR $lag $spx
./regress.sh SP15HOUS $lag $spx
./regress.sh SP15COMS $lag $spx
./regress.sh SP15FDBT $lag $spx
./regress.sh SP15TELS $lag $spx
./regress.sh SP15MATR $lag $spx
./regress.sh SP15SSEQ $lag $spx
./regress.sh SP15AIRL $lag $spx
./regress.sh SP15AUTO $lag $spx
./regress.sh SP15OILG $lag $spx
./regress.sh SP15OILP $lag $spx
./regress.sh SP15OILR $lag $spx
./regress.sh SP15OILE $lag $spx
./regress.sh SP15OGST $lag $spx
./regress.sh SP15OILD $lag $spx
./regress.sh SP15CHEM $lag $spx
./regress.sh SP15FDPR $lag $spx
./regress.sh SP15METL $lag $spx
./regress.sh SP15PHAR $lag $spx
./regress.sh SP15CSTE $lag $spx
./regress.sh SP15CSTM $lag $spx
./regress.sh SP15TIRE $lag $spx
./regress.sh SP15STEL $lag $spx
./regress.sh SP15AGRI $lag $spx

cat output/*.out > results_lag_"$lag"_"$spx"_daily.out
rm output/*.out
mv results_lag_"$lag"_"$spx"_daily.out output/
