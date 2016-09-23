#!/usr/bin/ksh

lag=$1
spx=$2

./regress.sh SP15THMF $lag $spx
./regress.sh SP15RBNK $lag $spx
./regress.sh SP15CFIN $lag $spx
./regress.sh SP15SPFN $lag $spx
./regress.sh SP15MSEC $lag $spx
./regress.sh SP15INBK $lag $spx
./regress.sh SP15AMGT $lag $spx
./regress.sh SP15INSB $lag $spx
#./regress.sh SP15DCAP $lag $spx
./regress.sh SP15MLIN $lag $spx
./regress.sh SP15REIN $lag $spx
./regress.sh SP15LIFE $lag $spx
./regress.sh SP15PROP $lag $spx
./regress.sh SP15CBNK $lag $spx
./regress.sh SP15DBNK $lag $spx

cat output/*.out > results_lag_"$lag"_"$spx"_daily.out
rm output/*.out
mv results_lag_"$lag"_"$spx"_daily.out output/
