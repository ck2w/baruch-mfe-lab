#!/usr/bin/ksh

lag=$1
sp500=$2

./regress.sh AIRL $lag $sp500
./regress.sh AUTO $lag $sp500
./regress.sh OILG $lag $sp500
./regress.sh CHEM $lag $sp500
./regress.sh FDPR $lag $sp500
./regress.sh METL $lag $sp500
./regress.sh PHAR $lag $sp500
./regress.sh CSTE $lag $sp500
./regress.sh TIRE $lag $sp500
./regress.sh STEL $lag $sp500
./regress.sh AGRI $lag $sp500
