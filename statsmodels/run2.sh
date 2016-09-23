#!/usr/bin/ksh

lag=$1
sp500=$2

./regress.sh BANK $lag $sp500
./regress.sh RETL $lag $sp500
./regress.sh AUCO $lag $sp500
./regress.sh MEDA $lag $sp500
./regress.sh INSU $lag $sp500
./regress.sh TRAN $lag $sp500
./regress.sh ENRS $lag $sp500
./regress.sh REAL $lag $sp500
./regress.sh SFTW $lag $sp500
./regress.sh TECH $lag $sp500
./regress.sh PHRM $lag $sp500
./regress.sh DIVF $lag $sp500
./regress.sh CPGS $lag $sp500
./regress.sh UTIL $lag $sp500
./regress.sh FDSR $lag $sp500
./regress.sh HCES $lag $sp500
./regress.sh CODU $lag $sp500
./regress.sh HOTR $lag $sp500
./regress.sh HOUS $lag $sp500
./regress.sh COMS $lag $sp500
./regress.sh FDBT $lag $sp500
./regress.sh TELS $lag $sp500
./regress.sh MATR $lag $sp500
./regress.sh SSEQ $lag $sp500
