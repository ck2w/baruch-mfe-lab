#!/usr/bin/ksh

industry_code=$1

DATA_FOLDER="./data"
OUTPUT_FOLDER="./output"

cl1_file=$DATA_FOLDER"/CL1T.csv"
cl3_file=$DATA_FOLDER"/CL3T.csv"
cl6_file=$DATA_FOLDER"/CL6T.csv"
co1_file=$DATA_FOLDER"/CO1T.csv"
co3_file=$DATA_FOLDER"/CO3T.csv"
co6_file=$DATA_FOLDER"/CO6T.csv"

if [[ -z $industry_code ]]; then
    echo "Industry code missing. Exit."
    exit
fi

industry_code_prefix=`echo $industry_code|cut -b 1-4`
if [[ $industry_code_prefix=="SP15" ]]; then
    industry_code=`echo $industry_code|cut -b 5-8`
    echo "Running crude oil correlation for S&P-1500" $industry_code
    industry_file=$industry_code".csv"
else
    echo "Industry code not supported!"
    exit
fi

industry_output=$OUTPUT_FOLDER/${industry_file%???}"corr"
industry_file=$DATA_FOLDER/$industry_file

./regress.py -y $industry_file -x $cl1_file -r 20140801-20160129 >> $industry_output
./regress.py -y $industry_file -x $cl3_file -r 20140801-20160129 >> $industry_output
./regress.py -y $industry_file -x $cl6_file -r 20140801-20160129 >> $industry_output
./regress.py -y $industry_file -x $co1_file -r 20140801-20160129 >> $industry_output
./regress.py -y $industry_file -x $co3_file -r 20140801-20160129 >> $industry_output
./regress.py -y $industry_file -x $co6_file -r 20140801-20160129 >> $industry_output

