#!/usr/bin/ksh

industry_code=$1
lag_tag=$2
spx_flag=$3
industry_file=""
sp500_option=""

DATA_FOLDER="./data"
OUTPUT_FOLDER="./output"

sp500_file=$DATA_FOLDER"/sp500.csv"
sp1500_file=$DATA_FOLDER"/sp1500.csv"
wti_spot_file=$DATA_FOLDER"/wti.csv"
brent_spot_file=$DATA_FOLDER"/brent.csv"
cl1_file=$DATA_FOLDER"/CL1T.csv"
cl3_file=$DATA_FOLDER"/CL3T.csv"
cl6_file=$DATA_FOLDER"/CL6T.csv"
co1_file=$DATA_FOLDER"/CO1T.csv"
co3_file=$DATA_FOLDER"/CO3T.csv"
co6_file=$DATA_FOLDER"/CO6T.csv"

AIRL=AIRL #airline
AMGT=AMGT #asset management
AGRI=AGRI #agriculture
AUCO=AUCO #automobil and components
AUTO=AUTO #automobil sub
BANK=BANK #bank
CBNK=CBNK #commercial banks
CFIN=CFIN #consumer finance
CHEM=CHEM #chemical
CODU=CODU #consumer durables
COMS=COMS #commercial professional services
COND=COND #consumer durables
CONS=CONS #consumer services
CPGS=CPGS #capital goods
CSTE=CSTE #construction engineering
CSTM=CSTM #construction materials
DBNK=DBNK #diversified banks
DCAP=DCAP #diversified capital markets
DIVF=DIVF #diversified financial
ENRS=ENRS #energy
FDBT=FDBT #food beverage & tobacco
FDPR=FDPR #food sub
FDSR=FDSR #food & staples
FINL=FINL #finacials
HCES=HCES #health care equipment
HLTH=HLTH #health care
HOTR=HOTR #consumer services
HOUS=HOUS #household and personal products
INBK=INBK #investment bank
INDU=INDU #industrials
INFT=INFT #technology
INSB=INSB #insurance broker
INSU=INSU #insurance
LIFE=LIFE #life and health insurance
MATR=MATR #materials
MEDA=MEDA #media
METL=METL #metals and mining
MLIN=MLIN #multiline insurance
MSEC=MSEC #multisector holdings
OILD=OILD #oil and gas drilling
OILG=OILG #oil and gas general
OILP=OILP #oil and gas exploration and production
OILR=OILR #oil and gas refining and marketing
OILE=OILE #oil and gas equipment and services
OGST=OGST #oil and gas storage and transportation
PHRM=PHRM #pharmaceutical biotech
PHAR=PHAR #pharmaceutical sub
PROP=PROP #property and casualty insurance
RBNK=RBNK #regional banks
REAL=REAL #real estate
REIN=REIN #reinsurance
RETL=RETL #retail
SFTW=SFTW #software
SPFN=SPFN #specialized finance
SSEQ=SSEQ #semiconductor
STEL=STEL #steel
TECH=TECH #technology hardware
TELS=TELS #telecommunication
THMF=THMF #thrifts and mortgage finance
TIRE=TIRE #tires and rubbers
TRAN=TRAN #transportation
UTIL=UTIL #utilities

if [[ -z $industry_code ]]; then
    echo "Industry code missing. Exit."
    exit
fi

industry_code_prefix=`echo $industry_code|cut -b 1-4`
if [[ $industry_code_prefix=="SP15" ]]; then
    industry_code=`echo $industry_code|cut -b 5-8`
    echo "Running regression for S&P-1500" $industry_code
    industry_file=$industry_code".csv"
else
    echo "Running regression for S&P-500" $industry_code
    if [[ $industry_code == $BANK ]]; then
        industry_file="bank.csv"
    elif [[ $industry_code == $RETL ]]; then
        industry_file="retail.csv"
    elif [[ $industry_code == $AUCO ]]; then
        industry_file="autocomponents.csv"
    elif [[ $industry_code == $MEDA ]]; then
        industry_file="media.csv"
    elif [[ $industry_code == $INSU ]]; then
        industry_file="insurance.csv"
    elif [[ $industry_code == $TRAN ]]; then
        industry_file="transportation.csv"
    elif [[ $industry_code == $ENRS ]]; then
        industry_file="energy.csv"
    elif [[ $industry_code == $REAL ]]; then
        industry_file="realestate.csv"
    elif [[ $industry_code == $SFTW ]]; then
        industry_file="software.csv"
    elif [[ $industry_code == $TECH ]]; then
        industry_file="hardware.csv"
    elif [[ $industry_code == $PHRM ]]; then
        industry_file="pharmbio.csv"
    elif [[ $industry_code == $DIVF ]]; then
        industry_file="financial.csv"
    elif [[ $industry_code == $CPGS ]]; then
        industry_file="capitalgoods.csv"
    elif [[ $industry_code == $UTIL ]]; then
        industry_file="utilities.csv"
    elif [[ $industry_code == $FDSR ]]; then
        industry_file="foodstaples.csv"
    elif [[ $industry_code == $HCES ]]; then
        industry_file="healthequip.csv"
    elif [[ $industry_code == $CODU ]]; then
        industry_file="consumerdurables.csv"
    elif [[ $industry_code == $HOTR ]]; then
        industry_file="consumerservices.csv"
    elif [[ $industry_code == $HOUS ]]; then
        industry_file="housepersonal.csv"
    elif [[ $industry_code == $COMS ]]; then
        industry_file="commercialprofessional.csv"
    elif [[ $industry_code == $FDBT ]]; then
        industry_file="foodbeverage.csv"
    elif [[ $industry_code == $TELS ]]; then
        industry_file="telecom.csv"
    elif [[ $industry_code == $MATR ]]; then
        industry_file="materials.csv"
    elif [[ $industry_code == $SSEQ ]]; then
        industry_file="semiconductor.csv"
    elif [[ $industry_code == $AIRL ]]; then
        industry_file="airline.csv"
    elif [[ $industry_code == $AUTO ]]; then
        industry_file="automobile.csv"
    elif [[ $industry_code == $OILG ]]; then
        industry_file="oilgas.csv"
    elif [[ $industry_code == $CHEM ]]; then
        industry_file="chemical.csv"
    elif [[ $industry_code == $FDPR ]]; then
        industry_file="food.csv"
    elif [[ $industry_code == $METL ]]; then
        industry_file="mining.csv"
    elif [[ $industry_code == $PHAR ]]; then
        industry_file="pharm.csv"
    elif [[ $industry_code == $CSTE ]]; then
        industry_file="construction.csv"
    elif [[ $industry_code == $TIRE ]]; then
        industry_file="rubber.csv"
    elif [[ $industry_code == $STEL ]]; then
        industry_file="steel.csv"
    elif [[ $industry_code == $AGRI ]]; then
        industry_file="agriculture.csv"
    else
        echo "Invalid industry code. exit."
        exit
    fi
fi

if [[ -z $spx_flag ]]; then
    sp500_option=""
elif [[ $spx_flag == "SP500" ]]; then
    sp500_option=" -x $sp500_file -l 0-0"
elif [[ $spx_flag == "SP1500" ]]; then
    sp500_option=" -x $sp1500_file -l 0-0"
else
    echo "Invalid command-line option. exit."
    exit
fi

industry_output=$OUTPUT_FOLDER/${industry_file%???}"out"
industry_file=$DATA_FOLDER/$industry_file

./regress.py -y $industry_file -x $wti_spot_file -l $lag_tag $sp500_option -r 20090201-20091231 \
                                                                           -r 20100101-20101231 \
                                                                           -r 20110101-20111231 \
                                                                           -r 20120101-20121231 \
                                                                           -r 20130101-20140731 \
                                                                           -r 20140801-20160127 > $industry_output

