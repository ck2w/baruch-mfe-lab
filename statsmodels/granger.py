#!/usr/bin/anaconda3/bin/python

import sys
from optparse import OptionParser
import pandas as pd
import numpy as np
import statsmodels.api as sm
import math
from datetime import datetime

from statsmodels.tsa.api import VAR
from scipy.signal import lfilter

from statsmodels.tsa.stattools import grangercausalitytests as granger

FREQUENCY = 1
LOG_PREFIX = "log_"
LOG_DIFF_PREFIX = "ret_"
LAG_PREFIX = "_lag_"

def usage():
    print("granger.py -y <dependent variable ts> -x <indepdent variable ts> -l <max lag> -r YYYYMMDD-YYYYMMDD")
    print("where: -l, --lag, default to 1")
    print("       -r, --range, YYYYMMDD-YYYYMMDD")
    sys.exit(1)

def get_commandline():

    if len(sys.argv) < 5: usage()

    parser = OptionParser()
    parser.add_option("-y", "--dependent")
    parser.add_option("-x", "--independent")
    parser.add_option("-l", "--lag")
    parser.add_option("-r", "--range", action="append")

    options, args = parser.parse_args()
    assert(options.dependent is not None)
    assert(options.independent is not None)

    config = vars(options)
    return config

if __name__ == "__main__":

    commandline = get_commandline()
    dependent_variable_file = commandline["dependent"]
    independent_variable_file = commandline["independent"]
    lag = int(commandline["lag"])

    dependent_variable_df = pd.read_csv(dependent_variable_file, index_col='Date', parse_dates=['Date'])
    dependent_variable_df.columns = ["y"]
    dependent_variable_df[LOG_PREFIX + "y"] = dependent_variable_df["y"].apply(math.log)
    dependent_variable_df[LOG_DIFF_PREFIX + "y"] = dependent_variable_df[LOG_PREFIX + "y"].diff(FREQUENCY)

    independent_variable_df = pd.read_csv(independent_variable_file, index_col='Date', parse_dates=['Date'])
    independent_variable_df.columns = ["x"]
    independent_variable_df[LOG_PREFIX + "x"] = independent_variable_df["x"].apply(math.log)
    independent_variable_df[LOG_DIFF_PREFIX + "x"] = independent_variable_df[LOG_PREFIX + "x"].diff(FREQUENCY)

    total_df = dependent_variable_df
    total_df = total_df.merge(independent_variable_df, how='inner', right_index=True, left_index=True)
    
    # max available date range
    first_date = total_df.index[0]
    last_date = total_df.index[-1]
    
    # set regression date range
    date_ranges = []
    date_range_tags = commandline["range"]
    if date_range_tags is None:
        date_ranges.append((first_date, last_date))
    else:
        for date_range in date_range_tags:
            start_end_date_pair = date_range.split("-")
            start_date = datetime.strptime(start_end_date_pair[0], '%Y%m%d')
            end_date = datetime.strptime(start_end_date_pair[1], '%Y%m%d')
            date_ranges.append((start_date, end_date))
  
    # slicing by date range
    for date_range in date_ranges:
        start_date = date_range[0]
        end_date = date_range[1]
        partial_df = total_df[start_date:end_date]
    
        mdata = partial_df[['ret_y', 'ret_x']] 

        results = granger(mdata, maxlag=lag)

        print("\n")
        print("\n")
        print("Date Range: ", start_date, " - ", end_date)
        print("###########################################################")

        for index in range(1, lag+1): 
            print("parameters ftest:     ", results[index][0]["params_ftest"])
            print("ssr chi2test:         ", results[index][0]["ssr_chi2test"])
            print("likelihood ratio test:", results[index][0]["lrtest"])
            print("ssr ftest:            ", results[index][0]["ssr_ftest"])
            print("--------------------------------------------------------------")
            print("restricted model:")
            reg = results[index][1][0]
            slopes = reg.params[1:]
            pvalues = reg.pvalues[1:]
            for slope, pvalue in zip(slopes, pvalues):
                print("{:1.6f}".format(slope), ", {:1.6f}".format(pvalue))
            
            print("--------------------------------------------------------------")
            
            print("full model:")
            reg = results[index][1][1]
            slopes = reg.params[1:]
            pvalues = reg.pvalues[1:]
            for slope, pvalue in zip(slopes, pvalues):
                print("{:1.6f}".format(slope), ", {:1.6f}".format(pvalue))
            print("==============================================================")




