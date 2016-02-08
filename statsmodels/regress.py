#!/home/quan/anaconda3/bin/python

import sys
from optparse import OptionParser
import pandas as pd
import numpy as np
from statsmodels import regression, stats
import statsmodels
import statsmodels.formula.api as smf
import math
from datetime import datetime

FREQUENCY = 1
LOG_PREFIX = "log_"
LOG_DIFF_PREFIX = "ret_"
LAG_PREFIX = "_lag_"

def usage():
    print("regress.py -y <dependent variable ts> -x <indepdent variable ts> -l #-# -r YYYYMMDD-YYYYMMDD")
    print("where: -l, --lag, #-#, default to 0-0")
    print("       -r, --range, YYYYMMDD-YYYYMMDD")
    sys.exit(1)

def get_commandline():

    if len(sys.argv) < 5: usage()

    parser = OptionParser()
    parser.add_option("-y", "--dependent")
    parser.add_option("-x", "--independent", action="append", default=[])
    parser.add_option("-l", "--lag", action="append", default=[])
    parser.add_option("-r", "--range", action="append")
    parser.add_option("-n", "--nwlag")

    options, args = parser.parse_args()
    assert(options.dependent is not None)
    assert(options.independent is not None)
    assert(len(options.lag) <= len(options.independent))

    config = vars(options)
    return config

if __name__ == "__main__":

    commandline = get_commandline()
    dependent_variable_file = commandline["dependent"]
    independent_variable_files = commandline["independent"]
    lags = commandline["lag"]
    number_of_specified_lags = len(lags)
    max_lag = 1

    dependent_variable_df = pd.read_csv(dependent_variable_file, index_col='Date', parse_dates=['Date'])
    dependent_variable_df.columns = ["y"]
    dependent_variable_df[LOG_PREFIX + "y"] = dependent_variable_df["y"].apply(math.log)
    dependent_variable_df[LOG_DIFF_PREFIX + "y"] = dependent_variable_df[LOG_PREFIX + "y"].diff(FREQUENCY)
    regression_formula = LOG_DIFF_PREFIX + "y ~ 1"
    datafile_summary = dependent_variable_file
    summary_header = "Date Range"
    
    independent_variable_dfs = []
    for index, independent_variable_file in enumerate(independent_variable_files):
       
        datafile_summary += "|" + independent_variable_file
        col = "x" + str(index)
        
        lag_head = 0
        lag_tail = 0
        if index < number_of_specified_lags:
            lag_head_tail_tag = lags[index]
            lag_head_tail_pair = lag_head_tail_tag.split("-")
            lag_head = int(lag_head_tail_pair[0])
            lag_tail = int(lag_head_tail_pair[1])
            assert(lag_head <= lag_tail)
            if lag_tail > max_lag: max_lag = lag_tail

        independent_variable_df = pd.read_csv(independent_variable_file, index_col='Date', parse_dates=['Date'])
        independent_variable_df.columns = [col]
        # log return
        independent_variable_df[LOG_PREFIX + col] = independent_variable_df[col].apply(math.log)
        independent_variable_df[LOG_DIFF_PREFIX + col + LAG_PREFIX + "0"] = independent_variable_df[LOG_PREFIX + col].diff(FREQUENCY)

        # lags
        if lag_head == 0:
            regression_formula += " + " + LOG_DIFF_PREFIX + col + LAG_PREFIX + "0"
            summary_header += ", " + LOG_DIFF_PREFIX + col + LAG_PREFIX + "0-beta , std-err , p-value"
            lag_head = 1

        for lag in range(lag_head, lag_tail+1):
            independent_variable_df[LOG_DIFF_PREFIX + col + LAG_PREFIX + str(lag)]\
                = independent_variable_df[LOG_DIFF_PREFIX + col + LAG_PREFIX + "0"].shift(lag)
            regression_formula += " + " + LOG_DIFF_PREFIX + col + LAG_PREFIX + str(lag)
            summary_header += ", " + LOG_DIFF_PREFIX + col + LAG_PREFIX + str(lag) + "-beta , std-err , p-value"
        
        independent_variable_dfs.append(independent_variable_df)

    summary_header += " , r-squared , r-squared-adjusted"

    # merge all data frames by date
    total_df = dependent_variable_df
    for independent_variable_df in independent_variable_dfs:
        total_df = total_df.merge(independent_variable_df, how='inner', right_index=True, left_index=True)
    
    total_df = total_df[max_lag:]

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
  
    if len(date_ranges) > 0:
        print("Data files: ", datafile_summary)
        #print("Regression formula: ", regression_formula)
        print(summary_header)

    # slicing by date range
    for date_range in date_ranges:
        start_date = date_range[0]
        end_date = date_range[1]
        partial_df = total_df[start_date:end_date]

        newey_west_lag = commandline["nwlag"]
        if newey_west_lag is not None: 
            newey_west_lag = int(newey_west_lag)
        else:
            newey_west_lag = 10

        reg = smf.ols(regression_formula, data=partial_df).fit(cov_type='HAC',cov_kwds={'maxlags': newey_west_lag})
        
        # print summary
        print(start_date.date(), end="~"),
        print(end_date.date(), end="")
        slopes = reg.params[1:]
        pvalues = reg.pvalues[1:]
        stderrs = reg.bse[1:]

        for slope, err, pvalue in zip(slopes, stderrs, pvalues):
            print(", {:1.6f}".format(slope), ", {:1.6f}".format(err), ", {:1.6f}".format(pvalue), end=" ")

        print(", {:1.6f}".format(reg.rsquared), ", {:1.6f}".format(reg.rsquared_adj))

        #print(reg.summary())
