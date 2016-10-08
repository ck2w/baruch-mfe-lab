#!/home/quan/anaconda3/bin/python

import pandas as pd
import numpy as np
from statsmodels import regression, stats
import statsmodels
import statsmodels.formula.api as smf
import math

# First, remove all the white spaces
# in the column header before processing

time_series_file = "financials2016-8.csv"
OFFSET=-1
LOG_PREFIX="log_"
LOG_RETURN="log_ret_"
PCT_RETURN="pct_ret_"


time_series_pd = pd.read_csv(time_series_file, index_col='Date', parse_dates=['Date'])

for col_name in time_series_pd.columns:
    #time_series_pd[LOG_PREFIX + col_name] = time_series_pd[col_name].diff(OFFSET)/time_series_pd[col_name].shift(-1)
    time_series_pd[LOG_PREFIX + col_name] = time_series_pd[col_name].apply(math.log)
    time_series_pd[LOG_RETURN + col_name] = time_series_pd[LOG_PREFIX + col_name].diff(OFFSET)

# drop the last row
time_series_pd = time_series_pd.ix[:-1]


bb = time_series_pd.drop('MS', 1)
bb = bb.drop('BAC', 1)
bb = bb.drop('BCS', 1)
bb = bb.drop('CS', 1)
bb = bb.drop('GS', 1)
bb = bb.drop('JPM', 1)
bb = bb.drop('RBS', 1)
bb = bb.drop('UBS', 1)
bb = bb.drop('log_MS', 1)
bb = bb.drop('log_BAC', 1)
bb = bb.drop('log_BCS', 1)
bb = bb.drop('log_CS', 1)
bb = bb.drop('log_GS', 1)
bb = bb.drop('log_JPM', 1)
bb = bb.drop('log_RBS', 1)
bb = bb.drop('log_UBS', 1)

print(bb.cov())


regression_formula = "log_ret_MS ~ log_ret_BAC + log_ret_BCS + log_ret_CS\
                      + log_ret_GS + log_ret_JPM + log_ret_RBS + log_ret_UBS" 

reg = smf.ols(regression_formula, data=time_series_pd).fit()

print(reg.summary())

residual_error=0
for resid in reg.resid:
    residual_error += resid*resid
print("Squared residual error = ", residual_error)
print("Residual error = ", math.sqrt(residual_error))

#print(time_series_pd.pct_change())
#print(time_series_pd.pct_change()["JPM"])
#print(time_series_pd.head())
#print(time_series_pd.tail())
#reg = smf.ols(regression_formula, data=partial_df).fit(cov_type='HAC',cov_kwds={'maxlags': newey_west_lag})
#time_series_df[LOG_RETURN + col_name] = time_series_df[col_name].apply(math.log)
