#!/home/quan/anaconda3/bin/python

import sys
from optparse import OptionParser

from kalman import *
#from pykalman import KalmanFilter
from config import *
import numpy as np
import scipy as sp
import pandas as pd
from scipy import optimize

YIELDS_DATA_FILE = "./monthly_yields.csv"
H15T = "H15T"

def printConfig(params):

    print("Transition matrix:")
    print(transitionMatrix(params))
    print("Transition covariance:")
    print(transitionCovariance(params))
    print("Transition offset:")
    print([params[I_MU1], params[I_MU2], params[I_MU3]])
    print("Lambda =", params[I_LAMBDA])
    print("Observation covariance:") 
    for n in range(NUM_OBSERVATIONS):
        print("H" + str(n) + str(n), "=", params[I_LAMBDA+n+1])

def printParams(params):

    print("params[I_A11] =", params[I_A11])
    print("params[I_A12] =", params[I_A12])
    print("params[I_A13] =", params[I_A13])
    print("params[I_A21] =", params[I_A21])
    print("params[I_A22] =", params[I_A22])
    print("params[I_A23] =", params[I_A23])
    print("params[I_A31] =", params[I_A31])
    print("params[I_A32] =", params[I_A32])
    print("params[I_A33] =", params[I_A33])

    print("params[I_MU1] =", params[I_MU1])
    print("params[I_MU2] =", params[I_MU2])
    print("params[I_MU3] =", params[I_MU3])

    print("params[I_Q11] =", params[I_Q11])
    print("params[I_Q12] =", params[I_Q12])
    print("params[I_Q13] =", params[I_Q13])
    print("params[I_Q22] =", params[I_Q22])
    print("params[I_Q23] =", params[I_Q23])
    print("params[I_Q33] =", params[I_Q33])

    print("params[I_LAMBDA] =", params[I_LAMBDA])
    for n in range(NUM_OBSERVATIONS):
        print("params[I_LAMBDA+" + str(n+1) + "]=", params[I_LAMBDA+n+1])

CALL_COUNT = 0
PRINT_PERIOD = 100
def loglike(params, data):

    global CALL_COUNT
    global PRINT_PERIOD

    tm = transitionMatrix(params)
    om = observationMatrix(params)
    tc = transitionCovariance(params)
    oc = observationCovariance(params)
    to = transitionOffset(params)
    oo = observationOffset(params) 

    ism = np.zeros(NUM_FACTORS)
    isc = np.identity(NUM_FACTORS) 
    
    if CALL_COUNT % (NUM_PARAMS * PRINT_PERIOD) == 0:
        printParams(params)
    CALL_COUNT += 1

    #printParams(params)
    try:
        kf = KalmanFilter(transition_matrices = tm,
                          observation_matrices = om,
                          transition_covariance = tc,
                          observation_covariance = oc,
                          transition_offsets = to,
                          observation_offsets = oo,
                          initial_state_mean = ism,
                          initial_state_covariance = isc)

        v = -kf.loglikelihood(data)
        #print("Normal!!!")
    except:
        #printParams(params)
        v = 0
        #print("Exception!!!")

    return v 

def get_commandline():

    parser = OptionParser()
    parser.add_option("-t", action='store_true', dest="test")

    options, args = parser.parse_args()
    config = vars(options)
    return config

if __name__ == "__main__":

    commandline = get_commandline()
    isTest = commandline["test"]

    # observation data
    yields_pd = pd.read_csv(YIELDS_DATA_FILE, index_col='Date', parse_dates=['Date'])
    maturityList = [ H15T + maturity for maturity in OBS ]
    yields_pd = yields_pd[maturityList]
    X = yields_pd.as_matrix()

    params = np.zeros(NUM_PARAMS)
    bounds = [None]*NUM_PARAMS
   
    # parameter initial estimates
    params[I_A11] = 0.911100048321
    params[I_A12] = 0.00697032966505
    params[I_A13] = 0.0102672049677
    params[I_A21] = 0.258084542543
    params[I_A22] = 0.949707579317
    params[I_A23] = -0.0113413735843
    params[I_A31] = 0.0366336489751
    params[I_A32] = -0.094132012834
    params[I_A33] = 0.918088154256
    params[I_MU1] = 0.142199062278
    params[I_MU2] = 3.75595015217
    params[I_MU3] = -3.18174502981
    params[I_Q11] = 0.0473923092582
    params[I_Q12] = -0.0752755426836
    params[I_Q13] = -0.170428624113
    params[I_Q22] = 0.319996395733
    params[I_Q23] = -0.245021941271
    params[I_Q33] = 0.350440707609
    params[I_LAMBDA] = 0.508823195391
    params[I_LAMBDA+1]= 0.00142827136831
    params[I_LAMBDA+2]= 0.000179379237594
    params[I_LAMBDA+3]= 0.00289372767108
    params[I_LAMBDA+4]= 0.0042746935139
    params[I_LAMBDA+5]= 0.000601294697819
    params[I_LAMBDA+6]= 0.00709819089612
    params[I_LAMBDA+7]= 0.0185728543035
    params[I_LAMBDA+8]= 0.014184454474
    params[I_LAMBDA+9]= 3.06642765843e-08
    params[I_LAMBDA+10]= 0.0240409764867
    params[I_LAMBDA+11]= 0.0749994608772
    
    #for n in range(NUM_OBSERVATIONS):
    #    params[I_LAMBDA+n+1] = 0.01

    if isTest:

        tm = transitionMatrix(params)
        om = observationMatrix(params)
        tc = transitionCovariance(params)
        oc = observationCovariance(params)
        to = transitionOffset(params)
        oo = observationOffset(params) 

        ism = np.zeros(NUM_FACTORS)
        isc = np.identity(NUM_FACTORS) 

        kf = KalmanFilter(transition_matrices = tm,
                          observation_matrices = om,
                          transition_covariance = tc,
                          observation_covariance = oc,
                          transition_offsets = to,
                          observation_offsets = oo,
                          initial_state_mean = ism,
                          initial_state_covariance = isc)

        res = kf.smooth(X)
        df = pd.DataFrame(res[0])
        df.to_csv("kalmanSmooth.csv")
    else:
        # parameter bounds
        bounds[I_A11] = (0.5,1)
        bounds[I_A12] = (-0.5,0.8)
        bounds[I_A13] = (-0.4,0.7)
        bounds[I_A21] = (-0.3,0.6)
        bounds[I_A22] = (0.5,1)
        bounds[I_A23] = (-0.6,0.5)
        bounds[I_A31] = (-0.7,0.4)
        bounds[I_A32] = (-0.8,0.3)
        bounds[I_A33] = (0.5,1)

        bounds[I_MU1] = (1e-6,10)
        bounds[I_MU2] = (1e-6,10) 
        bounds[I_MU3] = (-5,5)

        bounds[I_Q11] = (0.01,1)
        bounds[I_Q12] = (-1,1) 
        bounds[I_Q13] = (-1,1)
        bounds[I_Q22] = (0.01,1)
        bounds[I_Q23] = (-1,1)
        bounds[I_Q33] = (0.01,1)

        bounds[I_LAMBDA] = (0.2,0.9)
        
        for n in range(NUM_OBSERVATIONS):
            bounds[I_LAMBDA+n+1] = (0.0,0.1)

        res = sp.optimize.minimize(loglike, params, args=X, method='L-BFGS-B', 
                                   bounds=bounds, options={'disp': True})

        printParams(res.x)

