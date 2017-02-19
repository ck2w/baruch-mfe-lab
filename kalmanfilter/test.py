#!/home/quan/anaconda3/bin/python

from kalman import *
#from pykalman import KalmanFilter
from config import *
import numpy as np
import scipy as sp
import pandas as pd
from scipy import optimize

YIELDS_DATA_FILE = "./treasury_yields.csv"
H15T = "H15T"

def printParams(params):

    print("A11 =", params[I_A11])
    print("A12 =", params[I_A12])
    print("A13 =", params[I_A13])
    print("A21 =", params[I_A21])
    print("A22 =", params[I_A22])
    print("A23 =", params[I_A23])
    print("A31 =", params[I_A31])
    print("A32 =", params[I_A32])
    print("A33 =", params[I_A33])

    print("Mu1 =", params[I_MU1])
    print("Mu2 =", params[I_MU2])
    print("Mu3 =", params[I_MU3])

    print("Q11 =", params[I_Q11])
    print("Q12 =", params[I_Q12])
    print("Q13 =", params[I_Q13])
    print("Q22 =", params[I_Q22])
    print("Q23 =", params[I_Q23])
    print("Q33 =", params[I_Q33])

    print("Lambda =", params[I_LAMBDA])
    
    for n in range(NUM_OBSERVATIONS):
        print("H" + str(n) + str(n), "=", params[I_LAMBDA+n+1])


def loglike(params, data):

    tm = transitionMatrix(params)
    om = observationMatrix(params)
    tc = transitionCovariance(params)
    oc = observationCovariance(params)
    to = transitionOffset(params)
    oo = observationOffset(params) 

    ism = np.zeros(NUM_FACTORS)
    isc = np.identity(NUM_FACTORS) 

    #try:
    print("parameters:")
    printParams(params)
    print("transitionMatrix:")
    print(transitionMatrix(params))
    print("transitionCovariance:")
    print(transitionCovariance(params))
    print("transitionOffset:")
    print(transitionOffset(params))
    print("observationMatrix:")
    print(observationMatrix(params))
    print("observationCovariance:")
    print(observationCovariance(params))
    print("observationOffset:")
    print(observationOffset(params))
    kf = KalmanFilter(transition_matrices = tm,
                      observation_matrices = om,
                      transition_covariance = tc,
                      observation_covariance = oc,
                      transition_offsets = to,
                      observation_offsets = oo,
                      initial_state_mean = ism,
                      initial_state_covariance = isc)

    v = -kf.loglikelihood(data)
    print("print: ", v)
    return v 

    #except:
    #    printParams(params)


if __name__ == "__main__":

    # observation data
    yields_pd = pd.read_csv(YIELDS_DATA_FILE, index_col='Date', parse_dates=['Date'])
    maturityList = [ H15T + maturity for maturity in OBS ]
    yields_pd = yields_pd[maturityList]
    X = yields_pd.as_matrix()

    params = np.zeros(NUM_PARAMS)
    bounds = [None]*NUM_PARAMS
   
    # parameter initial estimates
    params[I_A11] = -2.0
    params[I_A12] = 1.0
    params[I_A13] = -1.0

    params[I_A21] = 1.0
    params[I_A22] = -2.0
    params[I_A23] = 1.0

    params[I_A31] = -1.0
    params[I_A32] = 1.0
    params[I_A33] = 2.0

    params[I_MU1] = -10.0
    params[I_MU2] = -5.0
    params[I_MU3] = 5.0

    params[I_Q11] = 0.1 
    params[I_Q12] = 0.0
    params[I_Q13] = 0.0
    params[I_Q22] = 0.1
    params[I_Q23] = 0.0
    params[I_Q33] = 1.0

    params[I_LAMBDA] = 0.5

    params[I_LAMBDA+1] = 1e-9
    params[I_LAMBDA+2] = 1.0
    params[I_LAMBDA+3] = 1.0
    params[I_LAMBDA+4] = 1.0
    params[I_LAMBDA+5] = 1.0
    params[I_LAMBDA+6] = 1e-9
    params[I_LAMBDA+7] = 1e-9
    params[I_LAMBDA+8] = 1.0

    printParams(params)
    loglike(params, X)


