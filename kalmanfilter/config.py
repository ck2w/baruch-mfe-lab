#!/home/quan/anaconda3/bin/python

import math
import numpy as np

# full list:
OBS = ['1M', '3M', '6M', '1Y', '2Y', '3Y', '5Y', '7Y', '10Y', '20Y', '30Y']
#OBS = ['3M', '6M', '1Y', '2Y', '3Y', '5Y', '7Y', '10Y']

# states
NUM_FACTORS = 3
NUM_STATE_PARAMS = int(3*NUM_FACTORS*(NUM_FACTORS+1)/2)

I_A11 = 0
I_A12 = 1
I_A13 = 2
I_A21 = 3
I_A22 = 4
I_A23 = 5
I_A31 = 6
I_A32 = 7
I_A33 = 8
I_MU1 = 9
I_MU2 = 10
I_MU3 = 11
I_Q11 = 12
I_Q12 = 13
I_Q13 = 14
I_Q22 = 15
I_Q23 = 16
I_Q33 = 17

# observations 

def getMaturityList(obs):

    maturityList = []
    for maturity in obs:
        timeUnit = maturity[-1:]
        numOfTimeUnits = float(maturity[:-1])
        maturityLength = MATURITY_MAP[timeUnit]*numOfTimeUnits
        maturityList.append(maturityLength)
    return maturityList

I_LAMBDA = NUM_STATE_PARAMS

MATURITY_MAP = { 'M': 1.0/12.0, 'Y': 1.0 }

MATURITY_LIST = getMaturityList(OBS)

NUM_OBSERVATIONS = len(MATURITY_LIST)
NUM_OBS_PARAMS = 1+NUM_OBSERVATIONS

NUM_PARAMS = NUM_STATE_PARAMS + NUM_OBS_PARAMS

def transitionMatrix(params):

    A = np.array([[params[I_A11], params[I_A12], params[I_A13]],
                  [params[I_A21], params[I_A22], params[I_A23]],
                  [params[I_A31], params[I_A32], params[I_A33]]])
    return A

def transitionCovariance(params):
    
    U = np.array([[params[I_Q11], params[I_Q12], params[I_Q13]],
                  [          0.0, params[I_Q22], params[I_Q23]],
                  [          0.0,           0.0, params[I_Q33]]])

    return np.dot(U.T,U) 

def transitionOffset(params):

    mu = np.array([params[I_MU1], params[I_MU2], params[I_MU3]])
    return np.dot(np.identity(NUM_FACTORS) - transitionMatrix(params), mu)

def observationMatrix(params):

    F = np.ones((NUM_OBSERVATIONS, NUM_FACTORS))

    lam = params[I_LAMBDA]
    if abs(lam) > 1e-6:
        for n in range(NUM_OBSERVATIONS):
            tau = MATURITY_LIST[n]
            F[n][1] = 1.0 - (1.0-math.exp(-tau*lam))/(tau*lam)
            F[n][2] = (1.0-math.exp(-tau*lam))/(tau*lam) - math.exp(-tau*lam)
    else:
        for n in range(NUM_OBSERVATIONS):
            tau = MATURITY_LIST[n]
            F[n][1] = 0.5*lam*tau
            F[n][2] = 1.0 - math.exp(-tau*lam)

    return F

def observationCovariance(params):

    V = np.zeros((NUM_OBSERVATIONS, NUM_OBSERVATIONS))
    for n in range(NUM_OBSERVATIONS):
        V[n][n] = params[I_LAMBDA+n+1]
    
    return V

def observationOffset(params):
    
    return np.zeros(NUM_OBSERVATIONS)

if __name__ == "__main__":

    print("num of parameters = ", NUM_PARAMS)
    params = np.zeros(NUM_PARAMS)
    params[I_A11] = 1.1
    params[I_A12] = 1.2
    params[I_A13] = 1.3
    params[I_A21] = 2.1
    params[I_A22] = 2.2
    params[I_A23] = 2.3
    params[I_A31] = 3.1
    params[I_A32] = 3.2
    params[I_A33] = 3.3
    params[I_MU1] = 1.0
    params[I_MU2] = 2.0
    params[I_MU3] = 3.0
    params[I_Q11] = 11.1
    params[I_Q12] = 12.2
    params[I_Q13] = 13.3
    params[I_Q22] = 22.2
    params[I_Q23] = 22.3
    params[I_Q33] = 33.3
    params[I_LAMBDA] = 0.456

    for n in range(NUM_OBSERVATIONS):
        params[I_LAMBDA+n+1] =0.333

    print("parameters:")
    print(params)
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


