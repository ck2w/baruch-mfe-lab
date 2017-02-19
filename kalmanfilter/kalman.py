#!/home/quan/anaconda3/bin/python

import math
import numpy as np
import scipy as sp
from scipy import linalg

class KalmanFilter:

    def __init__(self, transition_matrices, observation_matrices,
                 transition_covariance, observation_covariance,
                 transition_offsets, observation_offsets,
                 initial_state_mean, initial_state_covariance):

        self.G = transition_matrices 
        self.F = observation_matrices 
        self.W = transition_covariance 
        self.V = observation_covariance
        self.c = transition_offsets
        self.d = observation_offsets
        self.mean = initial_state_mean
        self.cov = initial_state_covariance
        self.FT = np.transpose(self.F)
        self.GT = np.transpose(self.G)

    def filter(self, X):

        nT = len(X)
        if nT > 0:
            logLikelihood = 0
            nContracts = len(X[0])

            n_dim_state = len(self.c)

            predicted_state_means = np.zeros((nT, n_dim_state))
            predicted_state_covariances = np.zeros((nT, n_dim_state, n_dim_state))
            filtered_state_means = np.zeros((nT, n_dim_state))
            filtered_state_covariances = np.zeros((nT, n_dim_state, n_dim_state))

            for i in range(nT):
                # Prediction step:
                # a: prediction mean
                # R: prediction covariance
                # Prediction is done starting from the second step
                # according to the implementation of pykalman
                #print("G:")
                #print(self.G)
                #print("cov:")
                #print(self.cov)
                #print("W:")
                #print(self.W)
                if i == 0:
                    a = self.mean
                    R = self.cov
                else:
                    a = self.c + np.dot(self.G, self.mean)
                    R = np.linalg.multi_dot([self.G,self.cov,self.GT]) + self.W
                    # symmetrization
                    R = 0.5*(R + R.T)

                predicted_state_means[i] = a
                predicted_state_covariances[i] = R

                # Correction step
                u = X[i] - (self.d + np.dot(self.F,a))
                Q = np.linalg.multi_dot([self.F,R,self.FT]) + self.V
                # symmetrization
                Q = 0.5*(Q + Q.T)
                FR = np.dot(self.F,R)
                RF = np.transpose(FR)
                #print("Q:")
                #print(Q)
                #print("F:")
                #print(self.F)
                #print("R:")
                #print(R)
                #QiFR = np.dot(sp.linalg.pinv(Q),FR)
                QiFR = sp.linalg.solve(Q,FR,sym_pos=True)
                AQA = np.dot(RF,QiFR)
                Qiu = sp.linalg.solve(Q,u,sym_pos=True)
                self.mean = a + np.dot(RF,Qiu)
                self.cov = R - AQA

                filtered_state_means[i] = self.mean
                filtered_state_covariances[i] = self.cov

                # Compute log-likelihood
                (sign, logDet) = np.linalg.slogdet(Q)
                logLikelihood += logDet
                uQiu = np.dot(np.transpose(u),Qiu)
                logLikelihood += uQiu

            logLikelihood = -0.5*logLikelihood - 0.5*nT*nContracts*math.log(2*math.pi)
            return (self.mean, self.cov, logLikelihood,
                    predicted_state_means, predicted_state_covariances,
                    filtered_state_means, filtered_state_covariances)

        return (self.mean, self.cov, None, None, None, None, None)

    def loglikelihood(self, X):
        return self.filter(X)[2]

    def smooth(self, X):
            #filtered_state_means, filtered_state_covariances,
            #predicted_state_means, predicted_state_covariances):

        (_,_,_,predicted_state_means, predicted_state_covariances,
         filtered_state_means, filtered_state_covariances ) = self.filter(X)

        n_timesteps, n_dim_state = filtered_state_means.shape

        smoothed_state_means = np.zeros((n_timesteps, n_dim_state))
        smoothed_state_covariances = np.zeros((n_timesteps, n_dim_state, n_dim_state))
        kalman_smoothing_gains = np.zeros((n_timesteps - 1, n_dim_state, n_dim_state))

        smoothed_state_means[-1] = filtered_state_means[-1]
        smoothed_state_covariances[-1] = filtered_state_covariances[-1]

        for t in reversed(range(n_timesteps - 1)):
            filtered_state_mean = filtered_state_means[t] 
            filtered_state_covariance = filtered_state_covariances[t]
            predicted_state_mean = predicted_state_means[t + 1]
            predicted_state_covariance = predicted_state_covariances[t + 1]
            next_smoothed_state_mean = smoothed_state_means[t + 1]
            next_smoothed_state_covariance = smoothed_state_covariances[t + 1]

            kalman_smoothing_gain = (
                np.dot(filtered_state_covariance,
                       np.dot(self.GT,
                              sp.linalg.pinv(predicted_state_covariance)))
            )

            smoothed_state_mean = (
                filtered_state_mean
                + np.dot(kalman_smoothing_gain,
                         next_smoothed_state_mean - predicted_state_mean)
            )

            smoothed_state_covariance = (
                filtered_state_covariance
                + np.dot(kalman_smoothing_gain,
                         np.dot(
                             (next_smoothed_state_covariance - predicted_state_covariance),
                             kalman_smoothing_gain.T
                         ))
            )

            smoothed_state_means[t] = smoothed_state_mean
            smoothed_state_covariances[t] = smoothed_state_covariance
            kalman_smoothing_gains[t] = kalman_smoothing_gain

        return (smoothed_state_means, smoothed_state_covariances, kalman_smoothing_gains)


if __name__ == "__main__":

    tm = [[1, 1], [0, 1]]
    om = [[0.1, 0.5], [-0.3, 0.0]]
    tc = [[1,0], [0,1]]
    oc = [[1,0], [0,1]]
    ism = [0,0]
    isc = [[1,0], [0,1]]
    to = [0,0]
    oo = [0,0]

    kf = KalmanFilter(transition_matrices = tm, #[tm,tm,tm], 
                      observation_matrices = om, #[om,om,om],
                      transition_covariance = tc,
                      observation_covariance = oc,
                      transition_offsets = to, #[to,to,to],
                      observation_offsets = oo, #[oo,oo,oo],
                      initial_state_mean = ism,
                      initial_state_covariance = isc)

    measurements = np.asarray([[1,0], [0,0], [0,1]])  # 3 observations

    #res = kf.filter(measurements)
    #print(res[0])
    #print(res[1])

    res = kf.loglikelihood(measurements)
    print(res)
