#!/home/quan/anaconda3/bin/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats
from data_accessor import Data_accessor
from sklearn import ensemble
from sklearn.ensemble.gradient_boosting import GradientBoostingRegressor

C_FEATURES = [ #'x6', # Not very good as C. But even worse as D.
               #'x30', # me # key. Not bad as C.
               #'x25', # slightly better as a C than as a D.
               #'x13', # OK as a C.
               #'x28',  # OK as a C.
               #'x29', # OK as a C.
               #'x51', # key. OK as a C.
             ]

D_FEATURES = [ #'x2', # Very strange. Highly correlated, but makes it worse. Toxic.
               #'x13', # Not bad as a D. But slightly better as a C.
               #'x28', # Roughly the same as as a C.
               #'x29', # Roughpy the same as as a C.
               'x30', # OK. Not bad as D. Roughly the same as C.
               'x42', # Better as D feature than as C feature. Improves.
               'x46', # key
               'x51', # OK as a D. Slightly better than as a C.
             ]

Q_FEATURES = [ 'x0',
               'x17', 
               'x22', 
               'x49', 
               'x53', 
               'x61',
               #'x5', 
               #'x9', 
               #'x16', 
               #'x21', 
               #'x37', 
               #'x39', 
               #'x62',
               #'x63', 
               #'x64', 
               #'x65', 
               #'x66', 
               #'weight'
             ]

S_FEATURES = [ 
               #'x9',
               #'x16',
               #'x49',
               #'x53',
             ]

TIME = 'timestamp'
ID = 'id'
WEIGHT = 'weight'
OUTPUT = 'y'

MAX_DEPTHS = [ 6 ]
LEARNING_RATES = [ 0.015, 0.016, 0.017, 0.018]
SUBSAMPLES = [ 1.0 ]

N_ESTIMATORS = 75
MIN_SAMPLE_LEAF = 500
OFFSET = 5.0/6.0
OUTLIER_REMOVAL = True
DATA_FILE='./train/ml_finalproj_train_vF.pkl'

print("Q features:", Q_FEATURES)
print("D features:", D_FEATURES)
print("C features:", C_FEATURES)
print("S features:", S_FEATURES)
print("Num Estimators:", N_ESTIMATORS)
print("Min Samples Leaf:", MIN_SAMPLE_LEAF)
print("Offset ratio:", round(OFFSET/(1-OFFSET),2))

def get_data(offset):
    file_name = DATA_FILE 
    da = Data_accessor(file_name)
    data = da.load()
    data.sort_values([TIME, ID], ascending=[1, 1], inplace=True)
    
    data['x49'] = data['x49'].apply(lambda x: abs(x))
    #data['x61'] = data['x61'].apply(lambda x: abs(x))
    #data['x53'] = data['x53'].apply(lambda x: abs(x))

    # weight and output
    w = data[WEIGHT]
    y = data[OUTPUT]
    
    # quantitative features
    X_Q = data[Q_FEATURES]
    
    # discrete features
    X_D = data[D_FEATURES]
    
    # categorical features
    X_C = data[C_FEATURES]
    for col in C_FEATURES:
        one_hot = pd.get_dummies(X_C[col])
        X_C = pd.concat([X_C, one_hot], axis=1)
    X_C = X_C.drop(C_FEATURES, axis=1)

    # complete feature matrix
    X = pd.concat([X_Q, X_D, X_C], axis=1)
    
    # apply sigmoid transformation to a few spiky features
    #for col in S_FEATURES:
    #    X[col] = X[col].apply(lambda x: 1/(1+np.exp(-x)))
    
    offset = int(len(data) * offset)
    print("Offset =", offset)
   
    # prepare test data (not touched)
    X_test, y_test, w_test = X[offset:], y[offset:], w[offset:]
    w_test /= np.max(np.abs(w_test), axis=0)

    # outlier removal
    X_train, y_train, w_train = X[:offset], y[:offset], w[:offset]

    if OUTLIER_REMOVAL:
        Xwy_train = pd.concat([X_train, w_train, y_train], axis=1)
        Xwy_train = Xwy_train[np.abs(stats.zscore(Xwy_train['y']))<3]
        Xwy_train = Xwy_train[np.abs(stats.zscore(Xwy_train['x49']))<3]
        Xwy_train = Xwy_train[np.abs(stats.zscore(Xwy_train['x53']))<3]
        #Xwy_train = Xwy_train[np.abs(Xwy_train['y'])<0.05]
        #Xwy_train = Xwy_train[np.abs(Xwy_train['x49'])<0.009]
        #Xwy_train = Xwy_train[np.abs(Xwy_train['x53'])<0.009]
        #Xwy_train = Xwy_train[(np.abs(stats.zscore(Xwy_train)) < 3).all(axis=1)]
        w_train = Xwy_train[WEIGHT]
        y_train = Xwy_train[OUTPUT]
        X_train = Xwy_train.drop([WEIGHT,OUTPUT], axis=1)

    w_train /= np.max(np.abs(w_train), axis=0)

    print("# of training rows = ", len(y_train))
    print("# of testing rows = ", len(y_test))
    
    data = { 'X': X, 'X_train': X_train, 'X_test': X_test,
             'y': y, 'y_train': y_train, 'y_test': y_test,
             'w': w, 'w_train': w_train, 'w_test': w_test }

    return data

def print_config(data, params):
    print(data['X_train'])
    print(data['X_test'])
    print(data['y_train'])
    print(data['y_test'])
    print(data['w_train'])
    print(data['w_test'])
    print(params)

def get_score(clf, X_train, y_train, w_train, X_test, y_test, w_test):
    clf.fit(X_train, y_train, w_train)
    r2_in_sample = clf.score(X_train, y_train, w_train)
    r2_out_of_sample = clf.score(X_test, y_test, w_test)
    return (r2_out_of_sample, r2_in_sample)

def get_cross_score(clf, X, y, w):
    scores = cross_val_score(clf, X, y, fit_params = {'sample_weight': w})
    r2 = scores.mean()
    return (r2, None)

data = get_data(OFFSET)

X = data['X']
y = data['y']
w = data['w']
X_train = data['X_train']
y_train = data['y_train']
w_train = data['w_train']
X_test = data['X_test']
y_test = data['y_test']
w_test = data['w_test']

for max_depth in MAX_DEPTHS:
    for learning_rate in LEARNING_RATES:
        for subsample in SUBSAMPLES:

            params = {'n_estimators': N_ESTIMATORS, 
                      'max_depth': max_depth, 
                      'learning_rate': learning_rate, 
                      'subsample': subsample,
                      'min_samples_leaf': MIN_SAMPLE_LEAF,
                      'min_samples_split': MIN_SAMPLE_LEAF,
                      'loss': 'ls'}

            clf = ensemble.GradientBoostingRegressor(**params)

            r2 = get_score( clf, 
                            X_train, y_train, w_train,
                            X_test, y_test, w_test )
            r2_out_of_sample = r2[0]
            r2_in_sample = r2[1]
            
            print("max depth:", max_depth, 
                  ", learning rate:", learning_rate,
                  ", subsample:", subsample,
                  ", R^2 in sample:", r2_in_sample,
                  ", R^2 out of sample:", r2_out_of_sample)

