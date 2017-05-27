#!/home/quan/anaconda3/bin/python

from RandomForestHW import *

# parameters
num_obs = 50000
depth = 2 
min_points_in_leaf = 100
num_cv_folds = 5
    
noise = 10
X_train, Y_train = generate_test_data(num_obs, noise)
X_test, Y_test = generate_test_data(num_obs, noise)
    
regression_tree = RegressionTree(depth, min_points_in_leaf, num_cv_folds)
regression_tree.fit(X_train, Y_train)
