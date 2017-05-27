def run_test(X_train, Y_train, X_test, Y_test, num_obs, max_depth, min_points_in_leaf, num_cv_folds=1):
    """
        Parameters
        ----------
        X_train: X of training data set
        Y_train: Y of training data set
        X_test: X of test data set
        Y_test: Y of test data set
        num_obs: number of instances of observations (data points);
        max_depth: maximum tree depths to test;
        min_points_in_leaf: minimum number of data points in leaf node;
        num_cv_folds: number of cross validation folds.
        
        Returns
        ----------
        depth_list: list of max tree depths
        r2_fit_list: list of R2 for in-sample fit
        r2_predict_list: list of R2 for out-of-sample prediction
    """

    depth_list = []
    r2_fit_list = []
    r2_predict_list = []
    for i in range(max_depth):
        # fit
        regression_tree = RegressionTree(i+1, min_points_in_leaf, num_cv_folds)
        regression_tree.fit(X_train, Y_train)
        
        # in-sample fit
        Y_fit = regression_tree.predict(X_train)
        # out-of-sample prediction
        Y_predict = regression_tree.predict(X_test)

        # R2 for in-sample fit
        var_y_train = np.var(Y_train)
        var_residual_fit = np.var(Y_train - Y_fit)
        r2_fit = 1 - var_residual_fit/var_y_train

        # R2 for out-of-sample prediction
        var_y_test = np.var(Y_test)
        var_residual_predict = np.var(Y_test - Y_predict)
        r2_predict = 1 - var_residual_predict/var_y_test
    
        depth_list.append(i+1)
        r2_fit_list.append(r2_fit)
        r2_predict_list.append(r2_predict)
        
    return (depth_list, r2_fit_list, r2_predict_list)
