def generate_test_data(N, noise=100):
    x = np.random.randn(N, 5)
    y = np.where(x[:, 0] > 0, 2, 5)
    y = y + np.where(x[:, 1] > 0, -3, 3)
    y = y + np.where(x[:, 2] > 0, 0, 0.5)
    y = y + np.random.randn(N)*noise
    return x,y

class TreeNode:
    def predict(x, y):
        assert False

    def depth(self):
        assert False

class BranchNode(TreeNode):
    def __init__(self, left, right, split_var_index, split_var_value):
        self.left = left
        self.right = right
        self.split_var_index = split_var_index
        self.split_var_value = split_var_value

    def predict(self, x):
        svar = x[:, self.split_var_index]
        is_left = svar < self.split_var_value
        leftx = x[is_left]
        rightx = x[~is_left]

        rv = np.zeros(x.shape[0])
        rv[is_left] = self.left.predict(leftx)
        rv[~is_left] = self.right.predict(rightx)

        return rv

    def depth(self):
        return 1 + max(self.left.depth(), self.right.depth())

class LeafNode(TreeNode):
    def __init__(self, mu):
        self.mu = mu

    def predict(self, x):
        return np.repeat(self.mu, x.shape[0])

    def depth(self):
        return 1

class RegressionTree:
    def __init__(self, max_depth, min_points_in_leaf, num_cv_folds=1):
        self.max_depth = max_depth
        self.min_points_in_leaf = min_points_in_leaf
        self.num_cv_folds = num_cv_folds

    def predict(self, x):
        assert self.fitted
        return self.root.predict(x)
    
    def compute_cv_variance(self, y):
        fold_size = int(len(y)/self.num_cv_folds)
        var = 0
        for k in range(self.num_cv_folds):
            test = y[k*fold_size : (k+1)*fold_size]
            train = np.concatenate([y[0:k*fold_size], y[(k+1)*fold_size:len(y)]])
            train_mean = np.mean(train)
            test_demean = test - train_mean
            test_var = np.sum(test_demean*test_demean)/len(test_demean)
            var += test_var
        var /= self.num_cv_folds
        return var

    def fit(self, x, y):
        self.fitted = True
        self.root = self.fit_internal(x, y, 1)

    def fit_internal(self, x, y, current_depth):
        num_features = x.shape[1]
        num_rows = x.shape[0]
        var_orig = np.var(y)

        if current_depth == self.max_depth:
            return LeafNode(np.mean(y))

        best_variable = None
        best_split = None
        best_leftx = None
        best_lefty = None
        best_rightx = None
        best_righty = None

        # Here, we have to loop over all features and figure out which one
        # might be splittable, and if it is, how to split it to maximize Variance Reduction
        max_var_reduction = 0
        for i in range(num_features):
            split_var_index = i
            split_var_values = np.percentile(x[:,i], [20,40,60,80])
            for split_var_value in split_var_values:
                svar = x[:, split_var_index]
                is_left = svar < split_var_value
                lefty = y[is_left]
                righty = y[~is_left]

                if len(lefty) <= self.min_points_in_leaf or len(righty) <= self.min_points_in_leaf:
                    continue

                var_lefty = 0
                var_righty = 0
                var_reduction = var_orig
                if self.num_cv_folds > 1:
                    var_lefty = self.compute_cv_variance(lefty)
                    var_righty = self.compute_cv_variance(righty)
                else:
                    var_lefty = np.var(lefty)
                    var_righty = np.var(righty)
                    
                var_reduction -= var_lefty*len(lefty)/len(y)
                var_reduction -= var_righty*len(righty)/len(y)

                if var_reduction > max_var_reduction:
                    best_variable = split_var_index
                    best_split = split_var_value
                    best_leftx = x[is_left]
                    best_lefty = y[is_left]
                    best_rightx = x[~is_left]
                    best_righty = y[~is_left]
                    max_var_reduction = var_reduction
        

        if best_variable is None:
            return LeafNode(np.mean(y))
        else:
            left_tree = self.fit_internal(best_leftx, best_lefty, current_depth+1)
            right_tree = self.fit_internal(best_rightx, best_righty, current_depth+1)
            return BranchNode(left_tree, right_tree, best_variable, best_split)


    def depth(self):
        return self.root.depth()
