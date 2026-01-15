import numpy as np

from skglm import GeneralizedLinearEstimator
from skglm.datafits import PoissonGroup
from skglm.penalties import WeightedGroupL2
from skglm.solvers import GroupProxNewton
from skglm.utils.data import grp_converter  # helper to build grp_indices / grp_ptr

# X: shape (n_samples, n_features)
# y: non-negative counts, shape (n_samples,)
# X, y = ...




import scipy
from scipy.io import loadmat

data = loadmat("glm_input_2_node_1.mat")


print(data.keys())      # inspect variables
X = data["X"]           # access a variable


n_features = X.shape[1]

# Example: contiguous groups of fixed size
grp_size = 10
grp_indices, grp_ptr = grp_converter(grp_size, n_features)
n_groups = len(grp_ptr) - 1

# Group weights (often all ones to start)
weights_g = np.ones(n_groups, dtype=np.float64)

alpha = 1.0  # regularization strength (tune this)

datafit = PoissonGroup(grp_ptr=grp_ptr, grp_indices=grp_indices)
penalty = WeightedGroupL2(alpha=alpha, weights=weights_g,
                         grp_ptr=grp_ptr, grp_indices=grp_indices)

solver = GroupProxNewton()  # group-aware solver

est = GeneralizedLinearEstimator(datafit=datafit, penalty=penalty, solver=solver)
est.fit(X, y)

# For Poisson, predict(...) returns the mean parameter (inverse-link applied; exp(·))
mu_hat = est.predict(X)
