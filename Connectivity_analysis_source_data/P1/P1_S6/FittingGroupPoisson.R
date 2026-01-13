library(grpreg)
library(R.matlab)


dat <- readMat("glm_input_2_node_1.mat")
X <- as.matrix(dat$X)
y <- as.vector(dat$y)
group <- as.vector(dat$group)


storage.mode(X) <- "double"
y <- as.double(y)

fit <- grpreg(X, y, group=group, penalty="grLasso", family="poisson")


Lambda <- as.vector(fit$lambda)
Beta   <- fit$beta

writeMat("grpreg_path.mat", X=X, y=y, Lambda=Lambda, Beta=Beta, group = group)



