# library(grpreg)
# library(R.matlab)
# 
# 
# dat <- readMat("glm_input_4_nodes_region2_factor20.mat")
# X <- as.matrix(dat$X)
# y <- as.vector(dat$y)
# group <- as.vector(dat$group)
# 
# 
# storage.mode(X) <- "double"
# y <- as.double(y)
# 
# fit <- grpreg(X, y, group=group, penalty="grLasso", family="poisson")
# 
# 
# Lambda <- as.vector(fit$lambda)
# Beta   <- fit$beta
# 
# writeMat("grpreg_4_nodes_region2_factor20.mat", X=X, y=y, Lambda=Lambda, Beta=Beta, group = group)





library(grpreg)
library(R.matlab)

n <- 4
factor_n <- 20

for (region in 1:n) {
  
  in_file  <- sprintf("glm_input_%d_nodes_region%d_factor%d.mat", n, region, factor_n)
  out_file <- sprintf("grpreg_%d_nodes_region%d_factor%d.mat",    n, region, factor_n)
  
  dat <- readMat(in_file)
  
  X <- as.matrix(dat$X); 
  y <- as.vector(dat$y);
  group <- as.vector(dat$group)
  
  
  storage.mode(X) <- "double"
  y <- as.double(y)
  
  
  stopifnot(ncol(X) == length(group), nrow(X) == length(y))
  
  fit <- grpreg(X, y, group = group, penalty = "grLasso", family = "poisson")
  
  Lambda <- as.vector(fit$lambda)
  Beta   <- fit$beta
  
  writeMat(out_file,
           X=X, y=y, Lambda=Lambda, Beta=Beta, group = group)
}

