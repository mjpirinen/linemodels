## code to prepare `linemodels.ex1` dataset goes here
scales = c(0.2, 0.2, 0.2)
slopes = matrix(c(0, 0,
                  0.5, 0,
                  1, 1), byrow = TRUE, ncol = 2)
cors = c(0.995, 0.995, 0.995)
ns = c(20, 40, 40)
models = rep(c(1,2,3), ns)
n = sum(ns)
K = length(scales)

set.seed(16)
X = c()
for(ii in 1:K){
  X = rbind(X, sample.line.model(n = ns[ii], scale = scales[ii],
                                 slope = slopes[ii,], cor = cors[ii]))}

N = 50000 # Total GWAS sample size
phi = 0.5 # proportion of cases
f = runif(n, 0.05, 0.5) # minor allele frequency
SE = 1/sqrt(N*phi*(1-phi)*2*f*(1-f))
summary(SE)
Y = X + matrix(rnorm(K*nrow(X), 0, rep(SE,each = K)), ncol = K, byrow = T)
ii = (Y[,1] < 0)
Y[ii,] = -Y[ii,]

linemodels.ex1 = data.frame(beta1 = Y[,1], beta2 = Y[,2], beta3 = Y[,3],
                            se1 = SE, se2 = SE, se3 = SE, maf = f)

usethis::use_data(linemodels.ex1, overwrite = TRUE)
