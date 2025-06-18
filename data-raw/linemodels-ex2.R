## code to prepare `linemodels.ex2` dataset goes here
slopes = matrix(c(0, -2,
                  3,  6), byrow = TRUE, ncol = 2)
set.seed(48)
ns = c(50,50) #Blue and Red points
scale = 0.1
cor = 0.999
se = 0.005
X = c()

for(ii in 1:length(ns)){
  X = rbind(X, sample.line.model(n = ns[ii], scale = scale,
                                 slope = slopes[ii,],
                                 cor = cor))}
# To filter out points that have small effect on all coordinates,
# as these will be uncertain in the classification,
# use the following code:
#ind = (sqrt(rowSums(X^2)) > 0.1*scale)
#X = X[ind,]
#ns = c(sum(ind[1:ns[1]]), sum(ind[(ns[1]+1):length(ind)]))

X = X + rnorm(nrow(X)*ncol(X), 0, se) #add noise
# Add outlier point as the last point.
# Find vector 'u' that bisects the angle between the lines 1 and 2
v1 = c(1,slopes[1,]); v1 = v1/sqrt(sum(v1^2))
v2 = c(1,slopes[2,]); v2 = v2/sqrt(sum(v2^2))
u = v1 + v2; u = u/sqrt(sum(u^2))
#outlier point has max |component| =1.5*scale:
X = rbind(X, 1.5*scale*u/max(abs(u)))

#swap x-value to be positive to simplify plotting:
#ind = (X[,1] < 0)
#X[ind,] = -X[ind,]

n = nrow(X)

SE = matrix(se, nrow = n, ncol = ncol(X))
SE [n,] = 0.05 #larger SE for the outlier

linemodels.ex2 = data.frame(beta1 = X[,1], beta2 = X[,2], beta3 = X[,3],
                            se1 = SE[,1], se2 = SE[,2], se3 = SE[,3],
                            annotation = rep(c(0,1),c(ns[1],ns[2]+1)))

usethis::use_data(linemodels.ex2, overwrite = TRUE)
