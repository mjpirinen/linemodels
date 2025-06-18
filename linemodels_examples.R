### Example analyses using linemodels R functions
### Matti Pirinen, updated on 17-June-2025.

### Contents
### 1. Install linemodels package
### 2. Example 1. Simulated data with 100 observations in 2-dimensions
### 3. Example 2. COVID-19 Host Genetic Initiative release 6 data
### 4. Example 3. Simulated data with 100 observations in 3-dimensions
### 5. Example 4. Annotations with 101 observations in 2-dimensions
### 6. Example 5. Annotations with 101 observations in 3-dimensions



# To install linemodels from GitHub,
#  you need to install 'devtools' package (in case you don't have it).
install.packages("devtools")
library(devtools)
install_github("mjpirinen/linemodels")
library(linemodels)

#  If you can't get the above working, all functions are in a single R file,
#  which you can also read in directly to R using command below
#source("https://raw.githubusercontent.com/mjpirinen/linemodels/main/R/linemodels.R")
#  And you can read in the example data distributed with the package
#filename = "https://github.com/mjpirinen/linemodels/raw/main/data/linemodels.ex1.rda"
#download.file(filename, dest = "./linemodels.ex1.rda") #put here suitable path
#load("./linemodels.ex1.rda") #the same path here
#  Now you have data.frame called 'linemodels.ex1'
head(linemodels.ex1)
#Should show 7 cols:
#beta1, beta2, beta3, se1, se2, se3, maf
# Similarly you can read 'linemodels.ex2'.

#
##
### Example 1. Simulated data with 100 observations on 2 effect variables.
##
#

# In this example, we will use only 2-dimensional subset of the data.
# Effects are in columns beta1, beta2 and standard errors (SE) in cols se1 and se2.
Y = linemodels.ex1[,c("beta1","beta2")]
SE = linemodels.ex1[,c("se1","se2")]

#Plot the raw data:
plot(Y, col = "black", pch = 5,
     xlab = expression(beta[1]), ylab = expression(beta[2]),
     xlim = c(-0.0,0.6), ylim = c(-0.03, 0.6),
     cex.lab = 1.3, cex.axis = 1.3)
grid()

#Come up with 3 line models:
scales = c(0.2, 0.2, 0.2)
slopes = c(0, 0.5, 1)
cors = c(0.995, 0.995, 0.995)
model.names = c("M0","M.5","M1")

#Plot the models, and add the data on top to see whether models make sense
colors = c("red","orange","dodgerblue")
visualize.line.models(scales, slopes, cors,
                      model.names = model.names, model.cols = colors,
                      legend.position = "topleft",
                      xlim = c(-0.0, 0.6), ylim = c(-0.03, 0.6),
                      xlab = expression(beta[1]), ylab = expression(beta[2]),
                      emphasize.axes = FALSE, cex.lab = 1.3, cex.axis = 1.3)
points(Y, pch = 1, col = "black", cex = 0.5)

# Assign the observations to the linemodels assuming fixed model priors = 1/3
# for each linemodel
# (Note: model.priors will be normalized automatically to sum to 1.)
res.lm = line.models(Y, SE,
                     scales = scales,
                     slopes = slopes,
                     cors = cors,
                     model.names = model.names,
                     model.priors = c(1,1,1),
                     r.lkhood = 0)

# Output includes probabilities for each observation (rows) and each linemodel (cols)
res.lm[1:2,]
# Color those points that have posterior prob > 0.95 in one of the linemodels
cols = colors[apply(res.lm, 1, which.max)] #color by max model
ind = (apply(res.lm, 1, max) > 0.95) #color only high probability cases
points(Y[ind,], pch = 19, col = cols[ind], cex = 0.5)

# If we don't want to fix the model.priors, we can
#  run a Gibbs sampler that estimates the proportion parameters.

res.prop = line.models.with.proportions(Y, SE,
                                        scales = scales,
                                        slopes = slopes,
                                        cors = cors,
                                        model.names = NULL,
                                        r.lkhood = 0,
                                        n.iter = 200, n.burnin = 20)

# Show the proportion estimates with uncertainty
res.prop$params

# Show individuals probabilities
res.prop$groups[1:2,]

# Circle those that have posterior > 95% for one model
cols = colors[apply(res.prop$groups, 1, which.max)] #color by max model
ind = (apply(res.prop$groups, 1, max) > 0.95) #color only high probability cases
points(Y[ind,], pch = 1, col = cols[ind], cex = 0.9, lwd = 0.5)

# Coloring is very similar between the two versions of linemodel
# functions since the estimated proportions in the 2nd analysis
# are fairly close to
# the fixed prior values of 1/3 used in the first analysis.

# If we don't want to fix the linemodels, we can estimate their parameters
# by maximum likelihood.
# Optimize all 3 parameters of the middle model "M.5"
# For this, 'par.include' matrix determines which parameters
# are being optimized.
# NOTE: par.include can also be given as a list (see example 3 below)
par.include = rbind(c(FALSE, FALSE, FALSE), #1st model kept fixed
                    c(TRUE, TRUE, TRUE), #2nd model optimized all params
                    c(FALSE, FALSE, FALSE)) #3rd model kept fixed

# Perturb the initial values for the 2nd model
# from their "true" values (0.2, 0.5 and 0.995) used in data simulation
line.models.optimize(Y, SE, par.include = par.include,
                     init.scales = c(scales[1], 0.05, scales[3]),
                     init.slopes = c(slopes[1], 0.2, slopes[3]),
                     init.cors = c(cors[1], 0.1, cors[3]),
                     model.priors = c(1,1,1),
                     model.names = model.names,
                     r.lkhood = 0, tol.loglk = 1e-2,
                     op.method = "BFGS",
                     print.steps = 2)

# Optimum found when scale = 0.23, slope = 0.51 and cor = 1.00.
# These are close to the values used in the data simulation (0.2, 0.5, 0.995).
# NOTE: EM-algorithm could converge to a local maximum so several
#  runs with varying initial values are recommended in real data analyses.
#  by parameter 'op.method' one can choose the method to be used in optimization.
# default is "BFGS", others such as "Nelder-Mead" can be used to confirm the result.


# Is the 2nd model statistically useful?
# Compute empirical P-value for adding the 2nd model compare to the
#  NULL where only 1st and 3rd model were given.
# NOTE: Here simulating only 10 logLR values under the null,
#       To get more reliability to say e.g. that P < 0.05,
#       would be good to do at least n.sims = 100 simulations.

set.seed(1655)
simulate.loglr(
  X = Y, SE = SE,
  n.sims = 10,
  par.include.null = rbind(c(FALSE,FALSE,FALSE),
                           c(FALSE,FALSE,FALSE)),
  init.scales.null = c(0.2, 0.2),
  init.slopes.null = c(0, 1),
  init.cors.null = c(0.995, 0.995),
  par.include.alt = rbind(c(FALSE,FALSE,FALSE),
                          c(TRUE,TRUE,TRUE),
                          c(FALSE,FALSE,FALSE)),
  init.scales.alt = c(0.2, 0.2, 0.2),
  init.slopes.alt = c(0, 0.5, 1),
  init.cors.alt = c(0.995, 0.995, 0.995),
  r.lkhood = 0, tol.loglk = 1e-2,
  op.method = "BFGS",
  print.steps = c(1,0))

# Results tell that for the observed data, the log likelihood ratio
# between the alternative and null models is
# $obs.loglr
# 233.6457
# And when simulating data sets that are similar in structure to the
# observed data and that follow the estimated null model,
# then 10 examples of logLR values are
# $null.loglr
# 0.6318554 0.7195902 0.4896049 1.3947264 0.6974900
# 0.4009584 1.2168922 0.3478202 1.7297213 1.9565947

# This suggests that the observed logLR is highly statistically significant
# and hence the 2nd model seems highly useful in explaining the data.
# In general, one should run such a simulation for at least 100 times to
# get a reliable idea how significant the result is.
# Note that the run time would decrease if less parameters were optimized,
# or if we can use option 'assume.constant.SE = TRUE' (see below).

# If we can assume that, for every outcome variable, the SEs of the effects
# across the data points are constant, then we can speed up the simulate.loglr()
# function.
# In genetics application, SEs may become approx. constant when the effects
# are scaled by values sqrt(2*f*(1-f)) where f is the variant-specific
# minor allele frequency.
# Let's check this here:

sc = sqrt(2 * linemodels.ex1$maf * (1 - linemodels.ex1$maf))
summary(sc*linemodels.ex1[,c("se1","se2","se3")])
# Scaled SEs are indeed constant for each of the 3 outcome variables.
# Thus, for the scaled effects, we can set 'assume.constant.SE = TRUE'
# in computation.

Y.sc = Y*sc
SE.sc = SE*sc

# Let's redo the previous simulation. First we need to re-estimate scale parameters
# because of the new scales in the data.

par.include = rbind(c(TRUE, FALSE, FALSE), #1st model optimizing only scale
                    c(TRUE, TRUE, TRUE), #2nd model optimized all params
                    c(TRUE, FALSE, FALSE)) #3rd model optimizing only scale

# Perturb the initial values for the 2nd model
# from their "true" values (0.2,0.5 and 0.995) used in data simulation
line.models.optimize(Y.sc, SE.sc, par.include = par.include,
                     init.scales = c(1,1,1),
                     init.slopes = c(slopes[1], 0.2, slopes[3]),
                     init.cors = c(cors[1], 0.1, cors[3]),
                     model.priors = c(1,1,1),
                     model.names = model.names,
                     r.lkhood = 0, tol.loglk = 1e-2,
                     assume.constant.SE = TRUE,
                     op.method = "BFGS",
                     print.steps = 2)

# Results become available instantly. New scales are
# 0.1515844 0.1418315 0.1350342
# So we will use 0.15 as a reasonable initial value for every scale.

# Now application of simulate.loglr() with assume.constant.SE = TRUE:

set.seed(1655)
simulate.loglr(
  X = Y.sc, SE = SE.sc,
  n.sims = 10,
  par.include.null = rbind(c(FALSE,FALSE,FALSE),
                           c(FALSE,FALSE,FALSE)),
  init.scales.null = c(0.15, 0.15),
  init.slopes.null = c(0, 1),
  init.cors.null = c(0.995, 0.995),
  par.include.alt = rbind(c(FALSE,FALSE,FALSE),
                          c(TRUE,TRUE,TRUE),
                          c(FALSE,FALSE,FALSE)),
  init.scales.alt = c(0.15, 0.15, 0.15),
  init.slopes.alt = c(0, 0.5, 1),
  init.cors.alt = c(0.995, 0.995, 0.995),
  r.lkhood = 0, tol.loglk = 1e-2,
  assume.constant.SE = TRUE,
  op.method = "BFGS",
  print.steps = c(1,0))

# Results lead to same conclusions as above but took only a fraction of time.
# NOTE: The numerical values are not the same as above since input data have
# changed due to scaling.


#
##
### Example 2. Analysis of COVID-19 Host Genetic Initiative release 6 data.
##
#

# First load in the linemodels package as instructed at the top of this file.

#Read in COVID-19 data
data.file = "https://raw.githubusercontent.com/mjpirinen/covid19-hgi_subtypes/main/covid_hgi_v6_B2_C2_common.tsv"
dat = read.csv(data.file,sep = "\t", header = TRUE, as.is = TRUE)
dat = dat[dat[,"SNP"] != "12:112914354:T:C",] #this is in LD with its neighbor -- remove from analysis
X = dat[,c("B2_beta","C2_beta")] #B2 is hospitalization GWAS, C2 is infection GWAS
SE = dat[,c("B2_sebeta","C2_sebeta")]
ii = X[,1] < 0 # Flip reference allele to make all B2 effects non-negative for nicer visualisation
X[ii,] = -X[ii,]

# First, consider Models severity of diseases (SEV) and
#  susceptibility for infection (SUC) with slopes 0.2 and 1, resp.
# For motivation to set SEV slope = 0.2 see
#  https://github.com/mjpirinen/covid19-hgi_subtypes/blob/main/README_Subtype_assignment.pdf
# Let's add a third model, called 'BOTH'.
# Its slope is chosen to halve the angle between SEV and SUC,
# and the model represents variants that have an effect on both severity
# and susceptibility.
slope.both = tan(atan(0.2) + (atan(1) - atan(0.2))/2)
slope.both

scales = c(0.15, 0.15, 0.15)
slopes = c(0.2, 1, slope.both)
cors = c(0.999, 0.999, 0.999)
model.names = c("SEVER.", "SUSCEP.", "BOTH")


# Two GWAS are correlated. The derivation of this 'r.lkhood' value can be
# found in the document 'README_Subtype_assignment.pdf' mentioned above.
r.lkhood = 0.4539485


# Run basic linemodels with equal priors on each model
res.1 = line.models(X, SE, scales, slopes, cors, model.names,
                    model.priors = rep(1/length(slopes), length(slopes)),
                    r.lkhood = r.lkhood, scale.weights = c(1))

# Run Gibbs sampler version of linemodels to estimate the proportion parameters
# of the three models (rather than assuming them equal as above in res.1)
set.seed(91)
res.2 = line.models.with.proportions(X, SE, scales, slopes, cors,
                                     model.names, r.lkhood = r.lkhood,
                                     n.iter = 10000, n.burnin = 50)


# Plot the results from res.2 version of the results

n = nrow(X)
K = ncol(res.2$groups)
p.thresh = 0.95
colors = c("red","dodgerblue","orange")
pchs.grp = c(19, 19, 19)
cols = rep("white",n)
pchs = rep(19, n)
for(ii in 1:K){
  ind = (res.2$groups[,ii] > p.thresh)
  cols[ind] = colors[ii]
  pchs[ind] = pchs.grp[ii]
}


#
###
#### Figure BEGINS
###
#



#pdf("Fig_COVID.pdf", width = 5.5, height = 3.2)
layout(matrix(c(1,2), ncol = 2), widths = c(2,1))
par(mar = c(4,4.2,1,1))
visualize.line.models(scales, slopes, cors,
                      model.names, model.cols = colors,
                      legend.position = NULL,
                      xlim = c(0, 0.45), ylim = c(0,0.45),
                      emphasize.axes = F,
                      xlab = expression(beta[HOS]), ylab = expression(beta[INF]))
arrows(X[,1]-1.96*SE[,1], X[,2], X[,1]+1.96*SE[,1], X[,2],
       code = 3, angle = 90, length = 0.0, col = "black", lty = 1, lwd = 0.5)
arrows(X[,1], X[,2]-1.96*SE[,2], X[,1], X[,2]+1.96*SE[,2],
       code = 3, angle = 90, length = 0.0, col = "black", lty = 1, lwd = 0.5)
points(X, col = cols, pch = pchs, cex = 0.8)
points(X, col = "black", pch = 1, cex = 0.9, lwd = 0.6)
legend("topleft", leg = c("SEVER.","SUSCEPT.","BOTH"),
       col = c(colors), cex = 0.8, lwd = 2.5, bg = "white")
text(-0.13,0.48, "A.", cex = 1.5,xpd = TRUE)

chosen.i = c(17,3,1)

x = X[chosen.i,1] + c(0.03, -0.03, 0.03)
y = X[chosen.i,2] + c(-0.03, 0.03, -0.03)
text(x,y, c("1","2","3"), cex = 1.2)
points(x, y, pch = 1, lwd = 1, cex = 2.5)

barplot(t(res.2$groups[chosen.i,]), col = colors, ylab = "probability",
        names.arg = c(1,2,3), cex.names = 1.1, cex.axis = 1.0, cex.lab = 1.0)
points(c(0.72,1.9,3.1),rep(-0.15,3), pch = 1, cex = 3, xpd = TRUE)
text(-2.8, 1.02, "B.", cex = 1.5, xpd = TRUE)

#dev.off()

#
###
#### Figure DONE.
###
#

#Check the results
res.2$params # proportions
colSums(res.2$groups > p.thresh) #How many confident assignments in each group
nrow(res.2$groups) # How many variants we had in total.


#Write the results to a file
Y = cbind(dat[,c(1:2,14,3:8,15:17)], res.2$groups)
names(Y)[c(3,7:15)] = c("RSID","HOS_beta","HOS_se","HOS_pval",
                        "INF_beta","INF_se","INF_pval",
                        "Pr_SEV","Pr_SUC","Pr_BOTH")
row.names(Y) = NULL
#write.table(Y,file = "COVID19HGI_linemodels_example_results.txt",
#            row.names = F, quote = F)


#
##
### 4. Example 3. Simulated data with 100 observations in 3-dimensions
##
#

# We work with data set linemodels.ex1 as in Example 1, but now with
# M = 3 dimensions.

# To make it efficiently, we assume constant SEs after scaling the effect sizes as
# in the last part of Example 2.
sc = sqrt(2 * linemodels.ex1$maf * (1 - linemodels.ex1$maf))
summary(sc * linemodels.ex1[,c("se1","se2","se3")])
# Scaled SEs are constant for each of the 3 variables.
# Thus, for scaled effects, we can use 'assume.constant.SE = TRUE' in computation.
Y.sc = linemodels.ex1[,c("beta1","beta2","beta3")]*sc
SE.sc = linemodels.ex1[,c("se1","se2","se3")]*sc

# Main difference compared to the 2 dimensional case is in optimization function, where
# par.include should be defined as a list with 3 components:
# K-vector 'scales',
# (K x (M-1))-matrix 'slopes' and
# K-vector 'cors',
# where K is the number of linemodels.
# Here, we will optimize all scales and slopes of 3 linemodels and assume
# that scales are the same for every linemodel (force.same.scale = TRUE).
# NOTE: We do not optimize cors because we want to keep linemodels identifiable
#       whereas a low cor parameter would make interpretation of the linemodel vague.
par.include = list(scales = c(TRUE, TRUE, TRUE), #optimizing all scales
                   slopes = matrix(TRUE, ncol = 2, nrow = 3), #optimizing all slopes
                   cors = c(FALSE, FALSE, FALSE)) #keeping all cors fixed to initial values

op.res = line.models.optimize(Y.sc, SE.sc, par.include = par.include,
                              init.scales = rep(0.15, 3),
                              init.slopes = matrix(c(0,0,
                                                     1,1,
                                                     0.5,0.5), byrow = TRUE, ncol = 2),
                              init.cors = c(0.995,0.995,0.995),
                              force.same.scales = TRUE,
                              r.lkhood = c(0,0,0),
                              tol.loglk = 1e-2,
                              op.method = "BFGS",
                              assume.constant.SE = TRUE,
                              print.steps = 2)
op.res

# The results correspond well to the true slope values used in generating the data,
# which were (0,0), (0.5, 0.2) and (1,1).

#Estimate membership proportions
res.prop = line.models.with.proportions(Y.sc, SE.sc,
                                        scales = op.res$scales,
                                        slopes = op.res$slopes,
                                        cors = op.res$cors,
                                        r.lkhood = c(0,0,0),
                                        n.iter = 200, n.burnin = 20)
res.prop$params

# Plot results using pairs of variables and including the data points with model
# assignments in the plot.

par(mfrow = c(1,3)) #plotting area must be defined outside the function
# plot.pairs tell which pairs of variables are plotted, here we plot all 3 pairs.
visualize.multidimensional.line.models(
  scales = op.res$scales,
  slopes = op.res$slopes,
  cors = op.res$cors,
  plot.pairs = rbind(c(1,2), c(1,3), c(2,3)),
  var.names = paste0("beta",1:3),
  model.cols = c("red","dodgerblue","orange"),
  X = Y.sc,
  model.prob = res.prop$groups,
  model.thresh = 0.8,
  undetermined.col = "gray",
  emphasize.axes = TRUE)

# Here the points with membership probability > 80% in one of the linemodels are colored.
# The remaining ("undetermined") points are left gray.


#
##
### 5. Example 4. Annotations with 101 observations in 2-dimensions
##
#
# We use data set linemodels.ex2
# Setting:
#  There are two sets of variants
#  50 "Blue" variants have slope = (0, -2)
#  50 "Red" variants have slope = (3, 6)
#  Last "Orange" variant (index 101) is an outlier
#  that is located in between the Blue and Red linemodels
#  and it has a large standard error, so assignment probability
#   is ~50%:50% (without annotations).

# Here we run in 2D and take only beta1 and beta2
X = linemodels.ex2[,c("beta1","beta2")]
SE = linemodels.ex2[,c("se1","se2")]
slopes = matrix(c(0, 3), ncol = 1)
scales = rep(0.1, nrow(slopes))
cors = rep(0.999, nrow(slopes))
colors = c("dodgerblue","red","orange")
n = nrow(X)
# visualize the data and the models
par(mfrow = c(1,1))
visualize.line.models(scales, slopes, cors,
                      model.cols = colors,
                      legend.position = NULL,
                      xlab = expression(beta[1]), ylab = expression(beta[2]),
                      emphasize.axes = FALSE, cex.lab = 1.3, cex.axis = 1.3)
points(X, pch = 1, col = "black", cex = 0.5)
arrows(X[n,1]-1.96*SE[n,1],X[n,2],X[n,1]+1.96*SE[n,1],X[n,2],
       code = 3, angle = 90, length = 0.01, col = "gray")
arrows(X[n,1],X[n,2]-1.96*SE[n,2],X[n,1],X[n,2]+1.96*SE[n,2],
       code = 3, angle = 90, length = 0.01, col = "gray")
points(X[n,1],X[n,2], col = "orange", cex = 1.5, lwd = 1.5)

# Run the data without annotations
res.p = line.models.with.proportions(X, SE,
                                     scales, slopes, cors,
                                     model.names = NULL,
                                     r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
                                     n.iter = 500, n.burnin = 20,
                                     verbose = TRUE)

# Run the data with annotations
res.f = line.models.with.features(X, SE,
                                  scales, slopes, cors,
                                  model.names = NULL,
                                  r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
                                  n.iter = 500, n.burnin = 20,
                                  features = linemodels.ex2$annotation,
                                  verbose = TRUE)

# Posterior probability of each point in each model is in res.f$groups
# Thus for the outlier point we have
res.p$groups[n,] #without annotations
res.f$groups[n,] #with annotations

#plot assignment probability of the outlier without and with annotations
par(mfrow = c(1,1))
barplot(t(rbind(res.p$groups[n,],
                res.f$groups[n,])), col = colors,
        main = "Posterior of the outlier",
        names.arg = c("No annot",
              "Annotation"))
# We see that 0.5 probability without annotations has turned to about
# 80% probability when annotations draw the outlier to the "Red" model.

# Plot the posteriors of the red model across all points between
# case with no annotation (x-axis) and annotation (y-axis).
# Here the outlier had annotation = 1 and clusters with
# the red points in the y-axis values.
plot(res.p$groups[,2], res.f1$groups[,2],
     col = c(rep(colors[1:2], times = c(50,50)), col = "orange"),
     pch = 19, ylab = "Pr(Red) with annotations",
     xlab = "Pr(Red) without annotations",
     main = "")
points(res.p$groups[n,2], res.f1$groups[n,2], col = "orange", cex = 1.5, lwd = 1.5)
abline(0,1)
par(mfrow = c(1,1))
legend("topleft",legend = c("True R (a = 1)","True B (a = 0)","Outlier (a = 1)"),
       pch = 19, col = colors[c(2,1,3)], cex = 0.7)



#The regression coefficients of the annotation model can be found in
res.f$params
#It has coeffs for the intercept and for the annotation value ("X1")
# Baseline model is Model1 and its coefficients are 0
# Coefficients are from multinomial regression and thus on log-odds scale.



#
##
### 6. Example 5. Annotations with 101 observations in 3-dimensions
##
#

# We continue from previous example but use 3D data.
X = linemodels.ex2[,c("beta1","beta2","beta3")]
SE = linemodels.ex2[,c("se1","se2","se3")]
slopes = matrix(c(0, -2,
                  3,  6), byrow = T, ncol = 2)
scales = rep(0.1, nrow(slopes))
cors = rep(0.999, nrow(slopes))
colors = c("dodgerblue","red","orange")
n = nrow(X)

# visualize the data and the models
par(mfrow = c(1,3))
visualize.multidimensional.line.models(scales, slopes, cors,
                                       plot.pairs = matrix(c(1,2, 1,3, 2,3),
                                                           ncol = 2, byrow = T),
                                       model.cols = colors[1:2],
                                       legend.position = NULL,
                                       X = X, SE = SE,
                                       model.prob = cbind(c(rep(c(1,0),c(50,50)),0.5),
                                                          c(rep(c(0,1),c(50,50)),0.5)),
                                       model.thresh = 0.9,
                                       undetermined.col = colors[3])

# Run the data without annotations
res.p = line.models.with.proportions(X, SE,
                                     scales, slopes, cors,
                                     model.names = NULL,
                                     r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
                                     n.iter = 500, n.burnin = 20,
                                     verbose = TRUE)

# Run the data with annotations
res.f = line.models.with.features(X, SE,
                                  scales, slopes, cors,
                                  model.names = NULL,
                                  r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
                                  n.iter = 500, n.burnin = 20,
                                  features = linemodels.ex2$annotation,
                                  verbose = TRUE)

# Posterior probability of each point in each model is in res.f$groups
# Thus for the outlier point we have
res.p$groups[n,] #without annotations
res.f$groups[n,] #with annotations

#plot assignment probability of the outlier without and with annotations
par(mfrow = c(1,1))
barplot(t(rbind(res.p$groups[n,],
                res.f$groups[n,])), col = colors,
        main = "Posterior of the outlier",
        names.arg = c("No annot",
                      "Annotation"))
# We see that 0.5 probability without annotations has turned to about
# 80% probability when annotations draw the outlier to the "Red" model.

# Plot the posteriors of the red model across all points between
# case with no annotation (x-axis) and annotation (y-axis).
# Here the outlier had annotation = 1 and clusters with
# the red points in the y-axis values.
plot(res.p$groups[,2], res.f1$groups[,2],
     col = c(rep(colors[1:2], times = c(50,50)), col = "orange"),
     pch = 19, ylab = "Pr(Red) with annotations",
     xlab = "Pr(Red) without annotations",
     main = "")
points(res.p$groups[n,2], res.f1$groups[n,2], col = "orange", cex = 1.5, lwd = 1.5)
abline(0,1)
par(mfrow = c(1,1))
legend("topleft",legend = c("True R (a = 1)","True B (a = 0)","Outlier (a = 1)"),
       pch = 19, col = colors[c(2,1,3)], cex = 0.7)


#The regression coefficients of the annotation model can be found in
res.f$params
#It has coeffs for the intercept and for the annotation value ("X1")
# Baseline model is Model1 and its coefficients are 0
# Coefficients are from multinomial regression and thus on log-odds scale.

