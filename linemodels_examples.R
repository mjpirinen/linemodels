### Example analyses using linemodels R functions
### Matti Pirinen, updated on 20-Nov-2023.

### Contents
### 1. Install linemodels package
### 2. Example 1. Simulated data with 100 variables.
### 3. Example 2. COVID-19 Host Genetic Initiative release 6 data



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

#
##
### Example 1. Simulated data with 100 variables.
##
#


# Effects are in columns 1,2 and SEs in cols 3,4
Y = linemodels.ex1[,1:2]
SE = linemodels.ex1[,3:4]

#Plot the raw data:
plot(Y, col = "black", pch = 5,
     xlab = expression(beta[1]), ylab = expression(beta[2]),
     xlim = c(-0.0,0.6), ylim = c(-0.03, 0.6), cex.lab = 1.3, cex.axis = 1.3)
grid()

#Come up with 3 line models:
scales = c(0.2,0.2,0.2)
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

# Use basic linemodels with fixed model priors = 1/3 for each model
# (Note: model.priors will be normalized automatically to sum to 1.)
res.lm = line.models(Y, SE,
                     scales = scales,
                     slopes = slopes,
                     cors = cors,
                     model.names = model.names,
                     model.priors = c(1,1,1),
                     r.lkhood = 0)

# Output is probabilities for each variable (rows) and each model (cols)
res.lm[1:2,]
# Color those that have posterior prob > 0.95 in one of the models
cols = colors[apply(res.lm, 1, which.max)] #color by max model
ind = (apply(res.lm, 1, max) > 0.95) #color only high probability cases
points(Y[ind,], pch = 19, col = cols[ind], cex = 0.5)

#Run line models with Gibbs sampler that estimates the proportion parameters

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


# Optimize all 3 parameters of the middle model "M.5"
# For this, 'par.include' matrix determines which parameters
# are being optimized:
par.include = rbind(c(FALSE, FALSE, FALSE), #1st model kept fixed
                    c(TRUE, TRUE, TRUE), #2nd model optimized all params
                    c(FALSE, FALSE, FALSE)) #3rd model kept fixed

# Perturb the initial values for the 2nd model
# from their "true" values (0.2,0.5 and 0.995) used in data simulation
line.models.optimize(Y, SE, par.include = par.include,
                     init.scales = c(scales[1], 0.05, scales[3]),
                     init.slopes = c(slopes[1], 0.2, slopes[3]),
                     init.cors = c(cors[1], 0.1, cors[3]),
                     model.priors = c(1,1,1),
                     model.names = model.names,
                     r.lkhood = 0, tol.loglk = 1e-2,
                     print.steps = 2)

# Optimum found when scale = 0.23, slope = 0.50 and cor = 1.00.
# These are close to the values used in the data simulation (0.2, 0.5, 0.995)
#  showing that optimization works.
# But remember that EM-algorithm could converge to a local maximum so several
#  runs with varying initial values are recommended in real data analyses


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
  print.steps = 0)

# Results tell that for the observed data, the log likelihood ratio
# between the alternative and null models is
# $obs.loglr
#  263.1592
# And when simulating data sets that are similar in structure to the
# observed data and that follow the estimated null model,
# then 10 examples of logLR values are
# $null.loglr
#  0.9242224 0.5098504 1.1606126 0.6065173 0.5418850
#  0.3789954 4.5949677 1.1753738 3.3070391 0.6493173
#
# This suggests that the observed logLR is highly statistically significant
# and hence the 2nd model seems highly useful in explaining the data.
# In general, one should run such a simulation for at least 100 times to
# get a reliable idea how significant the result is.
# Note that the run time would decrease if less parameters were optimized.


#
##
### Example 2. Analysis of COVID-19 Host Genetic Initiative release 6 data.
##
#

# First load in the linemodels package as instructed at the top of this file.

#Read in COVID-19 data
data.file = "https://raw.githubusercontent.com/mjpirinen/covid19-hgi_subtypes/main/covid_hgi_v6_B2_C2_common.tsv"
dat = read.csv(data.file,sep = "\t", header = T, as.is = T)
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




