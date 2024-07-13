# linemodels (0.4.0): An R package to cluster multi-dimensional effects into
#  groups defined by linear relationships.
# DATE: 10-July-2024
# Copyright (C) 2022-2024, Matti Pirinen
# Contact: matti.pirinen@helsinki.fi
# LICENSE:
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################

# Consider M effect sizes.
# The linemodels are defined by three parameters:
# scale = standard deviation of the (largest) effect,
# slope = set of M-1 slopes defining the line around which the effects are scattered,
# cor = non-negative pairwise correlation between each effect pair,
# where cor = 1 means the same effect and cor = 0 means independent effects.
#
# A line model has the following properties:
#
# 1) Effects are scattered around line defined by vector (1,slope_1,...,slope_(M-1)).
#    Each slope_i is the slope between effect variables i+1 and 1.
#    Slope can be any real number or +/-Inf.
#    If any slope is infinite, the effect 1 and effects with finite slopes are 0.
#    It is recommended to keep the slopes finite to avoid problems with interpretation.
#
# 2) The largest of the prior variances of the effects is scale^2.
#    The effect that has the largest scale is determined by the slopes,
#    and the scales of the remaining effects are determined by the slopes as well.
#
# 3) Distribution around the line is determined as follows.
#    Consider a distribution scattered around the direction of the vector (1,1,...,1)
#    with the given constant value of the pairwise correlations between the effects.
#    Rotate that distribution by an orthogonal rotation that rotates the direction of the
#    vector (1,1,...,1) to the direction of the vector (1,slope_1,...,slope_(M-1))
#    and use the corresponding distribution, scaled so that the maximum variance
#    among the effects is set to scale^2.
#    NOTE: This is not same as having the given correlation around the line defined by the slopes,
#          because the shape of that distribution depends on the slopes but the
#          definition we use here is independent of slopes up to an orthogonal transformation.
#
# Examples:
# Set scale = 0 to get the null model. (Values of slope and cor do not matter in this case.)
# Use scale > 0, slope = 1 and cor = 0 to get independent effects model for 2-dim case.
# Use scale > 0, slope = c(1,1,1) and cor = 1 to get fixed effects model for 4-dim case.
# Use scale > 0, slope = 1 and 0 < cor < 1  to get correlated effects model for 2-dim case.
#               note that cor value needs to be quite near 1 to get very similar values
#               between the two effects with high probability (e.g. cor = 0.99...0.999).
#
# To choose a reasonable scale, note that 95% of the effects are assumed to be < 2*scale.

#
# 12 functions for user to call directly are listed below.
#

# beta.cor.case.control(s1, r1, s2 ,r2, s1s2 = 0, r1r2 = 0, s1r2 = 0, r1s2 = 0)
# -- To compute correlation in effect estimators of two case-control studies

# slope.for.pair(ii, jj, slopes)
# -- To extract the slopes for a given pair of effects from the slopes matrix

# visualize.line.models(
#   scales, slopes, cors,
#   model.names = NULL, model.cols = NULL,
#   legend.position = "bottomright",
#   xlim = NULL, ylim = NULL,
#   xlab = "EFFECT1", ylab = "EFFECT2",
#   cex.lab = 1, cex.axis = 1,
#   line.lty = 1, region.lty = 2,
#   line.lwd = 1, region.lwd = 1,
#   plot.grid = TRUE,
#   plot.new = TRUE,
#   emphasize.axes = TRUE)
# -- To visualize 2-dim linemodels and their 95% highest probability regions.

# visualize.multidimensional.line.models(
#   scales, slopes, cors,
#   plot.pairs = NULL, var.names = NULL,
#   model.names = NULL, model.cols = NULL,
#   legend.position = "bottomright",
#   X = NULL, SE = NULL,
#   SE.mult = 1.96,
#   model.prob = NULL, model.thresh = 0.8,
#   undetermined.col = "gray",
#   pch = 1, cex = 1,
#   xlim = NULL, ylim = NULL,
#   cex.lab = 1, cex.axis = 1,
#   line.lty = 1, region.lty = 2,
#   line.lwd = 1, region.lwd = 1,
#   plot.grid = TRUE,
#   emphasize.axes = TRUE)
# -- To visualize multi-dim linemodels and their 95% highest probability regions.

# visualize.classification.regions(
#   scales, slopes, cors, SE,
#   r.lkhood = 0,
#   model.priors = rep(1,length(scales)),
#   col.palette = NULL,
#   add.legend = TRUE,
#   xlim = NULL, ylim = NULL,
#   breakpoints = 200,
#   xlab = "EFFECT1", ylab = "EFFECT2",
#   cex.lab = 1, cex.axis = 1,
#   line.lty = 1, region.lty = 2,
#   line.lwd = 1, region.lwd = 1,
#   emphasize.axes = TRUE)
# -- To visualize the classification regions of two competing 2-dim linemodels.

# visualize.scales(scales, scale.weights = c(1),
#                  model.names = NULL, model.cols = NULL,
#                  legend.position = "topleft")
# -- To visualize the univariate distributions of the effect sizes.

# line.models(X, SE,
#             scales, slopes, cors,
#             model.names = NULL,
#             model.priors = rep(1/length(scales), length(scales)),
#             r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
#             scale.weights = c(1))
# -- To evaluate the model probabilities for each observation separately.

# line.models.with.proportions(X, SE,
#                              scales, slopes, cors,
#                              model.names = NULL,
#                              r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
#                              n.iter = 200, n.burnin = 20,
#                              diri.prior = rep(1/length(scales),length(scales)),
#                              verbose = TRUE)
# -- To evaluate the model probabilities together with the proportions of
#    observations coming from each model. Uses Gibbs sampler for a joint
#    estimation of membership probabilities across all observations.

# line.models.optimize(X, SE,
#                      par.include = NULL,
#                      force.same.scales = FALSE,
#                      init.scales, init.slopes, init.cors,
#                      model.names = NULL,
#                      model.priors = rep(1,length(init.scales)),
#                      r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
#                      tol.loglk =  1e-3,
#                      tol.par = 0,
#                      op.method = "BFGS",
#                      assume.constant.SE = FALSE,
#                      print.steps = 1)
# -- To use EM-algorithm to optimize chosen parameters of the line models.

# sample.line.model(n = 1, scale, slope, cor, scale.weights = c(1))
# -- To generate samples from a line model.

# simulate.linemodels.for.observations(
#   X, SE, linemodel.prob = NULL,
#   scales, slopes, cors,
#   r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2))
# -- To generate realistic samples from a set of line models mimicking the observed data.

# simulate.loglr(
#   X, SE, n.sims = 100,
#   par.include.null = matrix(TRUE,
#                             nrow = length(init.scales.null),
#                             ncol = 3),
#   force.same.scales.null = FALSE,
#   init.scales.null, init.slopes.null, init.cors.null,
#   model.priors.null = rep(1, length(init.scales.null)),
#   par.include.alt = matrix(TRUE,
#                            nrow = length(init.scales.alt),
#                            ncol = 3),
#   force.same.scales.alt = FALSE,
#   init.scales.alt, init.slopes.alt, init.cors.alt,
#   model.priors.alt = rep(1,length(init.scales.alt)),
#   r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
#   tol.loglk = 1e-3,
#   tol.par = 0,
#   op.method = "BFGS",
#   assume.constant.SE = FALSE,
#   print.steps = c(1,0))
# -- To compute log-likelihood ratio between alternative and null and estimate its null distribution


################################################################################
############## FUNCTIONS #######################################################
################################################################################


rotate.diagonal.to.line <- function(slopes){
  # Determines the rotation that takes the line of vector (1,1,...,1) to
  # the line of the vector y = (1, slope[1], ... slope[M-1]).
  # Returns the M-dimensional rotation matrix.

  M = length(slopes) + 1
  x = rep(1, M)/sqrt(M) #unit vector of diagonal
  y = c(1, slopes)
  if(any(is.infinite(slopes))){ #allow use of Inf and -Inf in slopes
    y[is.finite(y)] = 0
    y[is.infinite(y)] = sign(y[is.infinite(y)])
  }
  y = y/sqrt(sum(y^2)) #unit vector of direction of y
  v = (y - sum(x*y) * x)
  v = v/sqrt(sum(v^2)) #unit vector on the plane (x,y) orthogonal to x
  cost = sum(x*y) #cosine between x and y
  if(cost > 1) cost = 1
  if(cost < (-1)) cost = -1
  sint = sqrt(1 - cost^2) #choose positive sign for the sine

  #Rotation of angle t in 2D space
  R = matrix(c(cost, -sint,
               sint, cost),
             ncol = 2, byrow = T)

  A = matrix(c(x,v), ncol = 2)
  # Q is the rotation of line x to line y in M-dim space.
  # Projection orthogonal to both u and v will stay as it is
  # but the projection on the plane defined by u and v will be rotated by matrix R.
  Q = diag(M) - x %*% t(x) - v %*% t(v) + A %*% R %*% t(A)
  return(Q)
}



prior.V <- function(scale = 0, slopes = 1, cor = 0){

 # Returns prior covariance matrix for a linemodel defined by
 #  scale, slopes and cor.
 # INPUT
 # scale > 0
 # slopes, vector with elements in [-Inf,Inf]
 #        element 'i' is the slope between effect variables 'i+1' and 1
 #        While +/-Inf allowed, better to use finite values if there is more than 1 slope.
 # cor >= 0, NOTE: negatively correlated effects are modeled by a negative slope,
 #                 not by a negative correlation

  if(cor > 1) stop("prior.V: cor > 1")
  if(cor < 0) stop("prior.V: cor < 0")
  if(scale < 0) stop("prior.V: scale < 0")
  M = length(slopes) + 1

  Q = rotate.diagonal.to.line(slopes)
  C = matrix(cor, M, M)
  diag(C) = 1

  #Consider distribution MVN(0,C) and rotate it by Q to get variance matrix S
  S = Q %*% C %*% t(Q)
  S = 0.5*(S + t(S)) #Make symmetric to account for possible rounding errors
  S = scale^2 * S / max(as.vector(S)) #largest component has the given scale
  return(S)
}



rdirichlet <- function(alpha){

  # Random sampling from Dirichlet distribution with parameter vector alpha.

  g = rgamma(length(alpha), shape = alpha, scale = 1)
  return(g / sum(g))
}



log.dmvnorm <- function(x, mu = NULL, S = NULL ){

  # Returns log of density of MV-Normal(mean = mu, var = S) at x
  # x can be either vector or matrix (evaluation is done for vectors in columns).
  # Outputs a vector of density values for each column of input x.

  if(is.vector(x)) x = matrix(as.numeric(x), ncol = 1)
  if(nrow(x) == 1) x = matrix(as.numeric(x), ncol = 1) #matrix has only one row, turn into a column vector
  K = nrow(x)
  if(is.null(mu)) mu = rep(0, K)
  if(is.null(S)) S = diag(1, K)
  if(length(mu) != K) stop("log.dmvnorm: the length of mu doesn't match the dimension of x")
  if(!all(dim(S) == K)) stop("log.dmvnorm: the dimension of S doesn't match the dimension of x")

  chol.S = chol(S) # Cholesky decomposition
  log.det = 2*sum(log(diag(chol.S))) # log of det(S)
  inv.chol.S = solve(t(chol.S)) # inverse of cholesky^T
  A = inv.chol.S %*% (x - mu)
  return(-K/2*log(2*pi) - 0.5*(log.det + colSums(A*A)))
}



check.input <- function(msg, X, SE,
                        scales, slopes, cors,
                        model.names = NULL,
                        model.priors = NULL,
                        r.lkhood = NULL,
                        scale.weights = NULL){
  # Checks common input variables for consistency and outputs
  # updated values of each of them.
  # Additionally, returns values of K (no. of linemodels) and M (dimension).
  # Input parameter 'msg' is used to specify from which function this was called
  # so that the error messages are informative.

  M = ncol(X) #number of dimensions
  if(M < 2) stop(paste(msg,"Input data 'X' must have at least two columns."))
  if(ncol(SE) != M) stop("Input data 'SE' must have same no. of columns as X.")
  if(nrow(X) != nrow(SE)) stop("Input data 'X' and SE' must have the same no. of rows")

  if(!is.vector(scales) || length(scales) < 1){
    stop(paste(msg,"scales should be vector of length >= 1"))}
  K = length(scales) #number of models

  if(!is.vector(slopes) && !is.matrix(slopes)){
    stop(paste(msg,"slopes should be a vector or a matrix"))}
  if(!is.matrix(slopes)) slopes = matrix(slopes, ncol = 1)
  if(nrow(slopes) != K) stop(paste(msg,"Length of 'scales' and number of rows of 'slopes' do not match."))
  if(ncol(slopes) != (M-1)) stop(paste(msg,"No. of columns of 'slopes' should be ncol(X) - 1."))

  if(length(cors) != K) stop(paste(msg,"Length of 'scales' and 'cors' do not match."))
  if(any(cors > 1) || any(cors < 0)) stop(paste(msg,"Some value of 'cors' is outside [0,1]."))
  if(!is.null(model.priors)){
    if(length(model.priors) != K) stop(paste(msg,"Length of 'model.priors' and 'scales' do not match."))
    if(any(model.priors < 0)) stop(paste(msg,"All 'model.priors' must be non-negative."))
    model.priors = model.priors/sum(model.priors)
  }

  if(!is.null(model.names) && length(model.names) != K) {
    stop(paste(msg,"Length of 'model.names' and 'scales' do not match."))}
  if(is.null(model.names)) model.names = paste0("Mo",1:K)

  if(is.vector(r.lkhood)){
    if(length(r.lkhood) != M*(M-1)/2) stop(paste(msg,"When 'r.lkhood' is vector, it must give the upper triangular part of the correlation matrix.
                                           Thus, its length must by M*(M-1)/2, where M is the number of columns of 'X'."))
    tmp = matrix(NA, M, M)
    tmp[upper.tri(tmp)] = tmp[lower.tri(tmp)] = r.lkhood
    diag(tmp) = 1
    r.lkhood = tmp
  }
  if(max(abs(r.lkhood - t(r.lkhood))) > 1e-10) stop(paste(msg,"Matrix 'r.lkhood' must be symmetric."))
  if(max(abs(diag(r.lkhood)-1)) > 1e-10) stop(paste(msg,"Diagonal of 'r.lkhood' must be 1."))
  if(any(r.lkhood < (-1)) || any(r.lkhood > 1)) stop(paste(msg,"Some value of 'r.lkhood' is outside [-1,1]."))

  if(!is.null(scale.weights)){
    if(is.vector(scale.weights)){
      scale.weights = matrix(scale.weights, byrow = T,
                           ncol = length(scale.weights), nrow = K)}
    if(!is.matrix(scale.weights)) stop("'scale.weights' should be vector or matrix.")
    if(any(scale.weights < 0)) stop(paste(msg,"No value of scale.weights can be negative."))
    if(nrow(scale.weights) != K) stop("Dimensions of 'scales' and 'scale.weights' do not match.")
  }

  return(list(K = K, M = M, slopes = slopes, model.names = model.names,
              model.priors = model.priors, r.lkhood = r.lkhood,
              scale.weights = scale.weights))
  }



#' Compute correlation in likelihood based on case and control overlaps between 2 studies
#' Formula is from Bhattacharjee et al. (2012)
#' Subset-Based Approach Improves Power and Interpretation for the Combined Analysis
#' of Genetic Association Studies of Heterogeneous Traits,
#' The American Journal of Human Genetics, 90(5): 821-835.
#'
#' @param s1 cases of study 1
#' @param r1 controls of study 1
#' @param s2 cases of study 2
#' @param r2 controls of study 2
#' @param s1s2 cases in 1 and cases in 2
#' @param r1r2 controls in 1 and controls in 2
#' @param s1r2 cases in 1 and controls in 2
#' @param r1s2 controls in 1 and cases in 2
#' @return correlation value
#' @examples
#' beta.cor.case.control(s1 = 1000, r1 = 1000, s2 = 1000, r2 = 2000, s1s2 = 0, r1r2 = 1000)
#' @export
beta.cor.case.control <- function(s1, r1, s2 ,r2, s1s2 = 0, r1r2 = 0, s1r2 = 0, r1s2 = 0){
  return(sqrt(s1*r1/(s1+r1))*sqrt(s2*r2/(s2+r2))*(s1s2/s1/s2 - s1r2/s1/r2 - r1s2/r1/s2 + r1r2/r1/r2))
}



#' Returns slopes between two effects ii and jj,
#' that is, the slopes 'b' of the line x_ii = b*x_jj
#'
#' @param ii index of "y-axis" variable
#' @param jj index of "x-axis" variable
#' @param slopes matrix of size Kx(M-1), where column 'i' gives the
#' slopes between effects 'i+1' and 1 or a vector of length M-1
#' M is the number of dimensions
#' K is the number of linemodels
#' @return slopes for the pair specified for all K models
#' @examples
#' slope.for.pair(3,2, c(2,3))
#' slope.for.pair(3,2, matrix(c(2,3, 0.5,-1), byrow = TRUE, ncol = 2))
#' @export
slope.for.pair <- function(ii, jj, slopes){

  if(is.vector(slopes)) slopes = matrix(slopes, nrow = 1)
  if(!is.matrix(slopes)) stop("slope.for.pair: input 'slopes' must be a vector or a matrix")
  M = ncol(slopes) + 1
  K = nrow(slopes)
  if(ii < 1 || ii > M) stop("slope.for.pair: 1st index not in [1,length(slopes)+1]")
  if(jj < 1 || jj > M) stop("slope.for.pair: 2nd index not in [1,length(slopes)+1]")
  if(ii == 1) {sl.ii = rep(1,K)} else {sl.ii = slopes[,ii-1]}
  if(jj == 1) {sl.jj = rep(1,K)} else {sl.jj = slopes[,jj-1]}
  sl = sl.ii/sl.jj
  sl[abs(sl.ii) < 1e-10 & abs(sl.jj) < 1e-10] = 0
  ind = is.infinite(sl.ii) & is.infinite(sl.jj)
  sl[ind] = sign(sl.ii[ind])*sign(sl.jj[ind])
  return(sl)
}



#' Visualize linemodels
#'
#' Plots lines and 95\% highest probability regions of the 2 dimensional models
#' specified by scales, slopes and correlations.
#' See help of 'line.models( )' for details about specifying the models.
#'
#' @param scales vector of standard deviations of larger effect of each linemodel
#' @param slopes vector of slopes of each linemodel
#' @param cors vector of correlation parameters of each linemodel
#' @param model.names vector of names of each linemodel
#' @param model.cols vector of colors for each linemodel
#' @param legend.position position of legend, value NULL removes legend
#' @param xlim,ylim,xlab,ylab,cex.lab,cex.axis standard plotting parameters
#' @param line.lty plotting type for lines (1 = solid, 2 = dashed ...)
#' @param line.lwd width of lines
#' @param region.lty plotting type for 95\% regions (1 = solid, 2 = dashed ...)
#' @param region.lwd width of line for 95\% regions
#' @param plot.grid if TRUE, plots a grid
#' @param plot.new if TRUE, makes a new plot, if FALSE adds to existing plot
#' @param emphasize.axes if TRUE, coordinate axes are marked with a black line
#' @return none
#' @examples
#' visualize.line.models(c(0.1,0.2), c(1,0), c(0.995,0))
#' @export
# Note: Does NOT allow specification of 'scale.weights'
# but uses the given scales as the single component for each distribution.
visualize.line.models <- function(scales, slopes, cors,
                                  model.names = NULL, model.cols = NULL,
                                  legend.position = "bottomright",
                                  xlim = NULL, ylim = NULL,
                                  xlab = "EFFECT1", ylab = "EFFECT2",
                                  cex.lab = 1, cex.axis = 1,
                                  line.lty = 1, region.lty = 2,
                                  line.lwd = 1, region.lwd = 1,
                                  plot.grid = TRUE,
                                  plot.new = TRUE,
                                  emphasize.axes = TRUE){

  K = length(scales) #number of models
  if(length(slopes) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(length(cors) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(any(cors > 1) | any(cors < 0)) stop("Some value of 'cors' is outside [0,1].")
  lim = 3*max(scales) #Default: show models within 3 SDs
  if(is.null(xlim)) xlim = c(-lim,lim)
  if(is.null(ylim)) ylim = c(-lim,lim)
  if(plot.new){
    plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
         cex.lab = cex.lab, cex.axis = cex.axis)}
  if(is.null(model.names)) model.names = paste0("M",1:K)
  if(is.null(model.cols)) model.cols = 1:K
  slopes[is.nan(slopes)] = 0 #NaN is when 2 slopes are 0 w.r.t third effect
  prob.level = 0.95
  b = qchisq(prob.level, df = 2)
  if(plot.grid) grid()
  if(emphasize.axes){
    abline(h = 0, lwd = 1)
    abline(v = 0, lwd = 1)
  }

  for(ii in 1:K){
    if(scales[ii] < 1e-16) {
      points(0, 0, col = model.cols[ii], pch = 19)
      next
    }

    if(!is.finite(slopes[ii])) {
      abline(v = 0, col = model.cols[ii], lty = line.lty, lwd = line.lwd)}
    else{
      abline(0, slopes[ii], col = model.cols[ii], lty = line.lty, lwd = line.lwd)}
    if(cors[ii] < 1){
      Sigma = prior.V(scale = scales[ii], slopes = slopes[ii], cor = cors[ii])
      a = as.numeric(t(solve(Sigma)))
      x.lim = sqrt(abs(4*a[4]*b/( (a[2]+a[3])^2 - 4*a[4]*a[1])))
      x = seq(-x.lim*0.999, x.lim*0.999, length = 500)
      y.upper = (-(a[2]+a[3])*x + sqrt((a[2]+a[3])^2*x^2 - 4*a[4]*(a[1]*x^2-b)))/2/a[4]
      y.lower = (-(a[2]+a[3])*x - sqrt((a[2]+a[3])^2*x^2 - 4*a[4]*(a[1]*x^2-b)))/2/a[4]
      lines(x, y.upper, col = model.cols[ii], lty = region.lty, lwd = region.lwd)
      lines(x, y.lower, col = model.cols[ii], lty = region.lty, lwd = region.lwd)
    }
  }
  if(!is.null(legend.position)){
    legend(legend.position, pch = 15, leg = model.names, col = model.cols)}
}



#' Visualize multidimensional linemodels
#'
#' Plots lines and 95\% highest probability regions of the multi-dimensional models
#' specified by scales, slopes and correlations.
#' See help of 'line.models( )' for details about specifying the models.
#'
#' @param scales vector of standard deviations of larger effect of each model
#' @param slopes matrix of slopes of the models (one row per model)
#' @param cors vector of correlation parameters of each model
#' @param plot.pairs matrix of pairs of effect variables to be plotted (2 cols and 1 row per pair)
#' @param var.names names used as axis labels for the effect variables
#' @param model.names vector of names of each linemodel
#' @param model.cols vector of colors for each linemodel
#' @param legend.position position of legend, value NULL removes legend
#' @param X points to be plotted, 1 col per effect variable, if NULL no points plotted
#' @param SE standard errors to be plotted around points, 1 col per effect, if NULL no SEs plotted
#' @param SE.mult multiplier of SE to be plotted around each point, default 1.96 i.e. 95\% conf. intervals
#' @param model.prob matrix of probability of each observation (row) in each model (col), can be NULL
#' @param model.thresh threshold probability to color a point by model.cols, default = 0.8
#' @param undetermined.col color for points where no model.prob is > model.thresh
#' @param xlim,ylim,cex,pch,cex.lab,cex.axis standard plotting parameters
#' @param line.lty plotting type for lines (1 = solid, 2 = dashed ...)
#' @param line.lwd width of lines
#' @param region.lty plotting type for 95\% regions (1 = solid, 2 = dashed ...)
#' @param region.lwd width of line for 95\% regions
#' @param plot.grid if TRUE, plots a grid
#' @param emphasize.axes if TRUE, coordinate axes are marked with a black line
#' @return none
#' @examples
#' par(mfrow = c(1,2))
#' visualize.multidimensional.line.models(c(0.1, 0.2), matrix(c(1,0, 0,1), ncol = 2), c(0.99, 0.95),
#' plot.pairs = matrix(c(1,2, 1,3),ncol = 2, byrow = TRUE))
#' @export
# Note: Does NOT allow specification of 'scale.weights'
# but uses the given scales as the single component for each distribution.
visualize.multidimensional.line.models <- function(scales, slopes, cors,
                                                   plot.pairs = NULL, var.names = NULL,
                                                   model.names = NULL, model.cols = NULL,
                                                   legend.position = "bottomright",
                                                   X = NULL, SE = NULL,
                                                   SE.mult = 1.96,
                                                   model.prob = NULL, model.thresh = 0.8,
                                                   undetermined.col = "gray",
                                                   pch = 1, cex = 1,
                                                   xlim = NULL, ylim = NULL,
                                                   cex.lab = 1, cex.axis = 1,
                                                   line.lty = 1, region.lty = 2,
                                                   line.lwd = 1, region.lwd = 1,
                                                   plot.grid = TRUE,
                                                   emphasize.axes = TRUE){
  M = ncol(slopes) + 1
  K = length(scales)
  if(nrow(slopes) != K) stop("No. of rows of slopes should equal to length of scales")
  if(is.vector(plot.pairs)) plot.pairs = matrix(plot.pairs, nrow = 1)
  if(ncol(plot.pairs) != 2) stop("plot.pairs should have 2 columns")
  if(any(as.numeric(plot.pairs) < 1) || any(as.numeric(plot.pairs) > M))
    stop("Values in plot.pairs should be between 1 and M, dimension of data.")
  if(is.null(var.names)) var.names = paste0("EFFECT",1:M)
  if(length(var.names) != M) stop("Length of var.names should equal to ncol(slopes) + 1")

  if(!is.null(X)){
    if(ncol(X) != M) stop("No. of cols of X should be 1 + no. of cols of slopes")
    cols = rep(undetermined.col, nrow(X))
    if(!is.null(model.prob)){
      if(nrow(X) != nrow(model.prob)) stop("No. of rows of X and model.prob must match.")
      if(ncol(model.prob) != K) stop("No. of cols of model.prob must equal to the length of scales.")
      for(jj in 1:K){cols[model.prob[,jj] > model.thresh] = model.cols[jj]}
    }
    if(!is.null(SE)){
      if(nrow(X) != nrow(SE)) stop("No. of rows of X and SE must match.")
      if(ncol(X) != ncol(SE)) stop("No. of cols of X and SE must match.")
      if(any(SE < 0)) stop("SEs must be non-negative")
    }
  }

  for(ii in 1:nrow(plot.pairs)){
    i1 = plot.pairs[ii,1]
    i2 = plot.pairs[ii,2]
    visualize.line.models(scales = scales,
                          slopes = slope.for.pair(i2, i1, slopes),
                          cors = cors,
                          model.names = model.names, model.cols = model.cols,
                          legend.position = legend.position,
                          xlim = xlim, ylim = ylim,
                          xlab = var.names[i1], ylab = var.names[i2],
                          cex.lab = cex.lab, cex.axis = cex.axis,
                          line.lty = line.lty, region.lty = region.lty,
                          line.lwd = line.lwd, region.lwd = region.lwd,
                          plot.grid = plot.grid,
                          plot.new = TRUE,
                          emphasize.axes = emphasize.axes)
    if(!is.null(X)) points(X[,i1], X[,i2],
                           col = cols, cex = cex, pch = pch)
    if(!is.null(SE)){
      arrows(X[,i1]-SE.mult*SE[,i1], X[,i2],
             X[,i1]+SE.mult*SE[,i1], X[,i2],
             col = cols, code = 3, angle = 90, length = 0)
      arrows(X[,i1], X[,i2]-SE.mult*SE[,i2],
             X[,i1], X[,i2]+SE.mult*SE[,i2],
             col = cols, code = 3, angle = 90, length = 0)
    }
  }
}



#' Visualize classification regions of 2-dimensional line models
#'
#' Colors the plane according to the classification probabilities
#' of two competing line models.
#' See help of 'line.models( )' for details about specifying the models.
#'
#' @param scales vector of standard deviations of larger effect of linemodels
#' @param slopes vector of slopes of linemodels
#' @param cors vector of correlation parameters of linemodels
#' @param SE vector of two standard errors for effect1 and effect2, respectively
#' @param r.lkhood correlation between estimators of the two values
#' @param model.priors vector of (unnormalized) prior probabilities of the two linemodels
#' @param col.palette vector of seven colors (with reasonable defaults)
#' @param add.legend if TRUE, adds legend on right-hand side of the plot
#' @param breakpoints how many breakpoints for defining the grid on each axis
#' @param xlim,ylim,xlab,ylab,cex.lab,cex.axis standard plotting parameters
#' @param line.lty plotting type for lines (1 = solid, 2 = dashed ...)
#' @param line.lwd width of lines
#' @param region.lty plotting type for 95\% regions (1 = solid, 2 = dashed ...)
#' @param region.lwd width of line for 95\% regions
#' @param emphasize.axes if TRUE, coordinate axes are marked with a black line
#' @return none
#' @examples
#' visualize.classification.regions(c(0.1,0.2), c(1,0), c(0.995,0.995), c(0.05,0.05))
#' @export
# Note: Does NOT allow specification of 'scale.weights'
# but uses the given scales as the single component for each distribution.
visualize.classification.regions <-
  function(scales, slopes, cors, SE,
           r.lkhood = 0,
           model.priors = rep(1,length(scales)),
           col.palette = NULL,
           add.legend = TRUE,
           xlim = NULL, ylim = NULL,
           breakpoints = 200,
           xlab = "EFFECT1", ylab = "EFFECT2",
           cex.lab = 1, cex.axis = 1,
           line.lty = 1, region.lty = 2,
           line.lwd = 1, region.lwd = 1,
           emphasize.axes = TRUE){

  K = length(scales) #number of models
  if(K != 2) stop("Should specify exactly two models.")
  if(length(slopes) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(length(cors) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(any(cors > 1) | any(cors < 0)) stop("Some value of 'cors' is outside [0,1].")
  if(length(SE) != 2 | any(SE <= 0)) stop("SE should have exactly two positive values.")
  if(!is.null(col.palette) & length(col.palette) !=7) stop("col.palette should have exactly 7 colors")
  if(length(breakpoints) != 1 | breakpoints[1]<2) stop("Value of brekpoints must be a single integer > 1")
  lim = 3*max(scales) #Default: show models within 3 SDs
  if(is.null(xlim)) xlim = c(-lim,lim)
  if(is.null(ylim)) ylim = c(-lim,lim)
  if(is.null(col.palette)) col.palette = c("blue","skyblue","gray","salmon","red","white","black")

  K = breakpoints
  x = cbind(
    rep(seq(xlim[1], xlim[2], length = K), K),
    rep(seq(ylim[1], ylim[2], length = K), each = K),
    rep(SE[1], K^2),
    rep(SE[2], K^2))

  res = line.models(
    X = x[,1:2],
    SE = x[,3:4],
    scales = scales,
    slopes = slopes,
    model.priors = model.priors,
    cors = cors,
    r.lkhood = r.lkhood
    )

  cols = rep(col.palette[1], nrow(res))
  cols[res[,1] < 0.9] = col.palette[2]
  cols[res[,1] < 0.6] = col.palette[3]
  cols[res[,1] < 0.4] = col.palette[4]
  cols[res[,1] < 0.1] = col.palette[5]

  if(add.legend) layout(matrix(c(1,2),nrow = 1), widths = c(3.5, 1))

  visualize.line.models(scales, slopes, cors,
                        xlim = xlim, ylim = ylim,
                        legend.position = NULL, model.cols = col.palette[c(6,7)],
                        xlab = xlab, ylab = ylab,
                        cex.lab = cex.lab, cex.axis = cex.axis)

  points(x[,1], x[,2], pch = 15, cex = 0.7, col = cols)

  visualize.line.models(scales, slopes, cors,
                        xlim = xlim, ylim = ylim,
                        legend.position = NULL, plot.new = FALSE,
                        model.cols = col.palette[c(6,7)],
                        xlab = xlab, ylab = ylab,
                        plot.grid = FALSE,
                        cex.lab = cex.lab, cex.axis = cex.axis,
                        line.lwd = line.lwd, region.lwd = region.lwd)

  if(add.legend){
    par(mar = c(0,0,0,0))
    plot.new()
    text(0.1, 0.80, expression(M[1]))
    text(0.5, 0.80, paste0("scale ",signif(scales[1],3),
                        "\nslope ",signif(slopes[1],3),
                        "\ncor ",signif(cors[1],3)), cex = 0.7)
    text(0.1, 0.62, expression(M[2]))
    text(0.5, 0.62, paste0("scale ",signif(scales[2],3),
                         "\nslope ",signif(slopes[2],3),
                         "\ncor ",signif(cors[2],3)), cex = 0.7)

    legend(0.0, 0.5, leg = c(expression(paste(M[1],epsilon,"[0.9,1.0]")),
                             expression(paste(M[1],epsilon,"[0.6,0.9)")),
                             expression(paste(M[1],epsilon,"[0.4,0.6)")),
                             expression(paste(M[1],epsilon,"[0.1,0.4)")),
                             expression(paste(M[1],epsilon,"[0.0,0.1)"))),
           pch = 15, col = col.palette[1:5], cex = 0.8)

    text(0.5, 0.15,paste0("SE = (",signif(SE[1],3),", ",signif(SE[2],3),")"), cex = 0.8)
  }
}



#' Visualize effect distributions of line models
#'
#' Plots density functions of distributions of effect sizes
#' of line models specified by scales and scale.weights.
#' See help of 'line.models( )' for details about specifying the models.
#'
#' @param scales vector of standard deviations of larger effect of each linemodel
#' @param scale.weights vector of weights of each component or a matrix
#'  with one row per linemodel and one column per component
#' @param model.names vector of names of each linemodel
#' @param model.cols vector of colors for each linemodel
#' @param legend.position position of legend, value NULL removes legend
#' @return none
#' @examples
#' visualize.scales(c(0.1,0.2), c(1,0,0,1))
#' visualize.scales(c(0.1,0.2),
#'                  rbind(c(1,0,0,1),
#'                        c(3,0,0,1)))
#' @export
visualize.scales <- function(scales, scale.weights = c(1),
                             model.names = NULL, model.cols = NULL,
                             legend.position = "topleft"){

  # Plots the univariate distributions of the effect sizes
  # Input parameters 'scales' and 'scale.weights' as in line.models()
  # Can take in user given colors 'model.cols' but there must be equal
  # number of colors as there are models specified by scales vector.
  # Set legend.position to NULL to remove the legend.

  K = length(scales)
  if(any(scales <= 0)) stop("'scales' must be positive.")
  if(any(scale.weights < 0)) stop("'scale.weights' must be non-negative.")
  x.ran = 2.5 * max(scales)

  if(is.vector(scale.weights)){
    scale.weights = matrix(scale.weights, byrow = T,
                              ncol = length(scale.weights), nrow = K)
  }
  if(!is.matrix(scale.weights)) stop("'scale.weights' should be vector or matrix.")
  if(nrow(scale.weights) != K) stop("Dimensions of 'scales' and 'scale.weights' do not match.")
  scale.weights = scale.weights / rowSums(scale.weights)
  n.comps = ncol(scale.weights)
  x = seq(-x.ran, x.ran, length = 1000)
  if(is.null(model.cols)) model.cols = c(2:(K+1))
  if(length(model.cols) != K) stop(paste("Length of vector of model.cols must equal the number of models, here",K))

  y.max = 0
  for(draw in c(FALSE, TRUE)){
    if(draw) {
      plot(NULL, xlim = c(-x.ran, x.ran), ylim = c(0, y.max),
           xlab = "effect", ylab = "density")
      if(is.null(model.names)) model.names = paste0("M",1:K)
      if(!is.null(legend.position)) legend(legend.position, col = model.cols,
                                           leg = model.names, lty = 1, lwd = 1.5)
      cat(paste("Plotted the following effect distributions:\n"))
      }
    for(kk in 1:K){
      y = rep(0, length(x))
      for(jj in 1:n.comps){
        y = y + scale.weights[kk,jj]*dnorm(x, 0, scales[kk]/(n.comps - jj + 1))}
      if(draw) {
        lines(x, y, col = model.cols[kk], lwd = 1.5)
        if(scales[kk] < 1e-300) points(0, 0, col = model.cols[kk], pch =  19)
        cat(paste0(model.names[kk],": Scale = ",scales[kk],
                  " weights = ",paste0(scale.weights[kk,],collapse = ","),"."),"\n")
        }
      y.max = max(c(y, y.max))
    }
  }
}



#' Membership probabilities in line models
#'
#' Consider M effect sizes.
#' The linemodels are defined by three parameters:
#' scale = standard deviation of the (largest) effect,
#' slope = set of M-1 slopes defining the line around which the effects are scattered,
#' cor = non-negative pairwise correlation between each effect pair,
#' where cor = 1 means the same effect and cor = 0 means independent effects.
#'
#' A linemodel has the following properties:
#'
#' 1) Effects are scattered around line defined by vector (1,slope_1,...,slope_(M-1)).
#'    Each slope_i is the slope between effect variables i+1 and 1.
#'    Slope can be any real number or +/-Inf.
#'    If any slope is infinite, the effect of variable 1 and variables with finite slopes are 0.
#'    It is recommended to keep the slopes finite to avoid problems with interpretation.
#'
#' 2) The largest of the prior variances of the effects is scale^2.
#'    The effect that has the largest scale is determined by the slopes,
#'    and the scales of the remaining effects are determined by the slopes as well.
#'
#' 3) Distribution around the line is determined as follows.
#'    Consider a distribution scattered around the direction of the vector (1,1,...,1)
#'    with the given constant value of the pairwise correlations between the variables.
#'    Rotate that distribution by an orthogonal rotation that rotates the direction of the
#'    vector (1,1,...,1) to the direction of the vector (1,slope_1,...,slope_(M-1))
#'    and use the corresponding distribution, scaled so that the maximum variance
#'    among the effects is set to scale^2.
#'    NOTE: This is not same as having the given correlation around line defined by the slopes,
#'          because the shape of that distribution depends on the slopes but the
#'         definition we use here is independent of slopes up to an orthogonal transformation.
#'
#' Examples:
#' Set scale = 0 to get the null model. (Values of slope and cor do not matter in this case.)
#' Use scale > 0, slope = 1 and cor = 0 to get independent effects model for 2-dim case.
#' Use scale > 0, slope = c(1,1,1) and cor = 1 to get fixed effects model for 4-dim case.
#' Use scale > 0, slope = 1 and 0 < cor < 1  to get correlated effects model for 2-dim case.
#'               note that cor value needs to be quite near 1 to get very similar values
#'               between the two effects with high probability (e.g. cor = 0.99...0.999).
#'
#' Effect distribution for the largest effect can be specified as a mixture
#' of Gaussians.
#' If 'scale.weights' is a vector of length C > 1, then model 'Mo' is a mixture
#' of C M-dimensional Gaussians where mixture component i has mean of 0,
#' the standard deviation for the largest effect of 'scales[Mo]/(C - i + 1)'
#' and the mixture weight proportional to scale.weights[Mo].
#' If 'scale.weights' is a matrix with C > 1 columns and K rows,
#' then model 'Mo' is a mixture of 'C' Gaussians specified by row 'Mo'
#' of 'scale.weights'.
#'
#' @param X matrix of estimates with M columns (one col per effect variable and one row per observation)
#' @param SE matrix of standard errors with M columns (one col per effect variable and one row per observation)
#' @param scales K-vector of standard deviations of the largest effect, one per each linemodel
#' @param slopes Kx(M-1)-matrix of slopes, one row per linemodel
#'  If M = 2, can be given also as a K-vector of slopes, one per each linemodel
#' @param cors K-vector of correlation parameters, one per each linemodel
#' @param model.names K-vector of names of linemodels
#' @param model.priors K-vector of prior probabilities of linemodels up to a normalization constant
#' @param r.lkhood correlation matrix of the estimators of the M effects
#'  can be given also as a vector of the upper triangular values in the row-major order
#' @param scale.weights vector of weights of each mixture component or a matrix
#'  with one row per linemodel and one column per mixture component
#' @return matrix with one column per linemodel and one row per observation giving
#' membership probabilities of each observation in each linemodel
#' @examples
#' #2D example
#' line.models(X = linemodels.ex1[,c("beta1","beta2")], SE = linemodels.ex1[,c("se1","se2")],
#'            scales = c(0.2, 0.2, 0.2),
#'            slopes = c(0, 0.5, 1),
#'            cors = c(0.995, 0.995, 0.995))
#' #3D example
#' line.models(X = linemodels.ex1[,c("beta1","beta2","beta3")],
#'            SE = linemodels.ex1[,c("se1","se2","se3")],
#'            scales = c(0.2, 0.2, 0.2),
#'            slopes = matrix(c(0,0, 0.5,0, 1,1), byrow = TRUE, ncol = 2),
#'            cors = c(0.995, 0.995, 0.995))
#' @export
line.models <- function(X, SE,
                        scales, slopes, cors,
                        model.names = NULL,
                        model.priors = rep(1/length(scales), length(scales)),
                        r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
                        scale.weights = c(1)){

  # Evaluates model probabilities for each data point separately.
  # INPUT
  # X, N x M matrix where rows are observations/effect variables
  #  and columns are the observed effect sizes.
  # SE, N x M matrix where rows are observations/effect variables
  #  and columns are the standard errors of effects in X.
  # scales, vector of positive standard deviations of the largest effect size
  #         of each linemodel. Length is K, the number of linemodels.
  # slopes, matrix of slopes of the linemodels. One row per model.
  #         M-1 columns where col 'i' specifies the slope between effects 'i' and 1 (i.e. x_i/x_1).
  #         If M = 2, 'slopes' can be given as a vector, one element per model.
  # cors, vector of correlation parameters of the linemodels. Length is K.
  # model.names, names to be used in output for the models. Length is K.
  # model.priors, prior probability of each of the K linemodels.
  #               Default is 1/K for each linemodel.
  # r.lkhood, determines correlation of the effect estimators.
  #           Default is 0, which assumes that estimators were uncorrelated.
  #           Can be given either as a symmetric correlation matrix or
  #           as a vector containing the upper triangle of the matrix in the row-major order.
  # scale.weights:
  # Effect distribution for the largest effect can be specified as a mixture
  # of Gaussians.
  # If 'scale.weights' is a vector of length C > 1, then linemodel 'Mo' is a mixture
  # of C Gaussians:
  # sum_{i=1}^C scale.weights[i]*N(0, prior.V(scales[Mo]/(C - i + 1), slopes[Mo,], cols[Mo]))
  # where vector 'scale.weights' has been normalized to sum to 1.
  # If 'scale.weights' is a matrix with C > 1 columns and K rows,
  # then linemodel 'Mo' is a mixture of C Gaussians:
  #  sum_{i=1}^k scale.weights[Mo,i]*N(0, scales[Mo]/(C - i + 1))
  # where rows of matrix 'scale.weights' has been normalized to sum to 1.
  # Default value scale.weights = c(1) means that each linemodel has only one component
  #  whose standard deviation is 'scales[Mo]'.

  # OUTPUT:
  # Returns a matrix with membership probabilities in each linemodel (in columns)
  # for each observations (in rows).

  chk = check.input("line.models:", X, SE, scales, slopes, cors, model.names,
                    model.priors, r.lkhood, scale.weights)
  K = chk$K #number of models
  M = chk$M #number of dimensions
  slopes = chk$slopes
  r.lkhood = chk$r.lkhood
  model.names = chk$model.names
  scale.weights = chk$scale.weights
  n.comps = ncol(scale.weights)
  pis = chk$model.priors

  pr.V = list() #list of lists of prior matrices of each model and each scale
  for(kk in 1:K) { #over models
    w = scale.weights[kk,]
    w = w/sum(w)
    model.Vs = list()
    for(jj in 1:n.comps){ #over scales
      model.Vs[[jj]] = prior.V(scale = scales[kk]/(n.comps + 1 - jj),
                               slopes = slopes[kk,], cor = cors[kk])}
    pr.V[[kk]] = model.Vs
  }

  n = nrow(X)
  logdnorm = matrix(NA, ncol = K, nrow = n)
  for(ii in 1:n){
    V = diag(SE[ii,]) %*% r.lkhood %*% diag(SE[ii,]) #var of likelihood
    for(kk in 1:K){
      tmp = sapply(1:n.comps, function(jj){
        log.dmvnorm(X[ii,], mu = rep(0, M), S = V + pr.V[[kk]][[jj]])})
      tmp = tmp + log(scale.weights[kk,])
      tmp.max = max(tmp)
      logdnorm[ii,kk] = log(sum(exp(tmp - tmp.max))) + tmp.max
    }
  }
  tmp = t(t(logdnorm) + log(pis))
  tmp.max = apply(tmp, 1, max)
  p = exp(tmp - tmp.max)
  res = p / rowSums(p)
  colnames(res) = model.names
  return(res)
}



#' Proportion parameters and membership probabilities of linemodels
#'
#' Consider M effect sizes.
#' The linemodels are defined by three parameters:
#' scale = standard deviation of the (largest) effect,
#' slope = set of M-1 slopes defining the line around which the effects are scattered,
#' cor = non-negative pairwise correlation between each effect pair,
#' where cor = 1 means the same effect and cor = 0 means independent effects.
#'
#' A linemodel has the following properties:
#'
#' 1) Effects are scattered around line defined by vector (1,slope_1,...,slope_(M-1)).
#'    Each slope_i is the slope between effect variables i+1 and 1.
#'    Slope can be any real number or +/-Inf.
#'    If any slope is infinite, the effect of variable 1 and variables with finite slopes are 0.
#'    It is recommended to keep the slopes finite to avoid problems with interpretation.
#'
#' 2) The largest of the prior variances of the effects is scale^2.
#'    The effect that has the largest scale is determined by the slopes,
#'    and the scales of the remaining effects are determined by the slopes as well.
#'
#' 3) Distribution around the line is determined as follows.
#'    Consider a distribution scattered around the direction of the vector (1,1,...,1)
#'    with the given constant value of the pairwise correlations between the effect variables.
#'    Rotate that distribution by an orthogonal rotation that rotates the direction of the
#'    vector (1,1,...,1) to the direction of the vector (1,slope_1,...,slope_(M-1))
#'    and use the corresponding distribution, scaled so that the maximum variance
#'    among the effects is set to scale^2.
#'    NOTE: This is not same as having the given correlation around line defined by the slopes,
#'          because the shape of that distribution depends on the slopes but the
#'         definition we use here is independent of slopes up to an orthogonal transformation.
#'
#' Examples:
#' Set scale = 0 to get the null model. (Values of slope and cor do not matter in this case.)
#' Use scale > 0, slope = 1 and cor = 0 to get independent effects model for 2-dim case.
#' Use scale > 0, slope = c(1,1,1) and cor = 1 to get fixed effects model for 4-dim case.
#' Use scale > 0, slope = 1 and 0 < cor < 1  to get correlated effects model for 2-dim case.
#'               note that cor value needs to be quite near 1 to get very similar values
#'               between the two effects with high probability (e.g. cor = 0.99...0.999).
#'
#'
#' The prior distribution of the mixture proportions is Dirichlet(diri.prior) and
#' a Gibbs sampler is used for estimating the posterior.
#'
#' @param X matrix of estimates with M columns (one col per effect variable and one row per observation)
#' @param SE matrix of standard errors with M columns (one col per effect variable and one row per observation)
#' @param scales K-vector of standard deviations of the largest effect, one per each linemodel
#' @param slopes Kx(M-1)-matrix of slopes, one row per linemodel
#'  If M = 2, can be given also as a K-vector of slopes, one per each linemodel
#' @param cors K-vector of correlation parameters, one per each linemodel
#' @param model.names K-vector of names of linemodels
#' @param r.lkhood correlation matrix of the estimators of the M effects
#'  can be given also as a vector of the upper triangular values in the row-major order
#' @param n.iter number of Gibbs sampler iterations (after burn-in)
#' @param n.burnin number of burn-in iterations that will be discarded
#' @param diri.prior parameters of the Dirichlet distribution used as prior for proportions
#' @param verbose if TRUE, prints the index of every 100th iteration
#' @return List with two components: (1) 'params', matrix of
#' posterior distribution of proportion parameters for each linemodel.
#' Columns are mean, lower and upper points of 95\% credible interval
#' and standard deviation.
#' (2)  'groups', matrix with one column per linemodel and one row per observation giving
#' membership probabilities of each observation in each model.
#' @examples
#' #2D example
#' line.models.with.proportions(
#'            X = linemodels.ex1[,c("beta1","beta2")], SE = linemodels.ex1[,c("se1","se2")],
#'            scales = c(0.2, 0.2, 0.2),
#'            slopes = c(0, 0.5, 1),
#'            cors = c(0.995, 0.995, 0.995))
#' #3D example
#' line.models.with.proportions(
#'            X = linemodels.ex1[,c("beta1","beta2","beta3")],
#'            SE = linemodels.ex1[,c("se1","se2","se3")],
#'            scales = c(0.2, 0.2, 0.2),
#'            slopes = matrix(c(0,0, 0.5,0, 1,1), byrow = TRUE, ncol = 2),
#'            cors = c(0.995, 0.995, 0.995))
#' @export
line.models.with.proportions <- function(X, SE,
                                         scales, slopes, cors,
                                         model.names = NULL,
                                         r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
                                         n.iter = 200, n.burnin = 20,
                                         diri.prior = rep(1/length(scales),length(scales)),
                                         verbose = TRUE){

  # Gibbs sampler to evaluate model probabilities for all data points jointly
  #  together with the posterior of the mixture proportions 'pi'.
  # INPUT
  # X, N x M matrix where rows are observations/effect variables
  #  and columns are the observed effect sizes.
  # SE, N x M matrix where rows are observations/effect variables
  #  and columns are the standard errors of effects in X.
  # scales, vector of positive standard deviations of the largest effect size
  #         of each linemodel. Length is K, the number of linemodels.
  # slopes, matrix of slopes of the linemodels. One row per model.
  #         M-1 columns where col 'i' specifies the slope between effects 'i' and 1 (i.e. x_i/x_1).
  #         If M = 2, 'slopes' can be given as a vector, one element per linemodel.
  # cors, vector of correlation parameters of the linemodels. Length is K.
  # model.names, names to be used in output for the linemodels. Length is K.
  # r.lkhood, determines correlation of the effect estimators.
  #           Default is 0, which assumes that estimators were uncorrelated.
  #           Can be given either as a symmetric correlation matrix or
  #           as a vector containing the upper triangle of the matrix in the row-major order.
  # n.iter, number of Gibbs sampler iterations after burn-in (default 200)
  # n.burnin, number of burn-in iterations that are discarded from final results (default 20).
  # diri.prior, prior distribution for proportion parameters is
  #             Dirichlet(diri.prior)
  #             Thus, 'diri.prior' should be a positive valued vector of length K.
  #             Default is diri.prior = rep(1/K, K).
  # verbose, logical, if TRUE (default), prints counter of every 100th of iteration.

  # Does not allow 'scale.weights' from line.models( ).

  # OUTPUT:
  # Returns a list of 2 matrices labelled
  # 'groups' row per observation and membership probabilities in linemodels in columns
  # 'params' posterior distribution of the proportions.

  chk = check.input("line.models.with.proportions:", X = X, SE = SE,
                    scales = scales, slopes = slopes, cors = cors,
                    model.names = model.names,
                    r.lkhood = r.lkhood)
  K = chk$K #number of models
  M = chk$M #number of dimensions
  slopes = chk$slopes
  r.lkhood = chk$r.lkhood
  model.names = chk$model.names
  if(any(diri.prior <= 0)) stop("All values in 'diri.prior' must be positive.")
  if(length(diri.prior) != K) stop("Length of 'diri.prior' do not match with number of models.")
  if(n.burnin < 0) stop("'n.burnin' must be non-negative.")
  if(n.iter < 1) stop("'n.iter' must be positive.")

  pr.V = list()
  for(kk in 1:K) pr.V[[kk]] = prior.V(scale = scales[kk], slopes = slopes[kk,], cor = cors[kk])
  n = nrow(X)

  logdnorm = matrix(NA, ncol = K, nrow = n)
  for(ii in 1:n){
    V = diag(SE[ii,]) %*% r.lkhood %*% diag(SE[ii,]) #var of likelihood
     for(kk in 1:K){
        logdnorm[ii,kk] = log.dmvnorm(X[ii,], mu = rep(0, M), S = V + pr.V[[kk]])
     }
  }

  R.par = matrix(NA, ncol = K, nrow = n.iter) #results for parameters
  R.ind = matrix(0, ncol = K, nrow = n) #results for individual effects
  colnames(R.par) = colnames(R.ind) = model.names

  pis = rep(1, K) / K #initial probabilities, vector pi
  grs = sample(1:K, size = n, prob = pis, replace = TRUE) #initially random grouping among variants

  for(ii in 1:(n.burnin + n.iter)){ #Gibbs sampler

      #count how many are assigned to each group
      tmp = rle(sort(grs))
      gr.counts = rep(0, K)
      gr.counts[tmp$values] = tmp$length

      #sample pis
      pis = rdirichlet(diri.prior + gr.counts)
      log.pis = log(pis)

      #sample group indicators for variants
      tmp = t(t(logdnorm) + log.pis)
      p = exp(tmp - apply(tmp, 1, max))
      #p = p/rowSums(p) #no need to normalize for 'sample()'
      grs = apply(p, 1, function(pr){sample(1:K, size = 1, prob = pr)})
      if(ii > n.burnin){ #save results if burn-in is over
        R.ind[(1:n) + n*(grs-1)] = 1 + R.ind[(1:n) + n*(grs-1)]
        R.par[ii - n.burnin,] = pis
      }
      if(verbose && (ii %% 100 == 0)) print(paste(ii,"iterations done."))
  }

  res.par = cbind(apply(R.par, 2, mean),
                t(apply(R.par ,2 , function(x){quantile(x, c(0.025, 0.975))})),
                apply(R.par, 2, sd))
  colnames(res.par) = c("mean","95%low","95%up","sd")

  res.ind = R.ind / n.iter
  return(list(params = res.par, groups = res.ind))
}



line.models.loglkhood <- function(X, SE,
                                  scales, slopes, cors,
                                  model.names = NULL,
                                  model.priors = rep(1/length(scales), length(scales)),
                                  r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
                                  assume.constant.SE = FALSE,
                                  return.posteriors = FALSE){

  # Evaluates log-likelihood of the data given the parameters
  #  X, SE, scales, slopes, cors, model.priors, model.names, r.lkhood
  #  are like in line.models().
  # Does not allow scale.weights from line-models().
  #
  #  Returns the log-likelihood of the mixture model
  #  where the components have the mixture weights proportional to 'model.priors'
  #  and the weights are the same for all the observations.

  # If assume.constant.SE = TRUE, then sets SE of each observation
  #   to the median SE of all the SEs of the corresponding dimension
  #   This speeds up the evaluation of likelihood.
  # If return.posteriors = TRUE, then returns a list with two members
  #  'loglkhood' and
  #  'membership.prob' matrix with
  #  posterior probabilities of observations' (rows)
  #  membership in each linemodel (columns).

  chk = check.input("line.models.loglkhood:", X = X, SE = SE,
                    scales = scales, slopes = slopes, cors = cors,
                    model.names = model.names, model.priors = model.priors,
                    r.lkhood = r.lkhood)
  K = chk$K #number of models
  M = chk$M #number of dimensions
  slopes = chk$slopes
  r.lkhood = chk$r.lkhood
  model.names = chk$model.names
  pis = chk$model.priors

  pr.V = list() #list of prior matrices of each model
  for(kk in 1:K) pr.V[[kk]] = prior.V(scale = scales[kk], slopes = slopes[kk,], cor = cors[kk])

  n = nrow(X)
  logdnorm = matrix(NA, ncol = K, nrow = n)

  if(assume.constant.SE){
    se = diag(apply(SE,2,median))
    V = se %*% r.lkhood %*% se
    for(kk in 1:K){
      logdnorm[,kk] = log.dmvnorm(t(X), mu = rep(0, M), S = V + pr.V[[kk]])}
  }else{
    for(ii in 1:n){
      V = diag(SE[ii,]) %*% r.lkhood %*% diag(SE[ii,]) #var of likelihood
      for(kk in 1:K){
        logdnorm[ii,kk] = log.dmvnorm(X[ii,], mu = rep(0, M), S = V + pr.V[[kk]])}
    }
  }

  tmp = t(t(logdnorm) + log(pis))
  tmp.max = apply(tmp, 1, max)
  log.rowsums = log(rowSums(exp(tmp - tmp.max))) + tmp.max
  log.lkhood = sum(log.rowsums)

  if(return.posteriors){
    post = exp(tmp - log.rowsums)
    colnames(post) = model.names
    return(list(loglkhood = log.lkhood, membership.prob = post))}
  else{
    return(log.lkhood)}
}



line.models.expected.loglkhood <- function(X, SE,
                                           scales, slopes, cors,
                                           r.lkhood = r.lkhood,
                                           posteriors = posteriors,
                                           assume.constant.SE = FALSE){

  # Evaluates expected log-likelihood of the data given the parameters
  # and the posterior probability of membership of each variant in each linemodel.
  #  X, SE, scales, slopes, cors, r.lkhood
  #  are like in line.models().
  # 'posteriors' is a matrix with one row per observation and one column per linemodel.
  # NOTE: input values are not checked because this function is expected to be used
  #  only within other functions, not directly by the user.

  K = length(scales)
  M = ncol(X)
  if(nrow(posteriors) != nrow(X)) stop("Input data 'X' and 'posteriors' must have the same no. of rows.")
  if(ncol(posteriors) != K) stop("posteriors must have one column per each linemodel")

  pr.V = list() #list of prior matrices of each linemodel
  for(kk in 1:K) pr.V[[kk]] = prior.V(scale = scales[kk], slopes = slopes[kk,], cor = cors[kk])

  n = nrow(X)
  logdnorm = matrix(NA, ncol = K, nrow = n)

  if(assume.constant.SE){
    se = diag(apply(SE,2,median))
    V = se %*% r.lkhood %*% se
    for(kk in 1:K){
      logdnorm[,kk] = log.dmvnorm(t(X), mu = rep(0, M), S = V + pr.V[[kk]])}
  }else{
    for(ii in 1:n){
      V = diag(SE[ii,]) %*% r.lkhood %*% diag(SE[ii,]) #var of likelihood
      for(kk in 1:K){
        logdnorm[ii,kk] = log.dmvnorm(X[ii,], mu = rep(0, M), S = V + pr.V[[kk]])}
    }
  }
  return(sum(posteriors * logdnorm)) #expected log-likelihood
}



current.to.par <- function(current, current.ind, force.same.scales){
  # Transforms 'current' matrix to 'par' vector
  # that has only those parameters that are optimized over
  # while 'current' has the current value of all parameters.
  # Since output of this function will be used as initial
  # value in optim(), the output values must be finite.
  # current.ind is a list with 3 components: scales, slopes and cors
  # these components are either NULL or then tell the indexes of the
  # target parameters in the current matrix

  par = c()
  #log transform scales with maximum scale 1e300
  if(length(current.ind$scales)) par = pmin(log(as.numeric(current[current.ind$scales])), log(1e300))
  if(force.same.scales && length(current.ind$scales)) par = par[1]
  #turn slopes to angles
  if(length(current.ind$slopes)) par = c(par, atan(as.numeric(current[current.ind$slopes])))
  #logit transform correlations
  if(length(current.ind$cors)) {
    tmp = as.numeric(current[current.ind$cors])
    tmp = log(tmp) - log(1 - tmp)
    # logit > -40 --> cor > 4e-18
    ii = (abs(tmp) > 40)
    tmp[ii] = 40*(-1)^as.integer(tmp[ii] < 0)
    par = c(par, tmp)
  }
  return(par)
}



par.to.current <- function(par, current, current.ind, force.same.scales){
  # Transforms 'par' vector to 'current' matrix
  # that has the current values of all parameters
  # while 'par' has only those that are optimized over.
  # See also current.to.par()

  n.sc = length(current.ind$scales)
  if(n.sc && force.same.scales) n.sc = 1
  if(n.sc) current[current.ind$scales] = exp(par[1:n.sc])
  n.sl = length(current.ind$slopes)
  if(n.sl) current[current.ind$slopes] = tan(par[(n.sc+1):(n.sc+n.sl)])
  n.co = length(current.ind$cors)
  if(n.co) current[current.ind$cors] = 1/(1 + exp(-par[(n.sc+n.sl+1):(n.sc+n.sl+n.co)]))
  return(current)
}



optim.fn <- function(par, current, current.ind, force.same.scales, X, SE,
                     r.lkhood, posteriors, assume.constant.SE){

  # Defines expectation of the log-likelihood function to be optimized by
  # the EM-algorithm.
  # 'current' has current values of the parameters (1 row, M+1 cols per linemodel)
  #       Matrix has 1 row and M+1 cols per model and columns correspond to
  #       1 = scale
  #       2...M = slope
  #       M+1 = correlation
  # 'current.ind' list defines which parameters are to be optimized.
  # 'current.ind' is a list with 3 components: scales, slopes and cors
  #  these components are either NULL or then tell the indexes of the
  #  target parameters in the current matrix.

  #NOTE: scales have been log-transformed and correlations logit-transformed before
  #      calling this function. Hence back-transforming them first.
  # If 'force.same.scales = TRUE, then sets all optimized scale parameters equal.
  #  Assumes that force.same.scales = TRUE only if at least one scale parameter is optimized
  # posteriors is matrix with one row per observation and one col per linemodel
  #   and gives the posterior probability of the membership of observation in each linemodel

  current = par.to.current(par, current, current.ind, force.same.scales)
  M = ncol(current) - 1
  -line.models.expected.loglkhood(X, SE,
                         scales = current[,1], slopes = matrix(current[,2:M], ncol = M-1),
                         cors = current[,M+1], r.lkhood = r.lkhood,
                         posteriors = posteriors,
                         assume.constant.SE = assume.constant.SE)
}



#' Optimize parameters of line models
#'
#' Consider M effect sizes.
#' The linemodels are defined by three parameters:
#' scale = standard deviation of the (largest) effect,
#' slope = set of M-1 slopes defining the line around which the effects are scattered,
#' cor = non-negative pairwise correlation between each effect pair,
#' where cor = 1 means the same effect and cor = 0 means independent effects.
#'
#' A linemodel has the following properties:
#'
#' 1) Effects are scattered around line defined by vector (1,slope_1,...,slope_(M-1)).
#'    Each slope_i is the slope between effect variables i+1 and 1.
#'    Slope can be any real number or +/-Inf.
#'    If any slope is infinite, the effect of variable 1 and variables with finite slopes are 0.
#'    It is recommended to keep the slopes finite to avoid problems with interpretation.
#'
#' 2) The largest of the prior variances of the effects is scale^2.
#'    The effect that has the largest scale is determined by the slopes,
#'    and the scales of the remaining effects are determined by the slopes as well.
#'
#' 3) Distribution around the line is determined as follows.
#'    Consider a distribution scattered around the direction of the vector (1,1,...,1)
#'    with the given constant value of the pairwise correlations between the effect variables.
#'    Rotate that distribution by an orthogonal rotation that rotates the direction of the
#'    vector (1,1,...,1) to the direction of the vector (1,slope_1,...,slope_(M-1))
#'    and use the corresponding distribution, scaled so that the maximum variance
#'    among the effects is set to scale^2.
#'    NOTE: This is not same as having the given correlation around line defined by the slopes,
#'          because the shape of that distribution depends on the slopes but the
#'         definition we use here is independent of slopes up to an orthogonal transformation.
#'
#' This function optimizes over any set of scales, slopes and cors.
#' The parameters not included in optimization are kept fixed to their initial values.
#' The proportion parameters are always optimized.
#'
#' @param X matrix of estimates with M columns (one col per effect variable and one row per observation)
#' @param SE matrix of standard errors with M columns (one col per effect variable and one row per observation)
#' @param par.include list with three components:
#' 'scales' K-vector, slopes' Kx(M-1) matrix and 'cors' K-vector,
#' where K is the no. of linemodels and M is the dimension of data.
#' Each element in scales, slopes and cors is a TRUE/FALSE value
#' where TRUE indicates the parameters that are to be optimized.
#' Alternatively, can be given as Kx3 matrix, where each row is for one linemodel and
#' column 1 is for scales, column 2 for slopes and column 3 for cors.
#'  If element (ii,2) is TRUE, then all slope parameters of model 'ii' will
#'  be optimized, and if it is FALSE, then none of the slope parameters of model 'ii'
#'  will be optimized.
#' @param force.same.scales logical; If TRUE, then forces all optimized scale parameters equal
#' @param init.scales K-vector of initial scale parameters, one per each linemodel
#' @param init.slopes Kx(M-1)-matrix of initial slopes, one row per linemodel
#'  If M = 2, can be given also as a K-vector of slopes, one per each linemodel
#' @param init.cors K-vector of initial correlation parameters, one per each linemodel
#' @param model.names K-vector of names of linemodels
#' @param model.priors K-vector of prior probabilities of linemodels up to a normalization constant
#' @param r.lkhood correlation matrix of the estimators of the M effects
#'  can be given also as a vector of the upper triangular values in the row-major order
#' @param tol.loglk tolerance for convergence in log-likelihoods of consecutive iterations
#' @param tol.par tolerance for convergence in maximum of absolute values of relative differences
#' across parameters between consecutive iterations. Can be set negative to determine
#' the convergence solely by log-likelihood.
#' @param op.method Optimization method in optim() fuction. By default "BFGS",
#' can also be "Nelder-Mead". If one method crashes, try the other.
#' Does not affect univariate optimization method that is always "Brent".
#' @param assume.constant.SE If TRUE, assumes that SEs of observations in any one dimension is constant across observations.
#' SEs of different dimensions can still be different.
#' Value TRUE speeds up greatly the computation but ignores possible differences in SEs.
#' Default is FALSE.
#' @param print.steps numeric;
#' If 0, no optimization results are printed on screen.
#' If 1, prints optimization results at start and end of each optimization call (default).
#' If >1 prints status after every optimization iteration.
#' @return A list with four components
#' scales, K-vector for scale parameters of linemodels 1,...,K
#' slopes, Kx(M-1) matrix for slopes where each linemodel is specified by one row
#' cors, K-vector for correlation parameters of models 1,...,K
#' weights, giving the mixture weights of the linemodels
#' @examples
#' #2D example optimizing slope and cor of 2nd model (w/o assuming constant SEs)
#' line.models.optimize(
#'      X = linemodels.ex1[,c("beta1","beta2")],
#'      SE = linemodels.ex1[,c("se1","se2")],
#'      par.include = rbind(c(F,F,F), c(F,T,T), c(F,F,F)),
#'      force.same.scales = FALSE,
#'      init.scales = c(0.2, 0.2, 0.2),
#'      init.slopes = c(0, 0.2, 1),
#'      init.cors = c(0.995, 0.5, 0.995))
#' #2D example assuming constant SEs after scaling the effects by minor allele frequencies
#' # and optimizing scales that are assumed same across the linemodels
#' # and slope and cor of the 2nd model as above.
#' sc = sqrt(2 * linemodels.ex1$maf * (1 - linemodels.ex1$maf))
#' line.models.optimize(
#'      X = linemodels.ex1[,c("beta1","beta2")]*sc,
#'      SE = linemodels.ex1[,c("se1","se2")]*sc,
#'      par.include = rbind(c(T,F,F), c(T,T,T), c(T,F,F)),
#'      force.same.scales = TRUE,
#'      init.scales = c(0.2, 0.2, 0.2),
#'      init.slopes = c(0, 0.2, 1),
#'      init.cors = c(0.995, 0.5, 0.995),
#'      assume.constant.SE = TRUE)
#' #3D example optimizing all slopes (you can also scale X and SEs and assume constant SE as above to speed up)
#' line.models.optimize(
#'      X = linemodels.ex1[,c("beta1","beta2","beta3")],
#'      SE = linemodels.ex1[,c("se1","se2","se3")],
#'      par.include = list(scales = c(F,F,F), slopes = matrix(TRUE, ncol = 2, nrow = 3), cors = c(F,F,F)),
#'      force.same.scales = FALSE,
#'      init.scales = c(0.2, 0.2, 0.2),
#'      init.slopes = rbind(c(0, 0), c(0.2, 0), c(0.5, 1)),
#'      init.cors = c(0.995, 0.995, 0.995),
#'      assume.constant.SE = FALSE)
#' @export
line.models.optimize <- function(X, SE,
                                 par.include = NULL,
                                 force.same.scales = FALSE,
                                 init.scales, init.slopes, init.cors,
                                 model.names = NULL,
                                 model.priors = rep(1,length(init.scales)),
                                 r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
                                 tol.loglk =  1e-3,
                                 tol.par = 0,
                                 op.method = "BFGS",
                                 assume.constant.SE = FALSE,
                                 print.steps = 1){

  # EM algorithm to optimize chosen parameters of the line models.
  #
  # X, SE, model.names, r.lkhood are as in line.models().
  #
  # par.include is a list that defines which of the parameters will be optimized
  #  TRUE = optimized, FALSE = not optimized.
  # par.include can have three components
  # scales, K-vector for scale parameters of linemodels 1,...,K
  # slopes, Kx(M-1) matrix for slopes where each linemodel is specified by one row
  #         if slopes is K-vector, then assumes M = 2
  # cors, K-vector for correlation parameters of linemodels 1,...,K
  #
  # init.scales, init.slopes, init.cors, model.priors, give initial values
  # for scales, slopes, cors and proportion parameters, respectively.
  # Note that those scales/slopes/cors that are not being optimized over
  #  will remain fixed to their initial values. Proportion parameters
  #  are always optimized by the algorithm.
  #
  # tol.loglk, defines tolerance for maximum difference in log-likelihood
  #  between consecutive iterations that is still considered to show convergence.
  # tol.par, defines tolerance for absolute value of maximum relative difference
  # in optimized parameters between consecutive iterations that is still
  # considered to show convergence.
  # Convergence happens when either of the tolerances above are achieved.
  # NOTE: tol.par can be set negative to determine convergence by only loglkhood.
  # op.method, method used for optimization as available in optim( ) function
  #            default is "BFGS", another option is "Nelder-Mead"
  #            If only one parameter is being optimized, then uses always "Brent"
  #            no matter what is value of 'op.method'

  # OUTPUT:
  # list with 4 components:
  # scales, K-vector for scale parameters of linemodels 1,...,K
  # slopes, Kx(M-1) matrix for slopes where each linemodel is specified by one row
  # cors, K-vector for correlation parameters of linemodels 1,...,K
  # weights, giving the mixture weights of the linemodels

  chk = check.input("line.models.optimize:", X, SE,
                    init.scales, init.slopes, init.cors,
                    model.names, model.priors, r.lkhood)

  K = chk$K #number of models
  M = chk$M #number of dimensions
  init.slopes = chk$slopes
  r.lkhood = chk$r.lkhood
  model.names = chk$model.names

  if(tol.loglk <= 0) stop("Tolerance 'tol.loglk' should be positive.")
  if(is.matrix(par.include)){
    if(nrow(par.include) != K || ncol(par.include) != 3){
      stop("If par.include is a matrix, it should have exactly 3 columns and one row per each model.")
    }
    par.include = list(scales = par.include[,1],
                       slopes = matrix(rep(par.include[,2], M-1), nrow = K),
                       cors = par.include[,3])
  }
  if(! ("scales" %in% names(par.include))) par.include$scales = rep(FALSE, K)
  if(! ("cors" %in% names(par.include))) par.include$cors = rep(FALSE, K)
  if(! ("slopes" %in% names(par.include))) par.include$slopes = matrix(FALSE, nrow = K, ncol = M-1)

  if(!(print.steps %in% c(0,1,2))) stop("print.steps should be one of values 0, 1 or 2.")

  current = cbind(init.scales, init.slopes, init.cors) #current values of parameters
  current.ind = list(scales = which(par.include$scales),
                     slopes = K + which(par.include$slopes),
                     cors = K*M + which(par.include$cors)) #indexes w.r.t current matrix of parameters to optimize

  if(force.same.scales){
    if(length(current.ind$scales) > 0) {
      current[current.ind$scales] = mean(current[current.ind$scales])} #initialize
    else{
      force.same.scales = FALSE} #no scales to update
  }
  n.optim = length(unlist(current.ind)) - as.numeric(force.same.scales)*(length(current.ind$scales)-1) #how many parameters to optimize
  w = chk$model.priors
  log.lkhood = line.models.loglkhood(X, SE,
                                     scales = current[,1], slopes = matrix(current[,2:M], ncol = M-1),
                                     cors = current[,M+1],
                                     model.priors = w,
                                     model.names = model.names, r.lkhood = r.lkhood,
                                     assume.constant.SE = assume.constant.SE,
                                     return.posteriors = FALSE)

  par.names = matrix(apply(matrix(0:(K*(M-1)-1),ncol = M - 1, nrow = K, byrow = TRUE), 2,
                    function(x){paste0("slope", 1+floor(x/(M-1)), "_", 1+(x%%(M-1)))}), ncol = M - 1)
  par.names = cbind(paste0("scale",1:K), par.names, paste0("cor",1:K))
  if(print.steps > 0){
    cat(paste0("\n\nInitial values\n"))
    cat(paste("scales:",paste(signif(current[,1],4),collapse=", ")),"\n")
    cat(paste("slopes:\n"))
    for(ii in 1:K){
      cat(paste0("\t", paste(signif(current[ii, 2:M],4), collapse = ", "),"\n"))}
    cat(paste("cors:",paste(signif(current[, M + 1],4),collapse=", ")),"\n")
    cat(paste("proportions:", paste(signif(w,3),collapse =", ")),"\n")
    cat(paste("Initial log-lkhood:",signif(log.lkhood,8),"\n"))
    if(n.optim > 0){
      cat(paste("Optimizing w.r.t:",
                paste(as.vector(c(
                  par.names[current.ind$scales],
                  par.names[current.ind$slopes],
                  par.names[current.ind$cors])), collapse = " ")))}
    else{
      cat(paste("All parameters are fixed. Optimizes only the proportions."))}
  }
  # For 1-dimensional optimization uses "Brent",
  # but Brent requires finite range of optimization so needs to set the range
  # depending on which type of the parameter is optimized.

  op.method = ifelse(n.optim > 1, op.method, "Brent")
  lower = -Inf # Leaves lower and upper as Infs if Brent is not used
  upper = Inf
  if(n.optim == 1){ # specify the range for the single parameter
    if(sum(par.include$scales)){ # log-transformed scale parameter (1e-6, 1e+6)
      lw.orig = 1e-6; lower = log(lw.orig)
      up.orig = 1e+6; upper = log(up.orig)
    }
    if(sum(par.include$slopes)) { # atan transformed slope parameter (-Inf, Inf)
      lw.orig = -Inf; lower = atan(lw.orig)
      up.orig = Inf; upper = atan(up.orig)
    }
    if(sum(par.include$cors)) { # log-odds transformed correlation parameter (1e-6, 1-1e-6)
      lw.orig = 1e-6; lower = log(lw.orig / (1 - lw.orig))
      up.orig = 1 - 1e-6; upper = log(up.orig / (1 - up.orig))
    }
    if(print.steps > 0){
      cat(paste0(", over the range: (",
                 signif(lw.orig,4),",",signif(up.orig,4),")"))}
  }
  if(print.steps > 0){
    if(force.same.scales) cat(paste("\nForcing all scales to be equal"))
    cat("\nConvergence criteria:\n")
    cat(paste(" Relative diff in parameters <=",tol.par," or\n"))
    cat(paste(" Difference in log-likelihood <",tol.loglk))
    cat("\n\n")
  }

  iter = 0
  converged = FALSE

  while(!converged){
    iter = iter + 1
    prev.log.lkhood = log.lkhood
    prev.vals = current

    #Maximizes over proportions:
    pr = line.models.loglkhood(X, SE,
                               scales = current[,1],
                               slopes = matrix(current[,2:M], ncol = M-1, nrow = K),
                               cors = current[,M+1],
                               model.names = model.names,
                               model.priors = w,
                               r.lkhood = r.lkhood,
                               assume.constant.SE = assume.constant.SE,
                               return.posteriors = TRUE)$membership.prob
    new.w = apply(pr, 2, mean)

    #Maximizes over other parameters than proportions:
    if(n.optim > 0){
      par = current.to.par(current, current.ind, force.same.scales)
      opt.out = optim(par, optim.fn, current = current, current.ind = current.ind,
                      force.same.scales = force.same.scales,
                      X = X, SE = SE, posteriors = pr,
                      r.lkhood = r.lkhood, method = op.method,
                      assume.constant.SE = assume.constant.SE,
                      lower = lower, upper = upper)
      new.par = par.to.current(opt.out$par, current, current.ind, force.same.scales)}
    else{ #no parameters to optimize
      new.par = current}

    if(force.same.scales) new.par[current.ind$scales] = new.par[current.ind$scales[1]]

    new.log.lkhood = line.models.loglkhood(X, SE,
                          scales = new.par[,1],
                          slopes = matrix(new.par[,2:M], ncol = M-1, nrow = K),
                          cors = new.par[,M+1],
                          model.names = model.names,
                          model.priors = new.w,
                          r.lkhood = r.lkhood,
                          assume.constant.SE = assume.constant.SE,
                          return.posteriors = FALSE)

    if(new.log.lkhood > log.lkhood){
      log.lkhood = new.log.lkhood
      w = new.w
      current = new.par}
    else{ #will stop because both convergence criteria below are fulfilled
      if(print.steps > 1){
        cat(paste("iter:",iter,"; Failed to increase log-likelihood further.",
                  "(previous:",signif(log.lkhood,8),"new:",signif(new.log.lkhood,8),")\n"))}}

    tmp = as.vector(abs( (current - prev.vals)/prev.vals))
    converged.par = all(tmp[!is.nan(tmp)] <= tol.par)
    converged.loglk = ((log.lkhood - prev.log.lkhood) < tol.loglk)
    converged = converged.par | converged.loglk
    if((print.steps > 1) | (converged & (print.steps > 0))){
      cat(paste("iter:",iter,"; log-lkhood:",signif(log.lkhood,8)),"\n")
      cat(paste("Relative diffs in optimized parameters:",
                paste(signif(as.vector(abs(((current - prev.vals)/prev.vals)[as.vector(unlist(current.ind))])), 3),
                      collapse =", ")),"\n")
      cat(paste("proportions:", paste(signif(w,3),collapse =", ")),"\n")
      cat(paste("scales:",paste(signif(current[,1],4),collapse=", ")),"\n")
      cat(paste("slopes:\n"))
      for(ii in 1:K){
        cat(paste0("\t", paste(signif(current[ii,2:M],4), collapse = ", "),"\n"))}
      cat(paste("cors:",paste(signif(current[,M+1],4),collapse=", ")),"\n")
    }
  }

  if(print.steps > 0){
    if(converged.par) cat("Parameter values converged.\n")
    if(converged.loglk) cat("Log-likelihood value converged.\n")}

  return(list(scales = as.numeric(current[,1]),
              slopes = matrix(current[,2:M], ncol = M-1),
              cors = as.numeric(current[,M+1]),
              weights = w))
}



#' Random sample from a line model
#'
#' See help of 'line.models( )' for details about specifying the models.
#'
#' @param n the number of samples
#' @param scale standard deviation of the largest effect
#' @param slope vector of slopes of the linemodel
#' @param cor correlation parameters of the linemodel
#' @param scale.weights vector of weights of each component when the effect size
#' distribution is a mixture of Gaussians
#' @return matrix with one row per sample and one column per effect on each outcome
#' @examples
#' sample.line.model(10, scale = 0.1, slope = 1, cor = 0.99)
#' sample.line.model(10, scale = 0.1, slope = 1, cor = 0.99, scale.weights = c(1,0,1))
#' sample.line.model(10, scale = 0.1, slope = c(1,-1), cor = 0.99)
#' @export
sample.line.model <- function(n = 1, scale, slope, cor, scale.weights = c(1)){

  # Outputs 'n' samples from a given linemodel defined by
  # scale, slope and cor and scale.weights.
  # See line.models() for definition of these parameters.
  # 'slope' must be a vector corresponding to the slopes of one line.
  # 'scale.weights' must be a vector and cannot be a matrix as in line.models( ).

  if(cor > 1) stop("cor > 1")
  if(cor < 0) stop("cor < 0")
  if(scale < 0) stop("scale < 0")
  if(!is.vector(scale.weights)) stop("scale.weights must be vector.")
  if(any(scale.weights < 0)) stop("scale.weights can have only positive values.")
  J = length(scale.weights)
  M = length(slope) + 1
  V = prior.V(scale = scale, slopes = slope, cor = cor)

  #Choose components for each sample
  ind = sample(1:J, size = n, prob = scale.weights, replace = T)

  A = t(chol(V)) #This may fail if determinant is ~0
  x = t(A %*% matrix(rnorm(M*n), nrow = M) / rep(J + 1 - ind, each = M))
  cat(paste0("Sampling effects with scale=", scale,
             " slope=(", paste(slope,collapse = ","),")"," cor=", cor),"\n")
  cat(paste0("   scale.weights=", paste(scale.weights, collapse = ","),"."),"\n")
  cat(paste0("Theoretical SD for the largest effect:",
             signif(sqrt(sum(scale^2*scale.weights/sum(scale.weights)/rev(1:J)^2)),5)))
  cat(paste0("; observed value:", signif(max(apply(x,2,sd)),5)),".\n")
  return(x)
}



#' Simulate a data set mimicking the observed multi-dimensional data given a set of linemodels
#'
#' See help of 'line.models( )' for details about specifying the models.
#'
#' Generates M-dimensional estimates for each observed data point as follows.
#' Pick up an underlying linemodel based on the given probabilities (default uniform probabilities).
#' Choose the primary coordinate to be the one with the largest expected magnitude based on the slopes.
#' For the primary coordinate, choose the true effect size from the univariate posterior
#'  distribution given the observed univariate estimate, its SE and the prior N(0, scale),
#'  where 'scale' is the scale parameter of the corresponding linemodel.
#' Project the point on the line based on the primary coordinate.
#' Rotate the point to be on the diagonal line.
#' Pick uniformly at random one coordinate to be kept fixed,
#'  sample the true effect of the remaining M-1 coordinates from the linemodel,
#'  conditional on the current value of the fixed coordinate.
#' Rotate the perturbed point back to correspond to the original slope.
#' Finally, sample the observed effect estimates from a Gaussian around the true values
#' by accounting for SEs and r.lkhood.
#'
#' @param X matrix of estimates with M columns (one col per effect variable and one row per observation)
#' @param SE matrix of standard errors with M columns (one col per effect variable and one row per observation)
#' @param scales K-vector of standard deviations of the largest effect, one per each linemodel
#' @param slopes Kx(M-1)-matrix of slopes, one row per linemodel
#'  If M = 2, can be given also as a K-vector of slopes, one per each linemodel
#' @param cors K-vector of correlation parameters, one per each linemodel
#' @param linemodel.prob matrix of probabilities of each observation (row) in each model (col) up to a normalization constant
#' @param r.lkhood correlation matrix of the estimators of the M effects
#'  can be given also as a vector of the upper triangular values in the row-major order
#' @return matrix of simulated observations with one column per effect variable and one row per observation in X
#' @examples
#' #2D example
#' simulate.linemodels.for.observations(
#'    X = linemodels.ex1[,c("beta1","beta2")], SE = linemodels.ex1[,c("se1","se2")],
#'    scales = c(0.25, 0.25),
#'    slopes = c(0, 1),
#'    cors = c(0.995, 0.995))
#' #3D example
#' simulate.linemodels.for.observations(
#'    X = linemodels.ex1[,c("beta1","beta2","beta3")], SE = linemodels.ex1[,c("se1","se2","se3")],
#'    scales = c(0.25, 0.25),
#'    slopes = rbind(c(0, 0), c(1,1)),
#'    cors = c(0.995, 0.95))
#' @export
simulate.linemodels.for.observations <- function(
    X, SE, linemodel.prob = NULL,
    scales, slopes, cors,
    r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2)){

  n = nrow(X)
  chk = check.input("simulate.linemodels.for.observations:",
                    X = X, SE =SE,
                    scales = scales, slopes = slopes, cors = cors,
                    r.lkhood = r.lkhood)
  K = chk$K #number of models
  M = chk$M #number of dimensions
  slopes = chk$slopes
  r.lkhood = chk$r.lkhood

  if(is.null(linemodel.prob)){ #set equal probabilities for the models if not otherwise specified
    linemodel.prob = matrix(1, ncol = K, nrow = n)}
  if(is.vector(linemodel.prob)){
    if(length(linemodel.prob) != K) stop("simulate.linemodels.for.observations: If linemodel.prob is vector, its length must be no. of linemodels.")
    linemodel.prob = matrix(linemodels.result, ncol = K, byrow = TRUE)
  }
  if(any(linemodel.prob < 0)) stop("All values in 'linemodel.prob' must be non-negative.")

  Y = matrix(NA, ncol = M, nrow = n) #simulated effect sizes
  for(ii in 1:n){
    k = sample(1:K, size = 1, prob = linemodel.prob[ii,])
    primary = which.max(abs(slopes[k,]))
    primary = ifelse(abs(slopes[k, primary]) > 1, 1 + primary, 1)
    v = 1/( 1/scales[k]^2 + 1/SE[ii, primary]^2) #posterior variance
    mu = v*X[ii, primary] / SE[ii, primary]^2 #posterior mean
    x = rnorm(1, mu, sqrt(v)) #sample a true effect for the primary coordinate
    #set the other coordinates to be exactly on the line
    if(primary == 1){x = c(1, as.numeric(slopes[k,]))*x}
    else{x = c(1,slopes[k,])/slopes[k,primary-1]*x}

    Q = rotate.diagonal.to.line(slopes[k,])
    # Rotate the point x to the diagonal line by t(Q)
    xx = t(Q) %*% x

    # Now xx[1] == xx[ii] for all ii and we will randomly choose one coordinate to be fixed
    # and sample other coordinates from the linemodel's distribution accounting for cor and scale
    # Conditional distribution for M-1 dimensions is
    #  mu = cor*xx[1]
    #  Sigma = scale^2*(1-cor)*B
    B = matrix(cors[k], M-1, M-1)
    diag(B) = 1 + diag(B)

    A = as.matrix(t(chol(B)))
    xx[sample(1:M, size = M-1)] = cors[k]*xx[1] + scales[k]*sqrt(1 - cors[k])*A %*% rnorm(M-1)
    # NOTE: the index that is not updated will keep its current value that was conditioned on here.
    # Rotate the point back to correspond to the original slope
    x = Q %*% xx

    # Sample the observation error around the point of true effects
    # using the observed SEs and r.lkhood
    V = diag(SE[ii,]) %*% r.lkhood %*% diag(SE[ii,])
    A = t(chol(V))
    Y[ii,] = x + A %*% rnorm(M)
  }
  return(Y)
}



#' Compute observed log of likelihood ratio (LR) and its null distribution for the given scenario.
#'
#' See help of 'line.models( )' for details about specifying the models.
#'
#' Compares two models (ALTernative and NULL) against each other by maximized likelihood ratio.
#' Both models are specified by input parameters that are passed to line.models.optimize().
#' The optimized log-likelihoods are computed and logLR of ALT vs. NULL is recorded.
#' Then a null distribution for logLR is estimated using user-specified number of simulated data sets.
#' Each data set is simulated to mimic the observed data set (w.r.t sample size, estimates, SEs)
#' and follows the optimized NULL model w.r.t the linemodel parameters and mixture weights of the linemodels.
#' Then both models (ALT and NULL) are optimized using the simulated data and logLRs are recorded.
#' Returns both the observed logLR and the simulated null distribution of logLR,
#' which can be used for deriving an empirical P-value.
#'
#' @param X matrix of estimates with M columns (one col per effect variable and one row per observation)
#' @param SE matrix of standard errors with M columns (one col per effect variable and one row per observation)
#' @param n.sims number of simulated data sets for estimating distribution of logLR
#' @param par.include.null list with three components for the NULL model:
#' 'scales' K-vector, slopes' Kx(M-1) matrix and 'cors' K-vector,
#' where K is the no. of linemodels and M is the dimension of data.
#' Each element in scales, slopes and cors is a TRUE/FALSE value
#' where TRUE indicates the parameters that are to be optimized.
#' Alternatively, can be given as Kx3 matrix, where each row is for one linemodel and
#' column 1 is for scales, column 2 for slopes and column 3 for cors.
#'  If element (ii,2) is TRUE, then all slope parameters of model 'ii' will
#'  be optimized, and if it is FALSE, then none of the slope parameters of model 'ii'
#'  will be optimized.
#' @param force.same.scales.null logical; If TRUE, then forces all optimized scale parameters equal in the NULL model
#' @param init.scales.null K-vector of initial standard deviations of the largest effect, one per each linemodel,
#' under the NULL
#' @param init.slopes.null Kx(M-1)-matrix of initial slopes, one row per linemodel, under the NULL
#'  If M = 2, can be given also as a K-vector of slopes, one per each linemodel
#' @param init.cors.null K-vector of initial correlation parameters, one per each linemodel, under the NULL
#' @param model.priors.null K-vector of initial prior probabilities of linemodels up to a normalization constant under the NULL
#' @param par.include.alt analogous to par.include.null but for the ALT model
#' @param force.same.scales.alt analogous to force.same.scales.null but for the ALT model
#' @param init.scales.alt analogous to init.scales.null but for the ALT model
#' under the ALT
#' @param init.slopes.alt analogous to init.slopes.null but for the ALT model
#'  If M = 2, can be given also as a K-vector of slopes, one per each model
#' @param init.cors.alt analogous to init.cors.null but for the ALT model
#' @param model.priors.alt analogous to model.priors.null but for the ALT model
#' @param r.lkhood correlation matrix of the estimators of the M effects
#'  can be given also as a vector of the upper triangular values in the row-major order
#' @param tol.loglk tolerance for convergence in log-likelihoods of consecutive iterations
#' @param tol.par tolerance for convergence in maximum of absolute values of relative differences
#' across parameters between consecutive iterations. Can be set negative to determine
#' the convergence solely by log-likelihood.
#' @param op.method Optimization method in optim() fuction. By default "BFGS",
#' can also be "Nelder-Mead". If one method crashes, try the other.
#' Does not affect univariate optimization method that is always "Brent".
#' @param assume.constant.SE If TRUE, assumes that SEs of observations in any one dimension is constant across observations.
#' SEs of different dimensions can still be different.
#' Value TRUE speeds up greatly the computation but ignores possible differences in SEs.
#' Default is FALSE.
#' @param print.steps vector with 2 elements;
#' If print.steps[1] is TRUE, then prints iteration number.
#' Value of print.steps[2] is passed to optimization and it must be 0, 1 or 2.
#' If 0, no optimization results are printed on screen.
#' If 1, prints optimization results at start and end of each optimization call (default).
#' If >1 prints status after every optimization iteration.
#' @return list with two components. 'obs.loglr' is the log-likelihood ratio for the observed data and
#' 'null.loglr' is a vector of 'n.sims' simulated loglr values under the NULL model.
#' @examples
#' simulate.loglr(
#'    X = linemodels.ex1[,c("beta1","beta2")], SE = linemodels.ex1[,c("se1","se2")],
#'    n.sims = 2,
#'    par.include.null = rbind(c(F,F,F), c(F,F,F)),
#'    init.scales.null = c(0.25, 0.25),
#'    init.slopes.null = c(0, 1),
#'    init.cors.null = c(0.995, 0.995),
#'    par.include.alt = rbind(c(F,F,F), c(F,F,F), c(F,T,F)),
#'    init.scales.alt = c(0.25, 0.25, 0.25),
#'    init.slopes.alt = c(0, 1, 0.5),
#'    init.cors.alt = c(0.995, 0.995, 0.995))
#' # Same as above but scaling effects and SEs to allow SEs to be treated as constant to speed up
#' # Also estimates the scales of all linemodels and forces same scales across linemodels.
#' sc = sqrt(2 * linemodels.ex1$maf * (1 - linemodels.ex1$maf))
#' simulate.loglr(
#'    X = linemodels.ex1[,c("beta1","beta2")]*sc, SE = linemodels.ex1[,c("se1","se2")]*sc,
#'    n.sims = 2,
#'    par.include.null = rbind(c(T,F,F), c(T,F,F)),
#'    init.scales.null = c(0.15, 0.15),
#'    init.slopes.null = c(0, 1),
#'    init.cors.null = c(0.995, 0.995),
#'    force.same.scales.null = TRUE,
#'    par.include.alt = rbind(c(T,F,F), c(T,F,F), c(T,T,F)),
#'    init.scales.alt = c(0.15, 0.15, 0.15),
#'    init.slopes.alt = c(0, 1, 0.5),
#'    init.cors.alt = c(0.995, 0.995, 0.995),
#'    force.same.scales.alt = TRUE,
#'    assume.constant.SE = TRUE)
#' @export
simulate.loglr <- function(
    X, SE, n.sims = 100,
    par.include.null = matrix(TRUE,
                              nrow = length(init.scales.null),
                              ncol = 3),
    force.same.scales.null = FALSE,
    init.scales.null, init.slopes.null, init.cors.null,
    model.priors.null = rep(1, length(init.scales.null)),
    par.include.alt = matrix(TRUE,
                             nrow = length(init.scales.alt),
                             ncol = 3),
    force.same.scales.alt = FALSE,
    init.scales.alt, init.slopes.alt, init.cors.alt,
    model.priors.alt = rep(1,length(init.scales.alt)),
    r.lkhood = rep(0, ncol(X)*(ncol(X)-1)/2),
    tol.loglk = 1e-3,
    tol.par = 0,
    op.method = "BFGS",
    assume.constant.SE = FALSE,
    print.steps = c(1,0)){

  #Check and initialize NULL model parameters
  chk = check.input("simulate.loglr null:", X, SE,
                    scales = init.scales.null, slopes = init.slopes.null, cors = init.cors.null,
                    model.priors = model.priors.null, r.lkhood = r.lkhood)
  M = chk$M #number of dimensions
  r.lkhood = chk$r.lkhood
  K.null = chk$K #number of models under the nULL
  init.slopes.null = chk$slopes
  model.priors.null = chk$model.priors

  if(is.matrix(par.include.null)){
    if(nrow(par.include.null) != K.null || ncol(par.include.null) != 3){
      stop("If par.include.null is a matrix, it should have exactly 3 columns and one row per each model.")
    }
    par.include.null = list(scales = par.include.null[,1],
                            slopes = matrix(rep(par.include.null[,2], M-1), nrow = K.null),
                            cors = par.include.null[,3])
  }
  if(! ("scales" %in% names(par.include.null))) par.include.null$scales = rep(FALSE, K.null)
  if(! ("cors" %in% names(par.include.null))) par.include.null$cors = rep(FALSE, K.null)
  if(! ("slopes" %in% names(par.include.null))) par.include.null$slopes = matrix(FALSE, nrow = K.null, ncol = M-1)

  #Check and initialize ALT model parameters
  chk = check.input("simulate.loglr alt:", X, SE,
                    scales = init.scales.alt, slopes = init.slopes.alt, cors = init.cors.alt,
                    model.priors = model.priors.alt, r.lkhood = r.lkhood)
  K.alt = chk$K
  init.slopes.alt = chk$slopes
  model.priors.alt = chk$model.priors

  if(is.matrix(par.include.alt)){
    if(nrow(par.include.alt) != K.alt || ncol(par.include.alt) != 3){
      stop("If par.include.alt is a matrix, it should have exactly 3 columns and one row per each model.")
    }
    par.include.alt = list(scales = par.include.alt[,1],
                           slopes = matrix(rep(par.include.alt[,2], M-1), nrow = K.alt),
                           cors = par.include.alt[,3])
  }
  if(! ("scales" %in% names(par.include.alt))) par.include.alt$scales = rep(FALSE, K.alt)
  if(! ("cors" %in% names(par.include.alt))) par.include.alt$cors = rep(FALSE, K.alt)
  if(! ("slopes" %in% names(par.include.alt))) par.include.alt$slopes = matrix(FALSE, nrow = K.alt, ncol = M-1)

  if(tol.loglk <= 0) stop("Tolerance 'tol.loglk' should be positive.")
  if(!(print.steps[2] %in% c(0,1,2))) stop("print.steps[2] should be one of values 0, 1 or 2.")

  loglk = matrix(NA, ncol = 2, nrow = n.sims + 1)

  for(ii in 1:(n.sims + 1)){
    if(ii == 1) {XX = X} #First, compute the log-likelihoods for the observed data.
    #Later, for simulations ii > 1, use parameters learned from the observed data at iteration ii == 1.
    else {
      XX = simulate.linemodels.for.observations(
        X, SE, linemodel.prob = post.null,
        scales = scales.null,
        slopes = slopes.null,
        cors = cors.null,
        r.lkhood = r.lkhood)}

    if(print.steps[1]) cat(paste("\nNULL model optimization for data set",ii-1,"\n"))
    lm.opt.null = line.models.optimize(
      XX, SE, par.include = par.include.null,
      force.same.scales = force.same.scales.null,
      init.scales = init.scales.null,
      init.slopes = init.slopes.null,
      init.cors = init.cors.null,
      model.priors = model.priors.null,
      r.lkhood = r.lkhood, tol.loglk = tol.loglk, tol.par = tol.par,
      op.method = op.method,
      assume.constant.SE = assume.constant.SE,
      print.steps = print.steps[2])

    loglk.null = line.models.loglkhood(
      XX, SE,
      scales = lm.opt.null$scales,
      slopes = lm.opt.null$slopes,
      cors = lm.opt.null$cors,
      model.priors = lm.opt.null$weights,
      r.lkhood = r.lkhood,
      return.posteriors = (ii == 1))

    if(ii == 1){ #save values for simulations from the null model to be used in later iterations ii > 1
      post.null = loglk.null$membership.prob
      loglk.null = loglk.null$loglkhood
      scales.null = lm.opt.null$scales
      slopes.null = lm.opt.null$slopes
      cors.null = lm.opt.null$cors}

    if(print.steps[1]) cat(paste("ALT model optimization for data set",ii-1,"\n"))
    lm.opt.alt = line.models.optimize(
      XX, SE, par.include = par.include.alt,
      force.same.scales = force.same.scales.alt,
      init.scales = init.scales.alt,
      init.slopes = init.slopes.alt,
      init.cors = init.cors.alt,
      model.priors = model.priors.alt,
      r.lkhood = r.lkhood, tol.loglk = tol.loglk, tol.par = tol.par,
      op.method = op.method,
      assume.constant.SE = assume.constant.SE,
      print.steps = print.steps[2])

    loglk.alt = line.models.loglkhood(
      XX, SE,
      scales = lm.opt.alt$scales,
      slopes = lm.opt.alt$slopes,
      cors = lm.opt.alt$cors,
      model.priors = lm.opt.alt$weights,
      r.lkhood = r.lkhood,
      return.posteriors = FALSE)

    loglk[ii,] = c(loglk.null, loglk.alt)
  }
  #Returns list where elements are
  # 1: observed logLR between alternative and null.
  # 2: vector of simulated logLR values under the null model.
  return(
    list( obs.loglr = loglk[1,2] - loglk[1,1],
          null.loglr = as.numeric(loglk[2:(n.sims + 1),2] - loglk[2:(n.sims + 1),1])))
}
