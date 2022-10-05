# linemodels (0.1.0): An R package to cluster 2-dimensional effects into
#  groups defined by linear relationships.
# DATE: 4-Oct-2022
# Copyright (C) 2022 Matti Pirinen
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

# Consider two effects:  x = effect1 and y = effect2.
# The line models are defined by three parameters
# scale = standard deviation of the (larger) effect
# slope = slope of line around which the effects are scattered
# cor = non-negative correlation between effects where cor = 1 means the same effect size
# A line model has the following properties:
#
# 1) effects are scattered around line y = slope*x.
#    'slope' can be any real number or Inf (in which case effect x is zero)
#
# 2) the larger of the prior variances of effects is scale^2
#    That is, if |slope| <= 1, then Var(x) = scale^2
#             if |slope| > 1, then Var(y) = scale^2
#
# 3) Distribution around the line is determined as follows:
#    Consider a distribution scattered around line y = x with correlation 'cor'.
#    Rotate that distribution by an orthogonal rotation defined by angle
#    theta = atan(slope) - pi/4
#    and use the corresponding distribution, scaled so that the maximum variance
#    of the two effects is scale^2.
#    NOTE: This is not same as "correlation 'cor' around line y = slope*x",
#          because the shape of that distribution depends on slope, but
#          we do not want such dependency in order to be consistent across slopes.
#
# Examples:
# Set scale = 0 to get the null model. (Values of slope and cor do not matter in this case.)
# Use scale > 0, slope = 1 and cor = 0 to get independent effects model.
# Use scale > 0, slope = 1 and cor = 1 to get fixed effects model.
# Use scale > 0, slope = 1 and 0 < cor < 1  to get correlated effects model.
#               note that values need to be quite near cor = 1 to get similar values
#               between the two effects with high probability (e.g. cor = 0.99 or 0.999).
#
# To choose a reasonable scale, note that 95% of the effects will be < 2*scale.

#
# There are 6 functions for user to call:
#

# visualize.line.models(scales, slopes, cors,
#                      model.names = NULL, model.cols = NULL,
#                      legend.position = "bottomright",
#                      xlim = NULL, ylim = NULL,
#                      xlab = "EFFECT1", ylab = "EFFECT2",
#                      cex.lab = 1, cex.axis = 1,
#                      line.lty = 1, region.lty = 2,
#                      line.lwd = 1, region.lwd = 1,
#                      emphasize.axes = TRUE)
# -- To visualize the lines and 95% highest probability regions of the models.

#visualize.scales(scales, scale.weights = c(1),
#                 model.names = NULL, model.cols = NULL,
#                 legend.position = "topleft")
# -- To visulaize the univariate distributions of the effect sizes.

# line.models(X, SE,
#             scales, slopes, cors,
#             model.names = NULL,
#             model.priors = rep(1/length(slopes), length(slopes)),
#             r.lkhood = 0, scale.weights = c(1))
# -- To evaluate the model probabilities for each opservation separately.

# line.models.with.proportions(X, SE,
#                             scales, slopes, cors,
#                             model.names = NULL,
#                             r.lkhood = 0,
#                             n.iter = 200, n.burnin = 20)
# -- To evaluate the model probabilities together with the proportions of
#    observations coming from each model. Uses Gibbs sampler for a joint
#    estimation of membership probabilities across all observations.

# line.models.optimize(X, SE,
#           par.include = matrix(TRUE, nrow = length(init.slopes), ncol = 3),
#           init.scales, init.slopes, init.cors,
#           model.priors = rep(1,length(init.slopes)),
#           model.names = NULL,
#           r.lkhood = 0, tol =  1e-3)
# -- To use EM-algorithm to optimize chosen parameters of the line models.

# sample.line.model(n = 1, scale, slope, cor, scale.weights = c(1))
# -- To generate samples from a line model.

################################################################################
############## FUNCTIONS #######################################################
################################################################################


prior.V <- function(scale = 0, slope = 1, cor = 0){

 # Returns prior covariance matrix for a line model defined by
 #  scale, slope and cor.
 # INPUT
 # scale > 0
 # slope in (-Inf,Inf]
 # cor >= 0, NOTE: negatively correlated effects are modeled by a negative slope

  if(cor > 1) stop("cor > 1")
  if(cor < 0) stop("cor < 0")
  if(scale < 0) stop("scale < 0")

  theta = atan(slope) - pi/4
  # R is rotation of angle theta
  R = matrix(c(cos(theta), -sin(theta),
               sin(theta), cos(theta)),
             ncol = 2, byrow = T)
  #Rotate diagonal case of correlation 'cor' by angle theta to S:
  S = R %*% matrix(c(1, cor, cor, 1), ncol = 2) %*% t(R)
  S = scale^2 * S / max(as.vector(S))
  return(S)
}



rdirichlet <- function(alpha){

  # Random sampling from Dirichlet distribution with parameter vector alpha.

  g = rgamma(length(alpha), shape = alpha, scale = 1)
  return(g/sum(g))
}



log.dmvnorm <- function(x, mu = rep(0, length(x)), S = diag(1, length(x)) ){

  # Returns log of density of MV-Normal(mean = mu, var = S) at x

  K = length(mu)
  stopifnot(all(dim(S) == K))
  stopifnot(length(x) == K)
  chol.S = chol(S) # Cholesky decomposition
  log.det = 2*sum(log(diag(chol.S))) # log of det(S)
  inv.chol.S = solve(t(chol.S)) # inverse of cholesky^T
  return(-K/2*log(2*pi) - 0.5*(log.det + crossprod(inv.chol.S %*% as.numeric(x-mu))))
}


#' Visualize line models
#'
#' Plots lines and 95\% highest probability regions of the models
#' specified by scales, slopes and correlations.
#' See help of 'line.models( )' for details about specifying the models.
#'
#' @param scales vector of standard deviations of larger effect of each model
#' @param slopes vector of slopes of each model
#' @param cors vector of correlation parameters of each model
#' @param model.names vector of names of each model
#' @param model.cols vector of colors for each model
#' @param legend.position position of legend, value NULL removes legend
#' @param xlim,ylim,xlab,ylab,cex.lab,cex.axis standard plotting parameters
#' @param line.lty plotting type for lines (1 = solid, 2 = dashed ...)
#' @param line.lwd width of lines
#' @param region.lty plotting type for 95\% regions (1 = solid, 2 = dashed ...)
#' @param region.lwd width of line for 95\% regions
#' @param emphasize.axes if TRUE, coordinate axes are marked with a black line
#' @return none
#' @examples
#' visualize.line.models(c(0.1,0.2), c(1,0), c(0.995,0))
#' @export
# Note: Does NOT allow specification of 'scale.weights'
# but uses the given scales as the single component for each distribution.
# (This is because HDR of a mixture of Gaussians is not simple to compute.)
visualize.line.models <- function(scales, slopes, cors,
                                  model.names = NULL, model.cols = NULL,
                                  legend.position = "bottomright",
                                  xlim = NULL, ylim = NULL,
                                  xlab = "EFFECT1", ylab = "EFFECT2",
                                  cex.lab = 1, cex.axis = 1,
                                  line.lty = 1, region.lty = 2,
                                  line.lwd = 1, region.lwd = 1,
                                  emphasize.axes = TRUE){

  K = length(slopes) #number of models
  if(length(scales) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(length(cors) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(any(cors > 1) | any(cors < 0)) stop("Some value of 'cors' is outside [0,1].")
  lim = 3*max(scales) #Default: show models within 3 SDs
  if(is.null(xlim)) xlim = c(-lim,lim)
  if(is.null(ylim)) ylim = c(-lim,lim)
  plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       cex.lab = cex.lab, cex.axis = cex.axis)
  if(is.null(model.names)) model.names = paste0("M",1:K)
  if(is.null(model.cols)) model.cols = 1:K
  prob.level = 0.95
  b = qchisq(prob.level, df = 2)
  grid()
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
      Sigma = prior.V(scale = scales[ii], slope = slopes[ii], cor = cors[ii])
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

#' Visualize effect distributions of line models
#'
#' Plots density functions of distributions of effect sizes
#' of line models specified by scales and scale.weights.
#' See help of 'line.models( )' for details about specifying the models.
#'
#' @param scales vector of standard deviations of larger effect of each model
#' @param scale.weights vector of weights of each component or a matrix
#'  with one row per model and one column per component
#' @param model.names vector of names of each model
#' @param model.cols vector of colors for each model
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
#' Consider two effects:  x = effect1 and y = effect2.
#' The line models are defined by three parameters:
#' scale = standard deviation of the (larger) effect,
#' slope = slope of line around which the effects are scattered,
#' cor = non-negative correlation between effects where cor = 1 means equal effect size.
#'
#' A line model has the following properties:
#'
#' 1) effects are scattered around line y = slope*x.
#'    'slope' can be any real number or Inf (in which case effect x is zero).
#'
#' 2) the larger of the prior variances of effects is scale^2.
#'    That is, if |slope| <= 1, then Var(x) = scale^2 and
#'             if |slope| > 1, then Var(y) = scale^2.
#'
#' 3) Distribution around the line is determined as follows.
#'    Consider a distribution scattered around line y = x with the given correlation.
#'    Rotate that distribution by an orthogonal rotation defined by angle
#'    theta = atan(slope) - pi/4
#'    and use the corresponding distribution, scaled so that the maximum variance
#'    of the two effects is scale^2.
#'    NOTE: This is not same as having the given correlation around line y = slope*x,
#'          because the shape of that distribution depends on the slope but
#'          definition we use here is independent of slope (up to an orthogonal transformation).
#'
#' Effect distribution for the larger effect can be specified as a mixture
#' of Gaussians.
#' If 'scale.weights' is a vector of length C > 1, then model 'M' is a mixture
#' of C two-dimensional Gaussians where component i has mean of 0,
#' the standard deviation for the larger effect of 'scales[M]/(C - i + 1)'
#' and the mixture weight is proportional to scale.weights[M].
#' If 'scale.weights' is a matrix with C > 1 columns and K rows,
#' then model 'M' is a mixture of 'C' Gaussians specified by row 'M'
#' of 'scale.weights'.
#'
#' @param X matrix of effect sizes with two columns (one col per trait)
#' @param SE matrix of standard errors with two columns (one col per trait)
#' @param scales vector of standard deviations of larger effect of each model
#' @param slopes vector of slopes of each model
#' @param cors vector of correlation parameters of each model
#' @param model.names vector of names of each model
#' @param model.priors vector of prior probabilities of models up to a normalization constant
#' @param r.lkhood correlation between estimators of the two effects
#' @param scale.weights vector of weights of each component or a matrix
#'  with one row per model and one column per component
#' @return matrix with one column per model and one row per variable giving
#' membership probabilities of each variable in each model
#' @examples
#' line.models(X = linemodels.ex1[,1:2], SE = linemodels.ex1[,3:4],
#'            scales = c(0.2, 0.2, 0.2),
#'            slopes = c(0, 0.5, 1),
#'            cors = c(0.995, 0.995, 0.995))
#' @export
line.models <- function(X, SE,
                        scales, slopes, cors,
                        model.names = NULL,
                        model.priors = rep(1/length(slopes), length(slopes)),
                        r.lkhood = 0, scale.weights = c(1)){

  # Evaluates model probabilities for each data point separately.
  # INPUT
  # X, N x 2 matrix where rows are observations/variables
  #  and columns are the two effect sizes (x and y coordinates).
  # SE, N x 2 matrix where rows are observations/variables
  #  and columns are the two standard errors of effects in X.
  # scales, vector of positive standard deviations of the larger effect size
  #         of the line models. Length is K, the number of models.
  # slopes, vector of slopes of the line models. Length is K.
  # cors, vector of correlation parameters of the line models. Length is K.
  # model.names, names to be used in output for the models. Length is K.
  # model.priors, prior probability of each of the K models.
  #               Default is 1/K for each model.
  # r.lkhood, determines correlation of the two effect estimators.
  #           Default is 0, which assumes that estimators were uncorrelated.
  # scale.weights:
  # Effect distribution for the larger effect can be specified as a mixture
  # of Gaussians.
  # If 'scale.weights' is a vector of length C > 1, then model 'M' is a mixture
  # of C Gaussians:
  # sum_{i=1}^C scale.weights[i]*N(0, prior.V(scales[M]/(C - i + 1), slopes[M], cols[M]))
  # where vector 'scale.weights' has been normalized to sum to 1.
  # If 'scale.weights' is a matrix with C > 1 columns and K rows,
  # then model 'M' is a mixture of C Gaussians:
  #  sum_{i=1}^k scale.weights[M,i]*N(0, scales[M]/(C - i + 1))
  # where rows of matrix 'scale.weights' has been normalized to sum to 1.
  # Default value scale.weights = c(1) means that each model has only one component
  #  whose standard deviation is 'scales[M]'.

  # OUTPUT:
  # Returns a matrix with membership probabilities (in columns)
  # for each observations (in rows).

  K = length(slopes) #number of models
  if(length(scales) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(length(cors) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(any(cors > 1) || any(cors < 0)) stop("Some value of 'cors' is outside [0,1].")
  if(r.lkhood < (-1) || r.lkhood > 1) stop("'r.lkhood' is outside [-1,1].")
  if(ncol(X) != 2) stop("Input data 'X' must have exactly two columns.")
  if(ncol(SE) != 2) stop("Input data 'SE' must have exactly two columns.")
  if(nrow(X) != nrow(SE)) stop("Input data 'X' and SE' must have the same no. of columns.")
  if(any(model.priors < 0)) stop("All 'model.priors' must be non-negative.")
  if(any(scale.weights < 0)) stop("All 'scale.weights' must be non-negative.")

  pis = model.priors/sum(model.priors)
  if(is.vector(scale.weights)){
    scale.weights = matrix(scale.weights, byrow = T,
                              ncol = length(scale.weights), nrow = K)
  }
  if(!is.matrix(scale.weights)) stop("'scale.weights' should be vector or matrix.")
  if(nrow(scale.weights) != K) stop("Dimensions of 'scales' and 'scale.weights' do not match.")
  n.comps = ncol(scale.weights)

  pr.V = list() #list of lists of prior matrices of each model and each scale
  for(kk in 1:K) { #over models
    w = scale.weights[kk,]
    w = w/sum(w)
    model.Vs = list()
    for(jj in 1:n.comps){ #over scales
      model.Vs[[jj]] = prior.V(scale = scales[kk]/(n.comps + 1 - jj),
                               slope = slopes[kk], cor = cors[kk])}
    pr.V[[kk]] = model.Vs
  }

  n = nrow(X)
  R.lkhood = matrix(c(1, r.lkhood, r.lkhood, 1), 2, 2, byrow = T)
  logdnorm = matrix(NA, ncol = K, nrow = n)
  for(ii in 1:n){
    V = diag(SE[ii,]) %*% R.lkhood %*% diag(SE[ii,]) #var of likelihood
    for(kk in 1:K){
      tmp = sapply(1:n.comps, function(jj){
        log.dmvnorm(X[ii,], mu = c(0,0), S = V + pr.V[[kk]][[jj]])})
      tmp = tmp + log(scale.weights[kk,])
      logdnorm[ii,kk] = log(sum(exp(tmp)))
    }
  }
  p = exp(t(t(logdnorm) + log(pis)))
  res = p/rowSums(p)
  colnames(res) = model.names
  return(res)
}


#' Proportion parameters and membership probabilities of line models
#'
#' Consider two effects:  x = effect1 and y = effect2.
#' The line models are defined by three parameters:
#' scale = standard deviation of the (larger) effect,
#' slope = slope of line around which the effects are scattered,
#' cor = non-negative correlation between effects where cor = 1 means equal effect size.
#'
#' A line model has the following properties:
#'
#' 1) effects are scattered around line y = slope*x.
#'    'slope' can be any real number or Inf (in which case effect x is zero).
#'
#' 2) the larger of the prior variances of effects is scale^2.
#'    That is, if |slope| <= 1, then Var(x) = scale^2 and
#'             if |slope| > 1, then Var(y) = scale^2.
#'
#' 3) Distribution around the line is determined as follows.
#'    Consider a distribution scattered around line y = x with the given correlation.
#'    Rotate that distribution by an orthogonal rotation defined by angle
#'    theta = atan(slope) - pi/4
#'    and use the corresponding distribution, scaled so that the maximum variance
#'    of the two effects is scale^2.
#'    NOTE: This is not same as having the given correlation around line y = slope*x,
#'          because the shape of that distribution depends on the slope but
#'          definition we use here is independent of slope (up to an orthogonal transformation).
#'
#' The prior distribution of the mixture proportions is Dirichlet(diri.prior) and
#' a Gibbs sampler is used to estimate the posterior.
#'
#' @param X matrix of effect sizes with two columns (one col per trait)
#' @param SE matrix of standard errors with two columns (one col per trait)
#' @param scales vector of standard deviations of larger effect of each model
#' @param slopes vector of slopes of each model
#' @param cors vector of correlation parameters of each model
#' @param model.names vector of names of each model
#' @param r.lkhood correlation between estimators of the two effects
#' @param n.iter number of Gibbs sampler iterations (after burn-in)
#' @param n.burnin number of burn-in iterations that will be discarded
#' @param diri.prior parameters of the Dirichlet distribution used as prior for proportions
#' @param verbose if TRUE, prints the index of every 100th iteration
#' @return List with two components: (1) 'params', matrix of
#' posterior distribution of proportion parameters for each model.
#' Columns are mean, lower and upper points of 95\% credible interval
#' and standard deviation. (2)  'groups',
#' matrix with one column per model and one row per variable giving
#' membership probabilities of each variable in each model.
#' @examples
#' line.models.with.proportions(
#'            X = linemodels.ex1[,1:2], SE = linemodels.ex1[,3:4],
#'            scales = c(0.2, 0.2, 0.2),
#'            slopes = c(0, 0.5, 1),
#'            cors = c(0.995, 0.995, 0.995))
#' @export
line.models.with.proportions <- function(X, SE,
                                         scales, slopes, cors,
                                         model.names = NULL,
                                         r.lkhood = 0,
                                         n.iter = 200, n.burnin = 20,
                                         diri.prior = rep(1/length(scales),length(scales)),
                                         verbose = TRUE){

 # Gibbs sampler to evaluates model probabilities for all data points jointly
 #  together with the posterior of the mixture proportions 'pi'.
 # INPUT
 # X, N x 2 matrix where rows are observations/variables
 #  and columns are the two effect sizes (x and y coordinates).
 # SE, N x 2 matrix where rows are observations/variables
 #  and columns are the two standard errors of effects in X.
 # scales, vector of positive standard deviations of the larger effect size
 #         of the line models. Length is K, the number of models.
 # slopes, vector of slopes of the line models. Length is K.
 # cors, vector of correlation parameters of the line models. Length is K.
 # model.names, names to be used in output for the models. Length is K.
 # r.lkhood, determines correlation of the two effect estimators.
 #           Default is 0, which assumes that estimators were uncorrelated.
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
 # 'groups' row per observation and membership probabilities in columns
 # 'params' posterior distribution of the proportions.

  K = length(slopes) #number of models
  if(length(scales) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(length(cors) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(any(cors > 1) || any(cors < 0)) stop("Some value of 'cors' is outside [0,1].")
  if(r.lkhood < (-1) || r.lkhood > 1) stop("'r.lkhood' is outside [-1,1].")
  if(ncol(X) != 2) stop("Input data 'X' must have exactly two columns.")
  if(ncol(SE) != 2) stop("Input data 'SE' must have exactly two columns.")
  if(nrow(X) != nrow(SE)) stop("Input data 'X' and SE' must have the same no. of columns.")
  if(any(diri.prior <= 0)) stop("All values in 'diri.prior' must be positive.")
  if(length(diri.prior) != K) stop("Length of 'diri.prior' do not match with number of models.")
  if(n.burnin < 0) stop("'n.burnin' must be non-negative.")
  if(n.iter < 1) stop("'n.iter' must be positive.")

  pr.V = list()
  for(kk in 1:K) pr.V[[kk]] = prior.V(scale= scales[kk], slope = slopes[kk], cor = cors[kk])
  n = nrow(X)

  logdnorm = matrix(NA, ncol = K, nrow = n)
  R.lkhood = matrix(c(1, r.lkhood, r.lkhood, 1), 2, 2, byrow = T)
  for(ii in 1:n){
    V = diag(SE[ii,]) %*% R.lkhood %*% diag(SE[ii,]) #var of likelihood
     for(kk in 1:K){
        logdnorm[ii,kk] = log.dmvnorm(X[ii,], mu = c(0,0), S = V + pr.V[[kk]])
     }
  }

  R.par = matrix(NA, ncol = K, nrow = n.iter) #results for parameters
  R.ind = matrix(0, ncol = K, nrow = n) #results for individual effects
  colnames(R.par) = colnames(R.ind) = model.names

  pis = rep(1, K) / K #initial probabilities, vector pi
  grs = sample(1:K, size = n, prob = pis, replace = T) #initially random grouping among variants

  for(ii in 1:(n.burnin + n.iter)){ #Gibbs sampler

      #count how many are assigned to each group
      tmp = rle(sort(grs))
      gr.counts = rep(0, K)
      gr.counts[tmp$values] = tmp$length

      #sample pis
      pis = rdirichlet(diri.prior + gr.counts)
      log.pis = log(pis)

      #sample group indicators for variants
      p = exp(t(t(logdnorm) + log.pis)) #no need to normalize for 'sample()'
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
                                  model.priors = rep(1/length(slopes), length(slopes)),
                                  model.names = NULL,
                                  r.lkhood = 0, return.model.post = FALSE){

  # Evaluates log-likelihood of the data given the parameters
  #  X, SE, scales, slopes, cors, model.priors, model.names, r.lkhood
  #  are like in line.models().
  # Does not allow scale.weights from line-models().
  #
  # If return.model.post = FALSE, then returns
  #  the expected log likelihood over posterior of memberships
  #  as needed by an EM-algorithm.
  # If return.model.post = TRUE, then returns
  #  posterior probabilities of observations' (rows)
  #  membership in each model (columns).

  K = length(slopes) #number of models
  if(is.null(model.names)) model.names = paste0("M",1:K)
  if(length(scales) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(length(cors) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(length(model.names) != K) stop("Length of 'model.names' and 'slopes' do not match.")
  if(any(cors > 1) || any(cors < 0)) stop("Some value of 'cors' is outside [0,1].")
  if(r.lkhood < (-1) || r.lkhood > 1) stop("'r.lkhood' is outside [-1,1].")
  if(ncol(X) != 2) stop("Input data 'X' must have exactly two columns.")
  if(ncol(SE) != 2) stop("Input data 'SE' must have exactly two columns.")
  if(nrow(X) != nrow(SE)) stop("Input data 'X' and SE' must have the same no. of columns.")
  if(any(model.priors < 0)) stop("All 'model.priors' must be non-negative.")

  pis = model.priors / sum(model.priors)

  pr.V = list() #list of prior matrices of each model
  for(kk in 1:K) pr.V[[kk]] = prior.V(scale = scales[kk], slope = slopes[kk], cor = cors[kk])

  n = nrow(X)
  R.lkhood = matrix(c(1, r.lkhood, r.lkhood, 1), 2, 2)
  logdnorm = matrix(NA, ncol = K, nrow = n)
  for(ii in 1:n){
    V = diag(SE[ii,]) %*% R.lkhood %*% diag(SE[ii,]) #var of likelihood
    for(kk in 1:K){
      logdnorm[ii,kk] = log.dmvnorm(X[ii,], mu = c(0,0), S = V + pr.V[[kk]])
    }
  }

  p = exp(t(t(logdnorm) + log(pis)))
  post = p/rowSums(p)
  colnames(post) = model.names
  if(return.model.post){
    return(post)
  }else{
    #Compute actual log-lkhood for given parameters:
    # log.lkhood = sum(log(rowSums(exp(t(t(logdnorm) + log(pis))))))

    #Compute expectation of log-lkhood under posterior of membership.
    # This is what is needed for EM-algorithm.
    log.lkhood = sum(log(rowSums(exp(logdnorm + log(post)))))
    return(log.lkhood)
  }
}



current.to.par <- function(current, par.include){
  # Transforms 'current' matrix to 'par' vector
  # that has only those parameters that are optimized over
  # while 'current' has the current value of all parameters.
  # Since output of this function will be used as initial
  # value in optim(), the output values must be finite.

  par = c()
  ii = 1
  for(k in 1:nrow(par.include)){
    if(par.include[k,1]) { #log transform scales
      par[ii] = min(c(log(current[k,1]), log(1e300))) #max scale 1e+300
      ii = ii+1}
    if(par.include[k,2]) { #turn slopes to angles
      par[ii] = atan(current[k,2])
      ii = ii+1}
    if(par.include[k,3]) { #logit transform correlations
      a = log(current[k,3]) - log(1 - current[k,3])
      # logit <= 40 --> cor > 4e-18
      if(abs(a) > 40) a = 40*(-1)^as.integer(a < 0)
      par[ii] = a
      ii = ii+1}
  }
  return(par)
}



par.to.current <- function(par, current, par.include){
  # Transforms 'par' vector to 'current' matrix
  # that has the current value of all parameters
  # while 'par' has only those that are optimized over.

  ii = 1
  K = nrow(par.include)
  for(k in 1:K){
    if(par.include[k,1]) {current[k,1] = exp(par[ii]); ii = ii+1}
    if(par.include[k,2]) {current[k,2] = tan(par[ii]); ii = ii+1}
    if(par.include[k,3]) {current[k,3] = 1/(1 + exp(-par[ii])); ii = ii+1}
  }
  return(current)
}



optim.fn <- function(par, current, par.include, X, SE,
                     model.priors, model.names, r.lkhood){

  # Defines log-lkhood function to be optimized by the EM-algorithm.
  # 'current' has current values of the parameters (1 row, 3 cols per model)
  # 'par.include' matrix defines which parameters are to be optimized.
  #       Matrix has 1 row and 3 cols per model and columns correspond to
  #       1 = scale
  #       2 = slope
  #       3 = correlation
  #NOTE: scales have been log-transformed and correlations logit-transformed before
  #      calling this function. Hence back-transforming them first.

  current = par.to.current(par, current, par.include)
  -line.models.loglkhood(X, SE,
                         scales = current[,1], slopes = current[,2], cors = current[,3],
                         model.priors = model.priors,
                         model.names = model.names, r.lkhood = r.lkhood,
                         return.model.post = FALSE)
}


#' Optimize parameters of line models
#'
#' Consider two effects:  x = effect1 and y = effect2.
#' The line models are defined by three parameters:
#' scale = standard deviation of the (larger) effect,
#' slope = slope of line around which the effects are scattered,
#' cor = non-negative correlation between effects where cor = 1 means equal effect size.
#'
#' A line model has the following properties:
#'
#' 1) effects are scattered around line y = slope*x.
#'    'slope' can be any real number or Inf (in which case effect x is zero).
#'
#' 2) the larger of the prior variances of effects is scale^2.
#'    That is, if |slope| <= 1, then Var(x) = scale^2 and
#'             if |slope| > 1, then Var(y) = scale^2.
#'
#' 3) Distribution around the line is determined as follows.
#'    Consider a distribution scattered around line y = x with the given correlation.
#'    Rotate that distribution by an orthogonal rotation defined by angle
#'    theta = atan(slope) - pi/4
#'    and use the corresponding distribution, scaled so that the maximum variance
#'    of the two effects is scale^2.
#'    NOTE: This is not same as having the given correlation around line y = slope*x,
#'          because the shape of that distribution depends on the slope but
#'          definition we use here is independent of slope (up to an orthogonal transformation).
#'
#' This function optimizes over any set of scales, slopes and cors.
#' The parameters not included in optimization are fixed to their initial values.
#' Proportion parameters are always optimized.
#'
#' @param X matrix of effect sizes with two columns (one col per trait)
#' @param SE matrix of standard errors with two columns (one col per trait)
#' @param par.include matrix of TRUE/FALSE values TRUE indicating which parameters
#' are optimized. One row per model and 3 columns corresponding
#' to 'scales', 'slopes' and 'cors', respectively.
#' @param init.scales vector of initial standard deviations of larger effect of each model
#' @param init.slopes vector of initial slopes of each model
#' @param init.cors vector of initial correlation parameters of each model
#' @param model.priors vector of initial prior probabilities of models up to a normalization constant
#' @param model.names vector of names of each model
#' @param r.lkhood correlation between estimators of the two effects
#' @param tol tolerance for convergence in log-likelihood between
#' adjacent iterations

#' @return matrix with models in rows and columns:
#' (1) scales, (2) slopes, (3) correlations.
#' @examples
#' line.models.optimize(
#'      X = linemodels.ex1[,1:2],
#'      SE = linemodels.ex1[,3:4],
#'      par.include = rbind(c(F,F,F), c(F,T,T), c(F,F,F)),
#'      init.scales = c(0.2, 0.2, 0.2),
#'      init.slopes = c(0, 0.2, 1),
#'      init.cors = c(0.995, 0.5, 0.995))
#' @export
line.models.optimize <- function(X, SE,
                                 par.include = matrix(TRUE, nrow = length(init.slopes), ncol = 3),
                                 init.scales, init.slopes, init.cors,
                                 model.priors = rep(1,length(init.slopes)),
                                 model.names = NULL,
                                 r.lkhood = 0, tol =  1e-3){

  # EM algorithm to optimize chosen parameters of the line models.
  #
  # X, SE, model.names, r.lkhood are as in line.models().
  #
  # par.include defines which of the parameters will be optimized
  #  TRUE = optimized, FALSE = not optimized.
  # Parameters are specified in a matrix with 1 row per model and 3 columns are
  #  1 = scale
  #  2 = slope
  #  3 = correlation
  #
  # init.scales, init.slopes, init.cors, model.priors, give initial values
  # for scales, slopes, cors and proportion parameters, respectively.
  # Note that those scales/slopes/cors that are not being optimized over
  #  will remain fixed to their initial values. Proportion parameters
  #  are always optimized over by the algorithm.
  #
  # tol, defines tolerance for maximum difference in log-likelihood
  #  in adjacent iterations that is still considered to show convergence.

  # OUTPUT:
  # matrix with all model parameters.
  #  Rows per models, 3 columns for
  #  (1) scales, (2) slopes, (3) correlations.

  K = length(init.slopes)
  if(length(init.scales) != K) stop("Lengths of init.slopes and init.scales do not match.")
  if(length(init.cors) != K) stop("Lengths of init.slopes and init.cors do not match.")
  if(length(model.priors) != K) stop("Lengths of model.priors and init.slopes do not match.")
  if(any(init.cors < 0) || any(init.cors > 1)) stop("Correlation values should be in [0,1].")
  if(abs(r.lkhood) > 1) stop("Correlation of likelihood, 'r.lkhood', should be in [-1,1].")
  if(any(init.scales <= 0)) stop("Scales should be positive.")
  if(any(model.priors <= 0)) stop("'model.priors' should be positive.")
  if(!is.matrix(par.include) || nrow(par.include) != K || ncol(par.include) != 3)
    stop("par.include should be a matrix with 1 row and 3 cols per model")

  iter = 0
  current = cbind(init.scales, init.slopes, init.cors) #current values of parameters
  converged = FALSE
  log.lkhood = -Inf

  cat(paste0("Initial values\n"))
  cat(paste("scales:",paste(signif(current[,1],4),collapse=", ")),"\n")
  cat(paste("slopes:",paste(signif(current[,2],4),collapse=", ")),"\n")
  cat(paste("cors:",paste(signif(current[,3],4),collapse=", ")),"\n")
  par.names = cbind(paste0("scale",1:K), paste0("slope",1:K), paste0("cor",1:K))
  cat(paste("Optimizing w.r.t:",
              paste(par.names[par.include], collapse = " ")))

  #Nelder-Mead is not robust for 1-dimensional optimization --> uses "Brent" in 1-dim case
  # but Brent requires finite range of optimization so needs to set the range
  # depending on which type of the parameter is optimized.

  op.method = ifelse(sum(par.include) > 1, "Nelder-Mead", "Brent")
  lower = -Inf # Leaves lower and upper as Infs if Nelder-Mead is used
  upper = Inf
  if(op.method == "Brent"){ # otherwise specify the range for the single parameter
    if(sum(par.include[,1])){ # log-transformed scale parameter (1e-6, 1e+6)
      lw.orig = 1e-6; lower = log(lw.orig)
      up.orig = 1e+6; upper = log(up.orig)
    }
    if(sum(par.include[,2])) { # atan transformed slope parameter (-Inf, Inf)
      lw.orig = -Inf; lower = atan(lw.orig)
      up.orig = Inf; upper = atan(up.orig)
    }
    if(sum(par.include[,3])) { # log-odds transformed correlation parameter (1e-6, 1-1e-6)
      lw.orig = 1e-6; lower = log(lw.orig / (1 - lw.orig))
      up.orig = 1 - 1e-6; upper = log(up.orig / (1 - up.orig))
    }
    cat(paste0(", over the range: (",
               signif(lw.orig,4),",",signif(up.orig,4),")"))
  }
  cat("\n\n")

  while(!converged){
    iter = iter + 1
    pr = line.models.loglkhood(X, SE,
                               scales = current[,1], slopes = current[,2], cors = current[,3],
                               model.priors = model.priors,
                               model.names = model.names, r.lkhood = r.lkhood,
                               return.model.post = TRUE)

    model.priors = apply(pr, 2, mean) #these proportions maximize the lkhood

    #Maximizes over other parameters than proportions:
    par = current.to.par(current, par.include)
    opt.out = optim(par, optim.fn, current = current, par.include = par.include,
                    X = X, SE = SE, model.priors = model.priors, model.names = model.names,
                    r.lkhood = r.lkhood, method = op.method, lower = lower, upper = upper)

    converged = (-opt.out$value - log.lkhood  < tol)
    if(-opt.out$value > log.lkhood){
      log.lkhood = -opt.out$value
      current = par.to.current(opt.out$par, current, par.include)
      cat(paste("iter:",iter,"; log-lkhood:",signif(log.lkhood,8)),"\n")
      cat(paste("proportions:", paste(signif(model.priors,3),collapse =", ")),"\n")
      cat(paste("scales:", paste(signif(current[,1],4),collapse =", ")),"\n")
      cat(paste("slopes:", paste(signif(current[,2],4),collapse =", ")),"\n")
      cat(paste("cors:", paste(signif(current[,3],4),collapse =", ")),"\n\n")
    }
    else{cat(paste("iter:",iter,"; Failed to increase log-likelihood further.\n"))}
  }
  colnames(current) = c("scales", "slopes", "cors")
  return(current)
}


#' Random sample from a line model
#'
#' See help of 'line.models( )' for details about specifying the models.
#'
#' @param n the number of samples
#' @param scale standard deviation of larger effect
#' @param slope slopes of line model
#' @param cor correlation parameters of line model
#' @param scale.weights vector of weights of each component when effect size
#' distribution is a mixture model
#' @return matrix with one row per sample and one column per effect on each outcome
#' @examples
#' sample.line.model(10, scale = 0.1, slope = 1, cor = 0.95)
#' sample.line.model(10, scale = 0.1, slope = 1, cor = 0.95, scale.weights = c(1,0,1))
#' @export
sample.line.model <- function(n = 1, scale, slope, cor, scale.weights = c(1)){

  # Outputs 'n' samples from a given line model defined by
  # scale, slope and cor and scale.weights.
  # See line.models() for definition of these parameters.
  # 'scale.weights' must be a vector and cannot be a matrix as in line.models( ).

  if(cor > 1) stop("cor > 1")
  if(cor < 0) stop("cor < 0")
  if(scale < 0) stop("scale < 0")
  if(!is.vector(scale.weights)) stop("scale.weights must be vector.")
  if(any(scale.weights < 0)) stop("scale.weights can have only positive values.")
  J = length(scale.weights)

  V = prior.V(scale = scale, slope = slope, cor = cor)

  #Choose components for each sample
  ind = sample(1:J, size = n, prob = scale.weights, replace = T)

  if(abs(V[1] * V[4] - V[2] * V[3]) < 1e-16){
    A = matrix(c(sqrt(V[1]), sqrt(V[4]), 0, 0), ncol = 2)}
  else {
    A = t(chol(V))}
  x = t(A %*% matrix(rnorm(2*n), nrow = 2) / rep(J + 1 - ind, each = 2))
  cat(paste0("Sampling effects with scale=", scale," slope=", slope," cor=", cor),"\n")
  cat(paste0("   scale.weights=", paste(scale.weights, collapse = ","),"."),"\n")
  cat(paste0("Theoretical SD for larger effect:",
             signif(sqrt(sum(scale^2*scale.weights/sum(scale.weights)/rev(1:J)^2)),5)))
  cat(paste0("; observed value:", signif(max(apply(x,2,sd)),5)),".\n")
  return(x)
}

