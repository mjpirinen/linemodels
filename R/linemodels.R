# linemodels (0.2.1): An R package to cluster 2-dimensional effects into
#  groups defined by linear relationships.
# DATE: 7-Aug-2023
# Copyright (C) 2022-2023, Matti Pirinen
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
# cor = non-negative correlation between effects,
# where cor = 1 means the same effect size and cor = 0 means independent effects
# A line model has the following properties:
#
# 1) effects are scattered around line y = slope*x.
#    'slope' can be any real number or Inf (in which case effect x = 0)
#
# 2) the larger of the prior variances of effects is scale^2
#    That is, if |slope| <= 1, then Var(x) = scale^2
#             if |slope| > 1, then Var(y) = scale^2
#
# 3) Distribution of effects around the line is determined as follows:
#    Consider a distribution scattered around line y = x with correlation 'cor'.
#    Rotate that distribution by an orthogonal rotation defined by angle
#    theta = atan(slope) - pi/4
#    and use the corresponding distribution, scaled so that the maximum variance
#    of the two effects is scale^2.
#    NOTE: This is not same as "correlation 'cor' around line y = slope*x",
#          because the shape of that distribution depends on the slope, but
#          definition we use here is independent of slope (up to an orthogonal transformation).
#
# Examples:
# Set scale = 0 to get the null model. (Values of slope and cor do not matter in this case.)
# Use scale > 0, slope = 1 and cor = 0 to get independent effects model.
# Use scale > 0, slope = 1 and cor = 1 to get fixed effects model.
# Use scale > 0, slope = 1 and 0 < cor < 1  to get correlated effects model.
#               note that cor value needs to be quite near 1 to get very similar values
#               between the two effects with high probability (e.g. cor = 0.99 or 0.999).
#
# To choose a reasonable scale, note that 95% of the effects are assumed to be < 2*scale.

#
# Seven functions that are meant for user to call directly are listed below.
#

# visualize.line.models(scales, slopes, cors,
#                      model.names = NULL, model.cols = NULL,
#                      legend.position = "bottomright",
#                      xlim = NULL, ylim = NULL,
#                      xlab = "EFFECT1", ylab = "EFFECT2",
#                      cex.lab = 1, cex.axis = 1,
#                      line.lty = 1, region.lty = 2,
#                      line.lwd = 1, region.lwd = 1,
#                      plot.grid = TRUE,
#                      plot.new = TRUE,
#                      emphasize.axes = TRUE)
# -- To visualize the line models and their 95% highest probability regions.

# visualize.classification.regions(
#           scales, slopes, cors, SE,
#           r.lkhood = 0,
#           model.priors = rep(1,length(scales)),
#           col.palette = NULL,
#           add.legend = TRUE,
#           xlim = NULL, ylim = NULL,
#           breakpoints = 200,
#           xlab = "EFFECT1", ylab = "EFFECT2",
#           cex.lab = 1, cex.axis = 1,
#           line.lty = 1, region.lty = 2,
#           line.lwd = 1, region.lwd = 1,
#           emphasize.axes = TRUE)
# -- To visualize the classification regions of two competing line models.

# visualize.scales(scales, scale.weights = c(1),
#                  model.names = NULL, model.cols = NULL,
#                  legend.position = "topleft")
# -- To visualize the univariate distributions of the effect sizes.

# line.models(X, SE,
#             scales, slopes, cors,
#             model.names = NULL,
#             model.priors = rep(1/length(slopes), length(slopes)),
#             r.lkhood = 0, scale.weights = c(1))
# -- To evaluate the model probabilities for each observation separately.

# line.models.with.proportions(X, SE,
#                             scales, slopes, cors,
#                             model.names = NULL,
#                             r.lkhood = 0,
#                             n.iter = 200, n.burnin = 20,
#                             diri.prior = rep(1/length(scales),length(scales)),
#                             verbose = TRUE)
# -- To evaluate the model probabilities together with the proportions of
#    observations coming from each model. Uses Gibbs sampler for a joint
#    estimation of membership probabilities across all observations.

# line.models.optimize(X, SE,
#                      par.include = matrix(TRUE, nrow = length(init.slopes), ncol = 3),
#                      force.same.scales = FALSE,
#                      init.scales, init.slopes, init.cors,
#                      model.priors = rep(1,length(init.slopes)),
#                      model.names = NULL,
#                      r.lkhood = 0,
#                      tol.loglk =  1e-3,
#                      tol.par = 0,
#                      return.weights = FALSE){
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
 # cor >= 0, NOTE: negatively correlated effects are modeled by a negative slope,
 #                 not by a negative correlation

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
  return(g / sum(g))
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
#' @param plot.grid, if TRUE, plots a grid
#' @param plot.new if TRUE, makes a new plot, if FALSE adds to existing plot
#' @param emphasize.axes if TRUE, coordinate axes are marked with a black line
#' @return none
#' @examples
#' visualize.line.models(c(0.1,0.2), c(1,0), c(0.995,0))
#' @export
# Note: Does NOT allow specification of 'scale.weights'
# but uses the given scales as the single component for each distribution.
# This is because probability regions of a mixture of Gaussians is not simple to compute.
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

  K = length(slopes) #number of models
  if(length(scales) != K) stop("Length of 'scales' and 'slopes' do not match.")
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

#' Visualize classification regions of line models
#'
#' Colors the plane according to the classification probabilities
#' of two competing line models.
#' See help of 'line.models( )' for details about specifying the models.
#'
#' @param scales vector of standard deviations of larger effect of each model
#' @param slopes vector of slopes of each model
#' @param cors vector of correlation parameters of each model
#' @param SE vector of two standard errors for effect1 and effect2, respectively
#' @param r.lkhood correlation between estimators of the two effects
#' @param model.priors vector of (unnormalized) prior probabilities of the two models
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

  K = length(slopes) #number of models
  if(K != 2) stop("Should specify exactly two models.")
  if(length(scales) != K) stop("Length of 'scales' and 'slopes' do not match.")
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
#' cor = non-negative correlation between effects,
#' where cor = 1 means the same effect and cor = 0 means independent effects.
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


#' Proportion parameters and membership probabilities of line models
#'
#' Consider two effects:  x = effect1 and y = effect2.
#' The line models are defined by three parameters:
#' scale = standard deviation of the (larger) effect,
#' slope = slope of line around which the effects are scattered,
#' cor = non-negative correlation between effects,
#' where cor = 1 means the same effect and cor = 0 means independent effects.
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
#' a Gibbs sampler is used for estimating the posterior.
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

 # Gibbs sampler to evaluate model probabilities for all data points jointly
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
                                  model.priors = rep(1/length(slopes), length(slopes)),
                                  model.names = NULL,
                                  r.lkhood = 0,
                                  return.posteriors = FALSE){

  # Evaluates log-likelihood of the data given the parameters
  #  X, SE, scales, slopes, cors, model.priors, model.names, r.lkhood
  #  are like in line.models().
  # Does not allow scale.weights from line-models().
  #
  #  Returns the log-likelihood of the mixture model
  #  where the components have the mixture weights proportional to 'model.priors'
  #  and the weights are the same for all the observations.

  # If return.posteriors = TRUE, then returns a list with two members
  #  'loglkhood' and
  #  'membership.prob' matrix with
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
                                           posteriors = posteriors){

  # Evaluates expected log-likelihood of the data given the parameters
  # and the posterior probability of membership of each variant in each model.
  #  X, SE, scales, slopes, cors, r.lkhood
  #  are like in line.models().
  # 'posteriors' is a matrix with one row per variant and one column per model.
  # NOTE: values of posteriors are not checked for non-negativity or row sums = 1.
  # This function is being optimized by the EM-algorithm.

  K = length(slopes) #number of models
  if(length(scales) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(length(cors) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if(any(cors > 1) || any(cors < 0)) stop("Some value of 'cors' is outside [0,1].")
  if(r.lkhood < (-1) || r.lkhood > 1) stop("'r.lkhood' is outside [-1,1].")
  if(ncol(X) != 2) stop("Input data 'X' must have exactly two columns.")
  if(ncol(SE) != 2) stop("Input data 'SE' must have exactly two columns.")
  if(nrow(X) != nrow(SE)) stop("Input data 'X' and SE' must have the same no. of columns.")
  if(nrow(posteriors) != nrow(X)) stop("Input data 'X' and 'posteriors' must have the same no. of rows.")
  if(ncol(posteriors) != K) stop("posteriors must have one column per each model")

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
  return(sum(posteriors * logdnorm)) #expected log-likelihood
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



optim.fn <- function(par, current, par.include, force.same.scales, X, SE,
                     r.lkhood, posteriors){

  # Defines expectation of the log-lkhood function to be optimized by
  # the EM-algorithm.
  # 'current' has current values of the parameters (1 row, 3 cols per model)
  # 'par.include' matrix defines which parameters are to be optimized.
  #       Matrix has 1 row and 3 cols per model and columns correspond to
  #       1 = scale
  #       2 = slope
  #       3 = correlation
  #NOTE: scales have been log-transformed and correlations logit-transformed before
  #      calling this function. Hence back-transforming them first.
  # If 'force.same.scales = TRUE, then sets all scale parameters equal.
  #     in this case, 1st col of par.include is (TRUE, FALSE, FALSE, ..., FALSE)
  # posteriors is matrix with one row per variant and one col per model
  #   and gives the posterior probability of the membership of variant in model

  current = par.to.current(par, current, par.include)
  if(force.same.scales) current[,1] = current[1,1]
  -line.models.expected.loglkhood(X, SE,
                         scales = current[,1], slopes = current[,2], cors = current[,3],
                         r.lkhood = r.lkhood, posteriors = posteriors)
}



#' Optimize parameters of line models
#'
#' Consider two effects:  x = effect1 and y = effect2.
#' The line models are defined by three parameters:
#' scale = standard deviation of the (larger) effect,
#' slope = slope of line around which the effects are scattered,
#' cor = non-negative correlation between effects,
#' where cor = 1 means the same effect and cor = 0 means independent effects.
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
#' @param force.same.scales logical; If TRUE, then forces all scale parameters equal
#' @param init.scales vector of initial standard deviations of larger effect of each model
#' @param init.slopes vector of initial slopes of each model
#' @param init.cors vector of initial correlation parameters of each model
#' @param model.priors vector of initial prior probabilities of models up to a normalization constant
#' @param model.names vector of names of each model
#' @param r.lkhood correlation between estimators of the two effects
#' @param tol.loglk tolerance for convergence in adjacent log-likelihoods
#' @param tol.par tolerance for convergence in maximum of absolute values of relative differences
#' across parameters between adjacent iterations. Can be set negative to determine
#' the convergence solely by log-likelihood.
#' @param return.weights logical; If TRUE returns both optimized mixture weights and
#' optimized parameters. If FALSE returns only optimized parameters.

#' @return If 'return.weights' = FALSE, returns matrix with models in rows and columns:
#' (1) scales, (2) slopes, (3) correlations.
#' If 'return.weights = TRUE, then returns a list with two components
# 'weights' giving the mixture weights of the components and
# 'parameters' giving the matrix as above.

#' @examples
#' line.models.optimize(
#'      X = linemodels.ex1[,1:2],
#'      SE = linemodels.ex1[,3:4],
#'      par.include = rbind(c(F,F,F), c(F,T,T), c(F,F,F)),
#'      force.same.scales = FALSE,
#'      init.scales = c(0.2, 0.2, 0.2),
#'      init.slopes = c(0, 0.2, 1),
#'      init.cors = c(0.995, 0.5, 0.995))
#' @export
line.models.optimize <- function(X, SE,
                                 par.include = matrix(TRUE,
                                                      nrow = length(init.slopes),
                                                      ncol = 3),
                                 force.same.scales = FALSE,
                                 init.scales, init.slopes, init.cors,
                                 model.priors = rep(1,length(init.slopes)),
                                 model.names = NULL,
                                 r.lkhood = 0,
                                 tol.loglk =  1e-3,
                                 tol.par = 0,
                                 return.weights = FALSE){

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
  # tol.loglk, defines tolerance for maximum difference in log-likelihood
  #  between adjacent iterations that is still considered to show convergence.
  # tol.par, defines tolerance for abolute value of maximum relative difference
  # in optimized parameters between adjacent iterations that is still
  # considered to show convergence.
  # Convergence happens when either of the tolearances above are achieved.
  # NOTE: tol.par can be set negative to determine convergence by only loglkhood.

  # OUTPUT:
  # matrix with all model parameters.
  #  Rows per models, 3 columns for
  #  (1) scales, (2) slopes, (3) correlations.
  # If 'return.weights = TRUE, then list with two components
  # 'proportions' giving the mixture weights of the components
  # 'parameters' matrix as above.

  K = length(init.slopes)
  if(length(init.scales) != K) stop("Lengths of init.slopes and init.scales do not match.")
  if(length(init.cors) != K) stop("Lengths of init.slopes and init.cors do not match.")
  if(length(model.priors) != K) stop("Lengths of model.priors and init.slopes do not match.")
  if(any(init.cors < 0) || any(init.cors > 1)) stop("Correlation values should be in [0,1].")
  if(abs(r.lkhood) > 1) stop("Correlation of likelihood, 'r.lkhood', should be in [-1,1].")
  if(any(init.scales <= 0)) stop("Scales should be positive.")
  if(any(model.priors <= 0)) stop("'model.priors' should be positive.")
  if(tol.loglk <= 0) stop("Tolerance 'tol.loglk' should be positive.")
  if(!is.matrix(par.include) || nrow(par.include) != K || ncol(par.include) != 3)
    stop("par.include should be a matrix with 1 row and 3 cols per model")

  current = cbind(init.scales, init.slopes, init.cors) #current values of parameters
  if(force.same.scales) {
    current[,1] = mean(current[,1]) #initialize
    par.include[, 1] = c(TRUE, rep(FALSE, nrow(par.include)-1)) #only 1 parameter
  }
  w = model.priors/sum(model.priors)
  log.lkhood = line.models.loglkhood(X, SE,
                                     scales = current[,1], slopes = current[,2], cors = current[,3],
                                     model.priors = w,
                                     model.names = model.names, r.lkhood = r.lkhood,
                                     return.posteriors = FALSE)

  cat(paste0("Initial values\n"))
  cat(paste("scales:",paste(signif(current[,1],4),collapse=", ")),"\n")
  cat(paste("slopes:",paste(signif(current[,2],4),collapse=", ")),"\n")
  cat(paste("cors:",paste(signif(current[,3],4),collapse=", ")),"\n")
  cat(paste("proportions:", paste(signif(w,3),collapse =", ")),"\n")
  cat(paste("Initial log-lkhood:",signif(log.lkhood,8),"\n"))
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
  if(force.same.scales) cat(paste("\nForcing all scales to be equal"))
  cat("\nConvergence criteria:\n")
  cat(paste(" Relative diff in parameters <=",tol.par," or\n"))
  cat(paste(" Difference in log-likelihood <",tol.loglk))
  cat("\n\n")

  iter = 0
  converged = FALSE

  while(!converged){
    iter = iter + 1
    prev.log.lkhood = log.lkhood
    prev.vals = current

    #Maximizes over proportions:
    pr = line.models.loglkhood(X, SE,
                               scales = current[,1], slopes = current[,2], cors = current[,3],
                               model.priors = w,
                               model.names = model.names, r.lkhood = r.lkhood,
                               return.posteriors = TRUE)$membership.prob

    #Maximizes over other parameters than proportions:
    par = current.to.par(current, par.include)
    opt.out = optim(par, optim.fn, current = current, par.include = par.include,
                    force.same.scales = force.same.scales,
                    X = X, SE = SE, posteriors = pr,
                    r.lkhood = r.lkhood, method = op.method, lower = lower, upper = upper)

    #these proportions (w) and parameters (par) maximize the lkhood
    new.w = apply(pr, 2, mean)
    new.par = par.to.current(opt.out$par, current, par.include)

    if(force.same.scales) new.par[,1] = new.par[1,1]

    new.log.lkhood = line.models.loglkhood(X, SE,
                          scales = new.par[,1], slopes = new.par[,2], cors = new.par[,3],
                          model.priors = new.w,
                          model.names = model.names, r.lkhood = r.lkhood,
                          return.posteriors = FALSE)

    if(new.log.lkhood > log.lkhood){
      log.lkhood = new.log.lkhood
      w = new.w
      current = new.par
      cat(paste("iter:",iter,"; log-lkhood:",signif(log.lkhood,8)),"\n")
      cat(paste("Relative diffs in optimized parameters:",
                paste(signif(as.vector(abs(((current - prev.vals)/prev.vals)[par.include])), 3),
                             collapse =", ")),"\n")
      cat(paste("proportions:", paste(signif(w,3),collapse =", ")),"\n")
      cat(paste("scales:", paste(signif(current[,1],4),collapse =", ")),"\n")
      cat(paste("slopes:", paste(signif(current[,2],4),collapse =", ")),"\n")
      cat(paste("cors:", paste(signif(current[,3],4),collapse =", ")),"\n\n")
    }
    else{ #will stop because both convergence criteria below are fulfilled
      cat(paste("iter:",iter,"; Failed to increase log-likelihood further.",
                "(previous:",signif(log.lkhood,8),"new:",signif(new.log.lkhood,8),")\n"))}

    converged.par = all(as.vector(abs((current - prev.vals)/prev.vals)[par.include]) <= tol.par)
    converged.loglk = ((log.lkhood - prev.log.lkhood) < tol.loglk)
    converged = converged.par | converged.loglk
  }

  if(converged.par) cat("Parameter values converged.\n")
  if(converged.loglk) cat("Log-likelihood value converged.\n")
  colnames(current) = c("scales", "slopes", "cors")
  if(return.weights) return(list(weights = w, parameters = current))
  else return(current)
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

