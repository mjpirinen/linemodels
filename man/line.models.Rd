% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linemodels.R
\name{line.models}
\alias{line.models}
\title{Membership probabilities in line models}
\usage{
line.models(
  X,
  SE,
  scales,
  slopes,
  cors,
  model.names = NULL,
  model.priors = rep(1/length(scales), length(scales)),
  r.lkhood = rep(0, ncol(X) * (ncol(X) - 1)/2),
  scale.weights = c(1)
)
}
\arguments{
\item{X}{matrix of estimates with M columns (one col per effect variable and one row per observation)}

\item{SE}{matrix of standard errors with M columns (one col per effect variable and one row per observation)}

\item{scales}{K-vector of standard deviations of the largest effect, one per each linemodel}

\item{slopes}{Kx(M-1)-matrix of slopes, one row per linemodel
If M = 2, can be given also as a K-vector of slopes, one per each linemodel}

\item{cors}{K-vector of correlation parameters, one per each linemodel}

\item{model.names}{K-vector of names of linemodels}

\item{model.priors}{K-vector of prior probabilities of linemodels up to a normalization constant}

\item{r.lkhood}{correlation matrix of the estimators of the M effects
can be given also as a vector of the upper triangular values in the row-major order}

\item{scale.weights}{vector of weights of each mixture component or a matrix
with one row per linemodel and one column per mixture component}
}
\value{
matrix with one column per linemodel and one row per observation giving
membership probabilities of each observation in each linemodel
}
\description{
Consider M effect sizes.
The linemodels are defined by three parameters:
scale = standard deviation of the (largest) effect,
slope = set of M-1 slopes defining the line around which the effects are scattered,
cor = non-negative pairwise correlation between each effect pair,
where cor = 1 means the same effect and cor = 0 means independent effects.
}
\details{
A linemodel has the following properties:

1) Effects are scattered around line defined by vector (1,slope_1,...,slope_(M-1)).
   Each slope_i is the slope between effect variables i+1 and 1.
   Slope can be any real number or +/-Inf.
   If any slope is infinite, the effect of variable 1 and variables with finite slopes are 0.
   It is recommended to keep the slopes finite to avoid problems with interpretation.

2) The largest of the prior variances of the effects is scale^2.
   The effect that has the largest scale is determined by the slopes,
   and the scales of the remaining effects are determined by the slopes as well.

3) Distribution around the line is determined as follows.
   Consider a distribution scattered around the direction of the vector (1,1,...,1)
   with the given constant value of the pairwise correlations between the variables.
   Rotate that distribution by an orthogonal rotation that rotates the direction of the
   vector (1,1,...,1) to the direction of the vector (1,slope_1,...,slope_(M-1))
   and use the corresponding distribution, scaled so that the maximum variance
   among the effects is set to scale^2.
   NOTE: This is not same as having the given correlation around line defined by the slopes,
         because the shape of that distribution depends on the slopes but the
        definition we use here is independent of slopes up to an orthogonal transformation.

Examples:
Set scale = 0 to get the null model. (Values of slope and cor do not matter in this case.)
Use scale > 0, slope = 1 and cor = 0 to get independent effects model for 2-dim case.
Use scale > 0, slope = c(1,1,1) and cor = 1 to get fixed effects model for 4-dim case.
Use scale > 0, slope = 1 and 0 < cor < 1  to get correlated effects model for 2-dim case.
              note that cor value needs to be quite near 1 to get very similar values
              between the two effects with high probability (e.g. cor = 0.99...0.999).

Effect distribution for the largest effect can be specified as a mixture
of Gaussians.
If 'scale.weights' is a vector of length C > 1, then model 'Mo' is a mixture
of C M-dimensional Gaussians where mixture component i has mean of 0,
the standard deviation for the largest effect of 'scales[Mo]/(C - i + 1)'
and the mixture weight proportional to scale.weights[Mo].
If 'scale.weights' is a matrix with C > 1 columns and K rows,
then model 'Mo' is a mixture of 'C' Gaussians specified by row 'Mo'
of 'scale.weights'.
}
\examples{
#2D example
line.models(X = linemodels.ex1[,c("beta1","beta2")], SE = linemodels.ex1[,c("se1","se2")],
           scales = c(0.2, 0.2, 0.2),
           slopes = c(0, 0.5, 1),
           cors = c(0.995, 0.995, 0.995))
#3D example
line.models(X = linemodels.ex1[,c("beta1","beta2","beta3")],
           SE = linemodels.ex1[,c("se1","se2","se3")],
           scales = c(0.2, 0.2, 0.2),
           slopes = matrix(c(0,0, 0.5,0, 1,1), byrow = TRUE, ncol = 2),
           cors = c(0.995, 0.995, 0.995))
}
