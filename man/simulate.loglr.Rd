% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linemodels.R
\name{simulate.loglr}
\alias{simulate.loglr}
\title{Compute observed log of likelihood ratio (LR) and its null distribution for the given scenario.}
\usage{
\method{simulate}{loglr}(
  X,
  SE,
  n.sims = 100,
  par.include.null = matrix(TRUE, nrow = length(init.slopes.null), ncol = 3),
  force.same.scales.null = FALSE,
  init.scales.null,
  init.slopes.null,
  init.cors.null,
  model.priors.null = rep(1, length(init.slopes.null)),
  par.include.alt = matrix(TRUE, nrow = length(init.slopes.alt), ncol = 3),
  force.same.scales.alt = FALSE,
  init.scales.alt,
  init.slopes.alt,
  init.cors.alt,
  model.priors.alt = rep(1, length(init.slopes.alt)),
  r.lkhood = 0,
  tol.loglk = 0.001,
  tol.par = 0,
  print.steps = 0
)
}
\arguments{
\item{X}{matrix of effect sizes with two columns (one col per trait)}

\item{SE}{matrix of standard errors with two columns (one col per trait)}

\item{n.sims}{number of simulated data sets for estimating distribution of logLR}

\item{par.include.null}{matrix of TRUE/FALSE values TRUE indicating which parameters
are optimized in the null model. One row per line model and 3 columns corresponding
to 'scales', 'slopes' and 'cors' of the line models undel NULL}

\item{force.same.scales.null}{logical; If TRUE, then forces all scale parameters equal in the NULL model}

\item{init.scales.null}{vector of initial standard deviations of larger effect of each line model
under the NULL model}

\item{init.slopes.null}{vector of initial slopes of each line model
under the NULL model}

\item{init.cors.null}{vector of initial correlation parameters of each line model
under the NULL model}

\item{model.priors.null}{vector of initial prior probabilities of line models up to a normalization constant
under the NULL model}

\item{par.include.alt}{matrix of TRUE/FALSE values TRUE indicating which parameters
are optimized in the ALT model. One row per line model and 3 columns corresponding
to 'scales', 'slopes' and 'cors' of the line models under ALT}

\item{force.same.scales.alt}{logical; If TRUE, then forces all scale parameters equal in the ALT model}

\item{init.scales.alt}{vector of initial standard deviations of larger effect of each line model
under the ALT model}

\item{init.slopes.alt}{vector of initial slopes of each line model
under the ALT model}

\item{init.cors.alt}{vector of initial correlation parameters of each line model
under the ALT model}

\item{model.priors.alt}{vector of initial prior probabilities of line models up to a normalization constant
under the ALT model}

\item{r.lkhood}{correlation between estimators of the two effects}

\item{tol.loglk}{tolerance for convergence of optimization in adjacent log-likelihoods}

\item{tol.par}{tolerance for convergence of optimization in maximum of absolute values of relative differences
across parameters between adjacent iterations. Can be set negative to determine
the convergence solely by log-likelihood.}

\item{print.steps}{numeric;
If 0, no optimization results are printed on screen.
If 1, prints optimization results at start and end of each optimization call.
If >1 prints status after every optimization iteration.}
}
\value{
list with two components. 'obs.loglr' is the log likelihood ratio for observed data and
'null.loglr' is vector of 'n.sims' simulated loglr values under the null model.
}
\description{
See help of 'line.models( )' for details about specifying the models.
}
\details{
Compares two models (ALTernative and NULL) against each other by maximized likelihood ratio.
Both models are specified by input parameters that are passed to line.models.optimize().
The optimized log likelihoods are computed and logLR of ALT vs. NULL is recorded.
Then a null distribution for logLR is estimated using user-specified number of simulated data sets.
Each data set is simulated to mimick the observed data set (w.r.t sample size, effects, SEs)
and follows the optimized NULL model w.r.t line model parameters and model probabilities.
Then both models (ALT and NULL) are optimized using the simulated data and logLRs are recorded.
The function returns both the observed logLR and the simulated null distribution of logLR,
which can be used for deriving an empirical P-value.
}
\examples{
simulate.loglr(
   X = linemodels.ex1[,1:2], SE = linemodels.ex1[,3:4],
   n.sims = 2,
   par.include.null = rbind(c(F,F,F), c(F,F,F)),
   init.scales.null = c(0.25, 0.25),
   init.slopes.null = c(0, 1),
   init.cors.null = c(0.995, 0.995),
   par.include.alt = rbind(c(F,F,F), c(F,F,F), c(F,T,F)),
   init.scales.alt = c(0.25, 0.25, 0.25),
   init.slopes.alt = c(0, 1, 0.5),
   init.cors.alt = c(0.995, 0.995, 0.995))
}
