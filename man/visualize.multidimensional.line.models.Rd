% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linemodels.R
\name{visualize.multidimensional.line.models}
\alias{visualize.multidimensional.line.models}
\title{Visualize multidimensional linemodels}
\usage{
visualize.multidimensional.line.models(
  scales,
  slopes,
  cors,
  plot.pairs = NULL,
  var.names = NULL,
  model.names = NULL,
  model.cols = NULL,
  legend.position = "bottomright",
  X = NULL,
  SE = NULL,
  SE.mult = 1.96,
  model.prob = NULL,
  model.thresh = 0.8,
  undetermined.col = "gray",
  pch = 1,
  cex = 1,
  xlim = NULL,
  ylim = NULL,
  cex.lab = 1,
  cex.axis = 1,
  line.lty = 1,
  region.lty = 2,
  line.lwd = 1,
  region.lwd = 1,
  plot.grid = TRUE,
  emphasize.axes = TRUE
)
}
\arguments{
\item{scales}{vector of standard deviations of larger effect of each model}

\item{slopes}{matrix of slopes of the models (one row per model)}

\item{cors}{vector of correlation parameters of each model}

\item{plot.pairs}{matrix of pairs of effect variables to be plotted (2 cols and 1 row per pair)}

\item{var.names}{names used as axis labels for the effect variables}

\item{model.names}{vector of names of each linemodel}

\item{model.cols}{vector of colors for each linemodel}

\item{legend.position}{position of legend, value NULL removes legend}

\item{X}{points to be plotted, 1 col per effect variable, if NULL no points plotted}

\item{SE}{standard errors to be plotted around points, 1 col per effect, if NULL no SEs plotted}

\item{SE.mult}{multiplier of SE to be plotted around each point, default 1.96 i.e. 95\% conf. intervals}

\item{model.prob}{matrix of probability of each observation (row) in each model (col), can be NULL}

\item{model.thresh}{threshold probability to color a point by model.cols, default = 0.8}

\item{undetermined.col}{color for points where no model.prob is > model.thresh}

\item{xlim, ylim, cex, pch, cex.lab, cex.axis}{standard plotting parameters}

\item{line.lty}{plotting type for lines (1 = solid, 2 = dashed ...)}

\item{region.lty}{plotting type for 95\% regions (1 = solid, 2 = dashed ...)}

\item{line.lwd}{width of lines}

\item{region.lwd}{width of line for 95\% regions}

\item{plot.grid}{if TRUE, plots a grid}

\item{emphasize.axes}{if TRUE, coordinate axes are marked with a black line}
}
\value{
none
}
\description{
Plots lines and 95\% highest probability regions of the multi-dimensional models
specified by scales, slopes and correlations.
See help of 'line.models( )' for details about specifying the models.
}
\examples{
par(mfrow = c(1,2))
visualize.multidimensional.line.models(c(0.1, 0.2), matrix(c(1,0, 0,1), ncol = 2), c(0.99, 0.95),
plot.pairs = matrix(c(1,2, 1,3),ncol = 2, byrow = TRUE))
}
