% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linemodels.R
\name{beta.cor.case.control}
\alias{beta.cor.case.control}
\title{Compute correlation in likelihood based on case and control overlaps between 2 studies
Formula is from Bhattacharjee et al. (2012)
Subset-Based Approach Improves Power and Interpretation for the Combined Analysis
of Genetic Association Studies of Heterogeneous Traits,
The American Journal of Human Genetics, 90(5): 821-835.}
\usage{
beta.cor.case.control(s1, r1, s2, r2, s1s2 = 0, r1r2 = 0, s1r2 = 0, r1s2 = 0)
}
\arguments{
\item{s1}{cases of study 1}

\item{r1}{controls of study 1}

\item{s2}{cases of study 2}

\item{r2}{controls of study 2}

\item{s1s2}{cases in 1 and cases in 2}

\item{r1r2}{controls in 1 and controls in 2}

\item{s1r2}{cases in 1 and controls in 2}

\item{r1s2}{controls in 1 and cases in 2}
}
\value{
correlation value
}
\description{
Compute correlation in likelihood based on case and control overlaps between 2 studies
Formula is from Bhattacharjee et al. (2012)
Subset-Based Approach Improves Power and Interpretation for the Combined Analysis
of Genetic Association Studies of Heterogeneous Traits,
The American Journal of Human Genetics, 90(5): 821-835.
}
\examples{
beta.cor.case.control(s1 = 1000, r1 = 1000, s2 = 1000, r2 = 2000, s1s2 = 0, r1r2 = 1000)
}
