#' linemodels example data set 1
#'
#' 100 genetic variants with associations to three outcomes.
#' Variants 1..20 from Model1, 21..60 from Model2 and 61..100 from Model3.
#' All models have scale = 0.2 and cor = 0.995.
#' Slopes are: Model1 = (0,0), Model2 = (0.5,0.2), Model3 = (1,1).
#'
#' @format
#' \describe{
#'   \item{beta1, beta2, beta3}{effects on outcomes 1, 2, 3}
#'   \item{se1, se2, se3}{standard errors of beta estimates}
#'   \item{maf}{minor allele frequency of variant}
#' }
"linemodels.ex1"

#' linemodels example data set 2
#'
#' 101 genetic variants with associations to three outcomes.
#' Variants 1..50 from Model1, 51..100 from Model2 and 101 is an outlier.
#' All models have scale = 0.1 and cor = 0.999.
#' Slopes are: Model1 = (0, 2), Model2 = (3, 6).
#'
#' @format
#' \describe{
#'   \item{beta1, beta2, beta3}{effects on outcomes 1, 2, 3}
#'   \item{se1, se2, se3}{standard errors of beta estimates}
#'   \item{annotation}{binary annotation that can be used as a feature}
#' }
"linemodels.ex2"

