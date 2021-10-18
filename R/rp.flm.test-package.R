

#' @title rp.flm.test - Goodness-of-fit tests for the functional linear model based on randomly projected empirical processes
#'
#' @description Package companion for the paper "Goodness-of-fit tests for the functional linear model based on randomly projected empirical processes" (Cuesta-Albertos et al., 2017). The package implements goodness-of-fit tests for the functional linear model with scalar response.
#'
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}) and Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}).
#' @docType package
#' @name rp.flm.test-package
#' @import fda.usc
#' @importFrom grDevices gray
#' @importFrom graphics legend lines par plot title
#' @importFrom stats density rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @useDynLib rp.flm.test, .registration = TRUE
#' @references
#' Cuesta-Albertos, J.A., García-Portugués, E., Febrero-Bande, M. and González-Manteiga, W. (2019). Goodness-of-fit tests for the functional linear model based on randomly projected empirical processes. \emph{Annals of Statistics}, 47(1):439-467. \url{https://doi.org/10.1214/18-AOS1693}
NULL
