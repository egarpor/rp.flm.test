

#' @title Utilities for the simulation of functional linear models
#'
#' @description The functions \code{r.cfs.2003}, \code{r.hh.2006} and \code{r.bridge} sample the functional covariate and construct functional coefficients for their use in functional linear models:
#' \itemize{
#'   \item{\code{r.cfs.2003} implements examples (a) and (b) in Section 5 of Cardot et al. (2003).}
#'   \item{\code{r.hh.2006} gives models (i), (ii) and (iii) in Section 5 of Hall and Hosseini-Nasab (2006).}
#'   \item{\code{r.bridge} samples a Brownian motion and creates a functional coefficient constructed from the eigenfunctions \eqn{\sqrt(2) * \sin(k t \pi)}{\sqrt2 * sin(k*t*\pi)}.}
#' }
#' @param n number of functions to sample.
#' @inheritParams r.mod
#' @param b coefficients of the functional coefficient in the theoretical basis of principal components of the Brownian motion.
#' @param type either example \code{"a"} or \code{"b"} from Cardot et al. (2003).
#' @param imod either \code{1}, \code{2} or \code{3} for denoting models (i), (ii) and (iii) in Hall and Hosseini-Nasab (2006).
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{X.fdata}}{the sample of functional data, an \code{\link[fda.usc]{fdata}} object of length \code{n}.}
#'   \item{\code{beta.fdata}}{the functional coefficient, an \code{\link[fda.usc]{fdata}} object.}
#' }
#' @examples
#' # Cardot et al. (2003)
#' plot(r.cfs.2003(n = 100, type = "a")$X.fdata)
#' plot(r.cfs.2003(n = 100, type = "b")$X.fdata)
#'
#' # Hall and Hosseini-Nasab (2006)
#' plot(r.hh.2006(n = 100, imod = 1)$X.fdata)
#' plot(r.hh.2006(n = 100, imod = 2)$X.fdata)
#' plot(r.hh.2006(n = 100, imod = 3)$X.fdata)
#'
#' # Sample bridge
#' plot(r.bridge(n = 100)$X.fdata)
#' @author Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}) and Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
#' @references
#' Cardot, H., Ferraty, F., Sarda, P. (2003) Spline estimators for the functional linear model. Statistica Sinica, 13(3), 571--592. \url{http://www3.stat.sinica.edu.tw/statistica/oldpdf/a13n31.pdf}
#'
#' Hall, P. and Hosseini-Nasab, M. (2006) On properties of functional principal components analysis. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 68(1), 109--126. \url{http://dx.doi.org/10.1111/j.1467-9868.2005.00535.x}
#' @export
r.cfs.2003 <- function(n, t = seq(0, 1, len = 201), b = c(2, 4, 5) / sqrt(2),
                       type = c("a", "b")[1]) {

  # X.fdata
  X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")

  if (type == "a") {

    # beta = b1 * v1 + b2 * v2 + b3 * v3
    v <- fda.usc::fdata(mdata = sqrt(2) * sin(t * pi * 0.5), argvals = t)
    for (k in 2:length(b)) {

      v <- c(v, fda.usc::fdata(mdata = sqrt(2) * sin(t * pi * (k - 0.5)),
                               argvals = t))

    }
    beta.fdata <- fda.usc::fdata(mdata = matrix(b, nrow = 1) %*% v$data, argvals = t)

  } else if (type == "b") {

    beta.fdata <- fda.usc::fdata(mdata = log(15 * t^2 + 10) + cos(4 * pi * t),
                                 argvals = t)

  } else {

    stop("Wrong type")

  }

  return(list("X.fdata" = X.fdata, "beta.fdata" = beta.fdata))

}


#' @rdname r.cfs.2003
#' @export
r.hh.2006 <- function(n, t = seq(0, 1, len = 201), imod = 1) {

  # Basis
  v <- fda.usc::fdata(mdata = sqrt(2) * cos(pi * t), argvals = t)
  for (k in 2:20) {

    v <- c(v, fda.usc::fdata(mdata = sqrt(2) * cos(k * pi * t), argvals = t))

  }

  # beta
  j <- 1:20
  beta.fdata <- fda.usc::fdata(mdata = (2^(3/2) * (-1)^j * j^(-2)) %*% v$data,
                               argvals = t)

  # X.fdata
  sdcoef <- j^(-imod)
  coefsim <- matrix(NA, ncol = 20, nrow = n)
  for (i in j) {

    coefsim[, i] <- rnorm(n, mean = 0, sd = sdcoef[i])

  }
  X.fdata <- fda.usc::fdata(mdata = coefsim %*% v$data, argvals = t)

  return(list("X.fdata" = X.fdata, "beta.fdata" = beta.fdata))

}


#' @rdname r.cfs.2003
#' @export
r.bridge <- function(n, t = seq(0, 1, len = 201), b = c(2, 4, 5) / sqrt(2)) {

  # X.fdata
  X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")
  lt <- length(t)
  for (i in 1:n) {

    X.fdata$data[i, ] <- X.fdata$data[i, ] - t * X.fdata$data[i, lt]

  }

  # beta = b1 * v1 + b2 * v2 + b3 * v3
  v <- fda.usc::fdata(mdata = sqrt(2) * sin(t * pi), argvals = t)
  for (k in 2:length(b)) {

    v <- c(v, fda.usc::fdata(mdata = sqrt(2) * sin(t * pi * k), argvals = t))

  }
  beta.fdata <- fda.usc::fdata(mdata = matrix(b, nrow = 1) %*% v$data, argvals = t)

  return(list("X.fdata" = X.fdata, "beta.fdata" = beta.fdata))

}


#' @title Geometric Brownian motion
#'
#' @description Sampling of paths of the geometric Brownian motion.
#'
#' @inheritParams r.cfs.2003
#' @param mu,sigma mean and diffusion of the underlying Brownian motion.
#' @param s0 a number or a vector of length \code{n} giving the initial value(s) of the geometric Brownian motion.
#' @return Functional sample, an \code{\link[fda.usc]{fdata}} object of length \code{n}.
#' @examples
#' plot(r.gbm(n = 10, s0 = 1))
#' plot(r.gbm(n = 10, mu = -1, s0 = rnorm(10, mean = 10)))
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
#' @export
r.gbm <- function(n, t = seq(0, 1, len = 201), mu = 0, sigma = 1, s0 = 1) {

  # Time-varying covariances
  St <- sigma^2 * outer(t, t, function(s, t) pmin(s, t))

  # Sample N((mu - sigma^2 / 2) * t, St)
  mdata <- mvtnorm::rmvnorm(n = n, mean = (mu - sigma^2 / 2) * t, sigma = St)
  mdata <- s0 * exp(mdata)

  # As fdata object
  return(fda.usc::fdata(mdata = mdata, argvals = t))

}


#' @title Ornstein-Uhlenbeck process
#'
#' @description Sampling of paths of the Ornstein-Uhlenbeck process.
#'
#' @inheritParams r.cfs.2003
#' @param alpha strength of the drift.
#' @param mu mean of the process.
#' @param sigma diffusion coefficient.
#' @param x0 a number or a vector of length \code{n} giving the initial value(s) of the Ornstein-Uhlenbeck process. By default, \code{n} points are sampled from the stationary distribution.
#' @return Functional sample, an \code{\link[fda.usc]{fdata}} object of length \code{n}.
#' @examples
#' plot(r.ou(n = 100))
#' plot(r.ou(n = 100, alpha = 2, sigma = 4, x0 = 1:100))
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}).
#' @export
r.ou <- function(n, t = seq(0, 1, len = 201), mu = 0, alpha = 1, sigma = 1,
                 x0 = rnorm(n, mean = mu, sd = sigma/sqrt(2 * alpha))){

  # Time-varying covariances
  St <- sigma^2/(2 * alpha) * outer(t, t, function(s, t) {
    exp(alpha * (2 * pmin(s, t) - (s + t))) - exp(-alpha * (s + t))
    })

  # Sample N(0, St) and add time-varying mean
  mdata <- mvtnorm::rmvnorm(n = n, mean = rep(0, length(t)), sigma = St)
  mdata <- mdata + outer(x0, t, function(x0, t) mu + (x0 - mu) * exp(-alpha * t))

  # As fdata object
  return(fda.usc::fdata(mdata = mdata, argvals = t))

}
