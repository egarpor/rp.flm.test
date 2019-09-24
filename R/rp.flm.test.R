

#' @title Statistics for testing the functional linear model using random projections
#'
#' @description Computes the Cramer-von Mises (CvM) and Kolmogorov-Smirnov (KS) statistics on the projected process
#' \deqn{T_{n,h}(u)=\frac{1}{n}\sum_{i=1}^n (Y_i-\langle X_i,\hat \beta\rangle)1_{\{\langle X_i, h\rangle\leq u\}},}{T_{n, h}(u)=1/n\sum_{i = 1}^n (Y_i - <X_i, \hat \beta>)1_{<X_i, h> \le u},}
#' designed to test the goodness-of-fit of a functional linear model with scalar response.
#'
#' @param proj.X matrix of size \code{c(n, n.proj)} containing, for each column, the projections of the functional data \eqn{X_1,\ldots,X_n} into a random direction \eqn{h}. Not required if \code{proj.X.ord} is provided.
#' @param residuals the residuals of the fitted functional linear model, \eqn{Y_i-\langle X_i,\hat \beta\rangle}{Y_i - <X_i, \hat \beta, Y_i>}. Either a vector of length \code{n} (same residuals for all projections) or a matrix of size \code{c(n.proj, n)} (each projection has an associated set residuals).
#' @param proj.X.ord matrix containing the row permutations of \code{proj.X} which rearranges them increasingly, for each column. So, for example \code{proj.X[proj.X.ord[, 1], 1]} equals \code{sort(proj.X[, 1])}. If not provided, it is computed internally.
#' @param F.code whether to use faster \code{FORTRAN} code or \code{R} code.
#' @return A list containing:
#' \describe{
#'   \item{\code{statistic}}{a matrix of size \code{c(n.proj, 2)} with the CvM (first column) and KS (second) statistics, for the \code{n.proj} different projections.}
#'   \item{\code{proj.X.ord}}{the computed row permutations of \code{proj.X}, useful for recycling in subsequent calls to \code{rp.flm.statistic} with the same projections but different residuals.}
#' }
#' @details \code{NA}'s are not allowed neither in the functional covariate nor in the scalar response.
#' @examples
#' # Simulated example
#' set.seed(345678)
#' t <- seq(0, 1, l = 101)
#' n <- 100
#' X <- r.ou(n = n, t = t)
#' beta0 <- fdata(mdata = cos(2 * pi * t) - (t - 0.5)^2, argvals = t,
#'                rangeval = c(0,1))
#' Y <- inprod.fdata(X, beta0) + rnorm(n, sd = 0.1)
#'
#' # Linear model
#' mod <- fregre.pc(fdataobj = X, y = Y, l = 1:3)
#'
#' # Projections
#' proj.X1 <- inprod.fdata(X, r.ou(n = 1, t = t))
#' proj.X2 <- inprod.fdata(X, r.ou(n = 1, t = t))
#' proj.X12 <- cbind(proj.X1, proj.X2)
#'
#' # Statistics
#' t1 <- rp.flm.statistic(proj.X = proj.X1, residuals = mod$residuals)
#' t2 <- rp.flm.statistic(proj.X = proj.X2, residuals = mod$residuals)
#' t12 <- rp.flm.statistic(proj.X = proj.X12, residuals = mod$residuals)
#' t1$statistic
#' t2$statistic
#' t12$statistic
#'
#' # Recycling proj.X.ord
#' rp.flm.statistic(proj.X.ord = t1$proj.X.ord, residuals = mod$residuals)$statistic
#' t1$statistic
#'
#' # Sort in the columns
#' cbind(proj.X12[t12$proj.X.ord[, 1], 1], proj.X12[t12$proj.X.ord[, 2], 2]) -
#' apply(proj.X12, 2, sort)
#'
#' # FORTRAN and R code
#' rp.flm.statistic(proj.X = proj.X1, residuals = mod$residuals)$statistic -
#' rp.flm.statistic(proj.X = proj.X1, residuals = mod$residuals,
#'                  F.code = FALSE)$statistic
#'
#' # Matrix and vector residuals
#' rp.flm.statistic(proj.X = proj.X12, residuals = mod$residuals)$statistic
#' rp.flm.statistic(proj.X = proj.X12,
#'                  residuals = rbind(mod$residuals, mod$residuals * 2))$statistic
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}) and Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}).
#' @references
#' Cuesta-Albertos, J.A., García-Portugués, E., Febrero-Bande, M. and González-Manteiga, W. (2019). Goodness-of-fit tests for the functional linear model based on randomly projected empirical processes. \emph{Annals of Statistics}, 47(1):439-467. \url{https://doi.org/10.1214/18-AOS1693}
#' @useDynLib rp.flm.test, .registration = TRUE
#' @export
rp.flm.statistic <- function(proj.X, residuals, proj.X.ord = NULL, F.code = TRUE) {

  # Number of projections
  n.proj <- ifelse(is.null(proj.X.ord), ncol(proj.X), ncol(proj.X.ord))

  # Residuals as a matrix
  if (!is.matrix(residuals)) {

    residuals <- matrix(residuals, nrow = n.proj, ncol = length(residuals),
                        byrow = TRUE)

  }
  n <- ncol(residuals)
  if (nrow(residuals) != n.proj) {

    stop("The number of rows in residuals must be the number of projections")

  }

  # Matrix of statistics (columns) projected in n.proj projections (rows)
  rp.stat <- matrix(0, nrow = n.proj, ncol = 2)

  # Order projections if not provided
  if (is.null(proj.X.ord)) {

    proj.X.ord <- apply(proj.X, 2, order)

  }

  # Compute statistics
  if (F.code) {

    # Statistic
    rp.stat <- .Fortran("rp_stat", proj_X_ord = proj.X.ord, residuals = residuals,
                        n_proj = n.proj, n = n, rp_stat_proj = rp.stat,
                        PACKAGE = "rp.flm.test")$rp_stat_proj

  } else {

    # R implementation
    for (i in 1:n.proj) {

      # Empirical process
      y <- cumsum(residuals[i, proj.X.ord[, i]])

      # Statistics (CvM and KS, rows)
      CvM <- sum(y^2)
      KS <- max(abs(y))
      rp.stat[i, ] <- c(CvM, KS)

    }

    # Standardize
    rp.stat[, 1] <- rp.stat[, 1] / (n^2)
    rp.stat[, 2] <- rp.stat[, 2] / sqrt(n)

  }

  # Return both statistics
  colnames(rp.stat) <- c("CvM", "KS")
  return(list(statistic = rp.stat, proj.X.ord = proj.X.ord))

}


#' @title Goodness-of fit test for the functional linear model using random projections
#'
#' @description Tests the composite null hypothesis of a Functional Linear Model with scalar response (FLM),
#' \deqn{H_0:\,Y=\langle X,\beta\rangle+\epsilon\quad\mathrm{vs}\quad H_1:\,Y\neq\langle X,\beta\rangle+\epsilon.}{H_0: Y = <X, \beta> + \epsilon vs H_1: Y != <X, \beta> + \epsilon.}
#' If \eqn{\beta=\beta_0}{\beta=\beta_0} is provided, then the simple hypothesis \eqn{H_0:\,Y=\langle X,\beta_0\rangle+\epsilon}{H_0: Y = <X, \beta_0> + \epsilon} is tested. The way of testing the null hypothesis is via a norm (Cramer-von Mises or Kolmogorov-Smirnov) in the empirical process indexed by the projections.
#'
#' @param X.fdata functional observations in the class \code{\link[fda.usc]{fdata}}.
#' @param Y scalar responses for the FLM. Must be a vector with the same number of elements as functions are in \code{X.fdata}.
#' @param beta0.fdata functional parameter for the simple null hypothesis, in the \code{\link[fda.usc]{fdata}} class. The \code{argvals} and \code{rangeval} arguments of \code{beta0.fdata} must be the same of \code{X.fdata}. If \code{beta0.fdata=NULL} (default), the function will test for the composite null hypothesis.
#' @param B number of bootstrap replicates to calibrate the distribution of the test statistic.
#' @param n.proj vector with the number of projections to consider. The maximum number must be smaller than \code{100} (which is more than reasonable for the use of the test).
#' @param est.method estimation method for \eqn{\beta}{\beta}, only used in the composite case. There are three methods:
#' \describe{
#'   \item{\code{"pc"}}{if \code{p} is given, then \eqn{\beta}{\beta} is estimated by \code{\link[fda.usc]{fregre.pc}}. Otherwise, \code{p} is chosen using \code{\link[fda.usc]{fregre.pc.cv}} and the \code{p.criterion} criterion.}
#'   \item{\code{"pls"}}{if \code{p} is given, \eqn{\beta}{\beta} is estimated by \code{\link[fda.usc]{fregre.pls}}. Otherwise, \code{p} is chosen using \code{\link[fda.usc]{fregre.pls.cv}} and the \code{p.criterion} criterion.}
#'   \item{\code{"basis"}}{if \code{p} is given, \eqn{\beta}{\beta} is estimated by \code{\link[fda.usc]{fregre.basis}}. Otherwise, \code{p} is chosen using \code{\link[fda.usc]{fregre.basis.cv}} and the \code{p.criterion} criterion. Both in \code{\link[fda.usc]{fregre.basis}} and \code{\link[fda.usc]{fregre.basis.cv}}, the same basis for \code{basis.x} and \code{basis.b} is considered.}
#' }
#' @param p number of elements for the basis representation of \code{beta0.fdata} and \code{X.fdata} with the \code{est.method} (only composite hypothesis). If not supplied, it is estimated from the data.
#' @param p.criterion for \code{est.method} equal to \code{"pc"} or \code{"pls"}, either \code{"SIC"}, \code{"SICc"} or one of the criteria described in \code{\link[fda.usc]{fregre.pc.cv}}. For \code{"basis"} a value for \code{type.CV} in \code{\link[fda.usc]{fregre.basis.cv}} such as \code{GCV.S}.
#' @param pmax maximum size of the basis expansion to consider in when using \code{p.criterion}.
#' @param type.basis type of basis if \code{est.method = "basis"}.
#' @param verbose whether to show or not information about the testing progress.
#' @param projs a \code{\link[fda.usc]{fdata}} object containing the random directions employed to project \code{X.fdata}. If numeric, the convenient value for \code{ncomp} in \code{\link{rdir.pc}}.
#' @param ... further arguments passed to \code{\link[fda]{create.basis}} (not \code{rangeval} that is taken as the \code{rangeval} of \code{X.fdata}).
#' @return An object with class \code{"htest"} whose underlying structure is a list containing the following components:
#' \describe{
#'   \item{\code{p.values.fdr}}{a matrix of size \code{c(n.proj, 2)}, containing in each row the FDR p-values of the CvM and KS tests up to that projection.}
#'   \item{\code{proj.statistics}}{a matrix of size \code{c(max(n.proj), 2)} with the value of the test statistic on each projection.}
#'   \item{\code{boot.proj.statistics}}{an array of size \code{c(max(n.proj), 2, B)} with the values of the bootstrap test statistics for each projection.}
#'   \item{\code{proj.p.values}}{a matrix of size \code{c(max(n.proj), 2)}}
#'   \item{\code{method}}{information about the test performed and the kind of estimation performed.}
#'   \item{\code{B}}{number of bootstrap replicates used.}
#'   \item{\code{n.proj}}{number of projections specified}
#'   \item{\code{projs}}{random directions employed to project \code{X.fdata}.}
#'   \item{\code{type.basis}}{type of basis for \code{est.method = "basis"}.}
#'   \item{\code{beta.est}}{estimated functional parameter \eqn{\hat \beta}{\hat \beta} in the composite hypothesis. For the simple hypothesis, \code{beta0.fdata}.}
#'   \item{\code{p}}{number of basis elements considered for estimation of \eqn{\beta}{\beta}.}
#'   \item{\code{p.criterion}}{criterion employed for selecting \code{p}.}
#'   \item{\code{data.name}}{the character string "Y = <X, b> + e"}
#' }
#' @details
#' No NA's are allowed neither in the functional covariate nor in the scalar response.
#' @examples
#' # Simulated example
#'
#' set.seed(345678)
#' t <- seq(0, 1, l = 101)
#' n <- 100
#' X <- r.ou(n = n, t = t, alpha = 2, sigma = 0.5)
#' beta0 <- fdata(mdata = cos(2 * pi * t) - (t - 0.5)^2, argvals = t,
#'                rangeval = c(0,1))
#' Y <- inprod.fdata(X, beta0) + rnorm(n, sd = 0.1)
#'
#' \dontrun{
#' # Test all cases
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pc")
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pls")
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "basis",
#'             p.criterion = fda.usc::GCV.S)
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pc", p = 5)
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pls", p = 5)
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "basis", p = 5)
#' rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0)
#' }
#' 
#' # Composite hypothesis: do not reject FLM
#' rp.test <- rp.flm.test(X.fdata = X, Y = Y, est.method = "pc")
#' rp.test$p.values.fdr
#' pcvm.test <- flm.test(X.fdata = X, Y = Y, est.method = "pc", B = 1e3,
#'                       plot.it = FALSE)
#' pcvm.test
#'
#' # Estimation of beta
#' par(mfrow = c(1, 3))
#' plot(X, main = "X")
#' plot(beta0, main = "beta")
#' lines(rp.test$beta.est, col = 2)
#' lines(pcvm.test$beta.est, col = 3)
#' plot(density(Y), main = "Density of Y", xlab = "Y", ylab = "Density")
#' rug(Y)
#'
#' \dontrun{
#' # Simple hypothesis: do not reject beta = beta0
#' rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0)$p.values.fdr
#' flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0, B = 1e3, plot.it = FALSE)
#'
#' # Simple hypothesis: reject beta = beta0^2
#' rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0^2)$p.values.fdr
#' flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0^2, B = 1e3, plot.it = FALSE)
#' }
#' 
#' \dontrun{
#' # Increasing n.proj
#' rp.flm.test(X.fdata = X, Y = Y, n.proj = 3, verbose = FALSE)$p.values.fdr
#' rp.flm.test(X.fdata = X, Y = Y, n.proj = 5, verbose = FALSE)$p.values.fdr
#' 
#' # Increasing B
#' rp.flm.test(X.fdata = X, Y = Y, B = 1e3, verbose = FALSE)$p.values.fdr
#' rp.flm.test(X.fdata = X, Y = Y, B = 5e3, verbose = FALSE)$p.values.fdr
#' }
#' 
#' \dontrun{
#' # Tecator dataset
#'
#' # Load data
#' data(tecator)
#' absorp <- tecator$absorp.fdata
#' ind <- 1:215 # sometimes ind <- 1:129
#' x <- absorp[ind, ]
#' y <- tecator$y$Fat[ind]
#'
#' # Composite hypothesis
#' rp.tecat <- rp.flm.test(X.fdata = x, Y = y, verbose = FALSE, B = 1e4)
#' rp.tecat$p.values.fdr
#' 
#' # Compare with flm.test
#' pcvm.tecat <- flm.test(X.fdata = x, Y = y, verbose = FALSE, p = rp.tecat$p, 
#'                        B = 1e4, plot.it = FALSE)
#' pcvm.tecat
#'
#' # Simple hypothesis
#' zero <- fdata(mdata = rep(0, length(x$argvals)), argvals = x$argvals,
#'               rangeval = x$rangeval)
#' rp.flm.test(X.fdata = x, Y = y, beta0.fdata = zero, verbose = FALSE, B = 1e4)
#' flm.test(X.fdata = x, Y = y, beta0.fdata = zero, B = 1e3)
#'
#' # With derivatives
#' rp.tecat.1 <- rp.flm.test(X.fdata = fdata.deriv(x, 1), Y = y, verbose = FALSE,
#'                           B = 1e4)
#' rp.tecat.1$p.values.fdr
#' rp.tecat.2 <- rp.flm.test(X.fdata = fdata.deriv(x, 2), Y = y, verbose = FALSE,
#'                           B = 1e4)
#' rp.tecat.2$p.values.fdr
#' }
#' 
#' \dontrun{
#' # AEMET dataset
#'
#' # Load data
#' data(aemet)
#' wind.speed <- apply(aemet$wind.speed$data, 1, mean)
#' temp <- aemet$temp
#' 
#' # Remove the 5% of the curves with less depth (i.e. 4 curves)
#' par(mfrow = c(1, 1))
#' res.FM <- depth.FM(temp, draw = TRUE)
#' qu <- quantile(res.FM$dep, prob = 0.05)
#' l <- which(res.FM$dep <= qu)
#' lines(aemet$temp[l], col = 3)
#' 
#' # Data without outliers
#' wind.speed <- wind.speed[-l]
#' temp <- temp[-l]
#' 
#' # Composite hypothesis
#' rp.aemet <- rp.flm.test(X.fdata = temp, Y = wind.speed, verbose = FALSE, 
#'                         B = 1e4)
#' rp.aemet$p.values.fdr
#' 
#' # Compare with flm.test
#' pcvm.aemet <- flm.test(X.fdata = temp, Y = wind.speed, B = 1e4,
#'                        p = rp.aemet$p, verbose = FALSE, plot.it = FALSE)
#' pcvm.aemet                        
#' 
#' # Simple hypothesis
#' zero <- fdata(mdata = rep(0, length(temp$argvals)), argvals = temp$argvals,
#'               rangeval = temp$rangeval)
#' rp.flm.test(X.fdata = temp, Y = wind.speed, beta0.fdata = zero, 
#'             verbose = FALSE, B = 1e4)
#' flm.test(X.fdata = temp, Y = wind.speed, beta0.fdata = zero, B = 1e4,
#'          plot.it = FALSE, verbose = FALSE)
#' }
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}) and Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}).
#' @references
#' Cuesta-Albertos, J.A., García-Portugués, E., Febrero-Bande, M. and González-Manteiga, W. (2019). Goodness-of-fit tests for the functional linear model based on randomly projected empirical processes. \emph{Annals of Statistics}, 47(1):439-467. \url{https://doi.org/10.1214/18-AOS1693}
#'
#' García-Portugués, E., González-Manteiga, W. and Febrero-Bande, M. (2014). A goodness-of-fit test for the functional linear model with scalar response. Journal of Computational and Graphical Statistics, 23(3), 761--778. \url{http://dx.doi.org/10.1080/10618600.2013.812519}
#' @export
rp.flm.test <- function(X.fdata, Y, beta0.fdata = NULL, B = 1000, n.proj = 3,
                        est.method = "pc", p = NULL, p.criterion = "SICc",
                        pmax = 10, type.basis = "bspline", projs = 0.95,
                        verbose = TRUE, ...) {

  # Sample size
  n <- dim(X.fdata)[1]

  # p data driven flag
  p.data.driven <- is.null(p)

  # Display progress
  if (verbose) {

    cat("Computing estimation of beta... ")

  }

  # Truncate maximum basis expansion
  pmax <- min(pmax, n)

  # Check max(n.proj)
  if (max(n.proj) > 100) {
    
    stop("max(n.proj) must be smaller than 100")
    
  }
  
  ## Estimation of beta

  # Composite hypothesis: optimal estimation of beta and the basis expansion
  if (is.null(beta0.fdata)) {

    # Center the data first
    X.fdata <- fda.usc::fdata.cen(X.fdata)$Xcen
    Y <- Y - mean(Y)

    # Method
    meth <- "Random projection based test for the functional linear model using"

    # PC
    if (est.method == "pc") {

      # Optimal p by p.criterion
      if (p.data.driven) {

        # Method
        meth <- paste(meth, "optimal PC basis representation")

        # Choose the number of basis elements
        mod <- fda.usc::fregre.pc.cv(fdataobj = X.fdata, y = Y, kmax = 1:pmax,
                                     criteria = p.criterion)
        p.opt <- length(mod$pc.opt)
        ord.opt <- mod$pc.opt

        # Return the best model
        mod <- mod$fregre.pc
        pc.comp <- mod$fdata.comp

      # Fixed p
      } else {

        # Method
        meth <- paste(meth, " a representation in a PC basis of ", p, "elements")

        # Estimation of beta on the given fixed basis
        mod <- fda.usc::fregre.pc(fdataobj = X.fdata, y = Y, l = 1:p)
        pc.comp <- mod$fdata.comp
        p.opt <- p
        ord.opt <- mod$l

      }

    # PLS
    } else if (est.method == "pls") {

      # Optimal p by p.criterion
      if (p.data.driven) {

        # Method
        meth <- paste(meth, "optimal PLS basis representation")

        # Choose the number of the basis: SIC is probably the best criteria
        mod <- fda.usc::fregre.pls.cv(fdataobj = X.fdata, y = Y, kmax = pmax,
                                      criteria = p.criterion)
        p.opt <- length(mod$pls.opt)
        ord.opt <- mod$pls.opt

        # Return the best model
        mod <- mod$fregre.pls

      # Fixed p
      } else {

        # Method
        meth <- paste(meth, "a representation in a PLS basis of ", p, "elements")

        # Estimation of beta on the given fixed basis
        mod <- fda.usc::fregre.pls(fdataobj = X.fdata, y = Y, l = 1:p)
        p.opt <- p
        ord.opt <- mod$l

      }

    # Deterministic basis
    } else if (est.method == "basis") {

      # Optimal p by p.criterion
      if (p.data.driven) {

        # Method
        meth <- paste(meth, "optimal", type.basis, "basis representation")

        # Choose the number of the bspline basis with GCV.S
        if (type.basis == "fourier") {

          basis.x <- seq(1, pmax, by = 2)

        } else {

          basis.x <- 5:max(pmax, 5)

        }
        mod <- fda.usc::fregre.basis.cv(fdataobj = X.fdata, y = Y,
                                        basis.x = basis.x, basis.b = NULL,
                                        type.basis = type.basis,
                                        type.CV = p.criterion, verbose = FALSE,
                                        ...)
        p.opt <- mod$basis.x.opt$nbasis
        ord.opt <- 1:p.opt

      # Fixed p
      } else {

        # Method
        meth <- paste(meth, "a representation in a", type.basis, "basis of ",
                      p, "elements")

        # Estimation of beta on the given fixed basis
        basis.opt <- do.call(what = paste("create.", type.basis,
                                          ".basis", sep = ""),
                             args = list(rangeval = X.fdata$rangeval,
                                         nbasis = p, ...))
        mod <- fda.usc::fregre.basis(fdataobj = X.fdata, y = Y,
                                     basis.x = basis.opt, basis.b = basis.opt)
        p.opt <- p
        ord.opt <- 1:p.opt

      }

    } else {

      stop(paste("Estimation method", est.method, "not implemented."))

    }

    # Estimated beta
    beta.est <- mod$beta.est

    # Hat matrix
    H <- mod$H

    # Residuals
    e <- mod$residuals

  # Simple hypothesis
  } else {

    # Method
    meth <- "Random projection based test for the simple hypothesis in a functional linear model"

    # Do not need to estimate beta
    beta.est <- beta0.fdata
    p.opt <- NA

    # Compute the residuals
    e <- drop(Y - fda.usc::inprod.fdata(X.fdata, beta.est))

  }

  ## Computation of the statistic

  # Fix seed for projections
  if (exists(".Random.seed")) {
    old <- .Random.seed
    on.exit({.Random.seed <<- old})
  }
  set.seed(987654321)
  
  # Sample random directions
  if (verbose) {

    cat("Done.\nComputing projections... ")

  }
  if (is.numeric(projs)) {

    # Compute PCs if not done yet
    if (est.method != "pc" | !is.null(beta0.fdata)) {

      pc.comp <- fda.usc::fdata2pc(X.fdata, ncomp = min(length(X.fdata$argvals),
                                                        nrow(X.fdata)))

    }

    # Random directions
    if (length(n.proj) > 1) {

      vec.nproj <- sort(n.proj)
      n.proj <- max(n.proj)

    } else {

      vec.nproj <- 1:n.proj

    }
    projs <- rdir.pc(n = n.proj, X.fdata = X.fdata, ncomp = projs,
                     fdata2pc.obj = pc.comp, sd = 0)
    
  } else {

    n.proj <- length(projs)
    vec.nproj <- 1:n.proj

  }

  # Compute projections for the statistic and the bootstrap replicates
  proj.X <- fda.usc::inprod.fdata(X.fdata, projs) # A matrix n x n.proj

  # Statistic
  rp.stat <- rp.flm.statistic(proj.X = proj.X, residuals = e, F.code = TRUE)

  ## Bootstrap calibration
  
  # Fix different seeds adequately to ensure randomness and coherency when 
  # n.proj is increased (i.e., for a run of the test with n.proj + 1, the first
  # n.proj p.values are the same)
  primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 
              59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 
              127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 
              191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 
              257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 
              331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 
              401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 
              467, 479, 487, 491, 499, 503, 509, 521, 523, 541)[1:n.proj]
  seeds <- matrix(1:B, nrow = n.proj, ncol = B, byrow = TRUE) * primes
  
  # Golden bootstrap sampling
  gold.val <- c((1 - sqrt(5)) / 2, (1 + sqrt(5)) / 2)
  gold.prob <- c((5 + sqrt(5)) / 10, (5 - sqrt(5)) / 10)
  golden.boot <- function(i) {
    t(sapply(1:n.proj, function(k) {
      set.seed(seeds[k, i])
      e * sample(x = gold.val, size = n, prob = gold.prob, replace = TRUE)
      }))
    }
  
  # Define required objects
  rp.stat.star <- array(NA, dim = c(n.proj, 2, B))
  if (verbose) {

    cat("Done.\nBootstrap calibration...\n ")
    pb <- txtProgressBar(style = 3)

  }

  # Composite hypothesis
  if (is.null(beta0.fdata)) {

    # Calculate the matrix that gives the residuals of the linear model
    # from the observed response. This allows to resample efficiently the
    # residuals without re-estimating again the beta
    Y.to.residuals.matrix <- diag(rep(1, n)) - H

    # Bootstrap resampling
    for (i in 1:B) {

      # Generate bootstrap errors
      e.star <- golden.boot(i)

      # Residuals from the bootstrap estimated model (implicit column recycling)
      Y.star <- t(t(e.star) + drop(Y - e))
      e.hat.star <- Y.star %*% Y.to.residuals.matrix

      # Calculate the bootstrap statistics
      rp.stat.star[, , i] <- rp.flm.statistic(residuals = e.hat.star,
                                              proj.X.ord = rp.stat$proj.X.ord,
                                              F.code = TRUE)$statistic

      # Display progress
      if (verbose) {

        setTxtProgressBar(pb, i / B)

      }

    }

  # Simple hypothesis
  } else {

    # Bootstrap resampling
    for (i in 1:B) {

      # Generate bootstrap errors
      e.hat.star <- golden.boot(i)

      # Calculate the bootstrap statistics
      rp.stat.star[, , i] <- rp.flm.statistic(residuals = e.hat.star,
                                              proj.X.ord = rp.stat$proj.X.ord,
                                              F.code = TRUE)$statistic

      # Display progress
      if (verbose) {

        setTxtProgressBar(pb, i/B)

      }

    }

  }

  # Compute the p-values of the projected tests
  positiveCorrection <- FALSE
  pval <- t(sapply(1:n.proj, function(i) {

    c(sum(rp.stat$statistic[i, 1] <= rp.stat.star[i, 1, ]),
      sum(rp.stat$statistic[i, 2] <= rp.stat.star[i, 2, ]))

    }) + positiveCorrection) / (B + positiveCorrection)

  # Compute p-values depending for the vector of projections
  rp.pvalue <- t(sapply(seq_along(vec.nproj), function(k) {

    apply(pval[1:vec.nproj[k], , drop = FALSE], 2, function(x) {

      l <- length(x)
      return(min(l / (1:l) * sort(x)))

    })

  }))
  colnames(rp.pvalue) <- colnames(pval) <- c("CvM", "KS")
  rownames(rp.pvalue) <- vec.nproj

  # Return result
  if (verbose) {

    cat("\nDone.\n")

  }
  options(warn = -1)
  mean.stats <- colMeans(rp.stat$statistic)
  names(mean.stats) <- paste("Mean", names(mean.stats))

  result <- structure(list(statistics.mean = mean.stats,
                           p.values.fdr = rp.pvalue,
                           proj.statistics = rp.stat$statistic,
                           boot.proj.statistics = rp.stat.star,
                           proj.p.values = pval, method = meth, B = B,
                           n.proj = vec.nproj, projs = projs,
                           type.basis = type.basis, beta.est = beta.est,
                           p = p.opt, p.criterion = p.criterion,
                           data.name = "Y = <X, b> + e"))
  class(result) <- "htest"
  return(result)

}
