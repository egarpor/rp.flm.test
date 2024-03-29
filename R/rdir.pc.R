

#' @title Data-driven sampling of random directions guided by a sample of functional data
#'
#' @description Generation of random directions based on the principal components \eqn{\hat e_1,\ldots,\hat e_k}{\hat e_1,...,\hat e_k} of a sample of functional data \eqn{X_1,\ldots,X_n}{X_1,...,X_n}. The random directions are sampled as
#' \deqn{h=\sum_{j=1}^kh_j\hat e_j,}{h=\sum_{j=1}^kh_j\hat e_j,}
#' with \eqn{h_j\sim\mathcal{N}(0, \sigma_j^2)}{h_j~N(0, \sigma_j^2)}, \eqn{j=1,\ldots,k}{j=1,...,k}. Useful for sampling non-orthogonal random directions \eqn{h}{h} such that they are non-orthogonal for the random sample.
#'
#' @param n number of curves to be generated.
#' @param X.fdata an \code{\link[fda.usc]{fdata}} object used to compute the functional principal components.
#' @param ncomp if an integer vector is provided, the index for the principal components to be considered. If a threshold between \code{0} and \code{1} is given, the number of components \eqn{k}{k} is determined automatically as the minimum number that explains at least the \code{ncomp} proportion of the total variance of \code{X.fdata}.
#' @param fdata2pc.obj output of \code{\link[fda.usc]{fdata2pc}} containing as many components as the ones to be selected by \code{ncomp}. Otherwise, it is computed internally.
#' @param sd if \code{0}, the standard deviations \eqn{\sigma_j} are estimated by the standard deviations of the scores for \eqn{e_j}. If not, the \eqn{\sigma_j}'s are set to \code{sd}.
#' @param zero.mean whether the projections should have zero mean. If not, the mean is set to the mean of \code{X.fdata}.
#' @param norm whether the samples should be L2-normalized or not.
#' @return A \code{\link[fda.usc]{fdata}} object with the sampled directions.
#' @examples
#' # Simulate some data
#' set.seed(345673)
#' X.fdata <- r.ou(n = 200, mu = 0, alpha = 1, sigma = 2, t = seq(0, 1, l = 201),
#'                 x0 = rep(0, 200))
#' pc <- fdata2pc(X.fdata, ncomp = 20)
#'
#' # Increasing n
#' set.seed(34567)
#' rdir.pc(n = 3, X.fdata = X.fdata)$data[, 1:10]
#' set.seed(34567)
#' rdir.pc(n = 5, X.fdata = X.fdata, fdata2pc.obj = pc)$data[, 1:10]
#' 
#' \donttest{
#' # Comparison for the variance type
#' set.seed(456732)
#' n.proj <- 100
#' set.seed(456732)
#' samp1 <- rdir.pc(n = n.proj, X.fdata = X.fdata, sd = 1, norm = FALSE, ncomp = 0.99)
#' set.seed(456732)
#' samp2 <- rdir.pc(n = n.proj, X.fdata = X.fdata, sd = 0, norm = FALSE, ncomp = 0.99)
#' set.seed(456732)
#' samp3 <- rdir.pc(n = n.proj, X.fdata = X.fdata, sd = 1, norm = TRUE, ncomp = 0.99)
#' set.seed(456732)
#' samp4 <- rdir.pc(n = n.proj, X.fdata = X.fdata, sd = 0, norm = TRUE, ncomp = 0.99)
#' par(mfrow = c(1, 2))
#' plot(X.fdata, col = gray(0.85), lty = 1)
#' lines(samp1[1:10], col = 2, lty = 1)
#' lines(samp2[1:10], col = 4, lty = 1)
#' legend("topleft", legend = c("Data", "Different variances", "Equal variances"),
#'        col = c(gray(0.85), 2, 4), lwd = 2)
#' plot(X.fdata, col = gray(0.85), lty = 1)
#' lines(samp3[1:10], col = 5, lty = 1)
#' lines(samp4[1:10], col = 6, lty = 1)
#' legend("topleft", legend = c("Data", "Different variances, normalized",
#'        "Equal variances, normalized"), col = c(gray(0.85), 5:6), lwd = 2)
#'
#' # Correlations (stronger with different variances and unnormalized;
#' # stronger with lower ncomp)
#' ind <- lower.tri(matrix(nrow = n.proj, ncol = n.proj))
#' median(abs(cor(sapply(1:n.proj, function(i) inprod.fdata(X.fdata, samp1[i]))))[ind])
#' median(abs(cor(sapply(1:n.proj, function(i) inprod.fdata(X.fdata, samp2[i]))))[ind])
#' median(abs(cor(sapply(1:n.proj, function(i) inprod.fdata(X.fdata, samp3[i]))))[ind])
#' median(abs(cor(sapply(1:n.proj, function(i) inprod.fdata(X.fdata, samp4[i]))))[ind])
#'
#' # Comparison for the threshold
#' samp1 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.25, fdata2pc.obj = pc)
#' samp2 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.50, fdata2pc.obj = pc)
#' samp3 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.90, fdata2pc.obj = pc)
#' samp4 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.95, fdata2pc.obj = pc)
#' samp5 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.99, fdata2pc.obj = pc)
#' cols <- rainbow(5, alpha = 0.25)
#' par(mfrow = c(3, 2))
#' plot(X.fdata, col = gray(0.75), lty = 1, main = "Data")
#' plot(samp1, col = cols[1], lty = 1, main = "Threshold = 0.25")
#' plot(samp2, col = cols[2], lty = 1, main = "Threshold = 0.50")
#' plot(samp3, col = cols[3], lty = 1, main = "Threshold = 0.90")
#' plot(samp4, col = cols[4], lty = 1, main = "Threshold = 0.95")
#' plot(samp5, col = cols[5], lty = 1, main = "Threshold = 0.99")
#'
#' # Normalizing
#' samp1 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.50, fdata2pc.obj = pc,
#'                  norm = TRUE)
#' samp2 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.90, fdata2pc.obj = pc,
#'                  norm = TRUE)
#' samp3 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.95, fdata2pc.obj = pc,
#'                  norm = TRUE)
#' samp4 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.99, fdata2pc.obj = pc,
#'                  norm = TRUE)
#' samp5 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.999, fdata2pc.obj = pc,
#'                  norm = TRUE)
#' cols <- rainbow(5, alpha = 0.25)
#' par(mfrow = c(3, 2))
#' plot(X.fdata, col = gray(0.75), lty = 1, main = "Data")
#' plot(samp1, col = cols[1], lty = 1, main = "Threshold = 0.50")
#' plot(samp2, col = cols[2], lty = 1, main = "Threshold = 0.90")
#' plot(samp3, col = cols[3], lty = 1, main = "Threshold = 0.95")
#' plot(samp4, col = cols[4], lty = 1, main = "Threshold = 0.99")
#' plot(samp5, col = cols[5], lty = 1, main = "Threshold = 0.999")
#' }
#' @author Eduardo García-Portugués (\email{edgarcia@@est-econ.uc3m.es}) and Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}).
#' @export
rdir.pc <- function(n, X.fdata, ncomp = 0.95, 
                    fdata2pc.obj = 
                      fda.usc::fdata2pc(X.fdata, 
                                        ncomp = min(length(X.fdata$argvals), 
                                                    nrow(X.fdata))),
                    sd = 0, zero.mean = TRUE, norm = FALSE) {
  
  # Check fdata
  if (class(X.fdata) != "fdata") {
    
    stop("X.fdata must be of class fdata")
    
  }
  
  # Consider PCs up to a threshold of the explained variance or up to max(ncomp)
  ej <- fdata2pc.obj
  m <- switch((ncomp < 1) + 1,
              max(ncomp),
              max(2, min(which(cumsum(ej$d^2) / sum(ej$d^2) > ncomp))))
  if (ncomp < 1) {
    
    ncomp <- 1:m
    
  }
  
  # Compute PCs with fdata2pc if ej contains less eigenvectors than m
  # The problem is that fda.usc::fdata2pc computes all the PCs and then returns
  # the eigenvalues (d) for all the components but only the eigenvectors (rotation)
  # for the ncomp components.
  if (nrow(ej$rotation) < m) {
    
    ej <- fda.usc::fdata2pc(X.fdata, ncomp = m)
    
  }
  
  # Standard deviations of the normal coefficients
  if (sd == 0) {
    
    # Standard deviations of scores of X.fdata on the eigenvectors
    sdarg <- apply(ej$x[, ncomp], 2, sd)
    
  } else {
    
    # Constant standard deviation
    sdarg <- rep(sd, length(ncomp))
    
  }
  
  # Eigenvectors
  eigv <- ej$rotation[ncomp]
  
  # Compute linear combinations of the eigenvectors with coefficients sampled
  # from a centred normal with standard deviations sdarg
  x <- matrix(rnorm(n * m), nrow = n, ncol = m, byrow = TRUE)
  # byrow = TRUE ensures that, for a common seed, generating n + 1 or n samples
  # gives the same result in the first n
  x <- t(t(x) * sdarg)
  rprojs <- fda.usc::fdata(mdata = x %*% eigv$data,
                           argvals = fda.usc::argvals(X.fdata))
  
  # Normalize
  if (norm) {
    
    rprojs$data <- rprojs$data / drop(fda.usc::norm.fdata(rprojs))
    
  }
  
  # Add mean
  if (!zero.mean) {
    
    rprojs$data <- t(t(rprojs$data) + drop(ej$mean$data))
    
  }
  
  return(rprojs)
  
}

