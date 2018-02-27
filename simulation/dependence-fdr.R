
# Required packages
library(mvtnorm)
library(car)
library(fda.usc)
library(viridis)
library(rgl)

## Dependence creates underrejection in FDR 

# Data from a Gaussian copula
x <- rmvnorm(1e4, sigma = rbind(c(1, 0.75), c(0.75, 1)))
u1 <- pnorm(x[, 1])
v1 <- pnorm(x[, 2])

# Independent data
x <- rmvnorm(1e4, sigma = rbind(c(1, 0), c(0, 1)))
u2 <- pnorm(x[, 1])
v2 <- pnorm(x[, 2])

# Compute FDR pvalues
A1 <- cbind(u1, v1, apply(cbind(u1, v1), 1, pvalue.FDR))[1:1e3, ]
A2 <- cbind(u2, v2, apply(cbind(u2, v2), 1, pvalue.FDR))[1:1e3, ]

# Scatterplot matrices
scatterplotMatrix(A1, diagonal = 'histogram')
scatterplotMatrix(A2, diagonal = 'histogram')

# 3D representation
plot3d(x = A1[, 1], y = A1[, 2], z = A1[, 3], zlab = "FDR pval", ylab = "pval2", xlab = "pval1")
points3d(x = A2[, 1], y = A2[, 2], z = A2[, 3], col = 2)

# The FDR for just two p-values
pvalue.FDR.3D <- function(pval1, pval2) {
  
  ind <- (pval1 < pval2)
  pmin(ind * pmin(2 * pval1, pval2) + (1 - ind) * pmin(2 * pval2, pval1), 1)
  
}

# Contour representaion of the FDR transformation
x <- seq(0, 1, l = 100)
filled.contour(x, x, outer(x, x, pvalue.FDR.3D))

# Dependent case
filled.contour(x, x, outer(x, x, pvalue.FDR.3D), plot.title = {title("Dependence")}, 
               plot.axes = {points(u1, v1, pch = 1, cex = 0.35)})
hist(A1[, 3], main = "p-value FDR dependence", freq = FALSE)

# Independent case
filled.contour(x, x, outer(x, x, pvalue.FDR.3D), plot.title = {title("Independence")}, 
               plot.axes = {points(u2, v2, pch = 1, cex = 0.35)})
hist(A2[, 3], main = "p-value FDR independence", freq = FALSE)



