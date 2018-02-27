
## The distribution of p-value FDR depending on the discretization of the 
## p-values and the number of p-values considered in the FDR
library(manipulate)
M <- 1e3
manipulate({
  
  # B: precision of p-values
  # K: number of p-values (e.g. projections)
  
  # FDR for a uniform sample of K p-vals with precision 1/B
  fdr <- sapply(1:M, function(i) {
    pval <- (sample(x = B + 1, size = K, replace = TRUE) - 1)/B 
    min(sort(pval) / (1:K)) * K
  })
  
  # Histogram
  hist(fdr, main = paste("B = ", B, ", K = ", K, ", mean = ", 
                         round(mean(fdr), 4), ", 10% = ", 
                         round(quantile(fdr, probs = 0.10), 4), sep = ""),
    breaks = seq(0, 1, l = 10), freq = FALSE, ylim = c(0, 1.5))
  abline(h = 1, lwd = 2, col = 2)

}, B = slider(min = 500, max = 1e4, step = 500), 
K = slider(min = 1, max = 100, step = 1))

## The empirical powers as a function of the discretization

# 95% confidence interval for the proportion p from a sample of size M
ci <- function(p, M) {
  
  p + c(-1, 1) * qnorm(0.025) * sqrt(p * (1 - p) / M)
  
}

# # Size as a function of K
# M <- 1e4
# K <- 1:50
# mK <- max(K)
# positiveCorrection <- TRUE
# av.fdr.pval <- function(b, seed = 12453252) {
# 
#   # A sample of p-vals with precision 1/b
#   set.seed(seed)
#   samp <- matrix((sample(x = b + 1, size = mK * M, replace = TRUE) - 1) +
#                    positiveCorrection, nrow = M, ncol = mK) /
#     (b + positiveCorrection)
# 
#   # Extract p-values
#   pvals <- sapply(K, function(k) {
#     apply(samp[, 1:k, drop = FALSE], 1, function(x) min(sort(x) / (1:k)) * k)
#   })
# 
#   # Power
#   apply(pvals, 2, function(x) c(mean(x < 0.10), mean(x < 0.05), mean(x < 0.01)))
# 
# }
# A <- vector("list", 5)
# i <- 1
# for (b in c(5e2, 1e3, 5e3, 1e4, 5e4)) {
# 
#   for (j in 1:5) {
# 
#     A[[i]][[j]] <- t(av.fdr.pval(b = b, seed = 202345 + j * 123))
# 
#   }
#   cat(i, "\n")
#   i <- i + 1
# 
# }

# Matplot: 5 trajectories of the empirical size of the FDR (indexed by the 
# number of p-values) with respect to different p-value discretizations
library(viridis)
cols <- rev(viridis(6))
load("discretization-fdr.RData")
# load("discretization-fdr-positiveCorrection.RData")
matplot(A[[1]][[1]], type = "n", ylim = c(0, 0.2), xlab = "Number of p-values", 
        ylab = "Empirical size for FDR p-value", cex.lab = 1.5, cex.axis = 1.5)
for (j in 1:5) matlines(A[[1]][[j]], type = "l", lty = 1, pch = 1:3, cex = 0.5,
                        col = cols[1])
for (j in 1:5) matlines(A[[2]][[j]], type = "l", lty = 1, pch = 1:3, cex = 0.5, 
                        col = cols[2])
for (j in 1:5) matlines(A[[3]][[j]], type = "l", lty = 1, pch = 1:3, cex = 0.5, 
                         col = cols[3])
for (j in 1:5) matlines(A[[4]][[j]], type = "l", lty = 1, pch = 1:3, cex = 0.5, 
                         col = cols[4])
for (j in 1:5) matlines(A[[5]][[j]], type = "l", lty = 1, pch = 1:3, cex = 0.5, 
                         col = cols[5])
abline(h = ci(p = 0.10, M = M), lty = 2, lwd = 1, col = 1)
abline(h = ci(p = 0.05, M = M), lty = 2, lwd = 1, col = 1)
abline(h = ci(p = 0.01, M = M), lty = 2, lwd = 1, col = 1)
abline(h = c(0.10, 0.05, 0.01), lty = 3, lwd = 1, col = 1)
legend("topleft", legend = c("B = 5e2", "B = 1e3", "B = 5e3", "B = 1e4", 
                             "B = 5e4"), col = cols[-6], lwd = 2, cex = 1.25)

