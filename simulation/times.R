
# Load packages
library(fda.usc)
library(rp.flm.test)
library(parallel)
library(simTool)

# For general simulation
simulation.function.time <- 
  function(n = 100, scenario = 1, delta = 0, n.proj = 3, B = 1e3, p = 10, 
           seed = 2345678, projs = rproc2fdata(5, t = seq(0, 1, l = 201)), 
           ...) {
    
    # Fix different seeds
    set.seed(seed)
    
    # Sample
    samp <- tryCatch(r.mod(n = n, scenario = scenario, delta = delta, 
                           R2 = 0.95, composite = TRUE), 
                     error = function(e) NA)
    
    # Get the p
    if (is.null(p)) {
      
      aux <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y,
                                  beta0.fdata = NULL, B = 2, pmax = 10,
                                  n.proj = n.proj, verbose = FALSE),
                      error = function(e) e)
      p <- tryCatch(aux$p, error = function(e) NULL)
      
    }
    
    # Flm test
    t.flm <- proc.time()
    test.flm <- tryCatch(flm.test(X.fdata = samp$X.fdata, Y = samp$Y, B = B,
                                  p = p, verbose = FALSE, plot.it = FALSE),
                         error = function(e) NA)
    t.flm <- proc.time() - t.flm
    
    # Random projection test
    t.rp <- proc.time()
    test.rp.flm <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y,
                                        beta0.fdata = NULL, B = B, 
                                        n.proj = n.proj, p = p, 
                                        verbose = FALSE, projs = projs),
                            error = function(e) e)
    t.rp <- proc.time() - t.rp
    
    # beta0
    beta0 <- samp$beta.fdata
    
    # Flm test simple
    t.flm.s <- proc.time()
    test.flm <- tryCatch(flm.test(X.fdata = samp$X.fdata, Y = samp$Y, 
                                  beta0.fdata = beta0, B = B,
                                  p = p, verbose = FALSE, plot.it = FALSE),
                         error = function(e) NA)
    t.flm.s <- proc.time() - t.flm.s
    
    # Random projection test simple
    t.rp.s <- proc.time()
    test.rp.flm <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y,
                                        beta0.fdata = beta0, B = B, 
                                        n.proj = n.proj, p = p, 
                                        verbose = FALSE, projs = projs),
                            error = function(e) e)
    t.rp.s <- proc.time() - t.rp.s
    
    # Return results
    return(c("RP" = t.rp[3], "FLM" = t.flm[3], "RP.S" = t.rp.s[3], 
             "FLM.S" = t.flm.s[3]))
    
  }

# Simulation
M <- 1e2
nn <- 2^(11:4)

# Fake data generation
fg <- function(n) n

# Create a cluster
cl <- makeCluster(4, type = "SOCK")

# Simulation study
eg <- evalGrids(
  
  dataGrid = expandGrid("fg", n = nn),
  procGrid = expandGrid(proc = "simulation.function.time"),
  replications = M,
  progress = TRUE,
  discardGeneratedData = TRUE,
  cluster = cl,
  clusterLibraries = c("fda.usc", "rp.flm.test"),
  fallback = "times",
  clusterSeed = rep(12345678, 6)
  
)

# Stop cluster
stopCluster(cl)

# As data frame
times <- as.data.frame.evalGrid(eg, summary.fun = mean)[, -c(1:3, 5:6)]

# Save brute data
save(list = "times", file = "times.RData")

# Times
save <- TRUE
if (save) pdf("times.pdf")
plot(times[, c(1, 3)], type = "o", pch = 19, cex = 0.5, col = 3,
     xlab = "Sample size", ylab = "Time (seconds)", cex.lab = 1.5,
     axes = FALSE)
axis(1, at = 2^(11:4), cex.axis = 1.25); axis(2, cex.axis = 1.25); box()
lines(times[, c(1, 2)], type = "o", pch = 19, cex = 0.5, col = 4)
lines(times[, c(1, 5)], type = "o", pch = 19, cex = 0.5, lty = 2, col = 3)
lines(times[, c(1, 4)], type = "o", pch = 19, cex = 0.5, lty = 2, col = 4)
legend("topleft", lwd = 2, col = rep(4:3, each = 2), lty = rep(1:2, 2), 
       legend = c("CvM and KS, composite", "CvM and KS, simple", 
                  "PCvM, composite", "PCvM, simple"), cex = 1.25)
if (save) dev.off()

