
# Load packages
library(fda.usc)
library(rp.flm.test)
library(foreach)
library(doMC)

# Run full simulation study
run.simulations <- FALSE
type <- c("A", "B")[1]

# The number of replicates for each block of simulations = number of total
# replicates (M) / number of MS pieces. For example, if you want M = 1000
# replicates, you run 10 blocks (indexed by MS from 1 to 10) of 100
# simulations (overwrite M = 100)

# The script assumes the variable MS, indicating the chunk number, is
# available in the workspace

######################
## Simulation study ##
######################

if (run.simulations) {
  
  # RData
  user.RData <- paste(switch(type,
                             "A" = "simus.dep095.b1e4.m1e4_",
                             "B" = "simus.dep095.flm.b1e4.m1e4_"),
                      MS, ".RData", sep = "")
  
    
  # For general simulation
  simulation.function <- function(scenario = 1, n = 100, delta = 0, B = 1e4,
                                  est.method = "pc", p.fixed = FALSE,
                                  n.proj = 10, type.basis = "bspline", 
                                  composite = TRUE, flm = TRUE, pmax = 10, 
                                  p.criterion = "SICc", type.projs = "sd0", 
                                  seed, ...) {
    
    # Seed
    set.seed(seed)
    
    # Is p data driven?
    if (!p.fixed) {
      
      p.fixed <- NULL
      
    } else {
      
      # Set p.fixed to the true/sensible number of PCs or to NULL if FALSE
      # p.fixed <- switch(scenario, 3, 4, 7, 2, 2, 2, 4, 2, 4)
      p.fixed <- p.fixed
      
    }
    
    # Sample
    samp <- tryCatch(r.mod(n = n, scenario = scenario, delta = delta,
                           t = seq(0, 1, len = 201), R2 = 0.95,
                           composite = composite), error = function(e) NA)
    
    # Composite or simple hypothesis?
    if (composite) {
      
      beta0 <- NULL
      
    } else {
      
      argvals <- argvals(samp$X.fdata)
      beta0 <- fda.usc::fdata(mdata = rep(0, length(argvals)), 
                              argvals = argvals)
      
    }
    
    # Projections
    # set.seed(M - seed + 1)
    projs <- switch(type.projs,
                    "sd0" = rdir.pc(n = n.proj, X.fdata = samp$X.fdata,
                                    ncomp = 0.95, norm = FALSE, sd = 0,
                                    zero.mean = TRUE),
                    "sd1" = rdir.pc(n = n.proj, X.fdata = samp$X.fdata,
                                    ncomp = 0.95, norm = FALSE, sd = 1,
                                    zero.mean = TRUE),
                    "ind" = r.ou(n = n.proj, t = seq(0, 1, l = 201),
                                 mu = 0, alpha = 0.5, sigma = 1))
    
    # Random projection test
    test.rp.flm <- tryCatch(rp.flm.test(X.fdata = samp$X.fdata, Y = samp$Y,
                                        beta0.fdata = beta0, B = B, pmax = pmax,
                                        est.method = est.method, 
                                        n.proj = n.proj,
                                        p.criterion = p.criterion, p = p.fixed,
                                        F.code = TRUE, verbose = FALSE,
                                        type.basis = type.basis, projs = projs),
                            error = function(e) e)
    
    # Flm test
    if (flm & composite) {
      
      test.flm <- tryCatch(fda.usc::flm.test(X.fdata = samp$X.fdata, Y = samp$Y,
                                             B = B, est.method = est.method,
                                             p = test.rp.flm$p, verbose = FALSE,
                                             plot.it = FALSE),
                           error = function(e) NA)
      
    } else {
      
      test.flm <- list()
      test.flm$p.value <- NA
      
    }
    
    # Attach call
    result <- c("p" = test.rp.flm$p,
                "CvM" = test.rp.flm$p.values.fdr[, 1],
                "KS" = test.rp.flm$p.values.fdr[, 2],
                "PCvM" = test.flm$p.value)
    return(result)
    
  }
  
  # For a
  simulation.function.a <- function(...) {
    
    res.0 <- simulation.function(..., B = 1e4, n.proj = 50, 
                                 flm = FALSE, type.projs = "sd0")
    res.1 <- simulation.function(..., B = 1e4, n.proj = 50, p.fixed = res.0[1],
                                 flm = FALSE, type.projs = "sd1")
    res.i <- simulation.function(..., B = 1e4, n.proj = 50, p.fixed = res.0[1],
                                 flm = FALSE, type.projs = "ind")
    l <- length(res.0)
    names(res.0) <- paste(names(res.0), "sd0", sep = ".")
    names(res.0)[c(1, l)] <- c("p", "PCvM")
    names(res.1) <- paste(names(res.1), "sd1", sep = ".")
    names(res.i) <- paste(names(res.i), "ind", sep = ".")
    c(res.0[-l], res.1[-c(1, l)], res.i[-c(1, l)], res.0[l])
    
  }
  
  # For b
  simulation.function.b <- function(...) {
    
    simulation.function(..., B = 1e4, n.proj = 10, flm = TRUE, 
                        type.projs = "sd0")
    
  }
  
  # Monte Carlo replicates
  M <- 1e4
  
  # Models
  scenarios <- 1:9
  
  # Sample sizes
  n <- c(50, 100, 250)
  
  # Deviations
  delta <- 0:2
  
  # Test
  # system.time(res <- simulation.function(scenario = 9, n = 100, seed = 345))
  # system.time(res <- simulation.function.a(scenario = 9, n = 100, seed = 345))
  # system.time(res <- simulation.function.b(scenario = 9, n = 100, seed = 345))
  
  # Compute the seeds
  seeds <- 1:M
  
  ######################################################################
  # Simulation study for the dependence of the tests on the number and #
  # kind of projections (B = 1e4, M = 1e4, n.proj = 50, flm = FALSE)   #
  ######################################################################
  
  if (type == "A") {
    
    # Cores
    nslots <- 20 # as.numeric(Sys.getenv("NSLOTS"))
    registerDoMC(nslots)
    
    # Look for MS variable
    if (exists("MS")) {
      
      # The number of replicates for each block of simulations = number of total
      # replicates (M) / number of MS pieces. For example, if you want M = 1000
      # replicates, you run 10 blocks (indexed by MS from 1 to 10) of 100
      # simulations (overwrite M = 100)
      M <- 20
      seeds <- seeds[((MS - 1) * M + 1):(MS * M)]
      
    }
    
    # Build up the cases
    cases <- expand.grid(n.M = 1:M, n.sce = scenarios,
                         n.n = n, n.delta = delta) # Do it reversely
    length.cases <- nrow(cases)
    cases <- cbind(1:length.cases, cases, rep(seeds, length.cases / M));
    colnames(cases)[c(1, ncol(cases))] <- c("case", "seed")
    
    # Launch parallel loop
    cat("Launch parallel loop\n")
    toc <- proc.time()
    results <- foreach(i = 1:length.cases, .combine = rbind) %dopar% {
      
      tac <- proc.time()[3]
      res <- c(cases[i, 1:ncol(cases)],
               tryCatch(simulation.function.a(scenario = cases$n.sce[i],
                                              n = cases$n.n[i],
                                              delta = cases$n.delta[i],
                                              seed = cases$seed[i]),
                        error = function(e) NA))
      tac <- proc.time()[3] - tac
      cat(paste("Iteration: ", i, "/", length.cases, " finished in ",
                sprintf("%.2f", tac)," secs. \tEstimated completion: ",
                Sys.time() + (length.cases - i) * tac / nslots,
                ". \tTime consumed: ",
                sprintf("%.2f", (proc.time() - toc)[3] / 3600), " hours.\n",
                sep = ""))
      res
      
    }
    
    # Time
    tic <- proc.time() - toc
    
    # Plot timing
    cat("Time elapsed :", tic[3], "\n")
    
    # Save results
    save("results", file = user.RData)
    
  }
  
  ########################################################################
  # Simulation study for selected number and kind of projections,        #
  # with comparison with flm (B = 1e4, M = 1e4, n.proj = 10, flm = TRUE) #
  ########################################################################
  
  if (type == "B") {
    
    # Cores
    nslots <- 20 # as.numeric(Sys.getenv("NSLOTS"))
    registerDoMC(nslots)
    
    # Look for MS variable
    if (exists("MS")) {
      
      # The number of replicates for each block of simulations = number of total
      # replicates (M) / number of MS pieces. For example, if you want M = 1000
      # replicates, you run 10 blocks (indexed by MS from 1 to 10) of 100
      # simulations (overwrite M = 100)
      M <- 40
      seeds <- seeds[((MS - 1) * M + 1):(MS * M)]
      
    }
    
    # Build up the cases
    cases <- expand.grid(n.M = 1:M, n.sce = scenarios,
                         n.n = n, n.delta = delta) # Do it reversely
    length.cases <- nrow(cases)
    cases <- cbind(1:length.cases, cases, rep(seeds, length.cases / M));
    colnames(cases)[c(1, ncol(cases))] <- c("case", "seed")
    
    # Launch parallel loop
    cat("Launch parallel loop\n")
    toc <- proc.time()
    results <- foreach(i = 1:length.cases, .combine = rbind) %dopar% {
      
      tac <- proc.time()[3]
      res <- c(cases[i, 1:ncol(cases)],
               tryCatch(simulation.function.b(scenario = cases$n.sce[i],
                                              n = cases$n.n[i],
                                              delta = cases$n.delta[i],
                                              seed = cases$seed[i]),
                        error = function(e) NA))
      tac <- proc.time()[3] - tac
      cat(paste("Iteration: ", i, "/", length.cases, " finished in ",
                sprintf("%.2f", tac)," secs. \tEstimated completion: ",
                Sys.time() + (length.cases - i) * tac / nslots,
                ". \tTime consumed: ",
                sprintf("%.2f", (proc.time() - toc)[3] / 3600), " hours.\n",
                sep = ""))
      res
      
    }
    
    # Time
    tic <- proc.time() - toc
    
    # Plot timing
    cat("Time elapsed :", tic[3], "\n")
    
    # Save results
    save("results", file = user.RData)
    
  }
  
}

##############
## Analysis ##
##############

if (!run.simulations) {

  ##################
  ## Read results ##
  ##################
  
  preprocess.results <- TRUE
  if (preprocess.results) {
    
    load("A/simus.dep095.b1e4.m1e4_1.RData")
    nr <- nrow(results)
    resultsA <- matrix(nrow = nr * 50, ncol = ncol(results))
    for (j in 1:10) {
      
      for (k in 1:50) {
        
        jk <- 50 * (j - 1) + k
        try({
          load(paste("A/simus.dep095.b1e4.m1e4_", jk, ".RData", sep = ""))
          resultsA[(nr * (k - 1) + 1):(nr * k), ] <- 
            t(apply(results, 1, unlist))
        })
        cat(jk, "\n")
    
      }
      colnames(resultsA) <- colnames(results)
      resultsA <- as.data.frame(resultsA)
      save(list = "resultsA", file = paste("resA_", j, ".RData", sep = ""))
    
    }
    
    load("B/simus.dep095.flm.b1e4.m1e4_1.RData")
    nr <- nrow(results)
    resB <- matrix(nrow = nrow(results) * 250, ncol = ncol(results))
    for (k in 1:250) {
  
      load(paste("B/simus.dep095.flm.b1e4.m1e4_", k, ".RData", sep = ""))
      resB[(nr * (k - 1) + 1):(nr * k), ] <- t(apply(results, 1, unlist))
      cat(k, "\n")
  
    }
    colnames(resB) <- colnames(results)
    resB <- as.data.frame(resB)
    save(list = "resB", file = "resB.RData")
    
  }

  ######################
  # Required functions #
  ######################
  
  # 95% confidence interval for the proportion p from a sample of size M
  ci <- function(p, M) {
    
    p + c(-1, 1) * qnorm(0.025) * sqrt(p * (1 - p) / M)
    
  }
  
  # Powers of the FDR p-values indexed in the number of projections
  trajs.power <- function(n, delta = 0, alpha = 0.05, stat = "CvM") {
    
    # Which level?
    pow <- switch(as.character(alpha),
                  "0.1" = pow.10,
                  "0.05" = pow.05,
                  "0.01" = pow.01)
    
    # Which statistic?
    ind <- switch(stat,
                  "CvM" = ind.CvM,
                  "KS" = ind.KS,
                  "PCvM" = ind.PCvM)
    
    # Transposed trajectory for matplot
    t(pow[pow$delta == delta & pow$n == n, ind])
    
  }
  
  # Trajectories plot
  trajs.plot <- function(n, scenarios = 1:9, alpha = 0.05, deltas = 0,
                         stat = "CvM", type.projs = "sd0", n.proj = 1:50, 
                         M = 1e4, main, cex.main = 1, cex.lab = 1, 
                         cex.text = 1, cex.axis = 1) {
    
    # Colors
    col.05 <- "darkorchid4"
    col.10 <- "darkorange3"
    col.01 <- "firebrick3"
    
    if (identical(deltas, 0)) {
      
      par(mfrow = c(length(stat), length(n)))
      
    } else {
      
      par(mfrow = c(1, length(stat) * length(deltas)))
      
    }
    
    # Indices scenarios (aggregate returns them sorted)
    ind.sce <- match(x = scenarios, table = sort(unique(pow.01$scenario)))
    
    # Type of projections
    shift <- switch(type.projs, "sd0" = 0, "sd1" = 1, "ind" = 2) / 3
    
    # Plots
    show.main <- missing(main)
    for (tt in rev(stat)) {
      
      for (nn in n) {
        
        # Caption
        if (show.main) {
          
          main <- paste(stat, ", n = ", nn, sep = "")
          
        }
        
        for (delta in deltas) {
          
          if (delta == 0) {
            
            A1 <- trajs.power(n = nn, delta = delta, alpha = 0.05, stat = stat)
            matplot(n.proj, A1[n.proj + shift * nrow(A1), ind.sce], type = "o",
                    pch = 16, cex = 0.5, col = col.05, ylim = c(0, 0.15), 
                    lty = 1,  xlim = c(1, max(n.proj)), 
                    xlab = "Number of projections", ylab = "", main = main, 
                    cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)
            mtext("Empirical rejection rate", side = 2, line = 2.75,
                  cex = cex.lab)
            A2 <- trajs.power(n = nn, delta = delta, alpha = 0.10, stat = stat)
            matlines(n.proj, A2[n.proj + shift * nrow(A2), ind.sce], type = "o",
                     pch = 16, cex = 0.5, lty = 1, col = col.10)
            A3 <- trajs.power(n = nn, delta = delta, alpha = 0.01, stat = stat)
            matlines(n.proj, A3[n.proj + shift * nrow(A3), ind.sce], type = "o", 
                     pch = 16, cex = 0.5, lty = 1, col = col.01)
            abline(h = ci(p = 0.05, M = M), lwd = 3, col = col.05, lty = 2)
            abline(h = ci(p = 0.10, M = M), lwd = 3, col = col.10, lty = 2)
            abline(h = ci(p = 0.01, M = M), lwd = 3, col = col.01, lty = 2)
            abline(h = 0.05, lwd = 1, col = col.05, lty = 3)
            abline(h = 0.10, lwd = 1, col = col.10, lty = 3)
            abline(h = 0.01, lwd = 1, col = col.01, lty = 3)
            
          } else {
            
            A <- trajs.power(n = nn, delta = delta, alpha = alpha, stat = stat)
            matplot(n.proj, A[n.proj + shift * nrow(A), ind.sce], type = "o", 
                    pch = 16, cex = 0.5, ylim = c(0, 1), lty = 1, 
                    xlim = c(1, max(n.proj)), xlab = "Number of projections", 
                    ylab = "", main = main, cex.main = cex.main, 
                    cex.lab = cex.lab, cex.axis = cex.axis, 
                    col = rainbow(length(scenarios)))
            mtext("Empirical rejection rate", side = 2, line = 2.75, 
                  cex = cex.lab)
            for (k in seq_along(scenarios)) {
              
              text(x = n.proj[1], y = A[n.proj[1] + shift * nrow(A), 
                                        ind.sce[k]], 
                   labels = k, pos = 2, offset = 0.35, 
                   col = rainbow(length(scenarios))[k], cex = cex.text)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  # Average power loss of PCvM wrt to CvM3
  power.loss <- function(n, d, k = 3) {
    
    pcvm <- trajs.power(n = n, delta = d, stat = "PCvM")[1, ]
    cvm <- trajs.power(n = n, delta = d, stat = "CvM")[k, ]
    if (d == 0) {
      
      return(summary(abs(pcvm - 0.05) - abs(cvm - 0.05)))
      
    } else {
      
      return(summary((cvm - pcvm) / pcvm))
      
    }
    
  }
  
  # Create table for alpha = 0.05
  block.1 <- function(d = 0, a = 0.05, ind.k = c(1, 3, 5)) {
    
    t(rbind(
      trajs.power(n = 100, delta = d, alpha = a, stat = "CvM")[ind.k, ],
      trajs.power(n = 100, delta = d, alpha = a, stat = "KS")[ind.k, ],
      trajs.power(n = 100, delta = d, alpha = a, stat = "PCvM"),
      trajs.power(n = 250, delta = d, alpha = a, stat = "CvM")[ind.k, ],
      trajs.power(n = 250, delta = d, alpha = a, stat = "KS")[ind.k, ],
      trajs.power(n = 250, delta = d, alpha = a, stat = "PCvM")
    ))
    
  }
  block.2 <- function(d = 0, a = 0.05, ind.k = c(1, 3, 5)) {
    
    t(rbind(
      trajs.power(n = 50, delta = d, alpha = a, stat = "CvM")[ind.k, ],
      trajs.power(n = 50, delta = d, alpha = a, stat = "KS")[ind.k, ],
      trajs.power(n = 50, delta = d, alpha = a, stat = "PCvM")
    ))
    
  }
  block.3 <- function(d = 0, a = 0.05, ind.k = c(1, 3, 5)) {
    
    t(rbind(
      trajs.power(n = 50, delta = d, alpha = a, stat = "CvM")[ind.k, ],
      trajs.power(n = 50, delta = d, alpha = a, stat = "KS")[ind.k, ],
      trajs.power(n = 100, delta = d, alpha = a, stat = "CvM")[ind.k, ],
      trajs.power(n = 100, delta = d, alpha = a, stat = "KS")[ind.k, ],
      trajs.power(n = 250, delta = d, alpha = a, stat = "CvM")[ind.k, ],
      trajs.power(n = 250, delta = d, alpha = a, stat = "KS")[ind.k, ]
    ))
    
  }

  ######################################################################
  # Simulation study for the dependence of the tests on the number and #
  # kind of projections (B = 1e4, M = 1e4, n.proj = 50, flm = FALSE)   #
  ######################################################################
  
  if (type == "A") {
  
    # Read data
    load("resA_1.RData")
    nr <- nrow(resultsA)
    resA <- matrix(nrow = nrow(resultsA) * 10, ncol = ncol(resultsA))
    for (k in 1:10) {
      
      load(paste("resA_", k, ".RData", sep = ""))
      resA[(nr * (k - 1) + 1):(nr * k), ] <- t(apply(resultsA, 1, unlist))
      cat(k, "\n")
      
    }
    
    df <- resA
    
    # Powers
    power <- function(x, alpha = 0.05) mean(x < alpha, na.rm = TRUE)
    by <- list(df$n.delta, df$n.n, df$n.sce)
    pow.10 <- aggregate.data.frame(x = df[, -c(1:7)], by = by, 
                                   FUN = power, alpha = 0.10)
    pow.05 <- aggregate.data.frame(x = df[, -c(1:7)], by = by, 
                                   FUN = power, alpha = 0.05)
    pow.01 <- aggregate.data.frame(x = df[, -c(1:7)], by = by, 
                                   FUN = power, alpha = 0.01)
    colnames(pow.10)[1:3] <- colnames(pow.05)[1:3] <- colnames(pow.01)[1:3] <-
      c("delta", "n", "scenario")
    gc()
    
    # CvM and KS indexes
    ind.CvM <- grep(pattern = glob2rx("CvM*"), x = colnames(pow.05))
    ind.KS <- grep(pattern = glob2rx("KS*"), x = colnames(pow.05))
    ind.PCvM <- grep(pattern = "PCvM", x = colnames(pow.05))
    
    ## Graphs FDR p-value vs projections for level
    
    trajs.plot(n = c(50, 100, 250), stat = "CvM")
    trajs.plot(n = c(50, 100, 250), stat = "KS")
    
    ## Graphs FDR p-values vs projections for power
    
    trajs.plot(n = 50, deltas = 1:2, stat = "CvM")
    trajs.plot(n = 100, deltas = 1:2, stat = "CvM")
    trajs.plot(n = 250, deltas = 1:2, stat = "CvM")
    trajs.plot(n = 50, deltas = 1:2, stat = "KS")
    trajs.plot(n = 100, deltas = 1:2, stat = "KS")
    trajs.plot(n = 250, deltas = 1:2, stat = "KS")
    
    ## Individual plots for sizes and powers
    
    save <- TRUE
    pow <- TRUE
    suffix <- c("sd0", "sd1", "ind")[1]#[2]
    nn <- c(50, 100, 250)#[2]
    dd <- c(1:2)#[1]
    for (n in nn) {
      
      ## H0
      
      # CvM
      if (save) pdf(paste("sizes_cvm_", n, "_", suffix, ".pdf", sep = ""))
      trajs.plot(n = n, deltas = 0, stat = "CvM", type.projs = suffix,
                 main = substitute(list(H[0], CvM, n == nn), list(nn = n)),
                 cex.main = 2.5, cex.lab = 2, cex.axis = 1.9)
      if (save) dev.off()
      
      # KS
      if (save) pdf(paste("sizes_ks_", n, "_", suffix, ".pdf", sep = ""))
      trajs.plot(n = n, deltas = 0, stat = "KS", type.projs = suffix,
                 main = substitute(list(H[0], KS, n == nn), list(nn = n)),
                 cex.main = 2.5, cex.lab = 2, cex.axis = 1.9)
      if (save) dev.off()
      
      ## H1, H2
      
      if (pow) {
        
        for (delta in dd) {
          
          # CvM
          if (save) pdf(paste("powers_cvm_", n, "_", delta, "_", suffix, 
                              ".pdf", sep = ""))
          trajs.plot(n = n, deltas = delta, stat = "CvM", type.projs = suffix,
                     main = substitute(list(H[1], CvM, n == nn, d == dd), 
                                       list(nn = n, dd = delta)),
                     cex.main = 2.5, cex.lab = 2, cex.axis = 1.9,
                     cex.text = 1.5)
          if (save) dev.off()
          
          # KS
          if (save) pdf(paste("powers_ks_", n, "_", delta, "_", suffix, 
                              ".pdf", sep = ""))
          trajs.plot(n = n, deltas = delta, stat = "KS", type.projs = suffix,
                     main = substitute(list(H[1], KS, n == nn, d == dd), 
                                       list(nn = n, dd = delta)),
                     cex.main = 2.5, cex.lab = 2, cex.axis = 1.9,
                     cex.text = 1.5)
          if (save) dev.off()
          
        }
        
      }
      
    }
    
  }
  
  ########################################################################
  # Simulation study for selected number and kind of projections,        #
  # with comparison with flm (B = 1e4, M = 1e4, n.proj = 10, flm = TRUE) #
  ########################################################################
  
  if (type == "B") {
    
    # Read data
    load("resB.RData")
    df <- resB
    
    # Powers
    power <- function(x, alpha = 0.05) mean(x < alpha, na.rm = TRUE)
    by <- list(df$n.delta, df$n.n, df$n.sce)
    pow.10 <- aggregate.data.frame(x = df[, -c(1:7)], by = by, 
                                   FUN = power, alpha = 0.10)
    pow.05 <- aggregate.data.frame(x = df[, -c(1:7)], by = by, 
                                   FUN = power, alpha = 0.05)
    pow.01 <- aggregate.data.frame(x = df[, -c(1:7)], by = by, 
                                   FUN = power, alpha = 0.01)
    colnames(pow.10)[1:3] <- colnames(pow.05)[1:3] <- colnames(pow.01)[1:3] <-
      c("delta", "n", "scenario")
    gc()
    
    # CvM and KS indexes
    ind.CvM <- grep(pattern = glob2rx("CvM*"), x = colnames(pow.05))
    ind.KS <- grep(pattern = glob2rx("KS*"), x = colnames(pow.05))
    ind.PCvM <- grep(pattern = "PCvM", x = colnames(pow.05))
    
    ## Tables n = 50, 100, 250 for alpha = 0.05
    
    # n = 100, 250
    tab.1 <- rbind(block.1(d = 0), block.1(d = 1), block.1(d = 2)) 
    # n = 50
    tab.2 <- rbind(block.2(d = 0), block.2(d = 1), block.2(d = 2))
    # n = 50, 100, 250
    tab.3 <- rbind(block.3(d = 0), block.3(d = 1), block.3(d = 2))
    
    # Print table
    # tab <- tab.3
    # tab <- tab.2
    tab <- tab.1
    k <- 1
    for (delta in 0:2) {
      
      for (scenario in seq_along(scenarios)) {
        
        H <- paste("$H_{", scenario, ",", delta, "}$", sep = "")
        cat(paste(H, paste("$", paste(sprintf(fmt = "%.1f", tab[k, ] * 100),
                                      collapse = "$ & $"), "$", sep = ""),
                  sep = " & "), "\\\\\n")
        k <- k + 1
        
      }
      if (delta == 2) {
        
        cat("\\bottomrule\\bottomrule\n")
        
      } else {
        
        cat("\\midrule\n")
        
      }
      
    }
    
    # CvM is more powerful than KS
    n1 <- 1:3; tab.3[-c(1:9), n1] >= tab.3[-c(1:9), n1 + 3]
    n2 <- 7:9; tab.3[-c(1:9), n2] >= tab.3[-c(1:9), n2 + 3]
    n3 <- 13:15; tab.3[-c(1:9), n3] >= tab.3[-c(1:9), n3 + 3]
    
    # Average power loss of PCvM wrt to CvM3
    power.loss(n = 50, d = 1)
    power.loss(n = 50, d = 2)
    power.loss(n = 100, d = 1)
    power.loss(n = 100, d = 2)
    power.loss(n = 250, d = 1)
    power.loss(n = 250, d = 2)
    
    ## Table average dn
    
    ind.p <- match(x = "p", table = colnames(df))
    summary.p <- aggregate.data.frame(x = df[, ind.p], by = by,
                                      FUN = function(x) {
                                        c(mean(x, na.rm = TRUE), 
                                          sd(x, na.rm = TRUE),
                                          quantile(x, probs = c(0, 0.25, 0.5, 
                                                                0.75, 1),
                                                   na.rm = TRUE))
                                      })
    summary.p <- as.matrix(summary.p)
    colnames(summary.p) <- c("delta", "n", "scenario", "mean.p", "sd.p",
                             "min.p", "25%.p", "50%.p", "75%.p", "max.p")
    summary.p <- as.data.frame(summary.p)
    summary.p
    
    # Table
    scenarios <- 1:9
    for (k in seq_along(scenarios)) {
      
      aux <- summary.p$mean.p[subset = summary.p$scenario == scenarios[k]]
      numbers <- sprintf("%.1f", aux)
      aux <- summary.p$sd.p[subset = summary.p$scenario == scenarios[k]]
      numbers <- paste(numbers, " (", sprintf("%.1f", aux), ")", sep = "")
      cat(paste("$H_{", k, ",\\delta}$ & ", 
                paste(numbers, collapse = " & "),  sep = ""), "\\\\\n")
      
    }
    
  }
  
}
