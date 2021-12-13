
# Set to the supplementary file directory.
setwd('~/Desktop/Single Mutant DMS/Data/Supplementary files/')

source('Scripts/basic_functions.R')


# Functions.

# Calculate normalized ddF.
calc_norm_ddF <- function(d, ddF, eta, return_matrix = FALSE) {
  
  sigma_constant <- ml_constant_rate(d, ddF, eta)
  norm_ddF <- ddF / sqrt(outer(sigma_constant ^ 2L, d) + eta ^ 2L)
  if(return_matrix) norm_ddF else as.vector(norm_ddF[is.finite(norm_ddF)])
}

# Calculate the maximum log-likelihood of the generic-rate model.
ml_generic_rate <- function(d, ddF, eta) {
  
  sigma_squared <- (ddF ^ 2L - eta ^ 2L) / outer(rep(1, nrow(ddF)), d)
  sigma_squared[which(sigma_squared < 0, arr.ind = TRUE)] <- 0
  
  mll <- -0.5 * ddF ^ 2L / (sigma_squared * outer(rep(1, nrow(ddF)), d) + eta ^ 2L) - 
    0.5 * log(sigma_squared * outer(rep(1, nrow(ddF)), d) + eta ^ 2L) - 0.5 * log(2 * pi)
  
  apply(mll, 1L, sum, na.rm = TRUE)
}

# Simulate data under the best-fit constant-rate model.
simulate_ddF_constant_rate <- function(proteins, intervals, dF, ddF, eta, lambda = NULL) {
  
  d <- calc_d(intervals)
  if(is.null(lambda)) lambda <- rep(1, length(d))
  
  # Best-fit constant-rate model.
  sigma <- ml_constant_rate(d, ddF, eta)
  
  eta_mod <- eta
  eta_mod[!is.finite(eta)] <- 0
  
  t(vapply(1L:nrow(ddF), function(i) {
    
    redo <- TRUE
    
    while(redo) {
      
      sim_ddF <- rnorm(ncol(ddF), 0, 1) * sqrt(lambda * sigma[i] ^ 2L * d + eta_mod[i, ] ^ 2L)
      
      # Simulated dF; first the AncNR3-C. teleta trajectory and then the AncSR-HsGR.
      sim_dF <- cumsum(c(0, sim_ddF[1L:6L]))
      sim_dF <- c(sim_dF, cumsum(sim_ddF[7L:8L]) + sim_dF[2L])
      
      # Randomly choosing a beginning point.
      begin <- sample(which(is.finite(dF[i, ])), 1L)
      sim_dF <- sim_dF - sim_dF[begin] + dF[i, begin]
      
      # Simulated ddF, imposing measurement bounds.
      anc_dF <- sim_dF[match(intervals[, 1L], proteins)]
      der_dF <- sim_dF[match(intervals[, 2L], proteins)]
      sim_ddF <- der_dF - anc_dF
      sim_ddF[anc_dF < L - WT_ACTIVITY + 0.05 & der_dF < L - WT_ACTIVITY + 0.05] <- NA_real_
      sim_ddF[anc_dF > U - WT_ACTIVITY - 0.05 & der_dF > U - WT_ACTIVITY - 0.05] <- NA_real_
      
      # Imposing the pattern of missing data.
      sim_ddF[!is.finite(eta[i, ])] <- NA_real_
      
      if(sum(is.finite(sim_ddF)) == sum(is.finite(ddF[i, ]))) redo <- FALSE
    }
    
    sim_ddF
  }, d))
}

# Simulate the normalized ddF distribution under the best-fit constant-rate model.
simulate_norm_ddF_constant_rate <- function(proteins, intervals, dF, ddF, eta, nrep) {
  
  d <- calc_d(intervals)
  sigma <- ml_constant_rate(d, ddF, eta)
  distr <- list()
  
  for(rep in 1L:nrep) {
    
    sim_ddF <- simulate_ddF_constant_rate(proteins, intervals, dF, ddF, eta)
    sim_sigma <- ml_constant_rate(d, sim_ddF, eta)
    distr[[rep]] <- sim_ddF / sqrt(outer(sim_sigma ^ 2L, d) + eta ^ 2L)
    if(rep %% 10L == 0L) cat('Progress: ', rep, '/', nrep, '\n', sep = '')
  }
  
  distr <- unlist(distr)
  distr <- distr[is.finite(distr)]
  
  # Enforcing symmetry.
  c(-distr, distr)
}

# Perform a likelihood-ratio test of constant-rate (null) vs. generic-rate (alternative).
compare_constant_vs_generic_rate <- function(proteins, intervals, dF, ddF, eta, niter) {
  
  # Calculate the null distribution of LRS by parametric bootstrapping.
  calc_null_distr <- function(sigma) {
    
    null_distr <- matrix(NA_real_, length(sigma), niter)
    
    for(i in 1L:niter) {
      
      # parametric bootstrap
      sim_ddF <- simulate_ddF_constant_rate(proteins, intervals, dF, ddF, eta)
      
      null_distr[, i] <- 2 * (ml_generic_rate(d, sim_ddF, eta) - ml_constant_rate(d, sim_ddF, eta, TRUE))
      
      if(i %% 10L == 0L) cat('Progress: ', i, '/', niter, '\n', sep = '')
    }
    
    null_distr
  }
  
  d <- calc_d(intervals)
  sigma <- ml_constant_rate(d, ddF, eta)
  
  # Observed likelihood-ratio statistics.
  lrs <- 2 * (ml_generic_rate(d, ddF, eta) - ml_constant_rate(d, ddF, eta, TRUE))
  
  # Null distribution of likelihood-ratio statistics for each mutation.
  null_distr <- calc_null_distr(sigma)
  
  # p value for each mutation.
  vapply(1L:length(sigma), function(i) sum(null_distr[i, ] >= lrs[i]) / niter, 0)
}


# Analyses.

DP <- read.table('dF.txt', stringsAsFactors = FALSE, header = TRUE)
proteins <- c(1L, 10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L)
# proteins <- c(2L, 11L, 20L, 6L, 7L, 9L, 12L, 5L, 14L) # Alt-All
intervals <- cbind(c(1L, 10L, 19L, 6L, 7L, 8L, 10L, 4L), c(10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L))
# intervals <- matrix(c(c(2L, 11L, 20L, 6L, 7L, 9L, 11L, 5L), c(11L, 20L, 6L, 7L, 9L, 12L, 5L, 14L)), nrow = 8L, ncol = 2L) # Alt-All


# Length of each phylogenetic interval.
d <- calc_d(intervals)

# The effects of mutations in each background.
dF <- calc_dF(proteins)

# The amount of change in each mutation's effect across each unit interval and its measurement noise.
ddF <- calc_ddF(intervals)
eta <- calc_eta(intervals)

# Mutations with 4 or less data points are excluded from analysis.
# These include mutations always at the lower bound of measurement.
retain <- apply(ddF, 1L, function(x) if(sum(!is.na(x)) <= 4L) FALSE else TRUE)
dF <- dF[retain, ]
ddF <- ddF[retain, ]
eta <- eta[retain, ]


# Analyses.

# Test of constant vs. generic rate.

p <- compare_constant_vs_generic_rate(proteins, intervals, dF, ddF, eta, 250L)
# p <- as.numeric(readLines('figure_3/p_likelihood_ratio_test.txt'))
h <- hist(p, breaks = seq(0, 1, by = 0.05))
h_sig <- hist(p[p <= calc_BH_cutoff(p, 0.2)], breaks = seq(0, 1, by = 0.05))
h_sig$counts <- h_sig$counts / sum(h$counts)
h$counts <- h$counts / sum(h$counts)
plot(h, ylim = c(0, 0.15))
plot(h_sig, add = TRUE, col = 'grey')


# Histogram of normalized ddF for mutations with or without significant p-value.

# Distribution expected under constant-rate model.
sim_norm_ddF <- simulate_norm_ddF_constant_rate(proteins, intervals, dF, ddF, eta, 50L)
sim_norm_ddF[sim_norm_ddF < -3] <- -3.01
sim_norm_ddF[sim_norm_ddF > 3] <- 3.01
h_sim <- hist(sim_norm_ddF, breaks = seq(-10/3, 10/3, by = 1/3), plot = FALSE)

# Observed distributions.
norm_ddF_sig <- calc_norm_ddF(d, ddF[p <= calc_BH_cutoff(p, 0.2), ], eta[p <= calc_BH_cutoff(p, 0.2), ])
norm_ddF_sig[norm_ddF_sig < -3] <- -3.01
norm_ddF_sig[norm_ddF_sig > 3] <- 3.01
h_sig <- hist(norm_ddF_sig, breaks = seq(-10/3, 10/3, by = 1/3), plot = FALSE)
h_sig$counts <- h_sig$density
plot(h_sig, ylim = c(0, 0.8))
points(h_sim$mids, h_sim$density, type = 'l', col = 'red')

norm_ddF_insig <- calc_norm_ddF(d, ddF[p > calc_BH_cutoff(p, 0.2), ], eta[p > calc_BH_cutoff(p, 0.2), ])
norm_ddF_insig[norm_ddF_insig < -3] <- -3.01
norm_ddF_insig[norm_ddF_insig > 3] <- 3
h_insig <- hist(norm_ddF_insig, breaks = seq(-10/3, 10/3, by = 1/3), plot = FALSE)
h_insig$counts <- h_insig$density
plot(h_insig, ylim = c(0, 0.4))
points(h_sim$mids, h_sim$density, type = 'l', col = 'red')


# Plotting mutations with a particular p-value and rate.

mutant_info <- multiple_comparison_data_frame(proteins, FALSE, FALSE)[[1L]][, c('SITE', 'WTAA', 'MTAA')]
mutant_info <- mutant_info[retain, ]

index <- which(p < 0.01 & sigma > 0.15)
k <- sample(index, 1L)
mutant_info[k, ]
dotplot_mutation_effect_trajectory(proteins, c(60, 'K', 'V'))
p[k]
sigma[k]

# Normalized ddF plot.
plot(NA, xlim = c(-3, 3), ylim = c(0, 1))
abline(v = norm_ddF[k, ])

