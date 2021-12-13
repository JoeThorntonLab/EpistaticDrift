
# Set to the supplementary file directory.
setwd('~/Desktop/Single Mutant DMS/Data/Supplementary files/')

source('Scripts/basic_functions.R')
library(ggplot2)


# Functions.

# Approximate the rate of epistatic divergence by sampling a given number of proteins.
approximate_sigma_by_sampling_proteins <- function(proteins, dF, noise, n_prot) {
  
  prot_combn <- combn(length(proteins), n_prot, simplify = FALSE)
  
  res <- matrix(NA_real_, nrow(dF), length(prot_combn))
  
  for(i in seq_along(prot_combn)) {
    
    prot <- prot_combn[[i]]
    all_pairs <- combn(prot, 2L, simplify = FALSE)
    
    ddF_prot <- vapply(all_pairs, function(pair) dF[, pair[2L]] - dF[, pair[1L]], dF[, 1L])
    eta_prot <- sqrt(vapply(all_pairs, function(pair) noise[, pair[2L]] ^ 2L + noise[, pair[1L]] ^ 2L, dF[, 1L]))
    d_prot <- vapply(all_pairs, function(pair) calc_sequence_divergence(proteins[pair[1L]], proteins[pair[2L]]), 0)
    
    sigma_prot <- ml_constant_rate(d_prot, ddF_prot, eta_prot)
    sigma_prot[rowSums(is.finite(dF[, prot])) < n_prot] <- NA_real_
    
    res[ ,i] <- sigma_prot
  }
  
  res
}


# Analyses.

DP <- read.table('dF.txt', stringsAsFactors = FALSE, header = TRUE)
proteins <- c(1L, 10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L)
intervals <- cbind(c(1L, 10L, 19L, 6L, 7L, 8L, 10L, 4L), c(10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L))

# Each interval's sequence divergence and number of substitutions.
d <- calc_d(intervals)

# Mutation effects and measurement noise.
dF <- calc_dF(proteins)
noise <- calc_noise(proteins)

# Amount of change in mutation effect (ddF) and noise in ddF (eta).
ddF <- calc_ddF(intervals, remove_out_of_range = FALSE)
eta <- calc_eta(intervals)

# Mutations with 4 or less data points.
# They are not explicitly excluded at this point.
retain <- apply(ddF, 1L, function(x) if(sum(!is.na(x)) <= 4L) FALSE else TRUE)

# Rate of change in mutation effect.
sigma <- ml_constant_rate(d, ddF, eta)
sigma[!retain] <- NA_real_


# Measuring rate of epistatic divergence using a subset of genetic backgrounds.

n_prot <- 2L:5L
correlation <- list()
bias <- list()

for(k in seq_along(n_prot)) {
  
  cat(k, '\n')
  sigma_k <- approximate_sigma_by_sampling_proteins(proteins, dF, noise, n_prot[k])
  correlation[[k]] <- as.vector(cor(sigma, sigma_k, use = 'pairwise.complete.obs') ^ 2L)
  bias[[k]] <- apply(sigma_k, 2L, function(x) lm(x ~ sigma + 0)$coefficients)
}

# Violin plot.

data <- data.frame(cor = unlist(correlation), bias = unlist(bias),
                   n = as.factor(unlist(lapply(n_prot, function(n) rep(n, choose(length(proteins), n))))))

ggplot(data, aes(n, cor)) +
  geom_violin() +
  ylim(0, 1) +
  stat_summary(fun.y = mean, geom = "point", size = 1)

ggplot(data, aes(n, bias)) +
  geom_violin() +
  ylim(0, 1.5) +
  stat_summary(fun.y = mean, geom = "point", size = 1)
