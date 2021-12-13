
# Set to the supplementary file directory.
setwd('~/Desktop/Single Mutant DMS/Data/Supplementary files/')

source('Scripts/basic_functions.R')


# Functions.

# Identify the number of proteins in which each mutation's effect is confidently characterized.
identify_number_of_measurements <- function(proteins) {
  
  # proteins : integer vector
  #   A list of proteins.
  #
  # All mutations identified by `multiple_comparison_data_frame` are examined.
  
  data_all <- multiple_comparison_data_frame(proteins)
  
  interval <- 0L:(length(proteins) - 1L) * (nrow(data_all) / length(proteins))
  
  n_measurement <- rep(NA_real_, nrow(data_all) / length(proteins))
  
  for(k in seq_along(n_measurement)) n_measurement[k] <- sum(apply(data_all[k + interval, ], 1L, function(x) !any(is.na(x))))
  
  n_measurement
}

# Identify mutations whose effect varies significantly across a set of proteins.
identify_significant_variation <- function(proteins, exclude, FDR, min_n = 5L) {
  
  # proteins : integer vector
  #   A list of proteins.
  # exclude : boolean vector
  #   Indicates whether each mutation should be excluded from analysis.
  # FDR : numeric
  #   False discovery rate.
  # min_n : integer
  #   Mutations characterized in less than `min_n` proteins are excluded from analysis.
  #
  # Heteroscedastic one-way ANOVA (Welch's one-way ANOVA) is used with false discovery rate controlled
  # using the Benjamini-Hochberg procedure.
  #
  # Returns a boolean vector indicating whether a mutation's effect varies significantly across `proteins`.
  
  data_all <- multiple_comparison_data_frame(proteins) - WT_ACTIVITY
  
  interval <- 0L:(length(proteins) - 1L) * (nrow(data_all) / length(proteins))
  
  p <- rep(NA_real_, nrow(data_all) / length(proteins))
  
  for(k in 1L:(nrow(data_all) / length(proteins))) {
    
    if(is.na(exclude[k]) | exclude[k]) next
    
    Y <- unname(data_all[k + interval, ])
    Y <- Y[apply(Y, 1L, function(x) sum(is.finite(x)) > 1L), , drop = FALSE]
    if(nrow(Y) < min_n) next
    
    df <- data.frame(effect = as.vector(t(Y)), protein = as.factor(rep(1L:nrow(Y), each = NREP)))
    
    p[k] <- onewaytests::welch.test(effect ~ protein, df, na.rm = TRUE, verbose = FALSE)$p.value
  }
  
  p <= calc_BH_cutoff(p, FDR)
}

# Calculate the range of each mutation's effect across a set of proteins.
calc_mutation_effect_range <- function(proteins, min_n = 5L) {
  
  # proteins : integer vector
  #   A list of proteins.
  # min_n : integer
  #   Mutations characterized in less than `min_n` proteins are excluded from analysis.
  #
  # Returns a two-column matrix for the minimum and maximum effect of each mutation.
  
  data_all <- multiple_comparison_data_frame(proteins) - WT_ACTIVITY
  
  interval <- 0L:(length(proteins) - 1L) * (nrow(data_all) / length(proteins))
  
  range <- matrix(NA_real_, nrow(data_all) / length(proteins), 2L)
  
  for(k in 1L:(nrow(data_all) / length(proteins))) {
    
    Y <- apply(data_all[k + interval, ], 1L, mean, na.rm = TRUE)
    if(sum(!is.na(Y)) < min_n) next
    
    range[k, 1L] <- min(Y, na.rm = TRUE)
    range[k, 2L] <- max(Y, na.rm = TRUE)
  }
  
  range
}

# Parse data calculated for all mutations into subsets of data for mutations accessible from each protein and sum.
calc_per_protein_sum <- function(data, proteins) {
  vapply(parse_data_by_protein(data, proteins), sum, numeric(1L), na.rm = TRUE)
}


# Analyses.

DP <- read.table('dF.txt', stringsAsFactors = FALSE, header = TRUE)
proteins <- c(1L, 10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L)
# proteins <- c(2L, 11L, 20L, 6L, 7L, 9L, 12L, 5L, 14L) # Alt-All

# Unconditionally destructive mutations.
is_destructive <- identify_destructive_mutations_multiple_comparison(proteins, L + 0.05, 0.1, 5L, 250L)

# Mutations with significant change in effect.
is_significant <- identify_significant_variation(proteins, is_destructive, 0.1, 5L)

# Range of mutation effect.
range <- calc_mutation_effect_range(proteins)
is_sign_change <- (range[, 1L] < 0 & range[, 2L] > 0)


# Plot: Minimum and maximum effects.

plot(range[is_destructive, 1L], range[is_destructive, 2L], col = rgb(0, 0, 1, 0.5), cex = 0.5, pch = 16L,
     xlim = c(L - WT_ACTIVITY - 0.05, U - WT_ACTIVITY + 0.05), ylim = c(L - WT_ACTIVITY - 0.05, U - WT_ACTIVITY + 0.05))
points(range[is_significant & is_sign_change, 1L], range[is_significant & is_sign_change, 2L], col = rgb(1, 0, 1, 0.5), cex = 0.5, pch = 16L)
points(range[is_significant & !is_sign_change, 1L], range[is_significant & !is_sign_change, 2L], col = rgb(1, 0, 0, 0.5), cex = 0.5, pch = 16L)
points(range[!is_significant, 1L], range[!is_significant, 2L], col = rgb(0, 0, 0, 0.5), cex = 0.5, pch = 16L)
abline(v = c(L, U) - WT_ACTIVITY, h = c(L, U) - WT_ACTIVITY)
abline(v = 0, h = 0)
abline(a = 0, b = 1)


# Plot: Classification.

total_obs_by_prot <- calc_per_protein_sum(is.finite(is_destructive), proteins)
is_destructive_by_prot <- calc_per_protein_sum(is_destructive, proteins)
is_significant_by_prot <- calc_per_protein_sum(is_significant, proteins)
is_not_significant_by_prot <- calc_per_protein_sum(!is_significant, proteins)
is_sign_change_by_prot <- calc_per_protein_sum(is_significant & is_sign_change, proteins)

frac <- c(mean(is_destructive_by_prot / total_obs_by_prot),
          mean(is_not_significant_by_prot / total_obs_by_prot),
          mean((is_significant_by_prot - is_sign_change_by_prot) / total_obs_by_prot),
          mean(is_sign_change_by_prot / total_obs_by_prot))
barplot(frac)
