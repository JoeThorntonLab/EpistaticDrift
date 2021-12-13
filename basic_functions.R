# Global constants and basic functions required for specific analyses.

# Constants for the transformation relating mutant phenotypes measured at different WT activity levels.
# They are inferred in the script 'normalizing_WT_activities.R'.
L <- -2.1482 # lower bound of measurement
R <- 1.6081 # dynamic range of measurement
U <- L + R # upper bound of measurement
WT_ACTIVITY <- -0.79 # standard WT activity
N <- 0.8700 # additional transformation parameter
PLOT_RANGE <- c(L - WT_ACTIVITY - 0.05, U - WT_ACTIVITY + 0.05)

# Mutations with standard error of the mean greater than `SEM_CUTOFF` are excluded from analyses.
SEM_CUTOFF <- 0.1

# Number of measurement replicates.
NREP <- 3L

# Number of sites.
NSITE <- 76L

# Amino acids in alphabetical order.
AA <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

# `calc_BH_cutoff`: Calculate the p value threshold for Benjamini-Hochberg control of false discovery rate.
# `identify_sequence`: Identify the sequence of a given protein.
# `calc_sequence_divergence`: Calculate pairwise sequence divergence.
# `calc_sequence_divergence_multiple`: Calculate pairwise sequence divergence for a set of proteins.
# `generate_f`: Generate a function for normalizing the WT activity level.
# `change_ref_state`: Recalculate the phenotypes of mutants at a site by assigning one mutant to have the standard WT activity.
# `multiple_comparison_data_frame`: Calculate the effects of mutations relative to all WT states in a set of proteins.
# `identify_destructive_mutations_multiple_comparison`: Identify unconditionally destructive mutations.
# `identify_uncond_positive_mutations_multiple_comparison`: # Identify mutations with effect indistinguishable from the upper bound of measurement.
# `dotplot_mutation_effect_trajectory`: Generate a dot plot for a mutation's effect along trajectory.
# `dotplot_with_CI`: Generate a dotplot with empirical confidence interval or standard deviation.
# `dotplot_with_CI_for_the_mean`: Generate a dotplot with t distribution-based confidence interval for the mean or standard error of the mean.
# `barplot_with_CI`: Generate a barplot with empirical confidence interval or standard deviation.
# `parse_data_by_protein`: Parse data calculated for all mutations into subsets of data for mutations accessible from each protein.
# `calc_mutation_effects`: Calculate the effects of specified mutations on a specified protein.
# `pairwise_comparison_plot`: Plot a pairwise comparison of mutation effects.

# For temporal dynamics of epistasis.

# `calc_d`: Calculate the amount of sequence divergence for given intervals.
# `calc_dF`: Calculate the effect of each mutation on given proteins.
# `calc_noise`: Calculate the standard error of the mean for each mutation's effect on given proteins.
# `calc_ddF`: Calculate the amount of change in each mutation's effect across given intervals.
# `calc_eta`: Calculate the standard error of the mean for the amount of change in each mutation's effect across given intervals.
# `ml_constant_rate`: Maximum-likelihood inference for the rate of change in mutation effect using the constant-rate model.
# `ml_interval_rate`: Maximum-likelihood inference for systematic among-interval variation in the rate of change in mutation effect.


# Calculate the p value threshold for Benjamini-Hochberg control of false discovery rate.
calc_BH_cutoff <- function(p, FDR) {
  
  # p : numeric vector
  #   Vector of p values. Missing values are ignored.
  # FDR : numeric
  #   False discovery rate.
  #
  # Returns the p value threshold.
  
  p <- sort(p)
  k <- which(p < FDR * (1L:length(p)) / length(p))
  if(length(k) == 0) 0 else p[max(k)]
}

# Identify the sequence of a given protein.
identify_sequence <- function(i, as_vector = TRUE) {
  
  # i : integer
  #   Protein number.
  # as_vector : bool
  #   If `TRUE`, returns a vector of length `NSITE`; otherwise returns a string.
  
  seq <- DP$WTAA[DP$PROT == i]
  seq <- seq[20L * (1L:(length(seq) / 20L))]
  if(as_vector) seq else paste(seq, collapse = '')
}

# Calculate pairwise sequence divergence.
calc_sequence_divergence <- function(i, j, return_nsub = FALSE) {
  
  # i, j : integers
  #   Protein numbers.
  # return_nsub : bool
  #   If `TRUE`, returns the number of different sites.
  #   If `FALSE`, returns percent sequence divergence.
  
  WTi <- DP$WTAA[DP$PROT == i]
  WTi <- WTi[20L * (1L:(length(WTi) / 20L))]
  
  WTj <- DP$WTAA[DP$PROT == j]
  WTj <- WTj[20L * (1L:(length(WTj) / 20L))]
  
  if(return_nsub) {
    sum(WTi != WTj)
  } else {
    100 * sum(WTi != WTj) / length(WTi)
  }
}

# Calculate pairwise sequence divergence for a set of proteins.
calc_sequence_divergence_multiple <- function(proteins) {
  
  # proteins : integer vector
  #   A list of protein numbers.
  #
  # Returns a square, upper triangular matrix containing all pairwise sequence divergences. 
  
  seq_div_matrix <- matrix(NA_real_, length(proteins), length(proteins))
  
  for(k in 1L:(length(proteins) - 1L)) {
    for(l in (k + 1L):length(proteins)) {
      
      seq_div_matrix[k, l] <- calc_sequence_divergence(proteins[k], proteins[l])
    }
  }
  
  seq_div_matrix
}

# Generate a function for normalizing the WT activity level.
generate_f <- function(ref_activity) {
  
  # ref_activity : numeric
  #   The WT activity level at which mutant phenotypes are originally measured.
  #
  # The generated function takes mutant phenotype measured at `ref_activity` as argument and
  # returns the phenotype expected when measured at `WT_ACTIVITY`.
  
  if(ref_activity < L | ref_activity > (L + R)) {
    warning('WT activity out of range; NULL returned\n')
    return(NULL)
  }
  
  norm_ref <- (ref_activity - L) / R
  norm_WT <- (WT_ACTIVITY - L) / R
  
  if(norm_ref == norm_WT) return(function(x) x) # no transformation necessary
  
  if(norm_ref < norm_WT) { # concave function
    
    a <- (norm_ref ^ N) * (1 - norm_WT) / (norm_WT - norm_ref ^ N)
    
    f_atomic <- function(x) {
      if(!is.finite(x)) {
        NA_real_
      } else if(x < L | x > (L + R)) {
        x
      } else {
        R * ((a + 1) * ((x - L) / R) ^ N) / (a + ((x - L) / R) ^ N) + L
      }
    }
    
  } else { # convex function
    
    a <- (norm_WT ^ N) * (1 - norm_ref) / (norm_ref - norm_WT ^ N)
    
    f_atomic <- function(x) {
      if(!is.finite(x)) {
        NA_real_
      } else if(x < L | x > (L + R)) {
        x
      } else {
        R * ((a * ((x - L) / R)) / (a + 1 - ((x - L) / R))) ^ (1 / N) + L
      }
    }
  } 
  
  function(x) vapply(x, f_atomic, numeric(1L))
}

# Recalculate the phenotypes of mutants at a site by assigning one mutant to have the standard WT activity.
change_ref_state <- function(data, ref_state, phenotype_only = TRUE) {
  
  # data : data frame
  #   A segment of mutant information (`D`) for mutants at a site in a protein.
  # ref_state : character
  #   Mutant amino acid state to be assigned the standard WT activity level.
  # phenotype_only : bool
  #   If `TRUE`, only returns the phenotypes; otherwise return the full modified data frame.
  #
  # The mutant with the new reference state is marked with NA.
  
  # If `ref_state` equals the WT state in `data`, no modification is necessary.
  if(unique(data$WTAA) == ref_state) {
    
    data[data$MTAA == ref_state, grep('F', colnames(data))] <- NA_real_
    
    if(phenotype_only) return(as.matrix(data[, grep('F', colnames(data))])) else return(data)
  }
  
  MTAA <- DP$MTAA[1L:20L]
  original_P <- as.matrix(data[, grep('F', colnames(data))]) # input phenotype data
  new_P <- matrix(NA_real_, length(MTAA), NREP) # phenotype data based on the new reference state
  
  delta <- 0.15 # reference phenotype outside the interval [L + delta, U - delta] is considered highly sensitive to measurement error
  delta_L <- 0.5 # when the reference phenotype is too close to L, the effect of a mutant with phenotype < L + delta_L is considered indeterminable (marked by NA)
  #delta_U <- 0.2 # when the reference phenotype is too close to U, the effect of a mutant with phenotype > U - delta_U is considered indeterminable (marked by NA)
  delta_U <- 0.1 # when the reference phenotype is too close to U, the effect of a mutant with phenotype > U - delta_U is considered indeterminable (marked by NA)
  
  for(j in 1L:NREP) { # for each replicate
    
    ref_activity <- original_P[which(MTAA == ref_state), j]
    
    if(!is.finite(ref_activity)) next # if measurement is missing
    
    if(ref_activity < L + delta) { # when `ref_activity` is too close to L
      
      f <- generate_f(if(ref_activity > L + 0.01) ref_activity else L + 0.01) # to avoid numerical problems
      
      for(k in 1L:20L) {
        
        if(!is.finite(original_P[k, j])) next
        
        if(k == which(MTAA == ref_state)) next # this is not a mutation
        
        # when the mutant phenotype is also close to L, the correction procedure is highly sensitive to measurement error; it is therefore not made
        if(original_P[k, j] > L + delta_L) new_P[k, j] <- f(original_P[k, j])
      }
      
    } else if(ref_activity > U - delta) { # when `ref_activity` is close to U
      
      f <- generate_f(if(ref_activity < U - 0.01) ref_activity else U - 0.01) # to avoid numerical problems
      
      for(k in 1L:20L) {
        
        if(!is.finite(original_P[k, j])) next
        
        if(k == which(MTAA == ref_state)) next # this is not a mutation
        
        # when the mutant phenotype is also close to U, the correction procedure is highly sensitive to measurement error; it is therefore not made
        if(original_P[k, j] < U - delta_U) new_P[k, j] <- f(original_P[k, j])
      }
      
    } else { # when `ref_activity` not too close to L or U
      
      f <- generate_f(ref_activity)
      
      new_P[, j] <- f(original_P[, j])
      
      new_P[which(MTAA == ref_state), j] <- NA_real_ # this is not a mutation
    }
  } # end of the loop over replicates
  
  if(phenotype_only) {
    new_P
    
  } else {
    
    data[, grep('F', colnames(data))] <- new_P
    data$WTAA <- ref_state
    data
  }
}

# Calculate the effects of mutations relative to all WT states in a set of proteins.
multiple_comparison_data_frame <- function(proteins, phenotype_only = TRUE, stack = TRUE) {
  
  # proteins : integer vector
  #   A list of proteins.
  # phenotype_only : bool
  #   If `TRUE`, returns only the phenotype matrix. Otherwise returns the full data frame.
  # stack : bool
  #   If `TRUE`, data for every protein is stacked in a single data frame or matrix.
  #   Otherwise, a list of data frames or matrices is returned.
  #
  # Mutations with standard error for the effect greater than `SEM_CUTOFF` are marked with NA.
  
  LEN <- max(DP$SITE)
  MTAA <- DP$MTAA[1L:20L]
  
  DP_proteins <- list()
  WT_proteins <- matrix('', length(proteins), LEN)
  
  for(i in 1L:length(proteins)) {
    DP_proteins[[i]] <- DP[DP$PROT == proteins[i], ]
    WT_proteins[i, ] <- DP_proteins[[i]]$WTAA[20L * 1L:LEN]
  }
  
  DP_expanded <- rep(list(NULL), length(proteins))
  
  for(k in 1L:LEN) {
    
    ref_states <- setdiff(sort(unique(WT_proteins[, k])), '-')
    
    if(length(ref_states) > 0L) {
      
      for(i in 1L:length(proteins)) {
        
        data <- DP_proteins[[i]][20L * (k - 1L) + 1L:20L, ]
        
        for(j in 1L:length(ref_states)) DP_expanded[[i]] <- rbind(DP_expanded[[i]], change_ref_state(data, ref_states[j], phenotype_only))
      }
      
    } else { # this site is a deletion in all proteins
      
      for(i in 1L:length(proteins)) {
        
        data <- DP_proteins[[i]][20L * (k - 1L) + 1L:20L, ]
        if(phenotype_only) data <- as.matrix(data[, grepl('F', colnames(data))])
        
        DP_expanded[[i]] <- rbind(DP_expanded[[i]], data)
      }
    }
  }
  
  # Marking mutations with large standard error.
  for(i in 1L:length(proteins)) {
    
    if(phenotype_only) {
      
      SEM <- apply(DP_expanded[[i]], 1L, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))))
      DP_expanded[[i]][which(SEM > SEM_CUTOFF), ] <- NA_real_
      
    } else {
      
      SEM <- apply(DP_expanded[[i]][, grep('F', colnames(DP_expanded[[i]]))], 1L, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))))
      DP_expanded[[i]][which(SEM > SEM_CUTOFF), grep('F', colnames(DP_expanded[[i]]))] <- NA_real_
    }
  }
  
  # Stacking data frames.
  if(stack) {
    
    for(i in 2L:length(proteins)) DP_expanded[[1L]] <- rbind(DP_expanded[[1L]], DP_expanded[[i]])
    DP_expanded[[1L]]
    
  } else {
    DP_expanded 
  }
}

# Identify unconditionally destructive mutations.
identify_destructive_mutations_multiple_comparison <- function(proteins, LB, FDR, min_n = 5L, nbootstrap = 250L) {
  
  # proteins : integer vector
  #   A list of proteins.
  # LB : numeric
  #   Lower bound of measurement.
  # FDR : numeric
  #   False discovery rate.
  # min_n : integer
  #   Mutations measured in less than `min_n` proteins are excluded from analysis.
  # nbootstrap : integer
  #   Number of null likelihood-ratio statistic to be simulated.
  #
  # Mutants with phenotype not significantly greater than `LB` across the proteins are identified by
  # the following likelihood-ratio test:
  #   Null hypothesis: Mutant has phenotype `LB` in all proteins.
  #   Alternative hypothesis: Mutant has phenotype greater than `LB` in some protein(s).
  # Parametric bootstrap likelihood-ratio test is performed.
  # False discovery rate is controlled at `FDR` using the Benjamini-Hochberg procedure.
  # Mutations measured in less than `min_n` proteins are excluded.
  
  # Calculate the likelihood-ratio statistic.
  calc_LRS <- function(Y) {
    
    # The rows of `Y` correspond to `proteins` and columns to measurement replicates.
    # There should be no missing values.
    
    ncol(Y) * sum(apply(Y, 1L, function(x) {if(mean(x) < LB) 0 else log(sum((x - LB) ^ 2L) / sum((x - mean(x)) ^ 2L))}))
  }
  
  # Perform the likelihood-ratio test.
  perform_LRT <- function(Y, nbootstrap) {
    
    LRS <- calc_LRS(Y)
    
    # Simulating the null distribution of likelihood-ratio statistic.
    Y <- Y - apply(Y, 1L, function(x) if(mean(x) < LB) 0 else mean(x) - LB)
    Y_bootstrap <- t(apply(Y, 1L, function(x) rnorm(NREP * nbootstrap, mean(x), sd(x))))
    null_distr <- vapply(1L:nbootstrap, function(i) calc_LRS(Y_bootstrap[, 1L:NREP + NREP * (i - 1)]), numeric(1L))
    
    sum(null_distr >= LRS) / nbootstrap
  }
  
  data_all <- multiple_comparison_data_frame(proteins)
  
  interval <- 0L:(length(proteins) - 1L) * (nrow(data_all) / length(proteins))
  
  p <- rep(NA_real_, nrow(data_all) / length(proteins))
  
  for(k in 1L:(nrow(data_all) / length(proteins))) {
    
    Y <- unname(data_all[k + interval, ])
    Y <- Y[apply(Y, 1L, function(x) !any(is.na(x))), , drop = FALSE]
    if(nrow(Y) < min_n) next
    
    p[k] <- perform_LRT(Y, nbootstrap)
  }
  
  p > calc_BH_cutoff(p, FDR)
}

# Identify mutations with effect indistinguishable from the upper bound of measurement.
identify_uncond_positive_mutations_multiple_comparison <- function(proteins, UB, FDR, min_n = 5L, nbootstrap = 250L) {
  
  # proteins : integer vector
  #   A list of proteins.
  # UB : numeric
  #   Upper bound of measurement.
  # FDR : numeric
  #   False discovery rate.
  # min_n : integer
  #   Mutations measured in less than `min_n` proteins are excluded from analysis.
  # nbootstrap : integer
  #   Number of null likelihood-ratio statistic to be simulated.
  #
  # Uses the same machinery as `identify_destructive_mutations_multiple_comparison` by inverting the sign of phenotype.
  
  # Calculates the likelihood-ratio statistic.
  calc_LRS <- function(Y) {
    
    # The rows of `Y` correspond to `proteins` and columns to measurement replicates.
    # There should be no missing values.
    
    ncol(Y) * sum(apply(Y, 1L, function(x) {if(mean(x) < LB) 0 else log(sum((x - LB) ^ 2L) / sum((x - mean(x)) ^ 2L))}))
  }
  
  # Perform the likelihood-ratio test.
  perform_LRT <- function(Y, nbootstrap) {
    
    LRS <- calc_LRS(Y)
    
    # Simulating the null distribution of likelihood-ratio statistic.
    Y <- Y - apply(Y, 1L, function(x) if(mean(x) < LB) 0 else mean(x) - LB)
    Y_bootstrap <- t(apply(Y, 1L, function(x) rnorm(NREP * nbootstrap, mean(x), sd(x))))
    null_distr <- vapply(1L:nbootstrap, function(i) calc_LRS(Y_bootstrap[, 1L:NREP + NREP * (i - 1)]), numeric(1L))
    
    sum(null_distr >= LRS) / nbootstrap
  }
  
  # Inverting the sign of mutant phenotypes.
  data_all <- -multiple_comparison_data_frame(proteins)
  LB <- -UB
  
  interval <- 0L:(length(proteins) - 1L) * (nrow(data_all) / length(proteins))
  
  p <- rep(NA_real_, nrow(data_all) / length(proteins))
  
  for(k in 1L:(nrow(data_all) / length(proteins))) {
    
    Y <- unname(data_all[k + interval, ])
    Y <- Y[apply(Y, 1L, function(x) !any(is.na(x))), , drop = FALSE]
    if(nrow(Y) < min_n) next
    
    p[k] <- perform_LRT(Y, nbootstrap)
  }
  
  p > calc_BH_cutoff(p, FDR)
}

# Generate a dot plot for a mutation's effect along trajectory.
dotplot_mutation_effect_trajectory <- function(trajectory, mutation, return_dF = FALSE, x = NULL) {
  
  # trajectory : integer vector
  #   An ordered list of proteins.
  # mutation : character vector
  #   A vector of (site, WT AA, mutant AA)
  #
  # Standard error of the mean for the effect is also plotted.
  
  data <- multiple_comparison_data_frame(trajectory, phenotype_only = FALSE)
  data_info <- data[data$PROT == trajectory[1L], c('SITE', 'WTAA', 'MTAA')]
  
  mutation <- which(data_info$SITE == as.integer(mutation[1L]) & data_info$WTAA == mutation[2L] & data_info$MTAA == mutation[3L])
  
  data <- as.matrix(data[, grep('F', colnames(data))]) - WT_ACTIVITY
  interval <- 0L:(length(trajectory) - 1L) * (nrow(data) / length(trajectory))
  
  data <- data[interval + mutation, ]
  
  data_mean <- apply(data, 1L, mean, na.rm = TRUE)
  sem <- apply(data, 1L, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))))
  
  data_mean[is.na(sem)] <- NA_real_ # removing data for which confidence interval cannot be calculated
  
  if(is.null(x)) x <- 1L:length(data_mean)
  
  plot(x, data_mean, ylim = c(-1.4, 0.4), main = paste('mutation', mutation), cex = 1.6, yaxt = 'n')
  axis(side = 2L, at = c(-1.2, -0.8, -0.4, 0, 0.4))
  
  for(i in 1L:length(data_mean)) {
    arrows(x0 = x[i], y0 = data_mean[i] - sem[i], x1 = x[i], y1 = data_mean[i] + sem[i], angle = 90, code = 3, length = 0.1)
  }
  
  abline(h = c(L, U) - WT_ACTIVITY)
  
  if(return_dF) data
}

# Generate a dotplot with empirical confidence interval or standard deviation.
dotplot_with_CI <- function(X, Y, CI, plot_sd = FALSE, metric = mean, ...) {
  
  # X, Y : numeric vector, matrix, or list
  #   If a matrix, rows correspond to data points and columns to replicate values.
  #   If a vector, elements correspond to data points. No confidence interval is plotted.
  # CI : numeric
  #   Confidence interval (between 0 and 1).
  # plot_sd : bool
  #   If `TRUE`, standard deviation is plotted instead of confidence interval.
  # ... : optional additional parameters
  #   Passed to `plot` for managing graphic parameters.
  
  if(is.vector(X) & !is.list(X)) X <- vapply(X, as.list, list(1))
  if(is.vector(Y) & !is.list(Y)) Y <- vapply(Y, as.list, list(1))
  if(is.matrix(X)) X <- split(X, rep(1L:nrow(X), ncol(X)))
  if(is.matrix(Y)) Y <- split(Y, rep(1L:nrow(Y), ncol(Y)))
  
  X_mean <- vapply(X, metric, numeric(1L), na.rm = TRUE)
  Y_mean <- vapply(Y, metric, numeric(1L), na.rm = TRUE)
  
  if(plot_sd) {
    
    X_CI <- cbind(X_mean - vapply(X, sd, numeric(1L), na.rm = TRUE), X_mean + vapply(X, sd, numeric(1L), na.rm = TRUE))
    Y_CI <- cbind(Y_mean - vapply(Y, sd, numeric(1L), na.rm = TRUE), Y_mean + vapply(Y, sd, numeric(1L), na.rm = TRUE))
    
  } else {
    
    X_CI <- t(vapply(X, function(x) if(length(x) > 1L) quantile(x, c((1 - CI) / 2, 1 - (1 - CI) / 2), na.rm = TRUE) else c(NA_real_, NA_real_), numeric(2L)))
    Y_CI <- t(vapply(Y, function(x) if(length(x) > 1L) quantile(x, c((1 - CI) / 2, 1 - (1 - CI) / 2), na.rm = TRUE) else c(NA_real_, NA_real_), numeric(2L)))
  }
  
  plot(X_mean, Y_mean, ...)
  
  for(k in 1L:length(X)) {
    
    if(all(is.finite(Y_CI[k, ]))) arrows(x0 = X_mean[k], y0 = Y_CI[k, 1L], x1 = X_mean[k], y1 = Y_CI[k, 2L], angle = 90, code = 3, length = 0.04)
    if(all(is.finite(X_CI[k, ]))) arrows(x0 = X_CI[k, 1L], y0 = Y_mean[k], x1 = X_CI[k, 2L], y1 = Y_mean[k], angle = 90, code = 3, length = 0.04)
  }
}

# Generate a dotplot with t distribution-based confidence interval for the mean or standard error of the mean.
dotplot_with_CI_for_the_mean <- function(X, Y, CI, plot_sem = FALSE, metric = mean, ...) {
  
  # X, Y : numeric vector or matrix
  #   If a matrix, rows correspond to data points and columns to replicate values.
  #   If a vector, elements correspond to data points. No confidence interval is plotted.
  # CI : numeric
  #   Confidence interval (between 0 and 1).
  # plot_sem : bool
  #   If `TRUE`, standard error of the mean is plotted instead of confidence interval.
  # ... : optional additional parameters
  #   Passed to `plot` for managing graphic parameters.
  
  calc_CI <- function(x) if(sum(!is.na(x)) < 2) NA_real_ else qt(1 - (1 - CI) / 2, sum(!is.na(x)) - 1L) * sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  
  if(is.vector(X) & !is.list(X)) X <- vapply(X, as.list, list(1))
  if(is.vector(Y) & !is.list(Y)) Y <- vapply(Y, as.list, list(1))
  if(is.matrix(X)) X <- split(X, rep(1L:nrow(X), ncol(X)))
  if(is.matrix(Y)) Y <- split(Y, rep(1L:nrow(Y), ncol(Y)))
  
  X_mean <- vapply(X, metric, numeric(1L), na.rm = TRUE)
  Y_mean <- vapply(Y, metric, numeric(1L), na.rm = TRUE)
  
  if(plot_sem) {
    
    X_CI <- vapply(X, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))), 0)
    Y_CI <- vapply(Y, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))), 0)
    
  } else {
    
    X_CI <- vapply(X, calc_CI, 0)
    Y_CI <- vapply(Y, calc_CI, 0)
  }
  
  plot(X_mean, Y_mean, ...)
  
  for(k in 1L:length(X)) {
    
    if(is.finite(Y_CI[k])) arrows(x0 = X_mean[k], y0 = Y_mean[k] - Y_CI[k], x1 = X_mean[k], y1 = Y_mean[k] + Y_CI[k], angle = 90, code = 3, length = 0.03)
    if(is.finite(X_CI[k])) arrows(x0 = X_mean[k] - X_CI[k], y0 = Y_mean[k], x1 = X_mean[k] + X_CI[k], y1 = Y_mean[k], angle = 90, code = 3, length = 0.03)
  }
}

# Generate a barplot with empirical confidence interval or standard deviation.
barplot_with_CI <- function(Y, CI, plot_sd = FALSE, ...) {
  
  # Y : numeric matrix
  #   Rows correspond to data points and columns to replicate values.
  # CI : numeric
  #   Confidence interval (between 0 and 1).
  # plot_sd : bool
  #   If `TRUE`, standard deviation is plotted instead of confidence interval.
  # ... : optional additional parameters
  #   Passed to `plot` for managing graphic parameters.
  
  Y <- as.matrix(Y)
  
  Y_mean <- apply(Y, 1L, mean, na.rm = TRUE)
  
  if(plot_sd) {
   
    Y_CI <- cbind(Y_mean - apply(Y, 1L, sd, na.rm = TRUE), Y_mean + apply(Y, 1L, sd, na.rm = TRUE))
  
  } else {
    
    Y_CI <- t(apply(Y, 1L, quantile, probs = c((1 - CI) / 2, 1 - (1 - CI) / 2)))
  }
  
  bx <- barplot(Y_mean, ...)
  
  for(k in 1L:nrow(Y)) {
    
    if(all(is.finite(Y_CI[k, ]))) arrows(x0 = bx[k, 1L], y0 = Y_CI[k, 1L], x1 = bx[k, 1L], y1 = Y_CI[k, 2L], angle = 90, code = 3, length = 0.04)
  }
}

# For each mutation at a site, plot the range of its effect across a specified set of proteins.
plot_mutation_effect_range <- function(proteins, site, WTAA = NULL) {
  
  data_all <- multiple_comparison_data_frame(proteins, phenotype_only = FALSE, stack = FALSE)
  
  dF <- sapply(data_all, function(data_i) {
    
    apply(data_i[data_i$SITE == site, grep('F', colnames(data_i))], 1L, mean) - WT_ACTIVITY
  })
  
  mutant_info <- data_all[[1L]][data_all[[1L]]$SITE == site, c('WTAA', 'MTAA')]
  
  if(!is.null(WTAA)) {
    dF <- dF[mutant_info$WTAA %in% WTAA, ]
    mutant_info <- mutant_info[mutant_info$WTAA %in% WTAA, ]
  }
  
  dF_range <- apply(dF, 1L, function(x) if(all(is.na(x))) c(NA_real_, NA_real_) else range(x, na.rm = TRUE))
  
  plot(1, type = 'n', xlim = c(1L, ncol(dF_range)), ylim = c(L - WT_ACTIVITY - 0.05, U - WT_ACTIVITY + 0.05), xaxt = 'n')
  rect(1L:ncol(dF_range) - 0.25, dF_range[1L, ], 1L:ncol(dF_range) + 0.25, dF_range[2L, ])
  axis(side = 1L, at = 1L:ncol(dF_range), labels = paste0(mutant_info$WTAA, mutant_info$MTAA))
  abline(h = c(L - WT_ACTIVITY, 0, U - WT_ACTIVITY))
}

# Parse data calculated for all mutations into subsets of data for mutations accessible from each protein.
parse_data_by_protein <- function(data, proteins, retain = NULL) {
  
  # Mutations defined by `multiple_comparison_data_frame` are all mutations whose reference state
  # is found as a WT state in one of the `proteins`. Therefore, some of these mutations are not
  # accessible from a given protein. In this function, `data` is a vector of data calculated
  # for every mutation. This function identifies mutations that are accessible from a given
  # protein and identifies their `data`. 
  
  # data : numeric vector
  #   Data calculated for each mutation defined by `multiple_comparison_data_frame` on `proteins`.
  # proteins : integer vector
  #   Proteins for which `data` is calculated.
  # retain : bool
  #   An optional vector specifying the subset of mutations for which `data` is calculated.
  #   If provided, `data` must be calculated only for mutations marked as `TRUE` in `retain`.
  #
  # Returns a list of length equal to the number of proteins. Each element of the list contains
  # `data` items for all mutants accessible from that protein.
  
  # All (site, WT AA) pairs in `proteins`.
  site_WTAA <- multiple_comparison_data_frame(proteins, FALSE, FALSE)[[1L]][, c('SITE', 'WTAA')]
  site_WTAA <- paste0(site_WTAA$SITE, site_WTAA$WTAA)
  
  # Mutations to be retained.
  if(is.null(retain)) retain <- rep(TRUE, length(site_WTAA))
  
  data_by_prot <- rep(list(NULL), length(proteins))
  
  for(i in seq_along(proteins)) {
    
    # (site, wild type AA) pairs for protein `i`.
    site_WTAA_prot <- unique(paste0(DP$SITE[DP$PROT == proteins[i]], DP$WTAA[DP$PROT == proteins[i]]))
    
    # Identifying the indices of point mutations for protein `i`.
    is_in_prot <- site_WTAA %in% site_WTAA_prot
    
    data_by_prot[[i]] <- data[is_in_prot[retain]]
  }
  
  data_by_prot
}

# Calculate the effects of specified mutations on a specified protein.
calc_mutation_effects <- function(mutations, background, return_mean = FALSE) {
  
  # mutations : data frame
  #   Must contain columns named 'SITE', 'WTAA', and 'MTAA' recording corresponding information for each mutation.
  # background : integer
  #   The protein on which the mutation effects are calculated.
  # return_mean : bool
  #   If `TRUE`, the mean effect is returned. Otherwise replicate measurements are returned in a matrix.
  
  effects <- matrix(NA_real_, nrow(mutations), length(grep('F', colnames(DP))))
  
  for(i in 1L:nrow(mutations)) {
    
    F_i <- DP[DP$PROT == background & DP$SITE == mutations$SITE[i], ]
    F_i <- change_ref_state(F_i, mutations$WTAA[i], phenotype_only = FALSE)
    
    effects[i, ] <- as.numeric(F_i[F_i$MTAA == mutations$MTAA[i], grep('F', colnames(F_i))]) - WT_ACTIVITY
  }
  
  SEM <- apply(effects, 1L, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))))
  effects[SEM > SEM_CUTOFF, ] <- NA_real_
  
  if(return_mean) {
    apply(effects, 1L, function(x) if(sum(is.finite(x)) < 2L) NA_real_ else mean(x, na.rm = TRUE))
  } else {
    effects
  }
}

# Plot a pairwise comparison of mutation effects.
pairwise_comparison_plot <- function(i, j, FDR) {
  
  # i, j : integer
  #   Two proteins.
  # FDR : numeric
  #   False discovery rate for identifying destructive mutations and mutations with varible effect.
  
  data_all <- multiple_comparison_data_frame(c(i, j)) - WT_ACTIVITY
  Mi <- apply(data_all[1L:(nrow(data_all) / 2L), ], 1L, mean)
  Mj <- apply(data_all[-(1L:(nrow(data_all) / 2L)), ], 1L, mean)
  
  is_destructive <- identify_destructive_mutations_multiple_comparison(c(i, j), L + 0.1, FDR = FDR, min_n = 2L)
  
  is_significant <- identify_significant_variation(c(i, j), FDR = FDR, min_n = 2L)
  
  cat(c(sum(is_destructive, na.rm = TRUE), sum(!is_destructive & is_significant, na.rm = TRUE),
        sum(!is_destructive & !is_significant, na.rm = TRUE)) / sum(is.finite(is_destructive)))
  
  range <- c(L - 0.05, U + 0.05) - WT_ACTIVITY
  
  plot(Mi[!is_destructive & !is_significant], Mj[!is_destructive & !is_significant], xlim = range, ylim = range,
       col = rgb(0, 0, 0, 0.5), cex = 0.6, pch = 16L)
  points(Mi[is_destructive], Mj[is_destructive], col = rgb(0, 0, 1, 0.5), cex = 0.6, pch = 16L)
  points(Mi[!is_destructive & is_significant], Mj[!is_destructive & is_significant], col = rgb(1, 0, 0, 0.5), cex = 0.6, pch = 16L)
  
  abline(v = c(L, U) - WT_ACTIVITY, h = c(L, U) - WT_ACTIVITY)
  abline(v = 0, h = 0)
  abline(a = 0, b = 1)
  
  # For modifying in illustrator.
  #plot(Mi[!is_destructive & !is_significant], Mj[!is_destructive & !is_significant], xlim = range, ylim = range,
  #     col = rgb(0, 0, 0, 0.5), cex = 0.6, pch = 16L)
  #points(Mi[is_destructive], Mj[is_destructive], col = rgb(0, 0, 1, 0.5), cex = 0.6, pch = 16L)
  #
  #plot(Mi[!is_destructive & !is_significant], Mj[!is_destructive & !is_significant], xlim = range, ylim = range,
  #     col = rgb(0, 0, 0, 0.5), cex = 0.6, pch = 16L)
  #points(Mi[!is_destructive & is_significant], Mj[!is_destructive & is_significant], col = rgb(1, 0, 0, 0.5), cex = 0.6, pch = 16L)
}

# Plot a pairwise comparison of mutation effects for a specified subset of mutations.
pairwise_comparison_plot_subset <- function(i, j, proteins, FDR, retain = NULL, separate = FALSE) {
  
  # i, j : integer
  #   Two proteins.
  # proteins : integer vector
  #   A list of proteins that defines the total mutation set.
  # FDR : numeric
  #   False discovery rate for identifying destructive mutations and mutations with varible effect.
  # retain : boolean or integer vector
  #   Boolean: Indicates whether a mutation should be included in the plot.
  #   Integer: Indices of mutations to be included.
  # separate : boolean
  #   If TRUE, separately plots mutations with/without significant differences in effect.
  
  proteins <- unique(c(proteins, i, j))
  dF <- multiple_comparison_data_frame(proteins, TRUE, FALSE)
  dFi <- dF[[which(proteins == i)]][retain, ] - WT_ACTIVITY
  dFj <- dF[[which(proteins == j)]][retain, ] - WT_ACTIVITY
  
  # Test of significant difference.
  p <- rep(NA_real_, nrow(dFi))
  for(k in 1L:nrow(dFi)) {
    
    x <- dFi[k, ]
    y <- dFj[k, ]
    if(sum(is.finite(x)) < 2L | sum(is.finite(y)) < 2L) next
    
    p[k] <- t.test(x, y, alternative = 'two.sided')$p.value
  }
  p_threshold <- calc_BH_cutoff(p, FDR)
  
  Mi <- apply(dFi, 1L, mean, na.rm = TRUE)
  Mj <- apply(dFj, 1L, mean, na.rm = TRUE)
  
  cat('r2 of dF:', cor(Mi, Mj, use = 'pairwise.complete.obs') ^ 2L, '\n')
  
  range <- c(L - 0.05, U + 0.05) - WT_ACTIVITY
  
  plot(Mi[p <= p_threshold], Mj[p <= p_threshold],
       xlim = range, ylim = range, col = rgb(1, 0, 0, 0.5), cex = 0.6, pch = 16L)
  abline(v = c(L, U) - WT_ACTIVITY, h = c(L, U) - WT_ACTIVITY)
  abline(v = 0, h = 0)
  abline(a = 0, b = 1)
  
  if(separate) {
    plot(Mi[p > p_threshold], Mj[p > p_threshold],
         xlim = range, ylim = range, col = rgb(0, 0, 0, 0.5), cex = 0.6, pch = 16L)
  } else {
    points(Mi[p > p_threshold], Mj[p > p_threshold], col = rgb(0, 0, 0, 0.5), cex = 0.6, pch = 16L)
  }
}

# Generate a dotplot with control of density.
dotplot_with_density_control <- function(x, y, xlim, ylim, xbin, ybin, max_count, ...) {
  
  x <- as.vector(x)
  y <- as.vector(y)
  retain <- is.finite(x) & is.finite(y)
  x <- x[retain]
  y <- y[retain]
  
  cat('Total number of points =', length(x))
  
  xgrid <- seq(xlim[1L], xlim[2L], by = xbin)
  ygrid <- seq(ylim[1L], ylim[2L], by = ybin)
  
  xi <- findInterval(x, xgrid)
  yi <- findInterval(y, ygrid)
  
  retain <- list()
  count <- 1L
  
  for(i in 0L:length(xgrid)) {
    for(j in 0L:length(ygrid)) {
      
      index <- which(xi == i & yi == j)
      if(length(index) == 0L) next
      
      if(length(index) <= max_count)
        retain[[count]] <- index
      else
        retain[[count]] <- sample(index, max_count)
      
      count <- count + 1L
    }
  }
  
  x <- x[unlist(retain)]
  y <- y[unlist(retain)]
  
  cat(' reduced to', length(x), '\n')
  plot(x, y, xlim = xlim, ylim = ylim, ...)
}


# Calculate the amount of sequence divergence for given intervals.
calc_d <- function(intervals, return_nsub = FALSE) {
  
  # intervals : two-column integer matrix
  #   Each row specifies the ancestral and descendant protein flanking an interval.
  # return_nsub : bool
  #   If `TRUE`, returns the number of substitutions. Otherwise retruns percent sequence divergence.
  
  apply(intervals, 1L, function(x) calc_sequence_divergence(x[1L], x[2L], return_nsub))
}

# Calculate the effect of each mutation on given proteins.
calc_dF <- function(proteins) {
  
  # proteins : integer vector
  #   A list of proteins.
  #
  # Mutation set defined by `multiple_comparison_data_frame` are analyzed.
  # Returns a matrix whose rows correspond to mutations and columns to proteins.
  
  if(length(proteins) == 1L) {
    apply(DP[DP$PROT == proteins, grep('F', colnames(DP))], 1L, mean) - WT_ACTIVITY
  } else {
    data_all <- multiple_comparison_data_frame(proteins) - WT_ACTIVITY
    matrix(apply(data_all, 1L, mean), ncol = length(proteins))
  }
}

# Calculate the standard error of the mean for each mutation's effect on given proteins.
calc_noise <- function(proteins) {
  
  # proteins : integer vector
  #   A list of proteins.
  #
  # Mutation set defined by `multiple_comparison_data_frame` are analyzed.
  # Returns a matrix whose rows correspond to mutations and columns to proteins.
  
  data_all <- multiple_comparison_data_frame(proteins) - WT_ACTIVITY
  matrix(apply(data_all, 1L, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))),
         ncol = length(proteins))
}

# Calculate the amount of change in each mutation's effect across given intervals.
calc_ddF <- function(intervals, remove_out_of_range = TRUE) {
  
  # intervals : two-column integer matrix
  #   Each row specifies the ancestral and descendant protein flanking an interval.
  # `remove_out_of_range`: bool
  #   If `TRUE`, intervals along which a mutation's effect lies outside measurement bounds are excluded.
  #
  # Mutation set defined by `multiple_comparison_data_frame` on all proteins in the interval are analyzed.
  # Returns a matrix whose rows correspond to mutations and columns to intervals.
  
  proteins <- unique(as.vector(intervals))
  data_all <- multiple_comparison_data_frame(proteins) - WT_ACTIVITY
  
  dF <- matrix(apply(data_all, 1L, mean), ncol = length(proteins))
  ddF <- matrix(NA_real_, nrow(dF), ncol(dF) - 1L)
  
  for(k in 1L:nrow(intervals)) {
    anc_dF <- dF[, proteins == intervals[k, 1L]]
    der_dF <- dF[, proteins == intervals[k, 2L]]
    ddF_k <- der_dF - anc_dF
    
    if(remove_out_of_range) {
      ddF_k[(der_dF < L - WT_ACTIVITY + 0.05) & (anc_dF < L - WT_ACTIVITY + 0.05)] <- NA_real_
      ddF_k[(der_dF > U - WT_ACTIVITY - 0.05) & (anc_dF > U - WT_ACTIVITY - 0.05)] <- NA_real_
    }
    
    ddF[, k] <- ddF_k
  }
  
  ddF
}

# Calculate the standard error of the mean for the amount of change in each mutation's effect across given intervals.
calc_eta <- function(intervals) {
  
  # intervals : two-column integer matrix
  #   Each row specifies the ancestral and descendant protein flanking an interval.
  #
  # Mutation set defined by `multiple_comparison_data_frame` on all proteins in the interval are analyzed.
  # Returns a matrix whose rows correspond to mutations and columns to intervals.
  
  proteins <- unique(as.vector(intervals))
  data_all <- multiple_comparison_data_frame(proteins) - WT_ACTIVITY
  
  noise_var <- matrix(apply(data_all, 1L, function(x) var(x, na.rm = TRUE) / sum(is.finite(x))), ncol = length(proteins))
  eta <- matrix(NA_real_, nrow(noise_var), ncol(noise_var) - 1L)
  
  for(k in 1L:nrow(intervals)) {
    anc_noise_var <- noise_var[, proteins == intervals[k, 1L]]
    der_noise_var <- noise_var[, proteins == intervals[k, 2L]]
    eta[, k] <- sqrt(anc_noise_var + der_noise_var)
  }
  
  eta
}

# Maximum-likelihood inference using the constant-rate model
ml_constant_rate <- function(d, ddF, eta, return_mll = FALSE) {
  
  # obtains the mle and mll for mutation 'i'
  calc_mle <- function(i) {
    
    if(sum(!is.na(ddF[i, ])) == 0L) return(NA_real_)
    
    # calculates the log-likelihood for 'sigma'
    calc_ll <- function(sigma) sum(-0.5 * ddF[i, ] ^ 2L / (sigma ^ 2L * d + eta[i, ] ^ 2L) - 0.5 * log(sigma ^ 2L * d + eta[i, ] ^ 2L) - 0.5 * log(2 * pi), na.rm = TRUE)
    
    mle <- optimize(calc_ll, c(0, 1), maximum = TRUE) # no mutation has sigma larger than 1
    
    if(return_mll) mle$objective else mle$maximum
  }
  
  vapply(1L:nrow(ddF), calc_mle, numeric(1L))
}

# Maximum-likelihood inference for systematic among-interval variation in the rate of change in mutation effect.
ml_interval_rate <- function(d, ddF, eta, nbootstrap = NULL) {
  
  # d : numeric vector
  #   Sequence divergence of each interval.
  # ddF : numeric matrix
  #   Amount of change in mutation effect. Generated by `calc_ddF`.
  # eta : numeric matrix
  #   Standard error of the mean for the amount of change in mutation effect.
  #   Generated by `calc_eta`.
  # nbootstrap : integer
  #   Number of bootstrap replicates to sample.
  #
  # This function infers lambda, which quantifies the relative rate of change in mutation effect for each interval.
  # Mutations are bootstrap sampled to determine estimation noise.
  
  # Calculate the log-likelihoood of `lambda`.
  calc_ll_lambda <- function(lambda) {
    
    if(all(lambda == 0)) lambda <- rep(1, length(d)) else lambda <- lambda / mean(lambda)
    
    # Calculate the log-likelihood for mutation `k` under given `lambda`.
    calc_ll_per_mutation <- function(k) {
      
      if(sum(!is.na(ddF[k, ])) == 0L) return(NA_real_)
      
      # Calculate the log-likelihood for mutation `k` under given `sigma` and `lambda`.
      calc_ll <- function(sigma) {
        
        total_var <- lambda * (sigma ^ 2L) * d + eta[k, ] ^ 2L
        
        sum(-0.5 * ddF[k, ] ^ 2L / total_var - 0.5 * log(total_var) - 0.5 * log(2 * pi), na.rm = TRUE)
      }
      
      optimize(calc_ll, c(0, 1), maximum = TRUE)$objective # no mutation has sigma larger than 1
    }
    
    sum(vapply(1L:nrow(ddF), calc_ll_per_mutation, numeric(1L)), na.rm = TRUE)
  }
  
  res <- optim(rep(1, length(d)), calc_ll_lambda, method = 'L-BFGS-B', lower = 0, upper = Inf, control = list(fnscale = -1))
  
  # No bootstrapping.
  if(is.null(nbootstrap)) return(res$par / mean(res$par))
  
  lambda_bootstrap <- matrix(NA_real_, nbootstrap, length(d))
  ddF_ori <- ddF
  eta_ori <- eta
  
  for(i in 1L:nbootstrap) {
    
    mutations_to_sample <- sample(nrow(ddF_ori), nrow(ddF_ori), TRUE)
    ddF <- ddF_ori[mutations_to_sample, ]
    eta <- eta_ori[mutations_to_sample, ]
    
    res <- optim(rep(1, length(d)), calc_ll_lambda, method = 'L-BFGS-B', lower = 0, upper = Inf, control = list(fnscale = -1))
    
    lambda_bootstrap[i, ] <- res$par / mean(res$par)
    
    cat('Progress: ', i, '/', nbootstrap, '\n', sep = '')
  }
  
  lambda_bootstrap
}
