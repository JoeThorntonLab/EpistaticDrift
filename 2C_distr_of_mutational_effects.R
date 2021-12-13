
# Set to the supplementary file directory.
setwd('~/Desktop/Single Mutant DMS/Data/Supplementary files/')

source('Scripts/basic_functions.R')


# Functions.

# Calculate the fraction of mutations with effect within a specified range.
calc_df_cdf <- function(proteins, range) {
  
  # proteins : integer vector
  #   Set of proteins.
  # range : numeric(2L)
  #   Range of mutational effect.
  
  f <- rep(NA_real_, length(proteins))
  
  for(i in 1L:length(proteins)) {
    
    dFi <- apply(DP[DP$PROT == proteins[i], grep('F', colnames(DP))], 1L, 
                 function(x) if(sum(is.finite(x)) < 2L) NA_real_ else mean(x, na.rm = TRUE)) - WT_ACTIVITY
    f[i] <- sum(dFi >= range[1L] & dFi <= range[2L], na.rm = TRUE) / sum(is.finite(dFi))
  }
  
  f
}


# Analyses.

DP <- read.table('dF.txt', stringsAsFactors = FALSE, header = TRUE)
proteins <- c(1L, 10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L)

f_negative <- calc_df_cdf(proteins, c(-Inf, 0))
f_destructive <- calc_df_cdf(proteins, c(-Inf, L + 0.05 - WT_ACTIVITY))

plot(f_negative, ylim = c(0, 1))
points(f_destructive, pch = 2L)
