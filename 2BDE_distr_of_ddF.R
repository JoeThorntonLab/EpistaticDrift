
# Set to the supplementary file directory.
setwd('~/Desktop/Single Mutant DMS/Data/Supplementary files/')

source('Scripts/basic_functions.R')


# Functions.

# Calculate the amount of difference in each mutation's effect between two proteins.
calc_ddF_pairwise <- function(i, j, FDR, exclude = NULL) {
  
  # i, j : integer
  #   Two proteins.
  # FDR : numeric
  #   False discovery rate for defining mutations outside dynamic range.
  # exclude : boolean vector
  #   Indicates whether each mutation should be excluded from analysis.
  #
  # ddF = effect on `j` minus effect on `i`.
  #
  # Returns a two-column matrix, rows corresponding to mutations and the columns
  # recording ddF and p value for t-test of ddF = 0.
  
  data_all <- multiple_comparison_data_frame(c(i, j), TRUE, FALSE)
  dFi <- as.matrix(data_all[[1L]]) - WT_ACTIVITY
  dFj <- as.matrix(data_all[[2L]]) - WT_ACTIVITY
  
  if(is.null(exclude)) exclude <- rep(FALSE, nrow(dFi))
  ddF <- rep(NA_real_, nrow(dFi))
  p <- rep(NA_real_, nrow(dFi))
  
  # Checking whether mutation effect is within dynamic range.
  LB <- L - WT_ACTIVITY + 0.05
  is_destructive_i <- apply(dFi, 1L, function(x) if(sum(is.finite(x)) < 2L) NA_real_ else t.test(x - LB, alternative = 'greater')$p.value)
  is_destructive_i <- is_destructive_i > calc_BH_cutoff(is_destructive_i, FDR)
  is_destructive_j <- apply(dFj, 1L, function(x) if(sum(is.finite(x)) < 2L) NA_real_ else t.test(x - LB, alternative = 'greater')$p.value)
  is_destructive_j <- is_destructive_j > calc_BH_cutoff(is_destructive_j, FDR)

  for(k in 1L:nrow(dFi)) {
    
    x <- dFi[k, ]
    y <- dFj[k, ]
    
    if(sum(is.finite(x)) < 2L | sum(is.finite(y)) < 2L) next
    if(is_destructive_i[k] & is_destructive_j[k]) next
    
    ddF[k] <- mean(y - x, na.rm = TRUE)
    p[k] <- t.test(x, y, alternative = 'two.sided')$p.value
  }
  
  cbind(ddF, p)
}


# Analyses.

DP <- read.table('dF.txt', stringsAsFactors = FALSE, header = TRUE)
intervals <- matrix(c(c(1L, 10L, 19L, 6L, 7L, 8L, 10L, 4L), c(10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L)), nrow = 8L, ncol = 2L)
# intervals <- matrix(c(c(2L, 11L, 20L, 6L, 7L, 9L, 11L, 5L), c(11L, 20L, 6L, 7L, 9L, 12L, 5L, 14L)), nrow = 8L, ncol = 2L) # Alt-All

d <- calc_d(intervals)
ddF <- apply(intervals, 1L, function(x) calc_ddF_pairwise(x[1L], x[2L], 0.1))


# Plot: ddF distribution for all mutations and intervals (Fig. 2B).

ddF_all <- unlist(lapply(ddF, function(x) x[is.finite(x[, 1L]), 1L]))
p_all <- unlist(lapply(ddF, function(x) x[, 2L]))
p_cutoff <- calc_BH_cutoff(p_all, 0.1)
ddF_insig <- unlist(lapply(ddF, function(x) x[x[, 2L] > p_cutoff & is.finite(x[, 1L]), 1L]))

ddF_all[ddF_all > 1] <- 1.01
ddF_all[ddF_all < -1] <- -1.01
ddF_insig[ddF_insig > 1] <- 1.01
ddF_insig[ddF_insig < -1] <- -1.01

h_all <- hist(ddF_all, breaks = seq(-1.6, 1.6, by = 0.1), plot = FALSE)
h_insig <- hist(ddF_insig, breaks = seq(-1.6, 1.6, by = 0.1), plot = FALSE)
h_insig$counts <- h_insig$counts / sum(h_all$counts)
h_all$counts <- h_all$counts / sum(h_all$counts)

plot(h_all, xlim = c(-1.1, 1.1), ylim = c(0, 0.3))
plot(h_insig, col = 'grey', add = TRUE)


# Plot: ddF distribution for a specific interval (Fig. 2D).

k <- 8L
ddFk <- ddF[[k]][, 1L]
pk <- ddF[[k]][, 2L]
pk_cutoff <- calc_BH_cutoff(pk, 0.1)
ddFk_insig <- ddFk[pk > pk_cutoff]

ddFk[ddFk > 1] <- 1.01
ddFk[ddFk < -1] <- -1.01
ddFk_insig[ddFk_insig > 1] <- 1.01
ddFk_insig[ddFk_insig < -1] <- -1.01

h_all <- hist(ddFk, breaks = seq(-1.6, 1.6, by = 0.1), plot = FALSE)
h_insig <- hist(ddFk_insig, breaks = seq(-1.6, 1.6, by = 0.1), plot = FALSE)
h_insig$counts <- h_insig$counts / sum(h_all$counts)
h_all$counts <- h_all$counts / sum(h_all$counts)

plot(h_all, xlim = c(-1.1, 1.1), ylim = c(0, 0.4))
plot(h_insig, col = 'grey', add = TRUE)


# Epistatic divergence as a function of sequence divergence (Fig. 2E).

all_intervals <- rbind(intervals, matrix(c(1L, 1L, 1L, 1L, 1L, 10L, 10L, 10L, 10L, 19L, 19L, 19L, 6L, 6L, 7L, 1L, 1L, 10L,
                                           19L, 6L, 7L, 8L, 12L, 6L, 7L, 8L, 12L, 7L, 8L, 12L, 8L, 12L, 12L, 4L, 14L, 14L), ncol = 2L))
#all_intervals <- rbind(intervals, matrix(c(2L, 2L, 2L, 2L, 2L, 11L, 11L, 11L, 11L, 20L, 20L, 20L, 6L, 6L, 7L, 2L, 2L, 11L,
#                                           20L, 6L, 7L, 9L, 12L, 6L, 7L, 9L, 12L, 7L, 9L, 12L, 9L, 12L, 12L, 5L, 14L, 14L), ncol = 2L)) # Alt-All
d_all <- apply(all_intervals, 1L, function(x) calc_sequence_divergence(x[1L], x[2L]))
ddF_all <- apply(all_intervals, 1L, function(x) calc_ddF_pairwise(x[1L], x[2L], 0.1))
ddF_var <- vapply(ddF_all, function(x) var(x[, 1L], na.rm = TRUE), numeric(1L))

l_all <- lm(log(ddF_var) ~ log(d_all)) # Fitting power function; all intervals
l_unit <- lm(log(ddF_var[1L:8L]) ~ log(d_all[1L:8L])) # Fitting power function; unit intervals

# Plot

plot(d_all[-(1L:8L)], ddF_var[-(1L:8L)], xlim = c(0, 45), ylim = c(0, 0.18), pch = 16L, cex = 1.2)
points(d_all[1L:8L], ddF_var[1L:8L], pch = 16L, cex = 1.2, col = 'red')
t <- seq(0, 45, by = 0.01)
points(t, exp(l_all$coefficients[1L]) * t ^ l_all$coefficients[2L], type = 'l')
points(t, exp(l_unit$coefficients[1L]) * t ^ l_unit$coefficients[2L], type = 'l', col = 'red')


# Analysis: Test of unbiasedness for individual mutations.

ddF_matrix <- calc_ddF(intervals)

# Retaining only mutations with 4 or more valid effects.
# This removes mutations always at the lower bound of measurement.
ddF_matrix <- ddF_matrix[apply(ddF_matrix, 1L, function(x) if(sum(!is.na(x)) <= 4L) FALSE else TRUE), ]

# Wilcoxon signed-rank test.

p <- apply(ddF_matrix, 1L, function(x) wilcox.test(x)$p.value)
h <- hist(p, breaks = seq(0, 1, by = 0.1), plot = FALSE)
h$counts <- h$counts / sum(h$counts)
plot(h, ylim = c(0, 0.2))


