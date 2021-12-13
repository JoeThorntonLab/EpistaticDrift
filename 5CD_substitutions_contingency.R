
# Set to the supplementary file directory.
setwd('~/Desktop/Single Mutant DMS/Data/Supplementary files/')

source('Scripts/basic_functions.R')


# Analyses.

SEM_CUTOFF <- 0.05 # More stringency required for these analyses

DP <- read.table('dF.txt', stringsAsFactors = FALSE, header = TRUE)

proteins <- c(1L, 10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L)
# proteins <- c(2L, 11L, 20L, 6L, 7L, 9L, 12L, 5L, 14L) # Alt-All
intervals <- cbind(c(1L, 10L, 19L, 6L, 7L, 8L, 10L, 4L), c(10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L))
# intervals <- matrix(c(c(2L, 11L, 20L, 6L, 7L, 9L, 11L, 5L), c(11L, 20L, 6L, 7L, 9L, 12L, 5L, 14L)), nrow = 8L, ncol = 2L) # Alt-All
d <- calc_d(intervals)

# Substitutions from AncSR to extant SRs.
substitutions <- read.table('Scripts/figure_5/substitutions_SR.txt', stringsAsFactors = FALSE, header = TRUE)


# The initial effect and memory half-life of each substitution.

# The effects and measurement noise of substitutions on `proteins`.

dF <- matrix(NA_real_, nrow(substitutions), length(proteins))
noise <- matrix(NA_real_, nrow(substitutions), length(proteins))

for(i in seq_along(proteins)) {
  
  effects <- calc_mutation_effects(substitutions, proteins[i])
  dF[, i] <- apply(effects, 1L, function(x) if(sum(is.finite(x)) < 2L) NA_real_ else mean(x, na.rm = TRUE))
  noise[, i] <- apply(effects, 1L, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))))
} 

# Amount of change in substitution effect (ddF) and measurement noise of ddF along `intervals`.

ddF <- matrix(NA_real_, nrow(substitutions), nrow(intervals))
eta <- matrix(NA_real_, nrow(substitutions), nrow(intervals))

for(k in 1L:nrow(intervals)) {
  anc_dF <- dF[, proteins == intervals[k, 1L]]
  der_dF <- dF[, proteins == intervals[k, 2L]]
  ddF[, k] <- der_dF - anc_dF
  
  anc_noise <- noise[, proteins == intervals[k, 1L]]
  der_noise <- noise[, proteins == intervals[k, 2L]]
  eta[, k] <- sqrt(anc_noise ^ 2L + der_noise ^ 2L)
}

# Rate of change in substitutions' effects.

retain <- apply(ddF, 1L, function(x) if(sum(!is.na(x)) <= 3L) FALSE else TRUE)
sigma <- ml_constant_rate(d, ddF, eta)
sigma[!retain] <- NA_real_

# Memory half-life of substitutions.

sub_m <- exp(-2.566) * sigma ^ -2.155 # Relation between rate of epistatic divergence and memory half-life
sub_m[sub_m > 200] <- 201

# Effects on AncSR.

sub_effect <- calc_mutation_effects(substitutions, 10L, TRUE)

# Only keeping mutations with measured initial effect and memory length.

retain <- is.finite(sub_m + sub_effect)
substitutions <- substitutions[retain, ]
sub_m <- sub_m[retain]

# Significance test of less than dF = -0.2 is performed to identify delterious mutations.

sub_effect <- calc_mutation_effects(substitutions, 10L, FALSE)
sub_p <- apply(sub_effect, 1L, function(x) if(sum(is.finite(x)) < 2L) NA_real_ else t.test(x + 0.2, alternative = 'two.sided')$p.value)
sub_effect <- calc_mutation_effects(substitutions, 10L, TRUE)
sub_init_del <- (sub_p <= calc_BH_cutoff(sub_p, 0.1) & sub_effect < -0.2)
sub_init_nondel <- !sub_init_del


# Initial effect and memory half-life of mutations.

mut_effect <- apply(DP[DP$PROT == 10L, grep('F', colnames(DP))], 1L,
                    function(x) if(sum(is.finite(x)) < 2L) NA_real_ else mean(x, na.rm = TRUE)) - WT_ACTIVITY
mut_p <- apply(DP[DP$PROT == 10L, grep('F', colnames(DP))], 1L,
               function(x) if(sum(is.finite(x)) < 2L) NA_real_ else t.test(x - WT_ACTIVITY + 0.2, alternative = 'less')$p.value)
mut_init_del <- (mut_p <= calc_BH_cutoff(mut_p, 0.1))
mut_init_nondel <- (mut_p > calc_BH_cutoff(mut_p, 0.1))

mut_m <- as.numeric(readLines('Scripts/figure_4/memory_half_life_AncSR_mutations.txt'))
# mut_m <- as.numeric(readLines('Scripts/figure_4/memory_half_life_AncSR_mutations_AltAll.txt'))

# Plot: Distribution of initial effect for mutations vs. substitutions (Fig. 5C).

h_mut_long <- hist(mut_effect[mut_m > 200], breaks = seq(-1.4, 0.4, by = 0.2), plot = FALSE)
h_mut_medium <- hist(mut_effect[mut_m > 50 & mut_m < 200], breaks = seq(-1.4, 0.4, by = 0.2), plot = FALSE)
h_mut_short <- hist(mut_effect[mut_m < 50], breaks = seq(-1.4, 0.4, by = 0.2), plot = FALSE)
h_sub_long <- hist(sub_effect[sub_m > 200], breaks = seq(-1.4, 0.4, by = 0.2), plot = FALSE)
h_sub_medium <- hist(sub_effect[sub_m > 50 & sub_m < 200], breaks = seq(-1.4, 0.4, by = 0.2), plot = FALSE)
h_sub_short <- hist(sub_effect[sub_m < 50], breaks = seq(-1.4, 0.4, by = 0.2), plot = FALSE)

h_mut_total <- sum(h_mut_long$counts) + sum(h_mut_medium$counts) + sum(h_mut_short$counts)
h_mut_long$counts <- h_mut_long$counts / h_mut_total
h_mut_medium$counts <- h_mut_medium$counts / h_mut_total
h_mut_short$counts <- h_mut_short$counts / h_mut_total

h_sub_total <- sum(h_sub_long$counts) + sum(h_sub_medium$counts) + sum(h_sub_short$counts)
h_sub_long$counts <- h_sub_long$counts / h_sub_total
h_sub_medium$counts <- h_sub_medium$counts / h_sub_total
h_sub_short$counts <- h_sub_short$counts / h_sub_total

plot(h_mut_long, xlim = c(-1.4, 0.4), ylim = c(0, 0.4))
points(h_sub_long$mids, h_sub_long$counts, col = 'red', cex = 1.5)
points(h_sub_long$mids, h_sub_long$counts, col = 'red', type = 'l')
abline(v = -0.2)

plot(h_mut_medium, xlim = c(-1.4, 0.4), ylim = c(0, 0.4))
points(h_sub_medium$mids, h_sub_medium$counts, col = 'red', cex = 1.5)
points(h_sub_medium$mids, h_sub_medium$counts, col = 'red', type = 'l')
abline(v = -0.2)

plot(h_mut_short, xlim = c(-1.4, 0.4), ylim = c(0, 0.4))
points(h_sub_short$mids, h_sub_short$counts, col = 'red', cex = 1.5)
points(h_sub_short$mids, h_sub_short$counts, col = 'red', type = 'l')
abline(v = -0.2)


# Fold-enrichment of mutations with dF > -0.2 among substitutions.

x <- mut_effect[mut_m > 50 & mut_m < 200]
y <- sub_effect[sub_m > 50 & sub_m < 200]
x <- x[is.finite(x)]
y <- y[is.finite(y)]
(sum(y > -0.2) / sum(y < -0.2)) / (sum(x > -0.2) / sum(x < -0.2))


# Frequency of 6 types of mutations.

frac_mut <- c(sum(mut_m < 50 & mut_init_nondel, na.rm = TRUE),
              sum(mut_m < 50 & mut_init_del, na.rm = TRUE),
              sum(mut_m > 50 & mut_m < 200 & mut_init_nondel, na.rm = TRUE),
              sum(mut_m > 50 & mut_m < 200 & mut_init_del, na.rm = TRUE),
              sum(mut_m > 200 & mut_init_nondel, na.rm = TRUE),
              sum(mut_m > 200 & mut_init_del, na.rm = TRUE)) / sum(is.finite(mut_m + mut_effect))

# Frequency of 6 types of substitutions.

frac_sub <- c(sum(sub_m < 50 & sub_init_nondel, na.rm = TRUE),
              sum(sub_m < 50 & sub_init_del, na.rm = TRUE),
              sum(sub_m > 50 & sub_m < 200 & sub_init_nondel, na.rm = TRUE),
              sum(sub_m > 50 & sub_m < 200 & sub_init_del, na.rm = TRUE),
              sum(sub_m > 200 & sub_init_nondel, na.rm = TRUE),
              sum(sub_m > 200 & sub_init_del, na.rm = TRUE)) / sum(is.finite(sub_m + sub_effect))


# (Not used) Later effects of initially nondeleterious short-memory substitutions.

dF <- dF[retain, ]
noise <- noise[retain, ]
index <- which(sub_m < 50 & sub_init_nondel)

init_dF <- dF[index, 2L]
init_noise <- noise[index, 2L]
min_dF <- dF[cbind(index, apply(dF[index, -1L], 1L, function(x) which.min(x) + 1L))]
min_noise <- noise[cbind(index, apply(dF[index, -1L], 1L, function(x) which.min(x) + 1L))]

plot(init_dF, ylim = c(-1.4, 0.4), cex = 1.8, pch = 16L)
abline(h = c(-0.2, 0.4, -1.2, L - WT_ACTIVITY, U - WT_ACTIVITY))
for(i in seq_along(init_noise)) {
  arrows(x0 = i, y0 = init_dF[i] - init_noise[i], x1 = i, y1 = init_dF[i] + init_noise[i], angle = 90, code = 3, length = 0.1)
}
points(min_dF, cex = 1.8, pch = 16L, col = 'red')
for(i in seq_along(init_noise)) {
  arrows(x0 = i, y0 = min_dF[i] - min_noise[i], x1 = i, y1 = min_dF[i] + min_noise[i], angle = 90, code = 3, length = 0.1)
}
