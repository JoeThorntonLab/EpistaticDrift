
# Set to the supplementary file directory.
setwd('~/Desktop/Single Mutant DMS/Data/Supplementary files/')

source('Scripts/basic_functions.R')


# Functions.

# Identify the substitutions along a trajectory and their effects on the immediate ancestor and descendant.
multiple_comparison_substitutions <- function(trajectory) {
  
  subs <- data.frame(NULL)
  
  for(i in 2L:length(trajectory)) {
    
    anc <- identify_sequence(trajectory[i - 1L])
    der <- identify_sequence(trajectory[i])
    
    sites <- which(anc != der & anc != '-' & der != '-')
    l <- length(sites)
    
    subs <- rbind(subs, data.frame(ANC = rep(trajectory[i - 1L], l), DER = rep(trajectory[i], l), SITE = sites, WTAA = anc[sites], MTAA = der[sites]))
  }
  
  data_all <- multiple_comparison_data_frame(trajectory, phenotype_only = FALSE)
  
  effects <- matrix(NA_real_, nrow = nrow(subs), ncol = 6L)
  colnames(effects) <- c(paste0('EF_ANC', 1:3), paste0('EF_DER', 1:3))
  
  for(i in 1L:nrow(subs)) {
    
    mutation_ANC <- which(data_all$PROT == subs$ANC[i] & data_all$SITE == subs$SITE[i] & data_all$WTAA == subs$WTAA[i] & data_all$MTAA == subs$MTAA[i])
    mutation_DER <- which(data_all$PROT == subs$DER[i] & data_all$SITE == subs$SITE[i] & data_all$WTAA == subs$WTAA[i] & data_all$MTAA == subs$MTAA[i])
    
    effects[i, ] <- as.numeric(c(data_all[mutation_ANC, grep('F', colnames(data_all))], data_all[mutation_DER, grep('F', colnames(data_all))])) - WT_ACTIVITY
  }
  
  cbind(subs, effects)
}


# Analyses.

DP <- read.table('dF.txt', stringsAsFactors = FALSE, header = TRUE)


# The effects of substitutions.

# 79 substitutions along the trajectories.
subs <- rbind(multiple_comparison_substitutions(c(1L, 10L, 19L, 6L, 7L, 8L, 12L)),
              multiple_comparison_substitutions(c(10L, 4L, 14L)))
#subs <- rbind(multiple_comparison_substitutions(c(2L, 11L, 20L, 6L, 7L, 9L, 12L)),
#              multiple_comparison_substitutions(c(11L, 5L, 14L))) # Alt-All

# Effects on ancestral and descendant protein.
sub_effects_anc <- apply(subs[, grep('EF_ANC', colnames(subs))], 1L, mean, na.rm = TRUE)
sub_effects_der <- apply(subs[, grep('EF_DER', colnames(subs))], 1L, mean, na.rm = TRUE)

# Only substitutions for which the two effects differ by less than 0.2 are used.
substitutions_to_use <- which(abs(sub_effects_anc - sub_effects_der) < 0.2)
sub_effects <- unname(apply(subs[substitutions_to_use, grep('EF', colnames(subs))], 1L, mean, na.rm = TRUE))


# The effects of mutations.

mut_effects <- apply(DP[DP$PROT %in% c(1L, 10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L), grep('F', colnames(DP))], 1L, mean) - WT_ACTIVITY
# mut_effects <- apply(DP[DP$PROT %in% c(2L, 11L, 20L, 6L, 7L, 9L, 12L, 5L, 14L), grep('F', colnames(DP))], 1L, mean) - WT_ACTIVITY # Alt-All
mut_effects <- mut_effects[is.finite(mut_effects)]


# Plot: Mutation spectrum vs. substitution spectrum.

min_effect <- -1.4
max_effect <- 0.4
increment <- 0.2

# Forcing very few mutations/substitutions with effect outside the range [min_effect, max_effect] to be within.
sub_effects[sub_effects < min_effect] <- min_effect + increment / 2
sub_effects[sub_effects > max_effect] <- max_effect - increment / 2
mut_effects[mut_effects < min_effect] <- min_effect + increment / 2
mut_effects[mut_effects > max_effect] <- max_effect - increment / 2

h_sub <- hist(sub_effects, breaks = seq(min_effect, max_effect, by = increment), plot = FALSE)
h_sub$counts <- h_sub$counts / sum(h_sub$counts)
h_mut <- hist(mut_effects, breaks = seq(min_effect, max_effect, by = increment), plot = FALSE)
h_mut$counts <- h_mut$counts / sum(h_mut$counts)

plot(h_mut, ylim = c(0, 0.6))
points(h_sub$mids, h_sub$counts, col = 'red')
points(h_sub$mids, h_sub$counts, col = 'red', type = 'l')
abline(v = -0.2)

# Enrichment of mutations with dF > -0.2 over mutations with dF < -0.2 in substitutions
sum(sub_effects > -0.2) / sum(sub_effects < -0.2) / (sum(mut_effects > -0.2) / sum(mut_effects < -0.2))
