
# Set to the supplementary file directory.
setwd('~/Desktop/Single Mutant DMS/Data/Supplementary files/')

source('Scripts/basic_functions.R')


# Functions.

# Calculate ddF and eta for all 400 possible mutations in each site.
calc_ddF_for_all_mutations <- function(proteins, intervals) {
  
  mutant_list <- list()
  dF_list <- list()
  noise_list <- list()
  
  # Calculating dF and noise for all 400 mutations in each site.
  
  for(site in 1L:76L) {
    
    site_dF <- rep(list(NULL), length(proteins))
    site_noise <- rep(list(NULL), length(proteins))
    
    for(i in seq_along(proteins)) {
      
      data_block <- DP[DP$PROT == proteins[i] & DP$SITE == site, ]
      
      for(WTAA in AA) {
        
        data_block_modified <- change_ref_state(data_block, WTAA)
        site_dF[[i]] <- c(site_dF[[i]], apply(data_block_modified, 1L, mean))
        site_noise[[i]] <- c(site_noise[[i]], apply(data_block_modified, 1L, function(x) sd(x) / sqrt(length(x))))
      }
    }
    
    dF_list[[site]] <- do.call(cbind, site_dF)
    noise_list[[site]] <- do.call(cbind, site_noise)
    mutant_list[[site]] <- data.frame(SITE = rep(site, length(AA) ^ 2L), WTAA = rep(AA, each = length(AA)),
                                      MTAA = rep(AA, length(AA)))
  }
  
  mutant_list <- do.call(rbind, mutant_list)
  dF_list <- do.call(rbind, dF_list)
  noise_list <- do.call(rbind, noise_list)
  
  # Calculating ddF and eta.
  
  ddF_list <- apply(intervals, 1L, function(k) dF_list[, proteins == k[2L]] - dF_list[, proteins == k[1L]])
  eta_list <- apply(intervals, 1L, function(k) sqrt(noise_list[, proteins == k[2L]] ^ 2L + noise_list[, proteins == k[1L]] ^ 2L))
  
  list(mutant_list, ddF_list, eta_list)
}


# Importing data.

DP <- read.table('dF.txt', stringsAsFactors = FALSE, header = TRUE)
proteins <- c(1L, 10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L)
intervals <- cbind(c(1L, 10L, 19L, 6L, 7L, 8L, 10L, 4L), c(10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L))
mutant_info <- multiple_comparison_data_frame(proteins, FALSE, FALSE)[[1L]][, c('SITE', 'WTAA', 'MTAA')]

# Mutations' effects.
dF <- calc_dF(proteins)

# Memory half-life.
h <- as.numeric(readLines('Scripts/figure_4/memory_half_life.txt')); h[h > 200] <- 200

# RSA
ER_RSA <- as.numeric(readLines('Scripts/figure_6/ER1_RSA.txt'))
GR_RSA <- as.numeric(readLines('Scripts/figure_6/GR_RSA.txt'))  
RSA <- apply(cbind(ER_RSA, GR_RSA), 1L, mean, na.rm = TRUE)

# Residue-residue distance (side chain center).
pairwise_dist <- read.table('Scripts/figure_6/residue_residue_distance.txt')

# Relative rate of substitution.
rel_rate <- as.numeric(readLines('Scripts/figure_6/relative_substitution_rate.txt'))


# Distribution of memory half-life among and within sites (Fig. 6A).

h_by_site <- unname(split(h, mutant_info$SITE))
h_by_site_mean <- vapply(unname(split(h, mutant_info$SITE)), mean, numeric(1L), na.rm = TRUE)
h_by_site_median <- vapply(unname(split(h, mutant_info$SITE)), median, numeric(1L), na.rm = TRUE)
dotplot_with_CI(1L:76L, h_by_site[order(h_by_site_median, h_by_site_mean)], 1, ylim = c(0, 200), pch = 16L, metric = median)


# Using median memory half-life of site as a predictor (Fig. 6B).

plot(h_by_site_median[mutant_info$SITE], h, xlim = c(0, 200), ylim = c(0, 200),
     cex = 0.5, col = rgb(0, 0, 0, 0.3), pch = 16L)
cor(h_by_site_median[mutant_info$SITE], h, use = 'pairwise.complete.obs') ^ 2L
l <- lm(h ~ h_by_site_median[mutant_info$SITE])
abline(a = l$coefficients[1L], b = l$coefficients[2L])


# RSA as an explanation of median per-site memory half-life (Fig. S9C).

plot(RSA, h_by_site_median, xlim = c(0, 1), ylim = c(0, 200), pch = 16L, cex = 1.3, col = rgb(0, 0, 0, 0.75))
l <- lm(h_by_site_median ~ RSA)
abline(a = l$coefficients[1L], b = l$coefficients[2L])
cor(RSA, h_by_site_median, use = 'pairwise.complete.obs') ^ 2L


# Rate of substitution as an explanation of median per-site memory half-life (Fig. S9C).

plot(rel_rate, h_by_site_median, xlim = c(0, 5), ylim = c(0, 200), pch = 16L, cex = 1.3, col = rgb(0, 0, 0, 0.75))
l <- lm(h_by_site_median ~ rel_rate)
abline(a = l$coefficients[1L], b = l$coefficients[2L])
cor(rel_rate, h_by_site_median, use = 'pairwise.complete.obs') ^ 2L


# Rate of substitution at neighboring sites as an explanation of median per-site memory half-life (Fig. S9C).

dist_cutoff <- 5.89 # Value found to maximize the correlation.
neighborhood_rate <- vapply(1L:76L, function(i) {
  
  neighborhood <- which(pairwise_dist[i, ] <= dist_cutoff)
  if(length(neighborhood) == 0L) 0 else sum(rel_rate[neighborhood])
}, 0)

plot(neighborhood_rate, h_by_site_median, xlim = c(0, 8), ylim = c(0, 200),
     pch = 16L, cex = 1.3, col = rgb(0, 0, 0, 0.75))
l <- lm(h_by_site_median ~ neighborhood_rate)
abline(a = l$coefficients[1L], b = l$coefficients[2L])
cor(neighborhood_rate, h_by_site_median, use = 'pairwise.complete.obs') ^ 2L


# Distance to the recognition helix as an explanation of median per-site memory half-life (Fig. S9C).

dist_RH <- apply(pairwise_dist, 1L, function(x) min(x[19L:30L], na.rm = TRUE))
dist_RH[19L:30L] <- 0
dist_RH[76L] <- NA_real_

plot(dist_RH, h_by_site_median, xlim = c(0, 20), ylim = c(0, 200),
     pch = 16L, cex = 1.3, col = rgb(0, 0, 0, 0.75))
l <- lm(h_by_site_median ~ dist_RH)
abline(a = l$coefficients[1L], b = l$coefficients[2L])
cor(dist_RH, h_by_site_median, use = 'pairwise.complete.obs') ^ 2L


# Distance to the dimerization interface as an explanation of median per-site memory half-life (Fig. S9C).

interface <- c(36L, 38L, 39L, 40L, 42L, 43L, 44L, 48L, 49L, 52L, 53L)
dist_interface <- apply(pairwise_dist, 1L, function(x) min(x[interface], na.rm = TRUE))
dist_interface[interface] <- 0
dist_interface[76L] <- NA_real_

plot(dist_interface, h_by_site_median, xlim = c(0, 31), ylim = c(0, 200),
     pch = 16L, cex = 1.3, col = rgb(0, 0, 0, 0.75))
l <- lm(h_by_site_median ~ dist_interface)
abline(a = l$coefficients[1L], b = l$coefficients[2L])
cor(dist_interface, h_by_site_median, use = 'pairwise.complete.obs') ^ 2L


# Using median memory half-life of mutation type as a predictor (Fig. 6C).

mutant_type <- paste0(mutant_info$WTAA, mutant_info$MTAA)
mutant_type_count <- table(mutant_type[is.finite(h)])
mutant_type_valid <- names(mutant_type_count)[mutant_type_count >= 4L]

h_by_mutant_type <- split(h, mutant_type)
h_by_mutant_type <- h_by_mutant_type[names(h_by_mutant_type) %in% mutant_type_valid]

h_predicted_by_mutant_type <- lapply(h_by_mutant_type, function(x) rep(median(x, na.rm = TRUE), length(x)))

plot(unlist(h_predicted_by_mutant_type), unlist(h_by_mutant_type), xlim = c(0, 200), ylim = c(0, 200),
     cex = 0.5, col = rgb(0, 0, 0, 0.3), pch = 16L)
l <- lm(unlist(h_by_mutant_type) ~ unlist(h_predicted_by_mutant_type))
abline(a = l$coefficients[1L], b = l$coefficients[2L])

cor(unlist(h_by_mutant_type), unlist(h_predicted_by_mutant_type), use = 'pairwise.complete.obs') ^ 2L
