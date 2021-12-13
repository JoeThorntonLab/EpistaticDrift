
# Set to the supplementary file directory.
setwd('~/Desktop/Single Mutant DMS/Data/Supplementary files/')

source('Scripts/basic_functions.R')


# Analyses.

DP <- read.table('dF.txt', stringsAsFactors = FALSE, header = TRUE)
proteins <- c(1L, 10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L)
# proteins <- c(2L, 11L, 20L, 6L, 7L, 9L, 12L, 5L, 14L) # Alt-All

# The effects of mutations in each background.
dF <- calc_dF(proteins)

# Precalculated memory half-life.
m <- as.numeric(readLines('Scripts/figure_4/memory_half_life.txt'))
# m <- as.numeric(readLines('Scripts/figure_4/memory_half_life_AltAll.txt')) # Alt-All

# Pattern of change in accessibility.

data_all <- multiple_comparison_data_frame(proteins, TRUE, FALSE)
classification <- matrix(NA_integer_, nrow(data_all[[1L]]), length(proteins))

for(i in seq_along(proteins)) {
  
  data <- data_all[[i]] - WT_ACTIVITY
  p <- apply(data, 1L, function(x) if(any(is.na(x))) NA_real_ else t.test(x - (-0.2), alternative = 'two.sided')$p.value)
  is_significant <- p <= calc_BH_cutoff(p, 0.2)
  delta <- apply(data - (-0.2), 1L, mean)
  classification[!is_significant, i] <- 0L
  classification[is_significant & delta > 0, i] <- 1L
  classification[is_significant & delta < 0, i] <- -1L
}


# Number of genetic backgrounds where each mutation is accessible.

n_accessible <- apply(classification, 1L, function(x) sum(x %in% c(0L, 1L)))
total_n <- apply(classification, 1L, function(x) sum(is.finite(x)))

res <- vapply(0L:9L, function(x) sum(total_n == 9L & n_accessible == x & m < 50, na.rm = TRUE), 0)
res <- res / sum(res)
barplot(res, ylim = c(0, 0.6))
sum(res[c(-1L, -10L)])



# Fraction of initially nondeleterious mutations that later becomes deleterious.

init_nondel <- classification[, 2L] %in% c(0L, 1L)
later_del <- apply(classification[, c(-1L, -2L)], 1L, function(x) any(x == -1L, na.rm = TRUE))

# Short-memory
sum(m < 50 & init_nondel & later_del, na.rm = TRUE) / 
  sum(m < 50 & init_nondel & is.finite(later_del), na.rm = TRUE)

# Medium-memory
sum(m > 50 & m < 200 & init_nondel & later_del, na.rm = TRUE) / 
  sum(m > 50 & m < 200 & init_nondel & is.finite(later_del), na.rm = TRUE)

# Long-memory
sum(m > 200 & init_nondel & later_del, na.rm = TRUE) / 
  sum(m > 200 & init_nondel & is.finite(later_del), na.rm = TRUE)


# Fraction of initially deleterious mutations that later becomes nondeleterious.

init_del <- classification[, 2L] == -1L
later_nondel <- apply(classification[, c(-1L, -2L)], 1L, function(x) any(x %in% c(0L, 1L), na.rm = TRUE))

# Short-memory
sum(m < 50 & init_del & later_nondel, na.rm = TRUE) / 
  sum(m < 50 & init_del & is.finite(later_nondel), na.rm = TRUE)

# Medium-memory
sum(m > 50 & m < 200 & init_del & later_nondel, na.rm = TRUE) / 
  sum(m > 50 & m < 200 & init_del & is.finite(later_nondel), na.rm = TRUE)

# Long-memory
sum(m > 200 & init_del & later_nondel, na.rm = TRUE) / 
  sum(m > 200 & init_del & is.finite(later_nondel), na.rm = TRUE)


