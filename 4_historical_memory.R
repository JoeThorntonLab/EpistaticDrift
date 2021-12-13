
# Set to the supplementary file directory.
setwd('~/Desktop/Single Mutant DMS/Data/Supplementary files/')

source('Scripts/basic_functions.R')


# Functions.

# Calculate the squared Pearson's correlation coefficient for the effects of a specified set of mutations between proteins.
calc_dF_cor <- function(proteins, intervals, retain = NULL, data_all = NULL) {
  
  # ARGUMENTS
  #
  # proteins : integer vector
  #     A list of proteins. Defines the mutation set to be analyzed.
  # intervals : two-column integer matrix
  #   Intervals for which correlation of mutation effects is calculated.
  #   Each row specifies the ancestral and descendant protein for an interval.
  # retain : boolean vector
  #   An optional vector specifying the subset of mutations to analyze.
  # data_all : data frame
  #   Output of `multiple_comparison_data_frame`. Speeds up the calculation.
  #
  # This function calculates the squared Pearson's correlation coefficient of mutation effects
  # between each pair of proteins. Attenuation of correlation due to measurement noise is accounted for.
  
  # Calculate the correlation of mean mutation effects, taking attenuation due to measurement noise into account.
  calc_squared_pearson_cor <- function(X, Y) {
    
    # X, Y : numeric matrices
    #     Each row corresponds to mutations, each column to measurement replicates.
    
    MX <- apply(X, 1L, mean)
    MY <- apply(Y, 1L, mean)
    
    # Squared standard error of the mean.
    SSEMX <- apply(X, 1L, var) / sum(!is.na(MX))
    SSEMY <- apply(Y, 1L, var) / sum(!is.na(MY))
    
    # Separation indices.
    RX <- (var(MX, na.rm = TRUE) - mean(SSEMX, na.rm = TRUE)) / var(MX, na.rm = TRUE)
    RY <- (var(MY, na.rm = TRUE) - mean(SSEMY, na.rm = TRUE)) / var(MY, na.rm = TRUE)
    
    cor(MX, MY, method = 'pearson', use = 'pairwise.complete.obs') ^ 2L / (RX * RY)
  }
  
  if(is.null(data_all)) data_all <- multiple_comparison_data_frame(proteins, phenotype_only = TRUE, stack = FALSE)
  
  if(!is.null(retain)) {
    for(i in seq_along(proteins)) {
      data_all[[i]] <- data_all[[i]][retain, ]
    }
  }
  
  cor2 <- rep(NA_real_, nrow(intervals))

  for(k in 1L:nrow(intervals)) {
    
    anc <- match(intervals[k, 1L], proteins)
    der <- match(intervals[k, 2L], proteins)
    cor2[k] <- calc_squared_pearson_cor(data_all[[anc]], data_all[[der]])
  }
  
  cor2
}

# Perform regression for decay of various types.
decay_regression <- function(y, x, curve, y0 = NULL, return_RSS = TRUE, show_plot = TRUE, ...) {
  
  # One of the following curves is fit to data:
  #
  # 1. Segmented linear decay ('linear')
  #   y = y0 + k * x, if y0 + k * x > 0
  #   y = 0, otherwise
  #
  # 2. Exponential decay ('exponential')
  #   y = y0 * exp(-k * x)
  #
  # 3. Harmonic decay ('harmonic')
  #   y = y0 / (1 + x / k)
  #
  # ARGUMENTS
  #
  # y, x : numeric vectors of same length
  #   Data.
  # curve : string; one of 'linear', 'exponential' and 'harmonic'
  #   Type of decay curve to be fit.
  # y0 : numeric
  #     If provided, y0 is fixed to the provided value.
  # return_RSS : bool
  #     If `TRUE`, returns the residual sum of squares. Otherwise the estimated parameters.
  # show_plot : bool
  #     Should a plot of data and fitted curve be generated?
  # ... : optional arguments to be passed to `plot`
  #
  # If `return_RSS` is `TRUE`, returns the residual sum of squares.
  # Otherwise returns the estimated memory length.
  # Memory length is calculated as follows:
  # 'linear': -y0 / k
  # 'exponential': log(2) / k
  # 'harmonic': k
  
  # Calculate the residual sum of squares given c(y0, k).
  calc_residual_sum_squares <- function(param) {
    
    k <- param[2L]
    if(!is.null(y0)) a <- y0 else a <- param[1L]
    
    if(curve == 'linear') {
      predicted <- a + k * x
      predicted[predicted < 0] <- 0
    } else if(curve == 'exponential') {
      predicted <- a * exp(-k * x)
    } else if(curve == 'harmonic') {
      predicted <- a / (1 + x / k)
    }
    
    sum((y - predicted) ^ 2L, na.rm = TRUE)
  }
  
  y <- as.vector(y[is.finite(x)])
  x <- as.vector(x[is.finite(x)])
  
  res <- optim(c(1, 1), calc_residual_sum_squares)
  if(is.null(y0)) y0 <- res$par[1L]
  k <- res$par[2L]
  
  if(show_plot) {
    plot(x, y, ...)
    
    # Best-fit curve
    args <- list(...)
    if(!is.null(args$xlim)) {
      t <- seq(args$xlim[1L], args$xlim[2L], length.out = 1000L)
    } else {
      t <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 1000L)
    }
    
    if(curve == 'linear') {
      predicted <- y0 + k * t
      predicted[predicted < 0] <- 0
    } else if(curve == 'exponential') {
      predicted <- y0 * exp(-k * t)
    } else if(curve == 'harmonic') {
      predicted <- y0 / (1 + t / k)
    }
    
    points(t, predicted, type = 'l')
  }
  
  if(return_RSS) {
    res$value
  } else {
    if(curve == 'linear') {
      -y0 / k
    } else if(curve == 'exponential') {
      log(2) / k
    } else if(curve == 'harmonic') {
      k
    }
  }
}


# Analyses.

DP <- read.table('dF.txt', stringsAsFactors = FALSE, header = TRUE)
proteins <- c(1L, 10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L)
# proteins <- c(2L, 11L, 20L, 6L, 7L, 9L, 12L, 5L, 14L) # Alt-All
intervals <- cbind(c(1L, 10L, 19L, 6L, 7L, 8L, 10L, 4L), c(10L, 19L, 6L, 7L, 8L, 12L, 4L, 14L))
# intervals <- matrix(c(c(2L, 11L, 20L, 6L, 7L, 9L, 11L, 5L), c(11L, 20L, 6L, 7L, 9L, 12L, 5L, 14L)), nrow = 8L, ncol = 2L) # Alt-All
all_intervals <- rbind(intervals, matrix(c(1L, 1L, 1L, 1L, 1L, 10L, 10L, 10L, 10L, 19L, 19L, 19L, 6L, 6L, 7L, 1L, 1L, 10L,
                                           19L, 6L, 7L, 8L, 12L, 6L, 7L, 8L, 12L, 7L, 8L, 12L, 8L, 12L, 12L, 4L, 14L, 14L), ncol = 2L))
#all_intervals <- rbind(intervals, matrix(c(2L, 2L, 2L, 2L, 2L, 11L, 11L, 11L, 11L, 20L, 20L, 20L, 6L, 6L, 7L, 2L, 2L, 11L,
#                                           20L, 6L, 7L, 9L, 12L, 6L, 7L, 9L, 12L, 7L, 9L, 12L, 9L, 12L, 12L, 5L, 14L, 14L), ncol = 2L)) # Alt-All

# Length of each phylogenetic interval.
d <- calc_d(intervals)
d_all <- calc_d(all_intervals)

# The effects of mutations in each background.
dF <- calc_dF(proteins)

# The amount of change in each mutation's effect across each unit interval and its measurement noise.
ddF <- calc_ddF(intervals, remove_out_of_range = FALSE)
eta <- calc_eta(intervals)

# Mutations with 4 or less data points.
# They are not explicitly excluded at this point.
retain <- apply(ddF, 1L, function(x) if(sum(!is.na(x)) <= 4L) FALSE else TRUE)

# Rate of epistatic divergence for each mutation.
sigma <- ml_constant_rate(d, ddF, eta)


# Memory half-life of all mutations (Fig. 4B).

res <- calc_dF_cor(proteins, all_intervals)
decay_regression(res, d_all, 'exponential', 1, FALSE, TRUE, xlim = c(0, 45), ylim = c(0, 1))


# Relationship between rate of epistatic divergence and memory half-life (Fig. 4, C and D).

# Using a fraction `window` of mutations, moving each time by fraction `slide`.
data_all <- multiple_comparison_data_frame(proteins, phenotype_only = TRUE, stack = FALSE)
window <- 0.1
slide <- 0.05
s <- rep(NA_real_, floor((1 - window) / slide) + 1L) # Mean rate of mutation group
m <- rep(NA_real_, floor((1 - window) / slide) + 1L) # Memory half-life of mutation group

for(i in 1:(floor((1 - window) / slide) + 1L)) {
  
  ql <- quantile(sigma[retain], slide * (i - 1L))
  qu <- quantile(sigma[retain], slide * (i - 1L) + window)
  mutations <- which(retain & sigma >= ql & sigma <= qu)
  s[i] <- mean(sigma[mutations])
  res <- calc_dF_cor(proteins, all_intervals, mutations, data_all)
  m[i] <- decay_regression(res, d_all, 'exponential', 1, FALSE, FALSE)
  
  # Plotting some subsets of mutations.
  # if(i %in% c(1L, 13L, 19L)) decay_regression(res, d_all, 'exponential', 1, FALSE, TRUE, xlim = c(0, 45), ylim = c(0, 1))
}

# Regression of power function; groups with memory half-life > 200% are excluded.

plot(s[-c(1L:5L)] ^ 2L, m[-c(1L:5L)], xlim = c(0, 0.0125), ylim = c(0, 225), cex = 1.2, pch = 16L)
l <- lm(log(m[-c(1L:5L)]) ~ log(s[-c(1L:5L)]))$coefficients
points(seq(0.01, 0.15, by = 0.001) ^ 2L, exp(l[1L]) * seq(0.01, 0.15, by = 0.001) ^ l[2L], type = 'l', col = 'red')


# Memory half-life of each mutation (Fig. 4E).

m <- exp(-2.566) * sigma ^ -2.155 # Best-fit power function from above
m[!retain] <- NA_real_
m[m > 200] <- 201

# Memory length of AncSR mutations: m <- parse_data_by_protein(m, proteins)[[which(proteins == 10L)]]

h <- hist(unlist(parse_data_by_protein(m, proteins)), breaks = seq(0, 225, by = 25), plot = FALSE)
h$counts <- h$counts / sum(h$counts)
plot(h, ylim = c(0, 0.55))
abline(v = c(50, 200))


# Comparing the effects of mutations (partitioned by memory length) between AncSR and HsGR (Fig. 2F).
pairwise_comparison_plot_subset(10L, 14L, proteins, 0.1, which(m > 200), separate = TRUE)

