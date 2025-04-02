library(lmFScreen)

####################################### CES-D #########################################

# STEP 1: load in information about the data
# Given data
means <- c(9.33, 7.60, 8.79)
sds <- c(7.16, 5.84, 5.68)
ns <- c(150, 449, 227)
p <- 3  
alpha <- 0.05

# STEP 2: calculate necessary statistics from this data
# Overall sample size
n <- sum(ns)
# Overall mean
overall_mean <- sum(means * ns) / sum(ns)
# Between-group sum of squares (SSB)
ssb <- sum(ns * (means - overall_mean)^2)
# Within-group sum of squares (SSW)
ssw <- sum((ns - 1) * sds^2)
# Total sum of squares (SST)
sst <- ssb + ssw
# Degrees of freedom
df_between <- length(means) - 1
df_within <- sum(ns) - length(means)
# F-statistic for omnibus test
f_statistic <- (ssb / df_between) / (ssw / df_within)
# R-squared
r_squared <- ssb / sst
# Mean Squared Error (MSE) and Residual Squared Error (RSE)
mse <- ssw / df_within
rse <- sqrt(mse)
# Calculate p-value for ANOVA
p_value_anova <- pf(f_statistic, df_between, df_within, lower.tail = FALSE)

# STEP 3: Make sure we can reproduce their inferential results
# Function to calculate SE, t-statistic, and p-value
calculate_t_and_p <- function(mean1, mean2, n1, n2, mse, df_within) {
  se <- sqrt(mse * (1/n1 + 1/n2))
  t_stat <- (mean1 - mean2) / se
  p_val <- 2 * pt(-abs(t_stat), df_within)  # two-tailed p-value
  return(list(t_stat = t_stat, p_val = p_val))
}
# Pairwise comparisons
results_1_2 <- calculate_t_and_p(means[1], means[2], ns[1], ns[2], mse, df_within)  # <65 vs 65-84
results_1_3 <- calculate_t_and_p(means[1], means[3], ns[1], ns[3], mse, df_within)  # <65 vs >=85
results_2_3 <- calculate_t_and_p(means[2], means[3], ns[2], ns[3], mse, df_within)  # 65-84 vs >=85
# Print t-test results (the Sidak values are the same as the Keil et al paper!)
cat("Follow-Up Paired T-Tests:\n")
cat("<65 vs 65-84: t-statistic =", results_1_2$t_stat, ", p-value =", results_1_2$p_val, "sidak =", 1 - (1-results_1_2$p_val)^3, "\n")
cat("<65 vs >=85: t-statistic =", results_1_3$t_stat, ", p-value =", results_1_3$p_val, "sidak =",1 - (1-results_1_3$p_val)^3, "\n")
cat("65-84 vs >=85: t-statistic =", results_2_3$t_stat, ", p-value =", results_2_3$p_val, "sidak =",  1 - (1-results_2_3$p_val)^3, "\n")

# STEP 4: compute the selective p-values 
set.seed(826)
p_under65_6584 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                      Fstat = (results_1_2$t_stat)^2, alpha_ov = alpha)
p_under65_over85 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                      Fstat = (results_1_3$t_stat)^2, alpha_ov = alpha)
p_6584_over85 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                      Fstat = (results_2_3$t_stat)^2, alpha_ov = alpha)

# Output the adjusted p-values
cat("Selective p-values:\n")
cat("p_under65_6584:", p_under65_6584, "\n")
cat("p_under65_over85:", p_under65_over85, "\n")
cat("p_6584_over85:", p_6584_over85, "\n")


######################################### Sleep #######################################

# STEP 1: load in information about the data
# Given data
means <- c(6.65, 6.93, 7.07)
sds <- c(0.87, 0.87, 1.02)
ns <- c(150, 449, 227)
p <- 3  
alpha <- 0.05

# STEP 2: calculate necessary statistics from this data
# Overall sample size
n <- sum(ns)
# Overall mean
overall_mean <- sum(means * ns) / sum(ns)
# Between-group sum of squares (SSB)
ssb <- sum(ns * (means - overall_mean)^2)
# Within-group sum of squares (SSW)
ssw <- sum((ns - 1) * sds^2)
# Total sum of squares (SST)
sst <- ssb + ssw
# Degrees of freedom
df_between <- length(means) - 1
df_within <- sum(ns) - length(means)
# F-statistic for omnibus test
f_statistic <- (ssb / df_between) / (ssw / df_within)
# R-squared
r_squared <- ssb / sst
# Mean Squared Error (MSE) and Residual Squared Error (RSE)
mse <- ssw / df_within
rse <- sqrt(mse)
# Calculate p-value for ANOVA
p_value_anova <- pf(f_statistic, df_between, df_within, lower.tail = FALSE)

# STEP 3: Make sure we can reproduce their inferential results
# Function to calculate SE, t-statistic, and p-value
calculate_t_and_p <- function(mean1, mean2, n1, n2, mse, df_within) {
  se <- sqrt(mse * (1/n1 + 1/n2))
  t_stat <- (mean1 - mean2) / se
  p_val <- 2 * pt(-abs(t_stat), df_within)  # two-tailed p-value
  return(list(t_stat = t_stat, p_val = p_val))
}
# Pairwise comparisons
results_1_2 <- calculate_t_and_p(means[1], means[2], ns[1], ns[2], mse, df_within)  # <65 vs 65-84
results_1_3 <- calculate_t_and_p(means[1], means[3], ns[1], ns[3], mse, df_within)  # <65 vs >=85
results_2_3 <- calculate_t_and_p(means[2], means[3], ns[2], ns[3], mse, df_within)  # 65-84 vs >=85
# Print t-test results (the Sidak values are the same as the Keil et al paper!)
cat("Follow-Up Paired T-Tests:\n")
cat("<65 vs 65-84: t-statistic =", results_1_2$t_stat, ", p-value =", results_1_2$p_val, "sidak =", 1 - (1-results_1_2$p_val)^3, "\n")
cat("<65 vs >=85: t-statistic =", results_1_3$t_stat, ", p-value =", results_1_3$p_val, "sidak =",1 - (1-results_1_3$p_val)^3, "\n")
cat("65-84 vs >=85: t-statistic =", results_2_3$t_stat, ", p-value =", results_2_3$p_val, "sidak =",  1 - (1-results_2_3$p_val)^3, "\n")

# STEP 4: compute the selective p-values 
set.seed(826)
p_under65_6584 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                            Fstat = (results_1_2$t_stat)^2, alpha_ov = alpha)
p_under65_over85 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                            Fstat = (results_1_3$t_stat)^2, alpha_ov = alpha)
p_6584_over85 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                            Fstat = (results_2_3$t_stat)^2, alpha_ov = alpha)

# Output the adjusted p-values
cat("Selective p-values:\n")
cat("p_under65_6584:", p_under65_6584, "\n")
cat("p_under65_over85:", p_under65_over85, "\n")
cat("p_6584_over85:", p_6584_over85, "\n")



######################################### MMSES #######################################

# STEP 1: load in information about the data
# Given data
means <- c(29.31, 28.12, 26.14)
sds <- c(1.19, 2.42, 4.15)
ns <- c(150, 449, 227)
p <- 3  
alpha <- 0.05

# STEP 2: calculate necessary statistics from this data
# Overall sample size
n <- sum(ns)
# Overall mean
overall_mean <- sum(means * ns) / sum(ns)
# Between-group sum of squares (SSB)
ssb <- sum(ns * (means - overall_mean)^2)
# Within-group sum of squares (SSW)
ssw <- sum((ns - 1) * sds^2)
# Total sum of squares (SST)
sst <- ssb + ssw
# Degrees of freedom
df_between <- length(means) - 1
df_within <- sum(ns) - length(means)
# F-statistic for omnibus test
f_statistic <- (ssb / df_between) / (ssw / df_within)
# R-squared
r_squared <- ssb / sst
# Mean Squared Error (MSE) and Residual Squared Error (RSE)
mse <- ssw / df_within
rse <- sqrt(mse)
# Calculate p-value for ANOVA
p_value_anova <- pf(f_statistic, df_between, df_within, lower.tail = FALSE)

# STEP 3: Make sure we can reproduce their inferential results
# Function to calculate SE, t-statistic, and p-value
calculate_t_and_p <- function(mean1, mean2, n1, n2, mse, df_within) {
  se <- sqrt(mse * (1/n1 + 1/n2))
  t_stat <- (mean1 - mean2) / se
  p_val <- 2 * pt(-abs(t_stat), df_within)  # two-tailed p-value
  return(list(t_stat = t_stat, p_val = p_val))
}
# Pairwise comparisons
results_1_2 <- calculate_t_and_p(means[1], means[2], ns[1], ns[2], mse, df_within)  # <65 vs 65-84
results_1_3 <- calculate_t_and_p(means[1], means[3], ns[1], ns[3], mse, df_within)  # <65 vs >=85
results_2_3 <- calculate_t_and_p(means[2], means[3], ns[2], ns[3], mse, df_within)  # 65-84 vs >=85
# Print t-test results (the Sidak values are the same as the Keil et al paper!)
cat("Follow-Up Paired T-Tests:\n")
cat("<65 vs 65-84: t-statistic =", results_1_2$t_stat, ", p-value =", results_1_2$p_val, "sidak =", 1 - (1-results_1_2$p_val)^3, "\n")
cat("<65 vs >=85: t-statistic =", results_1_3$t_stat, ", p-value =", results_1_3$p_val, "sidak =",1 - (1-results_1_3$p_val)^3, "\n")
cat("65-84 vs >=85: t-statistic =", results_2_3$t_stat, ", p-value =", results_2_3$p_val, "sidak =",  1 - (1-results_2_3$p_val)^3, "\n")

# STEP 4: compute the selective p-values 
set.seed(826)
p_under65_6584 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                             Fstat = (results_1_2$t_stat)^2, alpha_ov = alpha)
p_under65_over85 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                               Fstat = (results_1_3$t_stat)^2, alpha_ov = alpha)
p_6584_over85 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                            Fstat = (results_2_3$t_stat)^2, alpha_ov = alpha)

# Output the adjusted p-values
cat("Selective p-values:\n")
cat("p_under65_6584:", p_under65_6584, "\n")
cat("p_under65_over85:", p_under65_over85, "\n")
cat("p_6584_over85:", p_6584_over85, "\n")


######################################### MDRS #######################################

# STEP 1: load in information about the data
# Given data
means <- c(141.5, 139.0, 131.6)
sds <- c(2.57, 5.85, 14.89)
ns <- c(150, 449, 227)
p <- 3  
alpha <- 0.05

# STEP 2: calculate necessary statistics from this data
# Overall sample size
n <- sum(ns)
# Overall mean
overall_mean <- sum(means * ns) / sum(ns)
# Between-group sum of squares (SSB)
ssb <- sum(ns * (means - overall_mean)^2)
# Within-group sum of squares (SSW)
ssw <- sum((ns - 1) * sds^2)
# Total sum of squares (SST)
sst <- ssb + ssw
# Degrees of freedom
df_between <- length(means) - 1
df_within <- sum(ns) - length(means)
# F-statistic for omnibus test
f_statistic <- (ssb / df_between) / (ssw / df_within)
# R-squared
r_squared <- ssb / sst
# Mean Squared Error (MSE) and Residual Squared Error (RSE)
mse <- ssw / df_within
rse <- sqrt(mse)
# Calculate p-value for ANOVA
p_value_anova <- pf(f_statistic, df_between, df_within, lower.tail = FALSE)

# STEP 3: Make sure we can reproduce their inferential results
# Function to calculate SE, t-statistic, and p-value
calculate_t_and_p <- function(mean1, mean2, n1, n2, mse, df_within) {
  se <- sqrt(mse * (1/n1 + 1/n2))
  t_stat <- (mean1 - mean2) / se
  p_val <- 2 * pt(-abs(t_stat), df_within)  # two-tailed p-value
  return(list(t_stat = t_stat, p_val = p_val))
}
# Pairwise comparisons
results_1_2 <- calculate_t_and_p(means[1], means[2], ns[1], ns[2], mse, df_within)  # <65 vs 65-84
results_1_3 <- calculate_t_and_p(means[1], means[3], ns[1], ns[3], mse, df_within)  # <65 vs >=85
results_2_3 <- calculate_t_and_p(means[2], means[3], ns[2], ns[3], mse, df_within)  # 65-84 vs >=85
# Print t-test results (the Sidak values are the same as the Keil et al paper!)
cat("Follow-Up Paired T-Tests:\n")
cat("<65 vs 65-84: t-statistic =", results_1_2$t_stat, ", p-value =", results_1_2$p_val, "sidak =", 1 - (1-results_1_2$p_val)^3, "\n")
cat("<65 vs >=85: t-statistic =", results_1_3$t_stat, ", p-value =", results_1_3$p_val, "sidak =",1 - (1-results_1_3$p_val)^3, "\n")
cat("65-84 vs >=85: t-statistic =", results_2_3$t_stat, ", p-value =", results_2_3$p_val, "sidak =",  1 - (1-results_2_3$p_val)^3, "\n")

# STEP 4: compute the selective p-values 
set.seed(826)
p_under65_6584 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                             Fstat = (results_1_2$t_stat)^2, alpha_ov = alpha)
p_under65_over85 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                               Fstat = (results_1_3$t_stat)^2, alpha_ov = alpha)
p_6584_over85 <- psel_retro(n = n, p = p, R_squared = r_squared, RSE = rse, 
                            Fstat = (results_2_3$t_stat)^2, alpha_ov = alpha)

# Output the adjusted p-values
cat("Selective p-values:\n")
cat("p_under65_6584:", p_under65_6584, "\n")
cat("p_under65_over85:", p_under65_over85, "\n")
cat("p_6584_over85:", p_6584_over85, "\n")


