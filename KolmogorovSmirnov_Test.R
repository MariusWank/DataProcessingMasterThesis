# Load necessary libraries
library(fitdistrplus)
library(stats4)

# Manually transcribed data from the image
data <- c(
  101.44, 86.46, 60.08, 45.80, 70.34, 52.29, 40.73, 52.14, 86.44, 80.21,
  86.18, 82.54, 50.36, 74.10, 52.31, 65.13, 88.96, 21.63, 82.64, 141.74,
  106.27, 112.69, 32.64, 98.11, 59.74, 101.41, 93.16, 136.16, 96.89, 118.68,
  102.70, 47.51, 72.82, 40.43, 106.70, 75.41, 101.56, 84.14, 114.89, 96.61,
  102.48, 45.73, 88.11, 70.15, 79.30, 48.31, 56.64, 75.94, 83.20, 119.86,
  86.62, 51.63, 92.45, 125.19, 90.21, 148.50, 109.92, 125.74, 68.25, 47.91,
  33.79, 64.48, 88.78, 53.55, 78.00, 60.87, 58.74, 56.35, 63.09, 80.67,
  34.27, 37.52, 51.49, 55.51, 36.00, 64.42, 82.33, 27.03, 67.53, 56.42,
  55.59, 32.11, 96.36, 57.90, 57.06, 40.28, 83.53, 66.04, 63.98, 54.25,
  92.39, 115.95, 110.18, 99.45, 100.33, 126.58, 116.82, 134.59, 93.34, 114.39,
  104.61, 113.28, 114.96, 149.31, 104.94, 88.72, 75.34, 86.11, 131.55, 103.64,
  118.44, 90.90, 42.27, 135.12, 24.98, 112.06, 119.71, 58.32, 110.21, 42.60,
  63.27, 63.09, 169.17, 96.36, 36.27, 71.38, 40.76, 102.96, 71.72, 162.61,
  122.42, 114.94, 77.95, 151.23, 94.73, 149.75, 136.95, 156.23, 102.91, 65.53,
  111.35, 87.42, 150.60, 101.86, 115.22, 82.19, 138.43, 94.06, 126.14, 64.21,
  70.59, 122.88, 124.40, 51.70, 46.99, 115.41, 95.63, 96.80, 138.20, 108.46,
  128.68, 129.96
)

# Given mean and standard deviation for normal distribution
mean_given <- mean(data)  # Using sample mean
sd_given <- sd(data)  # Using sample standard deviation

# Fit Weibull distribution
weibull_fit <- fitdist(data, "weibull")

# Fit lognormal distribution
lognormal_fit <- fitdist(data, "lnorm")

# Log-likelihood for given normal distribution
logLik_normal <- sum(dnorm(data, mean = mean_given, sd = sd_given, log = TRUE))

# AIC for normal distribution
aic_normal <- -2 * logLik_normal + 2 * 2  # 2 parameters: mean and sd

# AIC for Weibull distribution
aic_weibull <- 2 * length(weibull_fit$estimate) - 2 * weibull_fit$loglik

# AIC for lognormal distribution
aic_lognormal <- 2 * length(lognormal_fit$estimate) - 2 * lognormal_fit$loglik

# KS test for normal distribution
ks_normal <- ks.test(data, "pnorm", mean = mean_given, sd = sd_given)
p_value_normal <- ks_normal$p.value

# KS test for Weibull distribution
ks_weibull <- ks.test(data, "pweibull", shape = weibull_fit$estimate["shape"], scale = weibull_fit$estimate["scale"])
p_value_weibull <- ks_weibull$p.value

# KS test for lognormal distribution
ks_lognormal <- ks.test(data, "plnorm", meanlog = lognormal_fit$estimate["meanlog"], sdlog = lognormal_fit$estimate["sdlog"])
p_value_lognormal <- ks_lognormal$p.value

# Combine AIC values and p-values into a data frame
results <- data.frame(
  Model = c("Normal", "Weibull", "Lognormal"),
  AIC = c(aic_normal, aic_weibull, aic_lognormal),
  P_Value = c(p_value_normal, p_value_weibull, p_value_lognormal)
)

# Rank models based on AIC
results <- results[order(results$AIC), ]
results$Rank <- 1:nrow(results)

# Print results
print(results)
