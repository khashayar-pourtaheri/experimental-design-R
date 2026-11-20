## Power and sample size: paired t-test and the F-test

# Paired t-test on two measurement methods
typical <- c(-0.255, -0.213, -0.19, -0.185, -0.045, -0.025, -0.015, 0.003,
             0.015, 0.02, 0.023, 0.04, 0.04, 0.05, 0.055, 0.058)
odd     <- c(-0.324, -0.185, -0.299, -0.144, -0.027, -0.039, -0.264, -0.077,
             -0.017, -0.169, -0.096, -0.330, -0.346, -0.191, -0.128, -0.182)
d <- typical - odd
n <- length(d)
sd_d <- sd(d)

res <- t.test(typical, odd, paired = TRUE, alternative = "two.sided")
list(t_stat = unname(res$statistic), df = n - 1, p = res$p.value)

# Power against a mean difference of 0.05
mu_H1 <- 0.05
ncp <- mu_H1 / (sd_d / sqrt(n))
ncp

power_res <- power.t.test(n = n, delta = mu_H1, sd = sd_d,
                          sig.level = 0.05, type = "paired",
                          alternative = "two.sided")
power_res$power

# Sample size for 90% power
n_req <- power.t.test(power = 0.90, delta = mu_H1, sd = sd_d,
                      sig.level = 0.05, type = "paired",
                      alternative = "two.sided")$n
ceiling(n_req)

# Power of the F-test for 4 groups, tau = (1, -5, 3, 0)
tau_i <- c(1, -5, 3, 0)
k <- length(tau_i)
tau_bar <- mean(tau_i)
sigma <- 20
n <- 5
alpha <- 0.05

lambda <- n * sum((tau_i - tau_bar)^2) / sigma^2
df1 <- k - 1
df2 <- k * (n - 1)
critical_F <- qf(1 - alpha, df1, df2)
power <- 1 - pf(critical_F, df1, df2, ncp = lambda)

cat("Noncentrality parameter:", round(lambda, 3), "\n")
cat("Critical F-value:", round(critical_F, 3), "\n")
cat("Power of the F-test:", round(power, 4), "\n")

# Smallest n per group with power >= 0.8 (sigma = 20)
target_power <- 0.80

calculate_power <- function(n) {
  df2 <- k * (n - 1)
  lambda <- n * sum((tau_i - tau_bar)^2) / sigma^2
  critical_F <- qf(1 - alpha, df1, df2)
  1 - pf(critical_F, df1, df2, ncp = lambda)
}

n_values <- 2:200
powers <- sapply(n_values, calculate_power)
n_min <- n_values[which(powers >= target_power)[1]]
cat("Minimum per-group sample size for 80% power:", n_min, "\n")
cat("Total sample size:", k * n_min, "\n")

# Largest sigma in 1:20 that keeps power >= 0.8 when n = 10
n <- 10
df2 <- k * (n - 1)
critical_F <- qf(1 - alpha, df1, df2)

calculate_power_sigma <- function(sigma) {
  lambda <- n * sum((tau_i - tau_bar)^2) / sigma^2
  1 - pf(critical_F, df1, df2, ncp = lambda)
}

sigma_values <- 1:20
powers_sigma <- sapply(sigma_values, calculate_power_sigma)
max_sigma <- max(sigma_values[powers_sigma >= target_power])
cat("Largest sigma with power of at least 0.8:", max_sigma, "\n")
