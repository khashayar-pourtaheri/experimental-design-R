## Two-sample t-tests, randomization test, CI for a difference in means

library(ggplot2)

# Simulated data: equal means, unequal variances and sample sizes
set.seed(123)
nA <- 20; nB <- 10
samples_A <- rnorm(nA, mean = 0, sd = 1)
samples_B <- rnorm(nB, mean = 0, sd = 2)

pooled_t <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  s1 <- var(x); s2 <- var(y)
  sp2 <- ((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2)
  (mean(x) - mean(y)) / sqrt(sp2 * (1/n1 + 1/n2))
}

t_obs <- pooled_t(samples_A, samples_B)
t_obs
t.test(samples_A, samples_B, var.equal = TRUE)$statistic   # same value

# Randomization test: resample both groups from the pooled data
set.seed(123)
combined <- c(samples_A, samples_B)
B <- 10000

get_boot_t <- function(pool, nA, nB) {
  x <- sample(pool, nA, replace = TRUE)
  y <- sample(pool, nB, replace = TRUE)
  pooled_t(x, y)
}

t_boot <- replicate(B, get_boot_t(combined, nA, nB))
p_rand <- mean(abs(t_boot) >= abs(t_obs))
p_rand

# Pooled t-test: critical values and p-value
df <- nA + nB - 2
crit <- qt(0.975, df)
tt_pooled <- t.test(samples_A, samples_B, var.equal = TRUE, alternative = "two.sided")
list(t_stat = unname(tt_pooled$statistic),
     df     = df,
     pvalue = tt_pooled$p.value,
     crit   = c(-crit, crit))

# Welch test (does not assume equal variances)
tt_welch <- t.test(samples_A, samples_B, var.equal = FALSE, alternative = "two.sided")
t_welch  <- unname(tt_welch$statistic)
df_welch <- unname(tt_welch$parameter)
crit_welch <- qt(0.975, df = df_welch)
c(t = t_welch, df = df_welch, p = tt_welch$p.value)
c(-crit_welch, crit_welch)

df_plot <- data.frame(
  value = c(samples_A, samples_B),
  group = rep(c("A (sd=1)", "B (sd=2)"), times = c(nA, nB))
)

ggplot(df_plot, aes(x = group, y = value)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.08, alpha = 0.6) +
  labs(x = NULL, y = "Value", title = "Observed samples: Group A vs Group B") +
  theme_minimal()

# 95% CI for a difference in means by hand (sleep data), checked against t.test
d <- sleep
A <- d$extra[d$group == 1]
B <- d$extra[d$group == 2]
n1 <- length(A)
n2 <- length(B)

sp2 <- ((n1 - 1) * var(A) + (n2 - 1) * var(B)) / (n1 + n2 - 2)
se <- sqrt(sp2) * sqrt(1/n1 + 1/n2)
diff_hat <- mean(A) - mean(B)
tcrit <- qt(0.975, df = n1 + n2 - 2)

lower_end <- diff_hat - tcrit * se
upper_end <- diff_hat + tcrit * se
cat("95% CI for mu_A - mu_B (equal var, two-sample): (",
    round(lower_end, 3), ", ", round(upper_end, 3), ")\n", sep = "")

t.test(A, B, var.equal = TRUE, conf.level = 0.95)$conf.int
