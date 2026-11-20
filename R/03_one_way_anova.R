## One-way ANOVA, residual checks, permutation test

# Blood coagulation time under four diets
A <- c(62, 60, 63, 59, 64)
B <- c(65, 67, 73, 65, 66)
C <- c(69, 66, 71, 67, 67, 68, 62)
D <- c(66, 62, 65, 61, 64, 65, 63)

dat <- data.frame(
  Diet = factor(rep(c("A", "B", "C", "D"),
                    times = c(length(A), length(B), length(C), length(D)))),
  Coag = c(A, B, C, D)
)

boxplot(Coag ~ Diet, data = dat,
        main = "Blood Coagulation Times by Diet",
        ylab = "Coagulation Time", xlab = "Diet")

means <- tapply(dat$Coag, dat$Diet, mean)
grand <- mean(dat$Coag)
means; grand

# Pooled variance (MSE) and treatment mean square (MST)
vars <- tapply(dat$Coag, dat$Diet, var)
ns <- table(dat$Diet)
m <- length(ns)
N <- nrow(dat)

MSE <- sum((ns - 1) * vars) / (N - m)
SST <- sum(ns * (means - grand)^2)
MST <- SST / (m - 1)
c(SST = SST, MST = MST, MSE = MSE)

fit <- aov(Coag ~ Diet, data = dat)
summary(fit)

Fval <- MST / MSE
pval <- pf(Fval, m - 1, N - m, lower.tail = FALSE)
c(F = Fval, df1 = m - 1, df2 = N - m, p = pval)

# Residuals around the group means
mu_hat_i <- means[dat$Diet]
residuals_explicit <- dat$Coag - as.numeric(mu_hat_i)

qqnorm(residuals_explicit, main = "Q-Q Plot of Residuals")
qqline(residuals_explicit)

plot(as.numeric(mu_hat_i), residuals_explicit,
     xlab = "Fitted values (group means)", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2)

# Species diversity at four zinc levels
zinc <- readRDS(file.path("data", "zinc.RDS"))

boxplot(DIVERSITY ~ ZINC, data = zinc,
        main = "Biodiversity by Zinc Level",
        ylab = "Diversity", xlab = "Zinc Level",
        col = "lightblue")

group_means <- aggregate(DIVERSITY ~ ZINC, data = zinc, mean)
overall_mean <- mean(zinc$DIVERSITY)
group_means; overall_mean

group_stats <- aggregate(DIVERSITY ~ ZINC, data = zinc,
                         FUN = function(x) c(n = length(x), var = var(x)))
group_stats <- do.call(data.frame, group_stats)
colnames(group_stats) <- c("ZINC", "n_i", "s_i2")

MSE <- sum((group_stats$n_i - 1) * group_stats$s_i2) /
  (sum(group_stats$n_i) - nrow(group_stats))

k <- nrow(group_means)
SS_between <- sum(group_stats$n_i * (group_means$DIVERSITY - overall_mean)^2)
MST <- SS_between / (k - 1)
F_ratio <- MST / MSE
c(MST = MST, MSE = MSE, F = F_ratio)

# Permutation test: shuffle the outcome across zinc levels
set.seed(502)
n_samples <- 1000
random_F <- numeric(n_samples)

for (i in seq_len(n_samples)) {
  shuffled_data <- data.frame(ZINC = zinc$ZINC, DIVERSITY = sample(zinc$DIVERSITY))

  g_means <- aggregate(DIVERSITY ~ ZINC, data = shuffled_data, mean)
  overall <- mean(shuffled_data$DIVERSITY)
  SSb <- sum(table(shuffled_data$ZINC) * (g_means$DIVERSITY - overall)^2)
  MST_perm <- SSb / (k - 1)

  g_vars <- aggregate(DIVERSITY ~ ZINC, data = shuffled_data, var)
  n_i <- table(shuffled_data$ZINC)
  MSE_perm <- sum((n_i - 1) * g_vars$DIVERSITY) / (sum(n_i) - k)

  random_F[i] <- MST_perm / MSE_perm
}

p_value <- mean(random_F >= F_ratio)
cat("P-value from randomization test:", round(p_value, 4), "\n")

hist(random_F, breaks = 30, col = "gray",
     main = "Randomization Distribution of F-statistic",
     xlab = "F-statistic", border = "white")
abline(v = F_ratio, col = "red", lwd = 2, lty = 2)
legend("topright", legend = "Observed F", col = "red", lty = 2, bty = "n")

# What departures from normality look like in a Q-Q plot
set.seed(123)
x_long  <- rt(1000, df = 2)                          # heavy tails
x_short <- runif(1000, -3, 3)                        # light tails
x_left  <- -rlnorm(1000, meanlog = 0, sdlog = 1)     # left skew
x_right <-  rlnorm(1000, meanlog = 0, sdlog = 1)     # right skew

par(mfrow = c(2, 2))
qqnorm(x_long,  main = "Long-tailed (t df=2)");      qqline(x_long,  col = "red")
qqnorm(x_short, main = "Short-tailed (Uniform)");    qqline(x_short, col = "red")
qqnorm(x_left,  main = "Left-skewed (-Log-normal)"); qqline(x_left,  col = "red")
qqnorm(x_right, main = "Right-skewed (Log-normal)"); qqline(x_right, col = "red")
par(mfrow = c(1, 1))
