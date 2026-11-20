## ANOVA as a linear model: design matrix, hat matrix, OLS by hand

# Four treatments, two replicates each, constraint tau_4 = tau_1 - tau_2 - tau_3
X <- matrix(c(
  1, 1,  0,  0,
  1, 1,  0,  0,
  1, 0,  1,  0,
  1, 0,  1,  0,
  1, 0,  0,  1,
  1, 0,  0,  1,
  1, 1, -1, -1,
  1, 1, -1, -1
), nrow = 8, byrow = TRUE)

XtX_inv <- solve(t(X) %*% X)
H <- X %*% XtX_inv %*% t(X)   # hat matrix
round(H, 3)

# Weight gain of animals on three diets (cell-means coding)
weight_gains <- c(78, 86.8, 103.8, 83.7, 89, 99.2, 83.8, 81.5, 86.2)

X <- matrix(c(1, 0, 0,
              1, 0, 0,
              1, 0, 0,
              0, 1, 0,
              0, 1, 0,
              0, 1, 0,
              0, 0, 1,
              0, 0, 1,
              0, 0, 1),
            nrow = 9, byrow = TRUE)

beta_hat <- solve(t(X) %*% X) %*% t(X) %*% weight_gains
colnames(beta_hat) <- "Estimate"
rownames(beta_hat) <- c("mu_beef", "mu_pork", "mu_grain")
beta_hat

residuals <- weight_gains - X %*% beta_hat
sigma_squared <- sum(residuals^2) / (length(weight_gains) - ncol(X))

XTX_inv <- solve(t(X) %*% X)
SE_beta_1_hat <- sqrt(sigma_squared * XTX_inv[1, 1])
SE_beta_1_hat

# 90% CI for the beef mean
df <- length(weight_gains) - ncol(X)
t_crit <- qt(1 - 0.10/2, df)
beta_1_hat <- beta_hat[1]
cat("90% CI for mu_beef: (", beta_1_hat - t_crit * SE_beta_1_hat, ", ",
    beta_1_hat + t_crit * SE_beta_1_hat, ")\n")

# Same estimates from lm()
diet <- factor(rep(c("beef", "pork", "grain"), each = 3),
               levels = c("beef", "pork", "grain"))
coef(lm(weight_gains ~ diet - 1))
