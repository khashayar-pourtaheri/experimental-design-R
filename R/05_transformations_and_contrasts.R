## Box-Cox transformation, contrasts, Bonferroni and Holm

library(MASS)
library(multcomp)

# Responses that differ by orders of magnitude across three treatments
dat <- data.frame(
  responses = c(2600, 2900, 2000, 2200, 3200,
                35000, 23000, 20000, 30000, 27000,
                2900000, 2300000, 1700000, 2900000, 2000000),
  treatment = as.factor(c(rep("1", 5), rep("2", 5), rep("3", 5)))
)

boxplot(responses ~ treatment, data = dat,
        main = "Responses by Treatment", ylab = "Responses")
tapply(dat$responses, dat$treatment, var)

anova_raw <- aov(responses ~ treatment, data = dat)
summary(anova_raw)
par(mfrow = c(2, 2)); plot(anova_raw); par(mfrow = c(1, 1))

bc <- boxcox(anova_raw, lambda = seq(-2, 2, 0.1))
lambda_opt <- bc$x[which.max(bc$y)]
lambda_opt   # close to 0, so use the log

dat$log_responses <- log(dat$responses)
anova_log <- aov(log_responses ~ treatment, data = dat)
summary(anova_log)
par(mfrow = c(2, 2)); plot(anova_log); par(mfrow = c(1, 1))

# Lentil yields under seven weed-control and fertilizer treatments
d.len <- read.table(file.path("data", "lentil.dat"), header = TRUE)
d.len$TR <- factor(d.len$TR)

fit.len <- aov(Y ~ TR, data = d.len)
summary(fit.len)

hist(residuals(fit.len), main = "Histogram of Residuals", xlab = "Residuals")
plot(fit.len$fitted.values, residuals(fit.len),
     main = "Fitted Values vs Residuals",
     xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")
qqnorm(residuals(fit.len)); qqline(residuals(fit.len), col = "red")

C1 <- c(-6,  1,  1,  1, 1,  1,  1)   # control vs all treatments
C2 <- c( 0, -1, -1, -1, 1,  1,  1)   # fertilizer vs none
C3 <- c( 0,  2, -1, -1, 2, -1, -1)   # hand weeding vs herbicide
C4 <- c( 0,  0, -1,  1, 0, -1,  1)   # herbicide before vs after
C5 <- c( 0, -2,  1,  1, 2, -1, -1)
C6 <- c( 0,  0,  1, -1, 0, -1,  1)

mat.contr <- rbind(C1 = C1, C2 = C2, C3 = C3, C4 = C4, C5 = C5, C6 = C6)

# Off-diagonal entries that are not zero mean the contrasts are not orthogonal
crossprod(t(mat.contr))

fit.mc <- glht(fit.len, linfct = mcp(TR = mat.contr))
summary(fit.mc, test = adjusted("bonferroni"))

# Response times of three circuit types
y <- c(9, 12, 10, 8, 15, 20, 21, 23, 17, 30, 6, 5, 8, 16, 7)
circ <- data.frame(Type = factor(rep(1:3, each = 5)), Y = y)

circ.fit <- aov(Y ~ Type, data = circ)
summary(circ.fit)

qqnorm(residuals(circ.fit)); qqline(residuals(circ.fit), col = "red")
plot(circ.fit$fitted.values, residuals(circ.fit),
     main = "Fitted Values vs Residuals", xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Type 2 vs the average of 1 and 3; type 1 vs type 3 (orthogonal)
circ.contr <- rbind(C1 = c(-1, 2, -1),
                    C2 = c( 1, 0, -1))
circ.mc <- glht(circ.fit, linfct = mcp(Type = circ.contr))
summary(circ.mc, test = adjusted("holm"))
