## Block designs: RCBD vs paired t-test, incomplete block design

library(ggplot2)

# Platelet aggregation before and after smoking, 11 people
smoke_data <- readRDS(file.path("data", "smoke.RDS"))
str(smoke_data)

ggplot(smoke_data, aes(x = PERIOD, y = AGGREG)) +
  geom_boxplot() +
  labs(x = "Period (1 = Before, 2 = After)", y = "Percent Aggregated Platelets") +
  theme_minimal()

ggplot(smoke_data, aes(x = PERIOD, y = AGGREG, group = PERSON, color = PERSON)) +
  geom_line() +
  geom_point() +
  labs(x = "Period (1 = Before, 2 = After)", y = "Percent Aggregated Platelets",
       title = "Interaction plot: person by period") +
  theme_minimal()

# One observation per person and period, so no interaction term
anova_model <- aov(AGGREG ~ PERSON + PERIOD, data = smoke_data)
anova(anova_model)

qqnorm(residuals(anova_model)); qqline(residuals(anova_model), col = "red")
shapiro.test(residuals(anova_model))

plot(fitted(anova_model), residuals(anova_model),
     xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, col = "red")

# Paired t-test gives the same answer as the block design
before <- smoke_data$AGGREG[smoke_data$PERIOD == 1]
after  <- smoke_data$AGGREG[smoke_data$PERIOD == 2]
paired_test <- t.test(before, after, paired = TRUE)
paired_test

F_stat <- anova(anova_model)["PERIOD", "F value"]
t_stat <- unname(paired_test$statistic)
c(F = F_stat, t_squared = t_stat^2)

# Two-way layout with no interaction: the lines in the interaction plot are parallel
interaction_data <- data.frame(
  music = factor(rep(1:3, each = 2)),
  snack = rep(1:2, 3),
  time  = c(370, 376, 380, 386, 374, 380)
)

ggplot(interaction_data, aes(x = snack, y = time, group = music, color = music)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(x = "Snack", y = "Mean Time Worked", color = "Music Type") +
  theme_minimal()

# Is this incomplete block design balanced?
ibd <- readRDS(file.path("data", "ibd.RDS"))
tab <- xtabs(~ block + treatment, data = ibd)
tab

# Number of blocks in which each pair of treatments appears together
concurrence <- crossprod(tab)
concurrence

pairs <- concurrence[upper.tri(concurrence)]
table(pairs)   # pairs that never meet, or meet unequally often, mean it is not balanced

g <- nlevels(ibd$treatment)
k <- max(rowSums(tab))
r <- unique(colSums(tab))
r * (k - 1) / (g - 1)   # lambda would have to be a whole number
