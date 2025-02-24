library(lgpr2) # v0.0.3
library(lgpr)
library(tidyverse)
library(MASS)
library(ggpubr)
ggplot2::theme_set(ggplot2::theme_bw())
source("functions_nonsep.R")

X <- simulate_input()
K <- nonsep_kernel(X$x, X$x, X$z, X$z)


df <- simulate_nonsep(X)
df <- simulate_obs(df)
df <- create_dummy_x(df)
plt <- ggplot(df, aes(x = x, y = f, color = z)) +
  geom_line() +
  geom_point(mapping = aes(x = x, y = y, color = z)) +
  facet_grid(. ~ z)
a <- split_train_test(df)
df_train <- a$train
df_test <- a$test
df <- a$full

y_pred <- gppred(df_train, df, 0.3, 1e-8)
df$y_pred_true <- y_pred


# Fit
m1 <- lgpr2:::LonModel$new(formula = y ~ gp(x, z) + gp(x))
m2 <- lgpr2:::LonModel$new(formula = y ~ gp(x1) + gp(x2) + gp(x3))
f1 <- m1$fit(data = df_train, chains = 1, iter_sampling = 600)
f2 <- m2$fit(data = df_train, chains = 1, iter_sampling = 600)

r1 <- f1$predict(df)$function_draws()$get_output()
r1 <- FunctionDraws$new(df, r1, "Shared kernel parameters")
r2 <- f2$predict(df)$function_draws()$get_output()
r2 <- FunctionDraws$new(df, r2, "Group-specific kernel parameters")
ymax <- max(df$y) + 1.2 * sd(df$y)
ymin <- min(df$y) - 1.2 * sd(df$y)



p1 <- plot_fit(r1, df, df_train, df_test)
p2 <- plot_fit(r2, df, df_train, df_test)

plt <- ggarrange(p1, p2, nrow = 2)
err <- compute_accuracy(df, r1, r2)
