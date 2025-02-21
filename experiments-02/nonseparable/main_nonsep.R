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

# Fit
m1 <- lgpr2:::LonModel$new(formula = y ~ gp(x, z) + gp(x))
m2 <- lgpr2:::LonModel$new(formula = y ~ gp(x1) + gp(x2) + gp(x3))
f1 <- m1$fit(data = df_train, chains = 1, iter_sampling = 600)
f2 <- m2$fit(data = df_train, chains = 1, iter_sampling = 600)

r1 <- f1$predict(df)$function_draws()$get_output()
r1 <- FunctionDraws$new(df, r1, "m1")
r2 <- f2$predict(df)$function_draws()$get_output()
r2 <- FunctionDraws$new(df, r2, "m2")
ymax <- max(df$y) + 0.3 * sd(df$y)
ymin <- min(df$y) - 0.3 * sd(df$y)
p1 <- r1$plot(x_var = "x", color_by = NULL) +
  geom_point(
    data = df_train, mapping = aes(x = x, y = y),
    inherit.aes = FALSE
  ) +
  geom_point(
    data = df_test, mapping = aes(x = x, y = y),
    inherit.aes = FALSE, pch = 4
  ) + ylim(ymin, ymax)
p2 <- r2$plot(x_var = "x", color_by = NULL) +
  geom_point(
    data = df_train, mapping = aes(x = x, y = y),
    inherit.aes = FALSE
  ) +
  geom_point(
    data = df_test, mapping = aes(x = x, y = y),
    inherit.aes = FALSE, pch = 4
  ) + ylim(ymin, ymax)
plt <- ggarrange(p1, p2, nrow = 2)
