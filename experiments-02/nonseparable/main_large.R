library(lgpr2) # v0.0.3
library(lgpr)
library(tidyverse)
library(MASS)
library(ggpubr)
ggplot2::theme_set(ggplot2::theme_bw())
source("functions_e2.R")

sigma_true <- 0.1

set.seed(123)
X <- simulate_input_vary(100, 10)
df <- simulate_fast(X)
df <- simulate_obs(df, sigma_true)
df <- create_dummy_x(df)
plt <- ggplot(df, aes(x = x, y = f, color = z)) +
  geom_line() +
  geom_point(mapping = aes(x = x, y = y, color = z)) +
  facet_grid(. ~ z)
a <- split_train_test(df, test_categ = c(2, 3))
df_train <- a$train
df_test <- a$test
df <- a$full

# Fit
m2 <- lgpr2:::LonModel$new(formula = y ~ gp(x, z) + gp(x))
t_start <- now()
f2 <- m2$fit(
  data = df_train, chains = 1, iter_sampling = 30, iter_warmup = 30,
  refresh = 1
)
t_end <- now()
t_dur <- t_end - t_start

r2 <- f2$predict(df)$function_draws()$get_output()
r2 <- FunctionDraws$new(df, r2, "Kernel 1")

p2 <- plot_fit(r2, df, df_train, df_test, FALSE, TRUE)
