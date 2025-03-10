library(lgpr2) # v0.0.3
library(lgpr)
library(tidyverse)
library(MASS)
library(ggpubr)
ggplot2::theme_set(ggplot2::theme_bw())
source("functions_e2.R")

sigma_true <- 0.1

set.seed(123)
X <- simulate_input_vary(100, 8)
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
f2 <- m2$fit(data = df_train, chains = 1, iter_sampling = 30, iter_warmup = 30)

r2 <- f2$predict(df)$function_draws()$get_output()
r2 <- FunctionDraws$new(df, r2, "Kernel 1")

p2 <- r2$plot() + geom_point(
  data = df_train, mapping = aes(x = x, y = y),
  inherit.aes = FALSE, color = "black"
) +
  geom_point(
    data = df_test, mapping = aes(x = x, y = y),
    inherit.aes = FALSE, color = "firebrick", pch = 4
  )
