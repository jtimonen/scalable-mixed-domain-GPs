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

# Fit
m1 <- lgpr2:::LonModel$new(formula = y ~ gp(x, z) + gp(x))
m2 <- lgpr2:::LonModel$new(formula = y ~ gp(x1) + gp(x2) + gp(x3))
f1 <- m1$fit(data = df, chains = 1, iter_sampling = 600)
f2 <- m2$fit(data = df, chains = 1, iter_sampling = 600)

p1 <- f1$plot()
p2 <- f2$plot()
