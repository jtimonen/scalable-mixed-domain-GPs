library(lgpr2) # v0.0.3
library(lgpr)
library(tidyverse)
library(MASS)
library(ggpubr)
ggplot2::theme_set(ggplot2::theme_bw())
source("functions_e2.R")

ell_true <- c(1, 0.5, 0.2)
sigma_true <- 0.3

set.seed(412345)
X <- simulate_input()
df <- simulate_nonsep(X, ell = ell_true)
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

y_pred <- gppred(df_train, df, sigma_true, 1e-8, ell_true)
df$y_pred_true <- y_pred


# Fit
m1 <- lgpr2:::LonModel$new(formula = y ~ gp(x, z))
m2 <- lgpr2:::LonModel$new(formula = y ~ gp(x, z) + gp(x))
f1 <- m1$fit(data = df_train, chains = 3, iter_sampling = 1200)
f2 <- m2$fit(data = df_train, chains = 3, iter_sampling = 1200)

r1 <- f1$predict(df)$function_draws()$get_output()
r1 <- FunctionDraws$new(df, r1, "Kernel 1")
r2 <- f2$predict(df)$function_draws()$get_output()
r2 <- FunctionDraws$new(df, r2, "Kernel 2")
ymax <- max(df$y) + 1.5 * sd(df$y)
ymin <- min(df$y) - 1.5 * sd(df$y)

p1 <- plot_fit(r1, df, df_train, df_test)
p2 <- plot_fit(r2, df, df_train, df_test)

plt <- ggarrange(p1, p2, nrow = 2, labels = c("a)", "b)"))
err <- compute_accuracy(df, r1, r2)

ggsave(plt, file = "nonsep.pdf", width = 7, height = 4.5)
