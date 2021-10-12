#!/usr/bin/env Rscript

# Main R script for the weather data experiment
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  array_idx <- 0
} else {
  array_idx <- as.numeric(args[1])
}
cat("Currently in", getwd(), "\n")

# Startup
source(normalizePath(file.path("..", "common.R")))
outdir <- startup()
dat <- load_weather_data()

# Split data
N <- nrow(dat)
# inds <-  sample.int(100, size=101, replace = FALSE)
#i_test <- which((dat$day > 180) & (dat$day <= 250))
rn <- CanadianWeather$region
i1 <- which(rn == "Atlantic")[1]
i2 <- which(rn == "Arctic")[2]
i3 <- which(rn == "Continental")[3]
i4 <- which(rn == "Pacific")[4]
test_stations <- names(c(i1, i2, i3, i4))
i_test <- which((dat$day > 100) & (dat$station %in% test_stations))
split <- lgpr:::split_data(dat, i_test = i_test)
train_dat <- split$train
test_dat <- split$test
print(nrow(train_dat))
print(nrow(test_dat))

# Settings
chains <- 1
iter <- 60
refresh <- 1
confs <- list(create_configuration(num_bf = 12, scale_bf = 1.5))

# Create model using lgpr
model_formula <- temperature ~ day + day | region + day | station
exact_model <- lgpr::create_model(model_formula, train_dat)

# Needs a bit of editing
# n_stat <- 35
# z_stat <- exact_model@stan_input$x_cat[2, ]
# miss <- c(1:n_stat) %in% unique(z_stat)
# i_missing <- which(miss == FALSE)
# z_stat[which(z_stat == n_stat)] <- i_missing
# exact_model@stan_input$x_cat[2, ] <- z_stat

# Sample approximate model
afits <- sample_approx(exact_model, confs, NULL,
  refresh = refresh,
  chains = chains,
  adapt_delta = 0.95,
  iter_warmup = iter / 2,
  iter_sampling = iter / 2
)

fa <- afits[[1]]
pa <- pred_approx(fa, dat)

# Plot f or it's component
plot_f <- function(pred, idx = 0, aesth) {
  if (idx == 0) {
    f <- pred@h
    ylab <- "h"
  } else {
    f <- pred@f_comp[[idx]]
    ylab <- names(pred@f_comp)[idx]
  }
  f_mean <- colMeans(f)
  f_std <- apply(f, 2, stats::sd)
  df <- cbind(dat, f_mean, f_std)
  plt <- ggplot(df, aesth) +
    geom_line() +
    ylab(ylab)
}


aes1 <- aes(
  x = day, y = f_mean, ymin = f_mean - 2 * f_std,
  ymax = f_mean + 2 * f_std
)
aes2 <- aes(
  x = day, y = f_mean, ymin = f_mean - 2 * f_std,
  ymax = f_mean + 2 * f_std, group = region, color = region,
  fill = region
)
aes3 <- aes(
  x = day, y = f_mean, ymin = f_mean - 2 * f_std,
  ymax = f_mean + 2 * f_std, group = station, color = station,
  fill = station
)

pf1 <- plot_f(pa, 1, aes1) + geom_ribbon(alpha = 0.3)
pf2 <- plot_f(pa, 2, aes2) + geom_ribbon(alpha = 0.3)
pf3 <- plot_f(pa, 3, aes3) + facet_wrap(. ~ region)
ph <- plot_f(pa, 0, aes3) +
  geom_point(
    data = dat,
    inherit.aes = FALSE,
    aes(x = day, y = temperature, group = station),
    pch = ".", alpha = 0.6,
  ) + facet_wrap(. ~ station) + geom_ribbon()
