#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  array_idx <- 0
} else {
  array_idx <- as.numeric(args[1])
}
cat("Currently in", getwd(), "\n")

# Startup
r_dir <- normalizePath("../R")
options(stan_dir = normalizePath("../stan"))
for (f in dir(r_dir)) {
  source(file.path(r_dir, f))
}
source("plotting_02.R")
source("predict_02.R")
outdir <- startup()
cat("\n * Results will be saved to:", outdir, "\n")

library(fda)
library(lgpr)
temp <- CanadianWeather$dailyAv[, , 1]
precip <- CanadianWeather$dailyAv[, , 3]
n_days <- nrow(temp)
n_stations <- ncol(temp)
regions <- rep(as.factor(CanadianWeather$region), each = n_days)
stations <- rep(as.factor(CanadianWeather$place), each = n_days)
day <- rep(c(1:n_days), times = n_stations)
dat <- data.frame(day, as.vector(temp), as.vector(precip), stations, regions)
colnames(dat) <- c("day", "temperature", "precipitation", "station", "region")
dat$day <- as.numeric(dat$day)

# Plot
p1 <- plot_data(dat,
  y_name = "temperature", group_by = "station",
  color_by = "station",
  x_name = "day", facet_by = "region",
  main = "Temperature"
)

p2 <- plot_data(dat,
  y_name = "precipitation", group_by = "station",
  color_by = "station",
  x_name = "day", facet_by = "region",
  main = "Precipitation (log10)"
)


# Settings
chains <- 2
iter <- 60
refresh <- 5
confs <- list(create_configuration(num_bf = 20, scale_bf = 1.5))

# Create model using lgpr
model_formula <- temperature ~ day + day | region
exact_model <- lgpr::create_model(model_formula, dat)

# Fit approximate model
afits <- sample_approx(exact_model, confs,
  refresh = refresh,
  chains = chains,
  adapt_delta = 0.95,
  iter_warmup = iter / 2,
  iter_sampling = iter / 2
)

fa <- afits[[1]]
pa <- pred_approx(fa, dat)

# Plot f or it's component
plot_f <- function(pred, idx = 0) {
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
  plt <- ggplot(df, aes(
    x = day, y = f_mean, ymin = f_mean - 2 * f_std,
    ymax = f_mean + 2 * f_std, group = region, color = region,
    fill = region
  )) +
    geom_line() +
    geom_ribbon(alpha = 0.3) +
    ylab(ylab)
}

pf1 <- plot_f(pa, 1)
pf2 <- plot_f(pa, 2)
ph <- plot_f(pa, 0)
