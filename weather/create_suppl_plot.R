# Startup
source(normalizePath(file.path("..", "common.R")))
outdir <- startup(create_dir = FALSE)
num_bf <- 32

# Load sampling results
fn <- file.path("results_temperature", paste0("res_", num_bf, ".rds"))
dat <- load_weather_data()
dat$day <- as.Date(dat$day, origin = "1964-01-01")
results <- readRDS(file = fn)
pa <- results$pred
fit <- get_cmdstanfit(results$fit)
total_time <- fit$time()$total

# Need this because this stores information how to scale from
# inference scale (normalized) to original data scale
emod <- results$fit@model@exact_model
y_scl <- emod@var_scalings$y

# Create new fit
of <- fit$output_files()
get_fn <- function(x) {
  parts <- strsplit(x[1], "/")[[1]]
  file.path("results_temperature", parts[length(parts)])
}
of <- c(get_fn(of[1]), get_fn(of[2]), get_fn(of[3]), get_fn(of[4]))
fit_new <- cmdstanr::as_cmdstan_fit(files = of)

# Sample from predictive distribution
draws_f <- pa@f
P <- ncol(draws_f)
draws_sig <- as.vector(posterior::merge_chains(fit_new$draws("sigma")))
draws_y <- array(0.0, dim = dim(draws_f))
for (s in 1:length(draws_sig)) {
  y_draw <- rnorm(n = P, mean = draws_f[s, ], sd = draws_sig[s])
  draws_y[s, ] <- lgpr:::apply_scaling(y_scl, y_draw, inverse = TRUE)
}
yp_mean <- colMeans(draws_y)
yp_std <- apply(draws_y, 2, stats::sd)

# Plot predictive distribution
df <- cbind(dat, yp_mean, yp_std)
colnames(df) <- c(colnames(dat), "y_mean", "y_std")
plt <- ggplot(df, aes(
  x = day, y = y_mean,
  ymin = y_mean - 2 * y_std,
  ymax = y_mean + 2 * y_std,
  group = station
)) +
  geom_line(alpha = 0.7) +
  ylab("Temperature (C)") +
  theme_bw() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme(panel.grid.minor = element_blank(), legend.position = "top") +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("") +
  facet_wrap(. ~ station) +
  geom_ribbon(fill = "steelblue", alpha = 0.7, color = NA) +
  geom_point(
    mapping = aes(x = day, y = temperature, group = station),
    pch = 20, alpha = 0.55, size = 0.06
  )

fn_out <- file.path("figures", "weather_suppl.pdf")
ggsave(plt, file = fn_out, width = 11, height = 8.3)
