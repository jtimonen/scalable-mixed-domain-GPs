#!/usr/bin/env Rscript
library(dplyr)
library(readr)

# Main R script for the vote share data experiment
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  array_idx <- 0
  num_bf <- 24
} else {
  array_idx <- as.numeric(args[1])
  num_bf <- array_idx
}
available_cores <- as.integer(Sys.getenv(
  "SLURM_JOB_CPUS_PER_NODE",
  parallel::detectCores()
))
parallel_chains <- available_cores
cat("Currently in", getwd(), "\n")
cat("num_bf =", num_bf, "\n")


# Startup
source(normalizePath(file.path("..", "common.R")))
outdir <- startup()
fn_out <- file.path(outdir, paste0("res_", array_idx, ".rds"))
cat(" * results will be saved to file:", fn_out, "\n")

# Load and format data
dat <- load_election_data()
plt <- lgpr::plot_data(dat,
  x_name = "year", y_name = "rep_share",
  group_by = "state", facet_by = "region"
)

# Create model
form <- rep_share ~ year + year | region + year | state
exact_model <- lgpr::create_model(form, dat,
  sample_f = TRUE
)
# exact_model@stan_input$x_cont <- exact_model@stan_input$x_cont_unnorm

# Settings
chains <- 4
iter <- 2000
refresh <- 5
confs <- list(create_configuration(num_bf = num_bf, scale_bf = 2.0))

# Sample approximate model
afits <- sample_approx_beta(exact_model, confs, dat,
  refresh = refresh,
  chains = chains,
  adapt_delta = 0.95,
  iter_warmup = iter / 2,
  iter_sampling = iter / 2,
  parallel_chains = parallel_chains,
  show_messages = TRUE
)

fa <- afits[[1]]

# Create test points and predict
x_star <- create_test_x(dat, seq(1976, 2024, by = 0.5))

pa <- pred_approx(fa, x_star)
y_rng <- pa$y_rng
pa <- pa$pred

# Save to file
results <- list(fit = fa, pred = pa)
# saveRDS(results, file = fn_out)

# Create y-axis labels for plots
ylab1 <- expression(paste("", "f", phantom()^{
  paste("", "(", "1", ")", "")
}, "(", "", "year", ")", "", ""))
ylab2 <- expression(paste("", "f", phantom()^{
  paste("", "(", "2", ")", "")
}, "(", "", "year", ",", "region", ")", "", ""))
ylab3 <- expression(paste("", "f", phantom()^{
  paste("", "(", "3", ")", "")
}, "(", "", "year", ",", "state", ")", "", ""))

# Plot
scc <- scale_x_continuous(breaks = c(unique(dat$year), 2020, 2024))
pf1 <- plot_f(x_star, pa, 1, aes1()) + geom_ribbon(alpha = 0.3) + scc +
  ylab(ylab1)
pf2 <- plot_f(x_star, pa, 2, aes2()) + scc + ylab(ylab2)
pf3 <- plot_f(x_star, pa, 3, aes3()) + facet_wrap(. ~ region) + scc +
  ylab(ylab3)


y_mean <- colMeans(y_rng)
y_std <- apply(y_rng, 2, stats::sd)
df <- cbind(x_star, y_mean, y_std)

# Plot predictions
dat2020 <- read_data_2020()
py <- ggplot(df, aes(
  x = year, y = y_mean, ymin = y_mean - 2 * y_std,
  ymax = y_mean + 2 * y_std, group = state
)) +
  geom_line() +
  ylab("Republican vote share") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("") +
  geom_hline(yintercept = 0.5, lty = 1, col = "steelblue") +
  geom_point(
    data = dat,
    inherit.aes = FALSE,
    aes(x = year, y = rep_share, group = state),
    pch = 20, alpha = 1, color = "steelblue"
  ) +
  facet_wrap(. ~ state) +
  geom_ribbon(alpha = 0.4) +
  geom_point(
    data = dat2020,
    inherit.aes = FALSE,
    aes(x = year, y = rep_share, group = state),
    pch = 4, alpha = 1, color = "red"
  ) +
  theme(strip.text.x = element_text(size = 7))

ggsave("all_states_large.pdf", py, width = 12, height = 10)
fit <- get_cmdstanfit(fa)

# Save
total_time <- fit$time()$total
cat(total_time, file = "runtime.txt")
fit$save_object(fn_out)

# Full components plot
plt12 <- ggarrange(pf1, pf2,
  nrow = 1, ncol = 2, labels = c("a", "b"),
  widths = c(0.7, 1.0)
)

plt_full <- ggarrange(plt12, pf3,
  heights = c(0.53, 1.0),
  labels = c("", "c"), nrow = 2, ncol = 1
)

ggsave(file = "components.pdf", plt_full, width = 7.65, height = 7.88)
