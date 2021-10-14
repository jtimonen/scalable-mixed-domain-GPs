#!/usr/bin/env Rscript
library(dplyr)

# Main R script for the vote share data experiment
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  array_idx <- 0
  num_bf <- 12
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
fn_out <- file.path(outdir, paste0("res_", array_idx, ".rds"))
cat(" * results will be saved to file:", fn_out, "\n")

# Startup
source(normalizePath(file.path("..", "common.R")))
outdir <- startup()

# Load and format data
dat <- load_election_data()
plt <- lgpr::plot_data(dat, x_name="year", y_name="rep_share", 
                       group_by = "state", facet_by = "region")

# Create model
exact_model <- lgpr::create_model(rinc ~ year + year | region + year | state,
                                  dat, likelihood="binomial")


# Settings
chains <- 4
iter <- 2000
refresh <- 5
confs <- list(create_configuration(num_bf = num_bf, scale_bf = 1.5))

# Sample approximate model
afits <- sample_approx(exact_model, confs, NULL,
                       refresh = refresh,
                       chains = chains,
                       adapt_delta = 0.95,
                       iter_warmup = iter / 2,
                       iter_sampling = iter / 2,
                       parallel_chains = parallel_chains,
                       show_messages = TRUE
)

fa <- afits[[1]]
pa <- pred_approx(fa, dat)

# Save to file
results <- list(fit = fa, pred = pa)
#saveRDS(results, file = fn_out)

# Plot
pf1 <- plot_f(dat, pa, 1, aes1()) + geom_ribbon(alpha = 0.3)
pf2 <- plot_f(dat,pa, 2, aes2()) + geom_ribbon(alpha = 0.3)
pf3 <- plot_f(dat, pa, 3, aes3()) + facet_wrap(. ~ region)
ph <- plot_f(dat, pa, 0, aes3()) +
  geom_point(
    data = dat,
    inherit.aes = FALSE,
    aes(x = day, y = temperature, group = station),
    pch = ".", alpha = 0.6,
  ) + facet_wrap(. ~ station) + geom_ribbon()



