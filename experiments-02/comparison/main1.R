library(lgpr2) # v0.0.3
library(ggplot2)
source("simulate.R")
source("utils.R")

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  idx <- 101
  OM <- 1
  N_indiv <- 30
} else {
  idx <- as.numeric(args[1])
  N_indiv <- as.numeric(args[2])
  OM <- as.numeric(args[3])
}
message("idx = ", idx)
message("N_indiv = ", N_indiv)
message("OM = ", OM)

# Setup depending on idx
CHAINS <- 2
ITER <- 1000
f_var <- 16
n_unrel <- 16
if (idx <= 50) {
  snr <- 0.1
} else if (idx <= 100) {
  snr <- 0.25
} else if (idx <= 150) {
  snr <- 0.5
} else {
  stop("too large idx!")
}
message("n_unrel = ", n_unrel)

# Simulate data
t_seq <- seq(6, 96, by = 6)
data_new <- simulate_gaussian(N_indiv, n_unrel, snr, t_seq)
true_snr <- compute_snr(data_new$components)
message("true_snr = ", true_snr)


# Create and compile reference model
B <- 24
scale_bf <- 1.5
full_form <- create_full_formula(data_new$xn, data_new$zn)
message("compiling")
model <- lgpr2::LonModel$new(formula = full_form)

message("fitting")

# Fit reference model
fit <- model$fit(
  data = data_new$dat, num_bf = B, scale_bf = scale_bf,
  iter_sampling = ITER, iter_warmup = ITER, chains = CHAINS,
  adapt_delta = 0.99, init = 0.1
)
message("diagnosing")
diag <- fit$diagnose()

message("relevance")

# Run selection with projection predictive method
r <- fit$relevances() # id, age, x..., z...
r <- mean(r)
r_path <- cumsum(sort(r, decreasing = TRUE))
names(r) <- c("id", "age", data_new$xn, data_new$zn, "noise")
D <- length(r) - 1
rels <- sort(r[1:D], index.return = TRUE, decreasing = TRUE)
path <- rels$ix
notidage <- !(path %in% c(1, 2))
idage_first <- c(1, 2, path[which(notidage)])
just_idage <- c(1, 2)

# Run searches
start_time <- Sys.time()
search_pp_fs <- pp_forward_search(fit, path = just_idage, num_steps = 6, B = B)
t_fs <- Sys.time() - start_time
start_time <- Sys.time()
search_pp_dir <- pp_forward_search(fit, path = idage_first, num_steps = 6, B = B)
t_dir <- Sys.time() - start_time
search_times <- c(t_fs, t_dir)
mcmc_time <- fit$get_stan_fit()$time()$total

# Results list
res <- list(
  r_path = r_path,
  search_pp_fs = search_pp_fs,
  search_pp_dir = search_pp_dir,
  dat = data_new,
  snr = true_snr,
  N_indiv = N_indiv,
  relevances = r,
  term_names = model$term_names(),
  diag = diag,
  search_times = search_times,
  mcmc_time = mcmc_time
)

# Save results
res_dir <- paste0("res/res-", N_indiv, "-", snr, "-", OM)
fig_dir <- paste0("res/figs-", N_indiv, "-", snr, "-", OM)
if (!dir.exists("res")) {
  dir.create("res")
}
if (!dir.exists(res_dir)) {
  dir.create(res_dir)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir)
}
fn <- paste0(res_dir, "/res-", idx, ".rds")
fn_fig <- paste0(fig_dir, "/fit-", idx, ".pdf")
saveRDS(res, file = fn)
ggsave(fit$plot(), width = 10, height = 7, file = fn_fig)
