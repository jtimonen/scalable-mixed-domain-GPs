library(lgpr2)
source("simulate.R")
source("utils.R")

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  idx <- 0
  OM <- 1
} else {
  idx <- as.numeric(args[1])
  OM <- as.numeric(args[2])
}
message("idx = ", idx)
message("OM = ", OM)

# Setup depending on idx
CHAINS <- 1
ITER <- 1000
f_var <- 16
N_indiv <- 12
if (idx <= 40) {
  n_unrel <- 6
} else if (idx <= 80) {
  n_unrel <- 8
} else if (idx <= 120) {
  n_unrel <- 10
} else {
  stop("too large idx!")
}
message("n_unrel = ", n_unrel)

# Simulate data
t_seq <- seq(3, 96, by = 3)
snr <- 0.2
data_new <- simulate_gaussian(N_indiv, n_unrel, snr, t_seq)
true_snr <- compute_snr(data_new$components)
message("true_snr = ", true_snr)


# Create and compile reference model
B <- 20
scale_bf <- 1.5
full_form <- create_full_formula(data_new$xn, data_new$zn)
model <- lgpr2::LonModel$new(formula = full_form)

# Fit reference model
fit <- model$fit(
  data = data_new$dat, num_bf = B, scale_bf = scale_bf,
  iter_sampling = ITER, iter_warmup = ITER, chains = CHAINS,
  adapt_delta = 0.99, init = 0.1
)
diag <- diagnose(fit)

# Run selection with projection predictive method
r <- fit$relevances() # id, age, x..., z...
r_path <- cumsum(sort(r, decreasing = TRUE))
sel_95 <- get_selected(r_path, 0.95)
sel_90 <- get_selected(r_path, 0.90)
sel_85 <- get_selected(r_path, 0.85)
names(r) <- c("id", "age", xn, zn, "noise")
D <- length(r) - 1
rels <- sort(r[1:D], index.return = TRUE, decreasing = TRUE)
path <- names(rels$x)
path <- rm$covar_names_to_covar_inds(path)
search_pp_fs <- selection_pp(rm, path = NULL)
search_pp_dir <- selection_pp(rm, path = path)


# Results list
tr <- c(0.85, 0.9, 0.95)
res <- list(
  r_path = r_path,
  search_pp_fs = search_pp_fs,
  search_pp_dir = search_pp_dir,
  sel_dir = results_dir(r_path, tr),
  sel_pp_fs = results_pp(search_pp_fs, tr),
  sel_pp_dir = results_pp(search_pp_dir, tr),
  dat = dat,
  snr = true_snr
)

# Save results
res_dir <- paste0("res-", likelihood)
fig_dir <- paste0("figs-", likelihood)
if (!dir.exists(res_dir)) {
  dir.create(res_dir)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir)
}
fn <- paste0(res_dir, "/res-", idx, ".rds")
fn_fig <- paste0(fig_dir, "/fit-", idx, ".pdf")
saveRDS(res, file = fn)
ggsave(rm$plot_by_id(), width = 10, height = 7, file = fn_fig)
