library(lgpr2) # v0.0.3
library(ggplot2)
library(tidyverse)
source("utils.R")

ITER <- 600
CHAINS <- 1

# Read results
snr <- 0.5
OM <- 1
N_indiv <- 100
idx <- 131
res_dir <- paste0("res33/res-", N_indiv, "-", snr, "-", OM)
fn <- paste0(res_dir, "/res-", idx, ".rds")
res <- readRDS(file = fn)

# Subset data
data_refit <- res$dat$dat
sub_ids <- unique(data_refit$id)[1:50]
data_refit <- data_refit %>% filter(id %in% sub_ids)

# Create and compile reference model
B <- 24
scale_bf <- 1.5
full_form <- create_full_formula(res$dat$xn, res$dat$zn)
message("compiling")
model <- lgpr2::LonModel$new(formula = full_form)

message("fitting")

# Fit reference model
fit <- model$fit(
  data = data_refit, num_bf = B, scale_bf = scale_bf,
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
names(r) <- c("id", "age", res$dat$xn, res$dat$zn, "noise")
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
  dat = res$dat,
  snr = snr,
  N_indiv = N_indiv,
  relevances = r,
  term_names = model$term_names(),
  diag = diag,
  search_times = search_times,
  mcmc_time = mcmc_time
)
