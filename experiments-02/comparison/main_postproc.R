library(lgpr2)
library(ggplot2)
source("simulate.R")
source("utils.R")

args <- commandArgs(trailingOnly = TRUE)
if (interactive()) {
  idx <- 0
  OM <- 1
  N_indiv <- 40
  snr <- 0.25
} else {
  idx <- as.numeric(args[1])
  N_indiv <- as.numeric(args[2])
  snr <- as.numeric(args[3])
  OM <- as.numeric(args[4])
}
message("idx = ", idx)
message("N_indiv = ", N_indiv)
message("snr = ", snr)
message("OM = ", OM)

# Load results
res_dir <- paste0("res/res-", N_indiv, "-", snr, "-", OM)
idx <- 0
fn <- paste0(res_dir, "/res-", idx, ".rds")
a <- readRDS(fn, file = fn)

# extract
get_var_path <- function(a, fs = FALSE) {
  if (fs) {
    idx <- a$search_pp_fs$path
  } else {
    idx <- a$search_pp_dir$path
  }
  a$term_names[idx]
}

path_dir <- get_var_path(a)
path_fs <- get_var_path(a, TRUE)
a$search_pp_dir$path_chr <- path_dir
a$search_pp_fs$path_chr <- path_fs

rel <- cumsum(sort(a$relevances, decreasing = TRUE))
