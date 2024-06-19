library(lgpr2)
library(tidyverse)

ls <- dir("res")
ls <- ls[grepl(pattern = "res-", ls)]

# Get scalar results
get_scalar_res <- function(dir) {
  rdir <- file.path("res", dir)
  files <- dir(rdir)
  out <- NULL
  for (j in seq_len(length(files))) {
    fn <- file.path(rdir, files[j])
    out <- rbind(out, res_to_df(fn))
  }
  tibble::as_tibble(out)
}

# Scalar results to data frame
res_to_df <- function(fn) {
  r <- readRDS(fn)
  x <- unlist(r[c("N_indiv", "snr")])
  x_diag <- r[["diag"]]
  num_terms <- length(r$term_names)
  x <- c(x, x_diag, num_terms)
  names(x)[length(x)] <- "num_terms"
  x <- data.frame(t(c(basename(fn), x)))
  colnames(x)[1] <- "file"
  x
}


# Get scalar results
get_path_res <- function(dir) {
  rdir <- file.path("res", dir)
  files <- dir(rdir)
  out <- NULL
  for (j in seq_len(length(files))) {
    fn <- file.path(rdir, files[j])
    out <- rbind(out, res_to_df(fn))
  }
  tibble::as_tibble(out)
}

# Scalar results
df_s <- NULL
for (d in ls) {
  df_s <- rbind(df_s, get_scalar_res(d))
}
