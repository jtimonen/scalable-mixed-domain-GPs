library(lgpr2)
library(tidyverse)

parent_res_dir <- "res-tri"
ls <- dir(parent_res_dir)
ls <- ls[grepl(pattern = "res-", ls)]

# Get scalar results
get_scalar_res <- function(dir) {
  rdir <- file.path(parent_res_dir, dir)
  files <- dir(rdir)
  out <- NULL
  for (j in seq_len(length(files))) {
    fn <- file.path(rdir, files[j])
    out <- rbind(out, scalar_res_to_df(fn))
  }
  tibble::as_tibble(out)
}


# Get path results
get_path_res <- function(dir) {
  rdir <- file.path(parent_res_dir, dir)
  files <- dir(rdir)
  out <- NULL
  for (j in seq_len(length(files))) {
    fn <- file.path(rdir, files[j])
    out <- rbind(out, path_res_to_df(fn))
  }
  tibble::as_tibble(out)
}

# Scalar results to data frame
scalar_res_to_df <- function(fn) {
  r <- readRDS(fn)
  x <- unlist(r[c("N_indiv", "snr")])
  x_diag <- r[["diag"]]
  num_terms <- length(r$term_names)
  x <- c(x, x_diag, num_terms)
  names(x)[length(x)] <- "num_terms"
  x <- data.frame(t(c(fn, x)))
  colnames(x)[1] <- "file"
  x
}

get_elpd_df <- function(search) {
  df <- search$history[c("elpd_loo", "elpd_loo_rel_diff")]
  df$num_sub_terms <- 0:(nrow(df) - 1)
  df
}

get_both_elpd_dfs <- function(r) {
  df1 <- get_elpd_df(r$search_pp_fs)
  df2 <- get_elpd_df(r$search_pp_dir)
  df1$method <- "forward search"
  df2$method <- "predefined path"
  tibble::as_tibble(rbind(df1, df2))
}

# Scalar results to data frame
path_res_to_df <- function(fn) {
  r <- readRDS(fn)
  a <- get_both_elpd_dfs(r)
  a$file <- fn
  a
}


# Scalar results
df_s <- NULL
df_p <- NULL
for (d in ls) {
  df_s <- rbind(df_s, get_scalar_res(d))
  df_p <- rbind(df_p, get_path_res(d))
}

df <- df_p %>% left_join(df_s, by = "file")
df$experim <- paste(df$file, df$method)
df$setup <- paste(df$num_terms, "_", df$N_indiv, "_", df$snr)
plt <- ggplot(
  df,
  aes(x = num_sub_terms, y = elpd_loo_rel_diff, color = method, group = experim)
) +
  geom_line() +
  geom_point() +
  facet_wrap(. ~ setup) +
  geom_hline(yintercept = c(-1, 1), lty = 2)
