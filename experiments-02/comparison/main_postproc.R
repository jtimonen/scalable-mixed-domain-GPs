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

get_elpd_df <- function(search, term_names) {
  df <- search$history[c("elpd_loo", "elpd_loo_rel_diff")]
  path <- search$path
  path_char <- term_names[path]
  path <- c(0, path)
  path_char <- c(NA, path_char)
  df$num_sub_terms <- 0:(nrow(df) - 1)
  df$path <- path
  df$path_char <- path_char
  df
}

get_both_elpd_dfs <- function(r) {
  df1 <- get_elpd_df(r$search_pp_fs, r$term_names)
  df2 <- get_elpd_df(r$search_pp_dir, r$term_names)
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
df$setup <- paste0(df$num_terms, "_", df$N_indiv, "_", df$snr)
df$method_in_setup <- paste0(df$method, "-", df$setup)
plt <- ggplot(
  df,
  aes(
    x = num_sub_terms, y = elpd_loo_rel_diff, color = method,
    fill = method, group = method_in_setup
  )
) +
  geom_line() +
  geom_point() +
  facet_wrap(. ~ method_in_setup) +
  geom_hline(yintercept = c(-1, 1), lty = 2)


# Summary plot
df_sum <- df %>%
  group_by(method, setup, num_sub_terms) %>%
  summarize(
    med = quantile(elpd_loo_rel_diff, na.rm = TRUE, probs = 0.5),
    mean = mean(elpd_loo_rel_diff, na.rm = TRUE),
    lower = quantile(elpd_loo_rel_diff, na.rm = TRUE, probs = 0.05),
    upper = quantile(elpd_loo_rel_diff, na.rm = TRUE, probs = 0.95),
    .groups = "drop"
  )
plt_elp <- ggplot(df_sum, aes(x = num_sub_terms, y = mean, color = method)) +
  facet_wrap(. ~ setup) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = c(-1, 1), lty = 2) +
  theme_bw() +
  ylab("LOO-ELPD relative diff") +
  xlab("Number of terms") +
  scale_x_continuous(breaks = unique(df_sum$num_sub_terms))


# Mst common first selection
df_path_char_freq <- df %>%
  filter(num_sub_terms == 1) %>%
  group_by(setup, method) %>%
  summarize(frequency = n(), .groups = "drop")
most_common <- df_path_char_freq %>%
  arrange(setup, method, desc(frequency)) %>%
  group_by(setup, method) %>%
  slice(1)


# Proportion of correct
correct <- c("f_baseline_id", "f_gp_age", "f_gp_ageXz", "f_gp_x")

df_cor <- NULL
for (j in 1:6) {
  df_j <- df %>%
    filter(num_sub_terms <= j) %>%
    filter(path_char %in% correct) %>%
    group_by(setup, method) %>%
    summarize(frequency = n(), .groups = "drop")
  df_j$num_terms <- j
  df_j <- df_j %>% mutate(percentage = frequency / (50 * num_terms))
  df_cor <- rbind(df_cor, df_j)
}
plt_cor <- ggplot(df_cor, aes(x = num_terms, y = percentage, color = method)) +
  facet_wrap(. ~ setup) +
  geom_line() +
  geom_point() +
  xlab("Number of terms in model") +
  ylab("Percentage of correct terms") +
  theme_bw() +
  scale_x_continuous(breaks = unique(df_cor$num_terms))
