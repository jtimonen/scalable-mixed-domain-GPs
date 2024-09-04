# Full data frame of search paths
create_result_df <- function(parent_res_dir) {
  ls <- dir(parent_res_dir)
  ls <- ls[grepl(pattern = "res-", ls)]
  df_s <- NULL
  df_p <- NULL
  for (d in ls) {
    df_s <- rbind(df_s, get_scalar_res(d))
    df_p <- rbind(df_p, get_path_res(d))
  }
  df <- df_p %>% left_join(df_s, by = "file")
  df$experim <- paste(df$file, df$method)
  if (length(unique(df$num_terms)) > 1) {
    df$setup <- paste0("N_terms = ", df$num_terms, ", SNR = ", df$snr)
  } else {
    df$setup <- paste0("SNR = ", df$snr)
  }
  df$method_in_setup <- paste0(df$method, "-", df$setup)
  df$term_desc <- Vectorize(term_name_to_desc)(df$term_char)
  df$method_x_file <- paste0(df$method, "-", df$file)
  df
}

# Summary
result_summary_df <- function(df) {
  df %>%
    group_by(method, setup, num_sub_terms) %>%
    summarize(
      med = quantile(elpd_loo_rel_diff, na.rm = TRUE, probs = 0.5),
      mean = mean(elpd_loo_rel_diff, na.rm = TRUE),
      lower = quantile(elpd_loo_rel_diff, na.rm = TRUE, probs = 0.05),
      upper = quantile(elpd_loo_rel_diff, na.rm = TRUE, probs = 0.95),
      .groups = "drop"
    )
}

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
  times <- r$search_times
  times <- data.frame(t(c(r$mcmc_time, 60 * as.numeric(times))))
  colnames(times) <- c("time_mcmc", "time_fs", "time_dir")
  x <- unlist(r[c("N_indiv", "snr")])
  x_diag <- r[["diag"]]
  num_terms <- length(r$term_names)
  x <- c(x, x_diag, num_terms)
  names(x)[length(x)] <- "num_terms"
  x <- data.frame(t(x))
  x$file <- fn
  cbind(x, times)
}

# Get main result df for one result
get_elpd_df <- function(search, term_names) {
  df <- search$history[c("elpd_loo", "elpd_loo_rel_diff")]
  path <- search$path
  path_char <- term_names[path]
  path <- c(0, path)
  path_char <- c(NA, path_char)
  df$num_sub_terms <- 0:(nrow(df) - 1)
  df$term <- path
  df$term_char <- path_char
  df
}

# Both
get_both_elpd_dfs <- function(r) {
  df1 <- get_elpd_df(r$search_pp_fs, r$term_names)
  df2 <- get_elpd_df(r$search_pp_dir, r$term_names)
  df1$cum_relevance <- 0
  df2$cum_relevance <- r$r_path[seq_len(nrow(df2))]
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


# Distribution of term selection at a given step
plot_sel_dist <- function(df_freq, k) {
  df_freq %>% ggplot(aes(x = term_desc, y = frequency, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(. ~ setup) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
    ggtitle(paste0("Variable selected at step k = ", k)) +
    xlab("Variable") +
    ylab("Selection rate")
}


# Term name to term description
term_name_to_desc <- function(term_name) {
  has_zu <- grepl(term_name, pattern = "z_u")
  has_xu <- grepl(term_name, pattern = "x_u")
  has_x <- grepl(term_name, pattern = "_x")
  has_z <- grepl(term_name, pattern = "_ageXz")
  if (is.na(term_name)) {
    return(NA)
  }
  if (term_name == "f_gp_w") {
    desc <- "TRUE w"
  } else if (term_name == "f_gp_ageXr") {
    desc <- "TRUE r"
  } else if (has_zu) {
    desc <- "FALSE z"
  } else if (has_xu) {
    desc <- "FALSE x"
  } else if (has_x) {
    desc <- "TRUE x"
  } else if (has_z) {
    desc <- "TRUE z"
  } else {
    desc <- term_name
  }
  desc
}



# ELPD plot
plot_elp <- function(df_sum, orange_line_x) {
  ggplot(df_sum, aes(x = num_sub_terms, y = mean, color = method)) +
    facet_wrap(. ~ setup) +
    geom_vline(xintercept = orange_line_x, color = "orange", lty = 2) +
    geom_hline(yintercept = c(-1, 1), lty = 2) +
    geom_line() +
    geom_point() +
    theme_bw() +
    ylab("LOO-ELPD rel. diff.") +
    xlab("Number of terms in model") +
    scale_x_continuous(breaks = unique(df_sum$num_sub_terms))
}


# Plot all experiments
plot_elp_full <- function(df, df_sum, orange_line_x) {
  ggplot(
    df,
    aes(x = num_sub_terms, y = elpd_loo_rel_diff, group = method_x_file)
  ) +
    geom_vline(xintercept = orange_line_x, color = "orange", lty = 1) +
    geom_hline(yintercept = c(-1, 1), lty = 1) +
    geom_line(color = "gray") +
    geom_line(
      data = df_sum, aes(x = num_sub_terms, y = mean, color = method),
      inherit.aes = FALSE
    ) +
    theme_bw() +
    ylab("LOO-ELPD rel. diff.") +
    xlab("Number of terms in model") +
    scale_x_continuous(breaks = unique(df_sum$num_sub_terms)) +
    facet_wrap(. ~ method + setup, labeller = "label_both")
}
