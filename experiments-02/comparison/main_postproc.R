library(lgpr2)
library(tidyverse)

parent_res_dir <- "res-100"
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

# Term name to term description
term_name_to_desc <- function(term_name) {
  has_zu <- grepl(term_name, pattern = "z_u")
  has_xu <- grepl(term_name, pattern = "x_u")
  if (has_zu) {
    desc <- "irrelevant categorical"
  } else if (has_xu) {
    desc <- "irrelevant continuous"
  } else {
    desc <- term_name
  }
  desc
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

# Full data frame of search paths
df <- df_p %>% left_join(df_s, by = "file")
df$experim <- paste(df$file, df$method)
df$setup <- paste0("N_terms = ", df$num_terms, ", SNR = ", df$snr)
df$method_in_setup <- paste0(df$method, "-", df$setup)
df$term_desc <- Vectorize(term_name_to_desc)(df$term_char)

# Plot each elp curve
plt <- ggplot(
  df,
  aes(
    x = num_sub_terms, y = elpd_loo_rel_diff, color = method,
    fill = method, group = file
  )
) +
  geom_line() +
  geom_point() +
  facet_wrap(. ~ method_in_setup) +
  geom_hline(yintercept = c(-1, 1), lty = 2)


# Summary data frame of elpd
df_sum <- df %>%
  group_by(method, setup, num_sub_terms) %>%
  summarize(
    med = quantile(elpd_loo_rel_diff, na.rm = TRUE, probs = 0.5),
    mean = mean(elpd_loo_rel_diff, na.rm = TRUE),
    lower = quantile(elpd_loo_rel_diff, na.rm = TRUE, probs = 0.05),
    upper = quantile(elpd_loo_rel_diff, na.rm = TRUE, probs = 0.95),
    .groups = "drop"
  )

# ELPD plot
plt_elp <- ggplot(df_sum, aes(x = num_sub_terms, y = mean, color = method)) +
  facet_wrap(. ~ setup) +
  geom_vline(xintercept = 4, color = "orange", lty = 2) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = c(-1, 1), lty = 2) +
  theme_bw() +
  ylab("LOO-ELPD rel. diff.") +
  xlab("Number of terms in model") +
  scale_x_continuous(breaks = unique(df_sum$num_sub_terms))


# Create a complete grid of all combinations
complete_grid <- expand.grid(
  setup = unique(df$setup),
  method = unique(df$method),
  term_desc = setdiff(unique(df$term_desc), c("f_baseline_id", NA))
)

# Compute total count for each setup-method combination
total_counts <- df %>%
  filter(num_sub_terms == 2) %>%
  group_by(setup, method) %>%
  summarize(total = n(), .groups = "drop")

# Most common first selection (excluding f_baseline_id)
df_freq <- df %>%
  filter(num_sub_terms == 2) %>%
  group_by(setup, method, term_desc) %>%
  summarize(count = n(), .groups = "drop") %>%
  left_join(total_counts, by = c("setup", "method")) %>%
  mutate(frequency = count / total) %>%
  right_join(complete_grid, by = c("setup", "method", "term_desc")) %>%
  replace_na(list(count = 0, total = 0, frequency = 0))

# Distribution of first selected term
plt_first <- df_freq %>% ggplot(aes(x = term_desc, y = frequency, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(. ~ setup) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
  ggtitle("First selected term") +
  xlab("Term") +
  ylab("Proportion of times selected first")


# Proportion of correct
correct <- c("f_baseline_id", "f_gp_age", "f_gp_ageXz", "f_gp_x")


df_cor <- NULL
df0 <- df %>% filter(num_sub_terms > 0)
for (j in 1:6) {
  # Compute total count for each setup-method combination
  tc_j <- df0 %>%
    filter(num_sub_terms <= j) %>%
    group_by(setup, method) %>%
    summarize(total = n(), .groups = "drop")

  df_j <- df0 %>%
    filter(num_sub_terms <= j) %>%
    filter(term_char %in% correct) %>%
    group_by(setup, method) %>%
    summarize(count = n(), .groups = "drop") %>%
    left_join(tc_j, by = c("setup", "method")) %>%
    mutate(percentage = count / total)

  df_j$num_terms <- j

  df_cor <- rbind(df_cor, df_j)
}
plt_cor <- ggplot(df_cor, aes(x = num_terms, y = percentage, color = method)) +
  facet_wrap(. ~ setup) +
  geom_vline(xintercept = 4, color = "orange", lty = 2) +
  geom_line() +
  geom_point() +
  xlab("Number of terms in model") +
  ylab("% of correct terms") +
  theme_bw() +
  scale_x_continuous(breaks = unique(df_cor$num_terms))

# Plt setup
refine_plot <- function(plt, pal = 6) {
  plt + facet_grid(. ~ setup) +
    scale_color_brewer(type = "qual", palette = pal) +
    scale_fill_brewer(type = "qual", palette = pal)
}

# Combined result plot
library(ggpubr)
plt_a <- refine_plot(plt_elp)
plt_b <- refine_plot(plt_cor)
plt_c <- refine_plot(plt_first + ggtitle(""))
plt_res <- ggarrange(plt_a, plt_b, plt_c,
  nrow = 3, ncol = 1, labels = "auto",
  heights = c(1, 1, 1.7)
)

# Save
ggsave(plt_res, filename = "resplot.pdf", width = 9.5, height = 8.2)
