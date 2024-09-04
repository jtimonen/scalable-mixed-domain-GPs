library(lgpr2)
library(tidyverse)
library(ggpubr)
source("postproc.R")

# File paths
parent_res_dir <- "res44" # "res22"
ls <- dir(parent_res_dir)
ls <- ls[grepl(pattern = "res-", ls)]

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
if (length(unique(df$num_terms)) > 1) {
  df$setup <- paste0("N_terms = ", df$num_terms, ", SNR = ", df$snr)
} else {
  df$setup <- paste0("SNR = ", df$snr)
}
df$method_in_setup <- paste0(df$method, "-", df$setup)
df$term_desc <- Vectorize(term_name_to_desc)(df$term_char)
df$method_x_file <- paste0(df$method, "-", df$file)

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
  geom_vline(xintercept = 6, color = "orange", lty = 2) +
  geom_hline(yintercept = c(-1, 1), lty = 2) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylab("LOO-ELPD rel. diff.") +
  xlab("Number of terms in model") +
  scale_x_continuous(breaks = unique(df_sum$num_sub_terms))

# Plot all experiments
plt_elp_better <- ggplot(
  df,
  aes(x = num_sub_terms, y = elpd_loo_rel_diff, group = method_x_file)
) +
  geom_vline(xintercept = 6, color = "orange", lty = 1) +
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


# Create a complete grid of all combinations
complete_grid <- expand.grid(
  setup = unique(df$setup),
  method = unique(df$method),
  term_desc = setdiff(unique(df$term_desc), c("f_baseline_id", "f_gp_age", NA))
)

# Compute total count for each setup-method combination
total_counts <- df %>%
  filter(num_sub_terms == 2) %>%
  group_by(setup, method) %>%
  summarize(total = n(), .groups = "drop")

# Most common first selection (excluding f_baseline_id and age)
selection_rate_at_step <- function(df_res, k) {
  df_res %>%
    filter(num_sub_terms == k) %>%
    group_by(setup, method, term_desc) %>%
    summarize(count = n(), .groups = "drop") %>%
    left_join(total_counts, by = c("setup", "method")) %>%
    mutate(frequency = count / total) %>%
    right_join(complete_grid, by = c("setup", "method", "term_desc")) %>%
    replace_na(list(count = 0, total = 0, frequency = 0))
}

df_freq3 <- selection_rate_at_step(df, 3)
df_freq4 <- selection_rate_at_step(df, 4)

# Distribution of first selected term
plt_sel3 <- plot_sel_dist(df_freq3, 3)
plt_sel4 <- plot_sel_dist(df_freq4, 4)

# Proportion of correct
correct <- c(
  "f_baseline_id", "f_gp_age", "f_gp_ageXz", "f_gp_x", "f_gp_w",
  "f_gp_ageXr"
)
df_cor <- NULL
df0 <- df %>% filter(num_sub_terms > 0)
for (j in 1:8) {
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

# Correctness of selection
plt_cor <- ggplot(df_cor, aes(x = num_terms, y = percentage, color = method)) +
  facet_wrap(. ~ setup) +
  geom_vline(xintercept = 6, color = "orange", lty = 2) +
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


# Summary data frame of cumulative relevance
df_sum_cr <- df %>%
  filter(method == "predefined path") %>%
  group_by(method, setup, num_sub_terms) %>%
  summarize(
    med = quantile(cum_relevance, na.rm = TRUE, probs = 0.5),
    mean = mean(cum_relevance, na.rm = TRUE),
    lower = quantile(cum_relevance, na.rm = TRUE, probs = 0.05),
    upper = quantile(cum_relevance, na.rm = TRUE, probs = 0.95),
    .groups = "drop"
  )

# Cumulative relevance plot
plt_cr <- ggplot(
  df_sum_cr,
  aes(x = num_sub_terms, y = med, ymin = lower, ymax = upper)
) +
  facet_wrap(. ~ setup) +
  geom_vline(xintercept = 6, color = "orange", lty = 2) +
  geom_hline(yintercept = c(0.95), lty = 2, color = "firebrick") +
  geom_errorbar(width = 0.4) +
  ylim(0.5, 1) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylab("Cumul. relevance") +
  xlab("Number of terms in model") +
  scale_x_continuous(breaks = unique(df_sum_cr$num_sub_terms))



# Combined result plot
plt_a <- refine_plot(plt_elp) + theme(legend.position = "top")
plt_b <- refine_plot(plt_cr)
plt_c <- refine_plot(plt_cor) + theme(legend.position = "none")
plt_d <- refine_plot(plt_sel3) + theme(legend.position = "none")
plt_e <- refine_plot(plt_sel4) + theme(legend.position = "none")
plt_res1 <- ggarrange(plt_a, plt_b, plt_c,
  nrow = 3, ncol = 1, labels = "auto",
  legend.grob = get_legend(plt_a)
)
plt_res2 <- ggarrange(plt_d, plt_e,
  nrow = 2, ncol = 1, labels = "auto",
  legend.grob = get_legend(plt_d)
)


# Diagnose
df_diagnose <- df %>%
  group_by(setup) %>%
  summarize(
    mean_mcmc_time = mean(time_mcmc),
    mean_forward_search_time = mean(time_fs),
    mean_predefined_path_time = mean(time_dir),
    mean_num_divergent = mean(num_divergent1 + num_divergent2),
    mean_num_max_treedepth = mean(num_max_treedepth1 + num_max_treedepth2),
    mean_max_rhat = mean(max_rhat),
    sd_max_rhat = sd(max_rhat),
    max_max_rhat = max(max_rhat)
  )

# Save
ggsave(plt_res1, filename = "fig_pred2.pdf", width = 8.5, height = 4.8)
ggsave(plt_res2, filename = "fig_sel2.pdf", width = 8.5, height = 6.8)

# Supplementary figures
plt_s1 <- plt_elp_better + scale_color_brewer(type = "qual", palette = 6) +
  scale_fill_brewer(type = "qual", palette = 6) +
  theme(legend.position = "top", legend.title = element_blank())
ggsave(plt_s1, filename = "fig_supp_pred2.pdf", width = 8, height = 5.2)

# Supplementary table
st <- xtable::xtable(t(df_diagnose))
