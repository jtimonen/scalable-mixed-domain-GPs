library(lgpr2)
library(tidyverse)
library(ggpubr)
source("postproc.R")

# File paths
parent_res_dir <- "res44" # "res22"

# Result dfs and plots
df <- create_result_df(parent_res_dir)
df_sum <- result_summary_df(df)
plt_elp <- plot_elp(df_sum, 6)
plt_elp_full <- plot_elp_full(df, df_sum, 6)

# Distribution of selected terms at different steps
df_freq3 <- selection_rate_at_step(df, 3)
df_freq4 <- selection_rate_at_step(df, 4)
df_freq5 <- selection_rate_at_step(df, 5)
df_freq6 <- selection_rate_at_step(df, 6)
plt_sel3 <- plot_sel_dist(df_freq3, 3)
plt_sel4 <- plot_sel_dist(df_freq4, 4)
plt_sel5 <- plot_sel_dist(df_freq3, 5)
plt_sel6 <- plot_sel_dist(df_freq4, 6)


# Proportion of correct
correct <- c(
  "f_baseline_id", "f_gp_age", "f_gp_ageXz", "f_gp_x", "f_gp_w",
  "f_gp_ageXr"
)
df_cor <- create_df_cor(df, 8, correct)
plt_cor <- plot_p_correct(df_cor, 6)


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

# Supplementary table
st <- xtable::xtable(t(df_diagnose))
