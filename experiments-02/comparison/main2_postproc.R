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
  "f_baseline_id", "f_gp_age", "f_gp_ageXz", "f_gp_x", "f_gp_w", "f_gp_ageXr"
)
df_cor <- create_df_cor(df, 8, correct)
plt_cor <- plot_p_correct(df_cor, 6)

# Cumulative relevance plot
df_sum_cr <- cumrel_summary_df(df)
plt_cr <- plot_cum_rel(df_sum_cr, 6)


# Combined result plot
plt_a <- refine_plot(plt_elp) + theme(legend.position = "top")
plt_b <- refine_plot(plt_cr)
plt_c <- refine_plot(plt_cor) + theme(legend.position = "none")
plt_d <- refine_plot(plt_sel3) + theme(legend.position = "none")
plt_e <- refine_plot(plt_sel4) + theme(legend.position = "none")
plt_f <- refine_plot(plt_sel5) + theme(legend.position = "none")
plt_g <- refine_plot(plt_sel6) + theme(legend.position = "none")
plt_res1 <- ggarrange(plt_a, plt_b, plt_c,
  nrow = 3, ncol = 1, labels = "auto",
  legend.grob = get_legend(plt_a)
)
plt_res2 <- ggarrange(plt_d, plt_e, plt_f, plt_g,
  nrow = 4, ncol = 1, labels = "auto",
  legend.grob = get_legend(plt_d)
)


# Diagnose
df_diagnose <- diagnosis(df)

# Save
ggsave(plt_res1, filename = "fig_pred2.pdf", width = 7.5, height = 6.5)
ggsave(plt_res2, filename = "fig_sel2.pdf", width = 7.5, height = 8)

# Supplementary table
st <- xtable::xtable(t(df_diagnose), digits = 3, caption = "Diagnostics 2.")
