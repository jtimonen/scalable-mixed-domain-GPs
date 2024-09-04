library(lgpr2)
library(tidyverse)
library(ggpubr)
source("postproc.R")


# File paths
parent_res_dir <- "res33" # "res22"

# Result dfs and plots
df <- create_result_df(parent_res_dir)
df_sum <- result_summary_df(df)
plt_elp <- plot_elp(df_sum, 4)

# Distribution of selected terms at different steps
df_freq3 <- selection_rate_at_step(df, 3)
df_freq4 <- selection_rate_at_step(df, 4)
plt_sel3 <- plot_sel_dist(df_freq3, 3)
plt_sel4 <- plot_sel_dist(df_freq4, 4)

# Proportion of correct
correct <- c("f_baseline_id", "f_gp_age", "f_gp_ageXz", "f_gp_x")
df_cor <- create_df_cor(df, 6, correct)
plt_cor <- plot_p_correct(df_cor, 4)

# Cumulative relevance
df_sum_cr <- cumrel_summary_df(df)
plt_cr <- plot_cum_rel(df_sum_cr, 4)

# Study wrong selections with snr=0.1
find_files_with_wrong_z <- function(df) {
  df24 <- df %>% filter(
    method == "forward search",
    num_sub_terms <= 4,
    num_sub_terms >= 2,
    snr == "0.1"
  )
  df24$afterX <- sapply(
    strsplit(df24$term_char, split = "X"),
    function(x) x[length(x)]
  )
  df24$has_zu <- grepl(df24$afterX, pattern = "z_u")
  df24 %>%
    filter(has_zu) %>%
    select(c("afterX", "file"))
}




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
df_diagnose <- diagnosis(df)

# Save
ggsave(plt_res1, filename = "fig_pred.pdf", width = 7.5, height = 6.5)
ggsave(plt_res2, filename = "fig_sel.pdf", width = 7.5, height = 5)

# Supplementary table
st <- xtable::xtable(t(df_diagnose), digits = 3, caption = "Diagnostics 1.")
