y_star <- test_dat[["y"]]
em <- compute_metrics(fits, preds, y_star)
rtables <- format_results(em, SCALES, NBFS)

# Plot denser predictions
x_dense <- create_x_dense(train_dat, test_dat)
preds_dense <- compute_predictions(fits, x_dense)

# Plots
plt1 <- plot_against_exact(1:4, train_dat, test_dat, preds_dense)
plt2 <- plot_against_exact(5:8, train_dat, test_dat, preds_dense)
plt3 <- plot_against_exact(9:12, train_dat, test_dat, preds_dense)

# Save plot
j <- 0
for (pp in pred_plots) {
  j <- j + 1
  nam <- names(pred_plots)[j]
  np <- gsub("=", "_", gsub("[.]", "-", gsub(" ", "_", nam)))
  fn <- file.path(outdir, paste0("preds_", N, "_", np, ".pdf"))
  cat("Saving", fn, "\n")
  ggsave(fn, plot = pp, width = 5.28, height = 4.75)
}

# Runtimes plot
# rt <- plot_runtimes_wrt_N(PRES, NUM_BF, N_sizes, scale_bf)
# ggsave("res/exp3/times3.pdf", plot = rt, width = 5.5, height = 4.3)

# Save results in Rdata
res_to_save <- results[c("runtimes", "num_div")]
res_to_save$metrics <- em

# Create latex table
library(xtable)
c1 <- paste0("exact = ", rtables$elpds_w1$exact)
c2 <- paste0("exact = ", rtables$elpds_w2$exact)
c3 <- paste0("exact = ", rtables$rmses$exact)
xtable(rtables$elpds_w1$approx, digits = 4, caption = c1)
xtable(rtables$elpds_w2$approx, digits = 4, caption = c2)
xtable(rtables$rmses$approx, digits = 4, caption = c3)
