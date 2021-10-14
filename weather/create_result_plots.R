# Startup
source(normalizePath(file.path("..", "common.R")))
outdir <- startup(create_dir = FALSE)
NUM_BF <- c(8, 12, 16, 24, 32)
TIMES <- rep(0, 5)
R_HATS <- list()
F_SUM <- list()
for (b in 1:length(NUM_BF)) {
  num_bf <- NUM_BF[b]
  fn <- file.path("results_temperature", paste0("res_", num_bf, ".rds"))
  dat <- load_weather_data()
  dat$day <- as.Date(dat$day, origin = "1964-01-01")
  results <- readRDS(file = fn)
  pa <- results$pred
  fit <- get_cmdstanfit(results$fit)
  total_time <- fit$time()$total
  TIMES[b] <- total_time
  # Create new fit and get summary
  of <- fit$output_files()
  get_fn <- function(x) {
    parts <- strsplit(x[1], "/")[[1]]
    file.path("results_temperature", parts[length(parts)])
  }
  of <- c(get_fn(of[1]), get_fn(of[2]), get_fn(of[3]), get_fn(of[4]))
  fit_new <- cmdstanr::as_cmdstan_fit(files = of)
  smr <- fit_new$summary()
  R_HATS[[b]] <- smr$rhat
  F_SUM[[b]] <- colMeans(pa@h)

  # Create y-axis labelss
  ylab1 <- expression(paste("", "f", phantom()^{
    paste("", "(", "1", ")", "")
  }, "(", "", "day", ")", "", ""))
  ylab2 <- expression(paste("", "f", phantom()^{
    paste("", "(", "2", ")", "")
  }, "(", "", "day", ",", "region", ")", "", ""))
  ylab3 <- expression(paste("", "f", phantom()^{
    paste("", "(", "3", ")", "")
  }, "(", "", "day", ",", "station", ")", "", ""))


  pf1 <- plot_f(dat, pa, 1, aes1()) + geom_ribbon(alpha = 0.3) + ylab(ylab1)
  pf2 <- plot_f(dat, pa, 2, aes2()) + geom_ribbon(alpha = 0.3) + ylab(ylab2)
  pf3 <- plot_f(dat, pa, 3, aes3()) + facet_wrap(. ~ region) + ylab(ylab3) +
    theme(legend.text = element_text(size = 8), legend.title = element_blank())
  ph <- plot_f(dat, pa, 0, aes3()) +
    geom_point(
      data = dat,
      inherit.aes = FALSE,
      aes(x = day, y = temperature, group = station),
      pch = ".", alpha = 0.7,
    ) + facet_wrap(. ~ station) + ylab("Temperature (C)")

  plt_12 <- ggarrange(pf1, pf2,
    nrow = 1, ncol = 2,
    labels = c("a", "b"), widths = c(0.8, 1)
  )
  full_plt <- ggarrange(plt_12, pf3,
    nrow = 2, ncol = 1, heights = c(0.6, 1),
    labels = c("", "c")
  )
  fn_out_comp <- file.path("figures", paste0("components_", num_bf, ".pdf"))
  fn_out_sum <- file.path("figures", paste0("sum_", num_bf, ".pdf"))
  ggsave(full_plt, file = fn_out_comp, width = 6.9, height = 6.7)
  ggsave(ph, file = fn_out_sum, width = 12, height = 16)
}
