library(ggplot2)
library(ggpubr)

read_res <- function(M) {
  fn <- paste0("res_rds/times_", M, ".rds")
  a <- readRDS(fn)
  n <- as.numeric(a$n)
  t1 <- as.numeric(a$t1)
  t2 <- as.numeric(a$t2)
  t3 <- as.numeric(a$t3)
  t4 <- as.numeric(a$t4)

  # Create data frame
  n <- rep(n, 4)
  time <- c(t1, t2, t3, t4)
  df <- data.frame(n, time)
}

create_plot <- function(df, M) {
  plt <- ggplot(df, aes(x = n, y = time)) +
    geom_point(pch = 4) +
    stat_smooth(method = "lm", col = "gray30") +
    theme_minimal() +
    ylab("Chain runtime (s)") +
    xlab("n") +
    ggtitle(paste0("Number of basis functions: ", M)) +
    ylim(0, 1000) +
    xlim(0, 1880)
  return(plt)
}

r1 <- read_res(25)
r2 <- read_res(30)
r3 <- read_res(35)
plt1 <- create_plot(r1, 25)
plt2 <- create_plot(r2, 30)
plt3 <- create_plot(r3, 35)

a <- list(plt1, plt2, plt3)
plt <- ggarrange(plt1, plt2, plt3, nrow = 1, ncol = 3)
