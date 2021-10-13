
# Create data frame for ggplot
create_pred_plot_df <- function(preds) {
  J <- length(preds)
  df <- NULL
  model <- c()
  for (j in 1:J) {
    p <- preds[[j]]
    if (is(p, "Prediction")) {
      h <- colMeans(p@h)
    } else {
      h <- colMeans(p@y_mean)
    }
    df_add <- cbind(p@x, h)
    df <- rbind(df, df_add)
    model <- c(model, rep(names(preds)[j], nrow(df_add)))
  }
  df$model <- as.factor(model)
  df$id_x_model <- paste(df$id, df$model)
  return(df)
}

# Create data frame for ggplot
create_f_plot_df <- function(pred) {
  df <- NULL
  if (is(pred, "Prediction")) {
    h <- pred@h
  } else {
    stop("pred must be a Prediction object")
  }
  f <- as.vector(h)
  S <- num_paramsets(pred)
  P <- num_evalpoints(pred)
  draw_idx <- as.factor(rep(1:P, each = S))
  xxx <- pred@x[rep(1:P, each = S), ]
  id <- xxx$id
  age <- xxx$age
  df <- data.frame(id, age, f, draw_idx)
  return(df)
}

# Modify id factor for labeling
modify_id_label <- function(df, levels, labels) {
  df$id <- factor(df$id,
    levels = levels(df$id),
    labels = levels(labels)
  )
  return(df)
}

# Function for plotting predictions
plot_preds <- function(train_dat, test_dat, preds) {
  pdat <- create_pred_plot_df(preds)
  # labels <- as.factor(paste(paste("id =", formatC(pdat$id, width = 2)),
  #  paste("z =", pdat$z),
  #  sep = ", "
  # ))
  # levels <- levels(labels)
  # pdat <- modify_id_label(pdat, levels, labels)
  pdat_exact <- pdat[which(pdat$model == "exact"), ]
  pdat_approx <- pdat[which(pdat$model != "exact"), ]
  # train_dat <- modify_id_label(train_dat, levels, labels)
  # test_dat <- modify_id_label(test_dat, levels, labels)
  num_fits <- length(levels(pdat_approx$model))
  my_colors <- RColorBrewer::brewer.pal(num_fits, "PuBu")[2:num_fits]
  plt <- ggplot2::ggplot(pdat_approx, aes(
    x = age, y = h, group = model,
    color = model
  )) +
    geom_line(lwd = 0.8) +
    facet_wrap(. ~ id) +
    scale_color_manual(values = my_colors)
  plt <- plt + geom_line(
    data = pdat_exact, aes(x = age, y = h),
    color = "orange", lty = 1, inherit.aes = FALSE,
    lwd = 0.3
  )
  plt <- plt + geom_point(
    data = train_dat, aes(x = age, y = y), inherit.aes = FALSE,
    col = "black", pch = 20
  )
  plt <- plt + geom_point(
    data = test_dat, aes(x = age, y = y), inherit.aes = FALSE,
    col = "black", pch = 4
  )
  return(plt)
}

# Plot f
plot_f_draws <- function(train_dat, test_dat, pred) {
  pdat <- create_f_plot_df(pred)
  plt <- ggplot2::ggplot(pdat, aes(x = age, y = f, group = draw_idx)) +
    geom_line(lwd = 0.8, alpha = 0.1) +
    facet_wrap(. ~ id)
  plt <- plt + geom_point(
    data = train_dat, aes(x = age, y = y), inherit.aes = FALSE,
    col = "black", pch = 20
  )
  plt <- plt + geom_point(
    data = test_dat, aes(x = age, y = y), inherit.aes = FALSE,
    col = "black", pch = 4
  )
  return(plt)
}

# Plot some approximate fits against exact
plot_against_exact <- function(inds, train_dat, test_dat, preds) {
  pd <- preds[c(inds, length(preds))]
  plot_preds(train_dat, test_dat, pd) + theme_bw() +
    theme(legend.position = "top") + ylab("")
}


# Plot results
plot_results_cb <- function(res, var_idx, scales, bfs) {
  ra <- res$approx[var_idx, , , ]
  ra_mean <- apply(ra, c(2, 3), mean)
  ra_std <- apply(ra, c(2, 3), stats::sd)
  re <- res$exact[var_idx, ]
  vname <- rownames(res$exact)[var_idx]
  re_mean <- mean(re)
  re_std <- mean(re)
  df_val_mean <- as.vector(ra_mean)
  df_val_std <- as.vector(ra_std)
  df_c <- rep(scales, times = length(bfs))
  df_b <- rep(bfs, each = length(scales))
  df <- data.frame(as.numeric(df_b), as.factor(df_c), df_val_mean, df_val_std)
  colnames(df) <- c("B", "c", "valm", "vals")
  aesth <- aes(
    x = B, y = valm, ymin = valm - vals, ymax = valm + vals,
    group = c, color = c, pch = c
  )
  plt <- ggplot(df, aesth) +
    geom_line() +
    geom_point() +
    ylab(vname) +
    scale_x_continuous(breaks = bfs) +
    geom_hline(yintercept = re_mean, lty = 2) +
    # geom_errorbar() +
    xlab("Number of basis functions (B)") +
    theme_bw() +
    labs(color = "Domain scale (c)", pch = "Domain scale (c)")
  return(plt)
}

# Plot MLPD results
plot_mlpd <- function(results, used_scales, used_nbfs, train = FALSE) {
  titles <- list(
    expression(paste("", "N", phantom()[{
      paste("train")
    }], " = 60", "")),
    expression(paste("", "N", phantom()[{
      paste("train")
    }], " = 90", "")),
    expression(paste("", "N", phantom()[{
      paste("train")
    }], " = 120", "")),
    expression(paste("", "N", phantom()[{
      paste("train")
    }], " = 150", "")),
    expression(paste("", "N", phantom()[{
      paste("train")
    }], " = 180", "")),
    expression(paste("", "N", phantom()[{
      paste("train")
    }], " = 210", ""))
  )
  if (train) {
    dname <- "train"
    var_idx <- 3
  } else {
    var_idx <- 4
    dname <- "test"
  }
  plots_mlpd <- list()
  var_idx <- 4
  for (j in 1:6) {
    plt_j <- plot_results_cb(results[[j]], var_idx, used_scales, used_nbfs)
    if (j %in% c(2, 3, 5, 6)) {
      plt_j <- plt_j + ylab(" ")
    } else {
      plt_j <- plt_j + ylab(paste0("MLPD (", dname, " data)"))
    }
    if (j != 5) {
      plt_j <- plt_j + xlab(" ")
    }
    plots_mlpd[[j]] <- plt_j + ylim(-5.25, -4.55) + ggtitle(titles[[j]])
  }

  # Return
  ggarrange(
    plotlist = plots_mlpd, nrow = 2, ncol = 3,
    common.legend = TRUE, label.y = 1.05
  )
}


# Plot runtime results
plot_runtimes_create_df <- function(res, scales, bfs, shown_c) {
  var_idx <- 1
  ra <- res$approx[var_idx, , , ]
  ra_mean <- apply(ra, c(2, 3), mean)
  ra_std <- apply(ra, c(2, 3), stats::sd)
  re <- res$exact[var_idx, ]
  vname <- rownames(res$exact)[var_idx]
  re_mean <- mean(re)
  re_std <- stats::sd(re)
  t_mean <- c(as.vector(ra_mean), re_mean)
  t_std <- c(as.vector(ra_std), re_std)
  df_c <- c(rep(scales, times = length(bfs)), shown_c)
  df_b <- c(rep(bfs, each = length(scales)), 0)
  ie <- as.factor(c(rep("Approximate", length(ra_mean)), "Exact"))
  df <- data.frame(as.factor(df_b), as.factor(df_c), t_mean, t_std, ie)
  colnames(df) <- c("B", "c", "t_mean", "t_std", "is_exact")
  inds_shown <- which(df$c == shown_c)
  return(df[inds_shown, ])
}


# Plot runtime results
plot_runtimes <- function(results, scales, bfs, N_TRAIN, shown_scale) {
  df <- NULL
  J <- length(N_TRAIN)
  for (j in 1:J) {
    df_j <- plot_runtimes_create_df(results[[j]], scales, bfs, shown_scale)
    n_train <- rep(N_TRAIN[j], nrow(df_j))
    df_j <- cbind(df_j, n_train)
    df <- rbind(df, df_j)
  }
  colnames(df)[ncol(df)] <- "n_train"
  plt <- ggplot(df, aes(
    x = n_train, y = t_mean, ymin = t_mean - t_std, ymax = t_mean + t_std,
    group = B, color = B, pch = B
  )) +
    scale_x_continuous(breaks = N_TRAIN) +
    geom_line() +
    geom_errorbar(alpha = 0.6) +
    geom_point() +
    ylab("Average runtime per chain (s)") +
    facet_wrap(. ~ is_exact, scales = "free_y") +
    xlab(expression(paste("", "N", phantom()[{
      paste("train")
    }], ""))) +
    theme_bw()
  return(plt)
}
