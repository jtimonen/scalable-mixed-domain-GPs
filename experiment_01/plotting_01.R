
# Create data frame for ggplot
create_pred_plot_df <- function(preds) {
  J <- length(preds)
  df <- NULL
  model <- c()
  for (j in 1:J) {
    p <- preds[[j]]
    if (isa(p, "Prediction")) {
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
  if (isa(pred, "Prediction")) {
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
