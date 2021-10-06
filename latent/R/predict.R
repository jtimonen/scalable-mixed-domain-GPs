# Create Stan input for pred
create_pred_stan_input <- function(amodel, x_star) {
  emodel <- amodel@exact_model
  xs_parsed <- lgpr:::kernelcomp.input_x(emodel, x_star)

  # Replace parts of  original Stan input according to test points
  si <- emodel@stan_input
  to_replace <- c(
    "num_obs", "x_cont", "x_cont_unnorm", "x_cont_mask",
    "x_cat", "idx_expand"
  )
  names(xs_parsed) <- to_replace
  si[to_replace] <- xs_parsed
  si_add <- amodel@added_stan_input
  si <- c(si, si_add)
}

# Like  lgpr:::posterior_f but with approximate model fit
posterior_f_approx <- function(fit, x_star, refresh = NULL) {
  amodel <- fit@model
  emodel <- amodel@exact_model
  fit <- fit@fit[[1]]
  xi <- posterior::merge_chains(fit$draws("xi"))
  alpha <- posterior::merge_chains(fit$draws("alpha"))
  ell <- posterior::merge_chains(fit$draws("ell"))
  S <- dim(alpha)[1]
  P <- nrow(x_star)
  si <- create_pred_stan_input(amodel, x_star)
  J <- length(component_names(emodel))
  F_PRED <- array(0.0, dim = c(S, J, P))
  tdata <- do_transformed_data(si)
  as <- alpha[, 1, , drop = T]
  es <- ell[, 1, , drop = T]
  xis <- xi[, 1, , drop = T]
  as <- matrix_to_list(as)
  es <- matrix_to_list(es)
  xis <- matrix_to_list(xis)
  if (is.null(refresh)) {
    refresh <- round(S / 4)
  }
  build_f_draws(si, tdata, as, es, xis, refresh)
}

# Extract draws of one component
get_approx_component <- function(fp, j) {
  S <- length(fp)
  P <- length(fp[[1]][[1]])
  fc <- matrix(0.0, S, P)
  for (s in seq_len(S)) {
    fc[s, ] <- fp[[s]][[j]]
  }
  return(fc)
}

# Predict with approximate model
pred_approx <- function(fit, x_star, refresh = NULL, c_hat_pred = NULL) {
  stopifnot(isa(fit, "ApproxModelFit"))
  om <- fit@model@obs_model
  fp <- posterior_f_approx(fit, x_star, refresh)
  emodel <- fit@model@exact_model
  S <- length(fp)
  J <- length(fp[[1]])
  P <- length(fp[[1]][[1]])
  f_comp <- list()
  f_sum <- matrix(0.0, S, P)
  for (j in seq_len(J)) {
    fj <- get_approx_component(fp, j)
    f_comp[[j]] <- fj
    f_sum <- f_sum + fj
  }
  names(f_comp) <- component_names(emodel)
  if (om == "gaussian") {
    yscl <- emodel@var_scalings$y
    h <- f_sum
    for (r in seq_len(nrow(h))) {
      h[r, ] <- lgpr:::apply_scaling(yscl, h[r, ], inverse = TRUE)
    }
  } else {
    c_hat_pred <- lgpr:::set_c_hat_pred(emodel, f_sum, c_hat_pred, TRUE)
    h <- lgpr:::map_f_to_h(emodel, f_sum, c_hat_pred, reduce = NULL)
  }

  # Return
  new("Prediction",
    f_comp = f_comp,
    f = f_sum,
    h = h,
    x = x_star,
    extrapolated = FALSE
  )
}


# Compute predictions using several fitted models
compute_predictions <- function(fits, x_star) {
  preds <- list()
  j <- 0
  for (f in fits) {
    j <- j + 1
    if (isa(f, "lgpfit")) {
      p <- lgpr::pred(f, x = x_star, reduce = NULL)
    } else {
      p <- pred_approx(f, x_star)
    }
    preds[[j]] <- p
  }
  names(preds) <- names(fits)
  return(preds)
}


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
  labels <- as.factor(paste(paste("id =", pdat$id), paste("z =", pdat$z),
    sep = ", "
  ))
  levels <- levels(labels)
  pdat <- modify_id_label(pdat, levels, labels)
  train_dat <- modify_id_label(train_dat, levels, labels)
  test_dat <- modify_id_label(test_dat, levels, labels)
  plt <- ggplot2::ggplot(pdat, aes(
    x = age, y = h, group = model,
    color = model
  )) +
    geom_line() +
    facet_wrap(. ~ id)
  plt <- plt + geom_point(
    data = train_dat, aes(x = age, y = y), inherit.aes = FALSE,
    col = "gray20", pch = 20
  )
  plt <- plt + geom_point(
    data = test_dat, aes(x = age, y = y), inherit.aes = FALSE,
    col = "red", pch = 4
  )
  return(plt)
}
