# Compute log predictive density
compute_lpd <- function(fit, df_star) {
  stopifnot(is(fit, "lgpfit"))
  p <- lgpr::pred(fit, x = df_star, reduce = NULL)
  y_name <- fit@model@var_names$y
  y_star <- df_star[[y_name]]
  gaussian_lpd(p, y_star)
}

# Gaussian log predictive density
gaussian_lpd <- function(pred, y_star) {
  stopifnot(is(pred, "GaussianPrediction"))
  P <- lgpr::num_evalpoints(pred)
  S <- lgpr::num_paramsets(pred)
  y_means <- pred@y_mean
  y_stds <- pred@y_std
  log_pds <- array(0.0, dim = dim(y_means))
  for (s in seq_len(S)) {
    mu <- y_means[s, ]
    sig <- y_stds[s, ]
    log_pds[s, ] <- stats::dnorm(y_star, mean = mu, sd = sig, log = TRUE)
  }
  return(rowMeans(log_pds))
}

# posterior_f but with approximate model fit
posterior_f_approx <- function(model, fit, df_star, num_bf, scale_bf,
                               refresh = NULL) {
  xi <- posterior::merge_chains(fit$draws("xi"))
  alpha <- posterior::merge_chains(fit$draws("alpha"))
  ell <- posterior::merge_chains(fit$draws("ell"))
  S <- dim(alpha)[1]
  P <- nrow(df_star)
  J <- length(component_names(model))
  pred_model <- lgpr::create_model(
    formula = formula(model@model_formula@call),
    data = df_star, prior = model@full_prior
  )
  decs <- categorical_kernel_decompositions(pred_model)
  si_add <- additional_stan_input(
    pred_model, num_bf, scale_bf,
    decs$decompositions
  )
  si <- c(pred_model@stan_input, si_add)
  F_PRED <- array(0.0, dim = c(S, J, P))
  tdata <- do_transformed_data(si)
  as <- alpha[, 1, , drop = T]
  es <- ell[, 1, , drop = T]
  xis <- xi[, 1, , drop = T]
  as <- matrix_to_list(as)
  es <- matrix_to_list(es)
  xis <- matrix_to_list(xis)
  if (is.null(refresh)) {
    refresh <- round(S / 10)
  }
  build_f_draws(si, tdata, as, es, xis, refresh)
}

# Predict with approximate model
pred_approx <- function(model, fit, df_star, num_bf, scale_bf,
                        refresh = NULL) {
  fp <- posterior_f_approx(model, fit, df_star, num_bf, scale_bf, refresh)
  return(fp)
}
