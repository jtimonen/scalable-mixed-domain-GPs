# Like  lgpr:::posterior_f but with approximate model fit
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
pred_approx <- function(model, fit, df_star, num_bf, scale_bf,
                        refresh = NULL, c_hat_pred = NULL) {
  fp <- posterior_f_approx(model, fit, df_star, num_bf, scale_bf, refresh)
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
  names(f_comp) <- component_names(model)
  c_hat_pred <- lgpr:::set_c_hat_pred(model, f_sum, c_hat_pred, TRUE)
  h <- lgpr:::map_f_to_h(model, f_sum, c_hat_pred, reduce = NULL)

  # Return
  new("Prediction",
    f_comp = f_comp,
    f = f_sum,
    h = h,
    x = df_star,
    extrapolated = FALSE
  )
}
