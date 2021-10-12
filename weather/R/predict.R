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
posterior_f_approx <- function(fit, x_star) {
  amodel <- fit@model
  emodel <- amodel@exact_model
  fit <- get_cmdstanfit(fit)
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
  build_f_draws(si, tdata, as, es, xis)
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
pred_approx <- function(fit, x_star, c_hat_pred = NULL) {
  stopifnot(is(fit, "ApproxModelFit"))
  om <- fit@model@obs_model
  fp <- posterior_f_approx(fit, x_star)
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

# Draw exact predictions
pred_exact <- function(fit, x_star) {
  x <- fit@model@data
  alpha <- lgpr::get_draws(fit, pars = "alpha")
  ell <- lgpr::get_draws(fit, pars = "ell")
  sig <- lgpr::get_draws(fit, pars = "sigma")
  S <- nrow(alpha)
  x_scl <- fit@model@var_scalings$x_cont$age
  y_scl <- fit@model@var_scalings$y
  age <- lgpr:::apply_scaling(x_scl, x$age)
  ages <- lgpr:::apply_scaling(x_scl, x_star$age)
  z <- x$z
  zs <- x_star$z
  y_norm <- lgpr:::apply_scaling(y_scl, x$y)
  compute_kernels <- function(alpha, ell, x1, x2, K_zs) {
    K1 <- alpha[1]^2 * lgpr:::kernel_eq(x1, x2, ell = ell[1])
    K2_b <- lgpr:::kernel_eq(x1, x2, ell = ell[2])
    K2 <- alpha[2]^2 * K_zs * K2_b
    return(list(K1, K2))
  }
  draw_f <- function(K, Ks, Kss, sigma, y_norm) {
    N <- nrow(K)
    K_y <- K + sigma**2 * diag(N)
    mu <- Ks %*% solve(K_y, y_norm)
    Sigma <- Kss - Ks %*% solve(K_y, t(Ks))
    MASS::mvrnorm(n = 1, mu, Sigma)
  }
  N <- length(age)
  P <- length(ages)
  f1_draws <- matrix(0.0, S, P)
  f2_draws <- matrix(0.0, S, P)
  f_draws <- matrix(0.0, S, P)
  h_draws <- f_draws
  K_zs <- lgpr:::kernel_zerosum(z, z, 3)
  K_zs_s <- lgpr:::kernel_zerosum(zs, z, 3)
  K_zs_ss <- lgpr:::kernel_zerosum(zs, zs, 3)
  for (s in 1:S) {
    alpha_s <- alpha[s, ]
    ell_s <- ell[s, ]
    K <- compute_kernels(alpha_s, ell_s, age, age, K_zs)
    Ks <- compute_kernels(alpha_s, ell_s, ages, age, K_zs_s)
    Kss <- compute_kernels(alpha_s, ell_s, ages, ages, K_zs_ss)
    sig_s <- sig[s, 1]
    K_sum <- K[[1]] + K[[2]]
    Ks_sum <- Ks[[1]] + Ks[[2]]
    Kss_sum <- Kss[[1]] + Kss[[2]]
    f1_draws[s, ] <- draw_f(K_sum, Ks[[1]], Kss[[1]], sig_s, y_norm)
    f2_draws[s, ] <- draw_f(K_sum, Ks[[2]], Kss[[2]], sig_s, y_norm)
    f_draws[s, ] <- draw_f(K_sum, Ks_sum, Kss_sum, sig_s, y_norm)
    h_draws[s, ] <- lgpr:::apply_scaling(y_scl, f_draws[s, ], inverse = TRUE)
    if (s %% 400 == 0) {
      cat(s, " ")
    }
  }
  cat("\n")
  # Return
  fc <- list(f_draws)
  names(fc) <- "one"
  new("Prediction",
    f_comp = fc,
    f = f_draws,
    h = h_draws,
    x = x_star,
    extrapolated = FALSE
  )
}

# Compute predictions using several fitted models
compute_predictions <- function(fits, x_star) {
  preds <- list()
  j <- 0
  nams <- names(fits)
  for (f in fits) {
    j <- j + 1
    msg <- paste0("computing predictions for: ", nams[j])
    message(msg)
    if (is(f, "lgpfit")) {
      p <- pred_exact(f, x_star)
    } else {
      p <- pred_approx(f, x_star)
    }
    preds[[j]] <- p
  }
  names(preds) <- nams
  return(preds)
}
