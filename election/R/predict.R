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
pred_approx <- function(fit, x_star) {
  stopifnot(is(fit, "ApproxModelFit"))
  om <- fit@model@obs_model
  cat(" * calling posterior_f_approx\n")
  fp <- posterior_f_approx(fit, x_star)
  emodel <- fit@model@exact_model
  S <- length(fp)
  J <- length(fp[[1]])
  P <- length(fp[[1]][[1]])
  f_comp <- list()
  f_sum <- matrix(0.0, S, P)

  # GP components
  for (j in seq_len(J)) {
    fj <- get_approx_component(fp, j)
    f_comp[[j]] <- fj
    f_sum <- f_sum + fj
  }

  # Other
  mu <- posterior::merge_chains(get_cmdstanfit(fit)$draws("mu"))
  f_comp_mu <- matrix(rep(mu[, 1, , drop = T], P), S, P, byrow = FALSE)
  f_comp[[J + 1]] <- f_comp_mu
  f_sum <- f_sum + f_comp_mu
  names(f_comp) <- c(component_names(emodel), "mu")
  c_hat_pred <- rep(0.0, P)
  nu <- as.vector(posterior::merge_chains(get_cmdstanfit(fit)$draws("nu")))
  h <- exp(f_sum) / (1 + exp(f_sum))

  # Draws of y
  cat(" * drawing y\n")
  y_rng <- matrix(0.0, S, P)
  for (s in 1:S) {
    mu_obs <- nu[s] * h[s, ]
    for (p in 1:P) {
      bb <- mu_obs[p]
      y_rng[s, p] <- stats::rbeta(n = 1, shape1 = bb, shape2 = nu[s] - bb)
    }
  }

  # Return
  pred <- new("Prediction",
    f_comp = f_comp,
    f = f_sum,
    h = h,
    x = x_star,
    extrapolated = FALSE
  )
  list(pred = pred, y_rng = y_rng)
}
