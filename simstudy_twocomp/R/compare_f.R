
# Wraps STAN_build_f_latent
build_f <- function(stan_data, tdata, alpha, ell, xi) {
  expose_stanfuns()
  N <- stan_data$num_obs
  seq_M <- STAN_seq_len(stan_data$num_bf)
  scale_bf <- stan_data$scale_bf
  num_xi <- as.array(stan_data$num_xi)
  comps <- matrix_to_list(stan_data$components)
  C_ranks <- as.array(stan_data$C_ranks)
  STAN_build_f(
    comps, num_xi, C_ranks, seq_M, tdata$L, tdata$PSI, alpha, ell, xi
  )
}

# Wraps STAN_build_f_draws
build_f_draws <- function(stan_data, tdata, ALPHA, ELL, XI) {
  expose_stanfuns()
  N <- stan_data$num_obs
  seq_M <- STAN_seq_len(stan_data$num_bf)
  scale_bf <- stan_data$scale_bf
  num_xi <- as.array(stan_data$num_xi)
  comps <- matrix_to_list(stan_data$components)
  C_ranks <- as.array(stan_data$C_ranks)
  STAN_build_f_draws(
    comps, num_xi, C_ranks, seq_M, tdata$L, tdata$PSI, ALPHA, ELL, XI
  )
}

# Draw xi and build latent
draw_f_latent <- function(stan_data, tdata, alpha = NULL, ell = NULL) {
  if (is.null(alpha)) {
    alpha <- rep(1.0, stan_data$num_comps)
  }
  if (is.null(ell)) {
    ell <- rep(1.0, stan_data$num_ell)
  }
  alpha <- as.array(alpha)
  ell <- as.array(ell)
  P <- sum(stan_data$num_xi)
  xi <- rnorm(n = P)
  build_f_latent(stan_data, tdata, alpha, ell, xi)
}

# Colsums for draws_rvar
colsums_rvar <- function(x) {
  D <- dim(x)[1]
  rv <- x[1, ]
  for (d in 2:D) {
    rv <- rv + x[d, ]
  }
  return(rv)
}
# Helper function
create_plotf_df <- function(data, fit, comp_idx) {
  if (is.null(comp_idx)) stop("comp_idx cannot be NULL")
  if (is(fit, "lgpfit")) {
    gpred <- pred(fit)
    if (comp_idx == 0) {
      f_mean <- as.vector(gpred@f_mean)
      f_sd <- as.vector(gpred@f_std)
    } else {
      f_mean <- as.vector(gpred@f_comp_mean[[comp_idx]])
      f_sd <- as.vector(gpred@f_comp_std[[comp_idx]])
    }
  } else if (is(fit, "stanfit")) {
    rv <- as_draws_rvars(fit)$f_latent
    if (comp_idx == 0) {
      rv <- colsums_rvar(rv)
    } else {
      rv <- rv[comp_idx, ]
    }
    f_mean <- as.vector(mean(rv))
    f_sd <- as.vector(sd(rv))
  } else if (is(fit, "CmdStanMCMC")) {
    rv <- as_draws_rvars(fit$draws())$f_latent
    if (comp_idx == 0) {
      rv <- colsums_rvar(rv)
    } else {
      rv <- rv[comp_idx, ]
    }
    f_mean <- as.vector(mean(rv))
    f_sd <- as.vector(sd(rv))
  } else {
    stop("fit should be a stanfit, lgpfit or CmdStanMCMC object!")
  }
  cbind(data, f_mean, f_sd)
}

# Plot
plot_f <- function(data, fit, aname, comp_idx) {
  df <- create_plotf_df(data, fit, comp_idx)
  map <- aes(x = age, ymin = f_mean - 2 * f_sd, ymax = f_mean + 2 * f_sd)
  plt <- ggplot(df, aes(x = age, y = f_mean)) +
    geom_ribbon(aes = map, alpha = 0.3) +
    geom_line() +
    ggtitle(aname) +
    ylab("") +
    theme_bw() +
    geom_point(
      data = data, aes(x = age, y = y), pch = 4,
      alpha = 0.3, color = "firebrick3"
    )
  return(plt)
}

# Plot comparison
plot_f_compare_with_exact <- function(data, fit, fit_approx, aname,
                                      comp_idx = 0,
                                      aes = NULL,
                                      plot_data = FALSE) {
  if (is.null(aes)) {
    aes <- aes(x = age, y = f_mean, group = id_x_fit, color = fit)
  }
  df <- create_plotf_df(data, fit, comp_idx)
  df_approx <- create_plotf_df(data, fit_approx, comp_idx)
  N1 <- nrow(df)
  N2 <- nrow(df_approx)
  fit <- as.factor(c(rep("exact", N1), rep(aname, N2)))
  df <- rbind(df, df_approx)
  df <- cbind(df, fit)
  df$z <- paste("z =", df$z)
  data$z <- paste("z =", data$z)
  df$id_x_fit <- paste(df$fit, df$id)
  df$z_x_fit <- paste(df$fit, df$z)
  plt <- ggplot(df, mapping = aes)
  ribbon <- TRUE
  if (ribbon) {
    plt <- plt + geom_ribbon(
      aes(x = age, ymin = f_mean - 2 * f_sd, ymax = f_mean + 2 * f_sd),
      alpha = 0.2, linetype = 2
    )
  }
  plt <- plt + geom_line() + theme_bw() + ylab("posterior f")
  if (plot_data) {
    plt <- plt + geom_point(
      data = data, aes(x = age, y = y), inherit.aes = FALSE,
      pch = 4, alpha = 0.7
    )
  }
  return(plt)
}

# Plot means comparison in same panel
plot_f_compare_same <- function(data, fits) {
  nams <- names(fits)
  J <- length(fits)
  DF <- NULL
  TYP <- NULL
  for (j in 1:J) {
    df <- create_plot_df(data, fits[[j]])
    DF <- rbind(DF, df)
    TYP <- c(TYP, c(rep(nams[j], nrow(df))))
  }
  fit <- as.factor(TYP)
  plt <- ggplot(DF, aes(x = age, y = f_mean, group = fit, color = fit))
  plt <- plt + geom_line()
  plt <- plt + geom_point(
    data = data, aes(x = age, y = y), inherit.aes = FALSE,
    pch = 4, alpha = 0.3
  )
  plt <- plt + theme_bw() + ylab("posterior mean f")
  return(plt)
}

# Plot mean and sd comparison in separate figures
plot_f_compare_separate <- function(data, fits, last_is_exact = FALSE, ...) {
  PLOTS <- list()
  J <- length(fits)
  nams <- names(fits)
  if (last_is_exact) {
    fit_exact <- fits[[J]]
    fits <- fits[1:(J - 1)]
    J <- length(fits)
  } else {
    fit_exact <- NULL
  }
  for (j in seq_len(J)) {
    if (!(is.null(fit_exact))) {
      plt <- plot_f_compare_with_exact(data, fit_exact, fits[[j]],
        aname = nams[j], ...
      )
    } else {
      plt <- plot_f(data, fits[[j]], nams[j], ...)
    }
    thm <- theme(legend.position = "top", legend.title = element_blank())
    PLOTS[[j]] <- plt + thm + ylab("")
  }
  return(PLOTS)
}
