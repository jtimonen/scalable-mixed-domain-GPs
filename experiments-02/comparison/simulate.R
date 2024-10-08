library(lgpr) # used only for data simulation


# Simulate data
simulate_gaussian <- function(N_indiv, n_unrel, snr, t_seq) {
  N <- N_indiv * length(t_seq)
  message("N = ", N)
  covariates <- c(1, 2)
  lengthscales <- rep(12, 2 + sum(covariates %in% c(0, 1, 2)))
  lengthscales[1] <- 10^5 # effectively make id*age effect be constant
  simdat <- lgpr::simulate_data(
    N = N_indiv,
    t_data = t_seq,
    covariates = covariates,
    relevances = c(1, 1, 1, 1),
    t_jitter = 1,
    noise_type = "gaussian",
    c_hat = 0,
    f_var = 16,
    snr = snr,
    lengthscales = lengthscales
  )
  n_flip_z <- round(N_indiv / 3)
  data_new <- add_unrelated(
    simdat@data,
    n_x = n_unrel,
    n_z = n_unrel,
    n_flip_z = n_flip_z
  )
  data_new$components <- simdat@components
  data_new
}

# Simulate data
simulate_gaussian_v2 <- function(N_indiv, n_unrel, snr, t_seq) {
  N <- N_indiv * length(t_seq)
  message("N = ", N)
  covariates <- c(1, 1, 2, 2)
  lengthscales <- rep(12, 2 + sum(covariates %in% c(0, 1, 2)))
  lengthscales[1] <- 10^5 # effectively make id*age effect be constant
  simdat <- lgpr::simulate_data(
    N = N_indiv,
    t_data = t_seq,
    covariates = covariates,
    relevances = c(1, 1, 1, 1, 1, 1),
    t_jitter = 1,
    noise_type = "gaussian",
    c_hat = 0,
    f_var = 16,
    snr = snr,
    lengthscales = lengthscales,
    n_categs = c(2, 3),
    names = c("x", "w", "z", "r")
  )
  n_flip_z <- round(N_indiv / 3)
  data_new <- add_unrelated(
    simdat@data,
    n_x = n_unrel,
    n_z = n_unrel,
    n_flip_z = n_flip_z
  )
  data_new$xn <- c(data_new$xn, "w")
  data_new$zn <- c(data_new$zn, "r")
  data_new$components <- simdat@components
  data_new
}


# Add unrelated variables
create_correlated_x <- function(x, p) {
  s <- stats::sd(x)
  rnorm(n = length(x), mean = x, sd = p * s)
}

# Flip categorical
flip_z_of_id <- function(z, id, id_to_flip) {
  levs <- as.numeric(levels(z))
  ids <- as.numeric(levels(id))
  rows_edit <- which(as.numeric(id) %in% id_to_flip)
  cat_to_flip <- unique(as.numeric(z)[rows_edit])
  hmm <- setdiff(levs, cat_to_flip)
  new_cat <- sample(as.character(hmm), size = 1)
  z_new <- as.numeric(z)
  z_new[rows_edit] <- as.numeric(new_cat)
  as.factor(z_new)
}

# Add unrelated variables
create_correlated_z <- function(id, z, n_edit) {
  levs <- as.numeric(levels(z))
  ids <- as.numeric(levels(id))
  ids_remaining <- ids
  for (n in seq_len(n_edit)) {
    id_to_flip <- sample(ids_remaining, size = 1)
    z <- flip_z_of_id(z, id, id_to_flip)
    ids_remaining <- setdiff(ids_remaining, id_to_flip)
  }
  z
}

# Add unrelated variables
create_new_x <- function(bn, dat, p) {
  for (j in 1:length(p)) {
    x_base <- dat[[bn]]
    xn <- paste0(bn, "_u", j)
    dat[[xn]] <- create_correlated_x(x_base, p = p[j])
  }
  dat
}

# Add unrelated variables
create_new_z <- function(bn, dat, n) {
  for (j in 1:length(n)) {
    z_base <- dat[[bn]]
    zn <- paste0(bn, "_u", j)
    dat[[zn]] <- create_correlated_z(dat$id, z_base, n[j])
  }
  dat
}

# Add unrelated variables
add_unrelated <- function(dat, n_x = 4, n_z = 4, n_flip_z = 4) {
  p_x <- rep(0.85, n_x)
  n_z_edit <- rep(n_flip_z, n_z)
  dat <- create_new_x("x", dat, p_x)
  dat <- create_new_z("z", dat, n_z_edit)
  cnd <- colnames(dat)
  list(
    dat = dat,
    xn = cnd[which(grepl(cnd, pattern = "x"))],
    zn = cnd[which(grepl(cnd, pattern = "z"))],
    px = p_x,
    nze = n_z_edit
  )
}
