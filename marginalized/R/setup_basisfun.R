# Create Stan input for bf model
setup_basisfun <- function(model, M_bf, L_bf) {
  stopifnot(is(model, "lgpmodel"))
  si <- model@stan_input
  c_bf <- 1.5 # not used
  si_bf <- list(M_bf = M_bf, c_bf = c_bf, L_bf = L_bf)
  c(si, si_bf)
}
