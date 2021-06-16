# Create Stan input for bf model
setup_basisfun <- function(model, M_bf, c_bf) {
  stopifnot(is(model, "lgpmodel"))
  si <- model@stan_input
  si_bf <- list(M_bf = M_bf, c_bf = c_bf) 
  c(si, si_bf)
}
