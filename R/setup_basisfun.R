# Create Stan input for bf model
setup_basisfun <- function(model, num_bf, width) {
  stopifnot(is(model, "lgpmodel"))
  si <- model@stan_input
}
