# Expose all Stan functions without creating a complete Stan model
expose_stanfuns <- function(STAN_HOME = "stan") {
  FILES <- c(
    "functions-utils.stan",
    "functions-kernels.stan",
    "functions-prior.stan",
    "functions-approx.stan"
  )

  # Create Stan model containing only a functions block with all the functions
  two_spaces <- "  "
  f_list <- lapply(file.path(STAN_HOME, FILES), FUN = readLines)
  functions <- paste(unlist(f_list), collapse = paste0("\n", two_spaces))
  functions <- paste0(two_spaces, functions)
  model_code <- paste(c("functions {", functions, "}"), collapse = "\n")
  header <- "// Automatically generated\n\n"
  model_code <- paste0(header, model_code, "\n")

  # Write Stan code to file and expose
  fn <- file.path(STAN_HOME, "all_functions.stan")
  cat(model_code, file = fn)
  rstan::expose_stan_functions(fn, verbose = TRUE)
  file.remove(fn)
}
