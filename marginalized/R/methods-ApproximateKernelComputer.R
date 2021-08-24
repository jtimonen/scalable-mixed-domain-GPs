#' @describeIn ApproximateKernelComputer Print a summary about the object.
setMethod("show", "ApproximateKernelComputer", function(object) {
  S <- num_paramsets(object)
  J <- num_components(object)
  P <- num_evalpoints(object)
  desc <- lgpr:::class_info("ApproximateKernelComputer")
  add_str <- paste0(
    "\n - ", S, " parameter sets",
    "\n - ", J, " components",
    "\n - ", P, " evaluation points "
  )
  add_str <- paste0(add_str)
  add_str <- paste0(
    add_str, "\n - full_covariance = FALSE\n"
  )
  desc <- paste0(desc, add_str)
  cat(desc)
})


#' @describeIn ApproximateKernelComputer Get number of components.
setMethod("num_components", "ApproximateKernelComputer", function(object) {
  K_const <- lgpr:::dollar(object@K_input, "K_const")
  length(K_const)
})

#' @describeIn ApproximateKernelComputer Get number of evaluation points.
setMethod("num_evalpoints", "ApproximateKernelComputer", function(object) {
  lgpr:::dollar(object@Ks_input, "n1")
})

# Get number of observations from ApproximateKernelComputer object
get_num_obs_kc <- function(object) {
  lgpr:::dollar(object@K_input, "n1")
}

#' @describeIn ApproximateKernelComputer Get number of parameter sets.
setMethod("num_paramsets", "ApproximateKernelComputer", function(object) {
  lgpr:::dollar(object@input, "num_paramsets")
})

#' @describeIn ApproximateKernelComputer Get component names.
setMethod("component_names", "ApproximateKernelComputer", function(object) {
  object@comp_names
})