#' Model class for approximate models
#'
#' @slot exact_model the corresponding exact model
#' @slot added_stan_input additional Stan input
#' @slot stan_file Stan file where the approximate model code is
#' @slot obs_model observation model name
setClass(
  "ApproxModel",
  representation(
    exact_model = "lgpmodel",
    added_stan_input = "list",
    stan_file = "character",
    obs_model = "character"
  )
)

#' Print info
setMethod("show", signature(object = "ApproxModel"), function(object) {
  num_bf <- object@added_stan_input$num_bf
  scale_bf <- object@added_stan_input$scale_bf
  desc <- paste0(
    "An approximate model with num_bf=",
    num_bf, ", scale_bf=", scale_bf, ". Observation model = ", object@obs_model,
    ". Stan model file = ", object@stan_file, "."
  )
  cat(desc, "\n")
})

#' Class for experiment configurations
#'
#' @slot num_bf number of basis functions
#' @slot scale_bf basis function domain scale
setClass(
  "ExperimentConfiguration",
  representation(
    num_bf = "integer",
    scale_bf = "numeric"
  )
)

# Get brief experiment info
num_bf_string <- function(object) {
  nbf <- formatC(object@num_bf, width = 3)
  paste0("B =", nbf)
}
scale_string <- function(object) {
  paste0("c = ", object@scale_bf)
}
experiment_info <- function(object) {
  paste0(num_bf_string(object), ", ", scale_string(object))
}

#' Print info
setMethod(
  "show", signature(object = "ExperimentConfiguration"),
  function(object) {
    ei <- experiment_info(object)
    desc <- paste0("<ExperimentConfiguration(", ei, ")>")
    cat(desc, "\n")
  }
)

#' Model class for approximate model fits
#'
#' @slot model the approximate model
#' @slot fit fit object
#' @slot backend used backend
#' @slot info additional info list
setClass(
  "ApproxModelFit",
  representation(
    model = "ApproxModel",
    fit = "list",
    backend = "character",
    info = "list"
  )
)

# Get the CmdStanFit of an ApproxModelFit
get_cmdstanfit <- function(fit) {
  stopifnot(is(fit, "ApproxModelFit"))
  fit@fit[[1]]
}

# Get full 'Stan' input of an ApproxModel
get_full_stan_input <- function(approx_model) {
  stopifnot(is(approx_model, "ApproxModel"))
  si <- approx_model@exact_model@stan_input
  si_add <- approx_model@added_stan_input
  c(si, si_add)
}
