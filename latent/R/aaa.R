
#' Model class for approximate models
#'
#' @slot exact_model the corresponding exact model
#' @slot add_stan_input additional Stan input
#' @slot stan_file Stan file where the approximate model code is
setClass(
  "ApproxModel",
  representation(
    exact_model = "lgpmodel",
    add_stan_input = "list",
    stan_file = "character"
  )
)

#' Print info
setMethod("show", signature(object = "ApproxModel"), function(object) {
  num_bf <- object@add_stan_input$num_bf
  scale_bf <- object@add_stan_input$scale_bf
  desc <- paste0(
    "An approximate model with num_bf=",
    num_bf, ", scale_bf=", scale_bf, "."
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
experiment_info <- function(object) {
  desc <- paste0(
    "num_bf=", object@num_bf, ", scale_bf=", object@scale_bf
  )
  return(desc)
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
setClass(
  "ApproxModelFit",
  representation(
    model = "ApproxModel",
    fit = "list",
    backend = "character"
  )
)
