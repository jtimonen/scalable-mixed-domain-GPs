#' An S4 class to represent input for approximate kernel matrix computations
#'
#' @slot input Common input (for example parameter values).
#' @slot K_input Input for computing kernel matrices between data points
#' (\code{N} x \code{N}). A list.
#' @slot Ks_input Input for computing kernel matrices between data and output
#' points (\code{P} x \code{N}). A list.
#' @slot Kss_input Input for computing kernel matrices between output
#' points (\code{P} x \code{P}). A list, empty if \code{full_covariance=FALSE}.
#' @slot comp_names Component names (character vector).
#' @slot full_covariance Boolean value determining if this can compute
#' full predictive covariance matrices (or just marginal variance at each point).
#' @slot no_separate_output_points Boolean value determining if
#' \code{Ks_input} and \code{Kss_input} are the same thing. Using this
#' knowledge can reduce unnecessary computations of kernel matrices.
#' @slot STREAM external pointer (for calling 'Stan' functions)
#' @param object The object for which to call a class method.
ApproximateKernelComputer <- setClass("ApproximateKernelComputer",
                           representation = representation(
                             input = "list",
                             K_input = "list",
                             Ks_input = "list",
                             Kss_input = "list",
                             comp_names = "character",
                             full_covariance = "logical",
                             STREAM = "externalptr",
                             no_separate_output_points = "logical"
                           )
)

