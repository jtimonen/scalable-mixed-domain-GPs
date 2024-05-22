# Create a configuration
create_configuration <- function(num_bf, scale_bf) {
  new(
    "ExperimentConfiguration",
    num_bf = as.integer(num_bf),
    scale_bf = scale_bf
  )
}
