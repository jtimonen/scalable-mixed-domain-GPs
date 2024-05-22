# Standardize variable to zero mean and unit variance
normalize_var <- function(x) (x - mean(x)) / stats::sd(x)

# Create reference model formula for lgpr2, out of covariate names
create_full_formula <- function(xn, zn) {
  part1 <- paste("gp(", xn, ")", collapse = " + ")
  part2 <- paste(paste0("gp(age, ", zn), ")", collapse = " + ")
  f <- paste0("y ~ offset(id) + gp(age) + ", part1, " + ", part2)
  as.formula(f)
}

# Compute true signal-to-noise ratio
compute_snr <- function(components) {
  y <- components$y
  h <- components$h
  var(h) / (var(y - h))
}
