# Aesthetic for component 1
aes1 <- function() {
  aes(
    x = year, y = f_mean, ymin = f_mean - 2 * f_std,
    ymax = f_mean + 2 * f_std
  )
}

# Aesthetic for component 2
aes2 <- function() {
  aes(
    x = year, y = f_mean, ymin = f_mean - 2 * f_std,
    ymax = f_mean + 2 * f_std, group = region, color = region,
    fill = region, lty = region
  )
}

# Aesthetic for component 3
aes3 <- function() {
  aes(
    x = year, y = f_mean, ymin = f_mean - 2 * f_std,
    ymax = f_mean + 2 * f_std, group = state, color = state,
    fill = state
  )
}

# Plot f or it's component
plot_f <- function(dat, pred, idx = 0, aesth) {
  if (idx == 0) {
    f <- pred@h
    ylab <- "h"
  } else {
    f <- pred@f_comp[[idx]]
    ylab <- names(pred@f_comp)[idx]
  }
  f_mean <- colMeans(f)
  f_std <- apply(f, 2, stats::sd)
  df <- cbind(dat, f_mean, f_std)
  plt <- ggplot(df, aesth) +
    geom_line() +
    ylab(ylab) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("")
}
