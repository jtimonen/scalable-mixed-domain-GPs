library(MASS)
library(lgpr)

# Melanoma (MASS)
dat <- Melanoma
dat$sex[dat$sex == 0] <- "Female"
dat$sex[dat$sex == 1] <- "Male"
dat$status[dat$status == 1] <- "Dead (melanoma)"
dat$status[dat$status == 2] <- "Alive"
dat$status[dat$status == 3] <- "Dead (other cause)"
dat$ulcer[dat$ulcer == 0] <- "No"
dat$ulcer[dat$ulcer == 1] <- "Yes"
dat$sex <- as.factor(dat$sex)
dat$status <- as.factor(dat$status)
dat$ulcer <- as.factor(dat$ulcer)
dat <- data.frame(dat)

# Plot
plot_data(dat,
  y_name = "thickness", group_by = "sex", color_by = "sex",
  facet_by = "status"
)
