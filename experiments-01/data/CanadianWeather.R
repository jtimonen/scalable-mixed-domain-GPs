library(fda)
library(lgpr)
temp <- CanadianWeather$dailyAv[, , 1]
precip <- CanadianWeather$dailyAv[, , 3]
n_days <- nrow(temp)
n_stations <- ncol(temp)
regions <- rep(as.factor(CanadianWeather$region), each = n_days)
stations <- rep(as.factor(CanadianWeather$place), each = n_days)
day <- rep(c(1:n_days), times = n_stations)
dat <- data.frame(day, as.vector(temp), as.vector(precip), stations, regions)
colnames(dat) <- c("day", "temperature", "precipitation", "station", "region")

# Plot
p1 <- plot_data(dat,
  y_name = "temperature", group_by = "station",
  color_by = "station",
  x_name = "day", facet_by = "region",
  main = "Temperature"
)

p2 <- plot_data(dat,
  y_name = "precipitation", group_by = "station",
  color_by = "station",
  x_name = "day", facet_by = "region",
  main = "Precipitation (log10)"
)
