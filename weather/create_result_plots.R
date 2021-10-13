# Startup
source(normalizePath(file.path("..", "common.R")))
outdir <- startup(create_dir = FALSE)
fn <- file.path("res_01", "res_16.rds")
dat <- load_weather_data()
dat$day <- as.Date(dat$day, origin = "1964-01-01")
results <- readRDS(file = fn)
pa <- results$pred

# Create y-axis labels
ylab1 <- expression(paste("", "f", phantom()^{
  paste("", "(", "1", ")", "")
}, "(", "", "day", ")", "", ""))
ylab2 <- expression(paste("", "f", phantom()^{
  paste("", "(", "2", ")", "")
}, "(", "", "day", ",", "region", ")", "", ""))
ylab3 <- expression(paste("", "f", phantom()^{
  paste("", "(", "3", ")", "")
}, "(", "", "day", ",", "station", ")", "", ""))


pf1 <- plot_f(pa, 1, aes1()) + geom_ribbon(alpha = 0.3) + ylab(ylab1)
pf2 <- plot_f(pa, 2, aes2()) + geom_ribbon(alpha = 0.3) + ylab(ylab2)
pf3 <- plot_f(pa, 3, aes3()) + facet_wrap(. ~ region) + ylab(ylab3) +
  theme(legend.text=element_text(size=8), legend.title = element_blank())
ph <- plot_f(pa, 0, aes3()) +
  geom_point(
    data = dat,
    inherit.aes = FALSE,
    aes(x = day, y = temperature, group = station),
    pch = ".", alpha = 0.6,
  ) + facet_wrap(. ~ station) + geom_ribbon() + ylab("Temperature (C)")

plt_12 <- ggarrange(pf1, pf2, nrow=1, ncol=2,
                    labels = c("a", "b"), widths = c(0.8, 1))
full_plt <- ggarrange(plt_12, pf3, nrow=2, ncol=1, heights = c(0.6, 1),
                      labels=c("", "c"))
ggsave(full_plt, file = "weather.pdf", width = 6.9, height = 6.7)
