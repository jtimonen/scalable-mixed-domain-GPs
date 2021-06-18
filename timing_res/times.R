# Runtime data
n1 <- 4806
n2 <- 2406
n3 <- 966
n4 <- 486
t1 <- c(2439.2, 2471.74, 2318.14, 2479.18)
t2 <- c(658.31, 674.21, 663.2, 649.32)
t3 <- c(461.99, 438.91, 441.18, 432.38)
t4 <- c(139.09, 142.48, 143.59, 143.28)
nodes <- c("csl40", "skl33", "csl39", "csl40")

# Create data frame
n <- rep(c(n1, n2, n3, n4), each = 4)
time <- c(t1, t2, t3, t4)
node <- rep(nodes, each = 4)
df <- data.frame(n, time, node)

# Plot
library(ggplot2)
plt <- ggplot(df, aes(x = n, y = time, color = node)) +
  geom_point(pch = 4) +
  stat_smooth(method = "lm", col = "gray30") +
  theme_minimal() +
  ylab("Runtime per chain (s)") +
  xlab("n") +
  ggtitle("Model: y ~ age + age | z", subtitle = "J=2, R=0+2=2, M=25, L=3.0") +
  ylim(0, 2500) +
  xlim(0, 5000)
