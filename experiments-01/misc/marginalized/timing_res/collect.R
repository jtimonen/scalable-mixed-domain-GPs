M <- 35
df <- list()
for (j in 0:29) {
  fn <- paste0("out/times_", M, "_", j, ".txt")
  a <- read.csv(fn)$x
  df <- rbind(df, a)
}
colnames(df) <- c("n", "t1", "t2", "t3", "t4")
df <- as.data.frame(df)
print(df)

fn_out <- paste0("times_", M, ".rds")
saveRDS(df, file = fn_out)
