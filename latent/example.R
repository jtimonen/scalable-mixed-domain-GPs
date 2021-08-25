# lgpr
library(lgpr)
if (packageVersion("lgpr") != "1.1.4") {
  stop("lgpr 1.1.4 is required!")
}

# Other requirements
library(rstan)
library(ggplot2)
library(ggpubr)
rstan_options(javascript = FALSE)
rstan_options(auto_write = TRUE)


# Simulate data using lgpr
n_per_N <- 14
sd <- simulate_data(
  N = 6, t_data = seq(1, 5, length.out = n_per_N),
  relevances = c(0, 1, 1),
  covariates = c(2),
  lengthscales = c(1, 0.5, 0.5), t_jitter = 0.5
)
dat <- sd@data

# Create model using lgpr
model <- create_model(y ~ age + age | z + age | id, dat, sample_f = TRUE)

# Source all R files
for (f in dir("R")) {
  path <- file.path("R", f)
  source(path)
}

# Create additional Stan input
num_bf <- 40
scale_bf <- 5.0
stan_data <- setup_approx(model, num_bf = num_bf, scale_bf = scale_bf)

# Create model and sample
sm1 <- stan_model("stan/lgp_latent_basisfun.stan")
f1 <- sampling(sm1, data = stan_data, cores = 4)

# Create model and sample
sm2 <- stan_model("stan/lgp_latent_covariance.stan")
f2 <- sampling(sm2, data = stan_data, cores = 4)

# Compare
N <- stan_data$num_obs
plt <- plot_params_comparison(f2, f1, ag_name = "approx", N = N)

# Compare functions
f1 <- extract(f1, pars="f_latent")$f_latent
f2 <- extract(f2, pars="f_latent")$f_latent
fl_1 <- apply(f1,c(2,3), mean)
fl_2 <- apply(f2,c(2,3), mean)
df1 <- data.frame(t(fl_1))
df2 <- data.frame(t(fl_2))
colnames(df1) <- c("f1", "f2", "f3")
colnames(df2) <- c("f1", "f2", "f3")
df1 <- cbind(model@data, df1)
df2 <- cbind(model@data, df2)

type <- as.factor(rep(c("approx", "exact"), each=N))
df <- rbind(df1, df2)
df$type <- type

p1 <-  ggplot(df, aes(x=age,y=f1,group=type,color=type)) +geom_line()
p2 <-  ggplot(df, aes(x=age,y=f2,group=z,color=type)) +geom_line()
p3 <-  ggplot(df, aes(x=age,y=f3,group=id,color=type)) +geom_line()

