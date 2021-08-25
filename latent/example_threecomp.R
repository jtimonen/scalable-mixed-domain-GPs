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
n_per_N <- 10
sd <- simulate_data(
  N = 10, t_data = seq(1, 5, length.out = n_per_N),
  relevances = c(0, 1, 1),
  covariates = c(2),
  lengthscales = c(1.5, 0.75, 0.75), t_jitter = 0.2
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
max_x <- max(abs(model@stan_input$x_cont))
num_bf <- 30
scale_bf <- 1.25*max_x
si_add <- setup_approx(model, num_bf = num_bf, scale_bf = scale_bf)
stan_data <- c(model@stan_input, si_add)

# Create model and sample
sm1 <- stan_model("stan/lgp_latent_basisfun_threecomp.stan")
f1 <- sampling(sm1, data = stan_data, cores = 4)

# Create model and sample
#sm2 <- stan_model("stan/lgp_latent_covariance.stan")
#f2 <- sampling(sm2, data = stan_data, cores = 4)

# Compare
N <- stan_data$num_obs
cat("N=",N,"\n",sep="")
#plt <- plot_params_comparison_onecomp(f2, f1, ag_name = "approx", N = N)

# Compare functions
f_draws1 <- extract(f1, pars = "f_latent")$f_latent
fl_1 <- apply(f_draws1, c(2, 3), mean)
df1 <- data.frame(t(fl_1))
colnames(df1) <- c("f1")
df1 <- cbind(dat, df1)

#f_draws2 <- extract(f2, pars = "f_latent")$f_latent
#fl_2 <- apply(f_draws2, c(2, 3), mean)
#df2 <- data.frame(t(fl_2))
#colnames(df2) <- c("f1")
#df2 <- cbind(dat, df2)

#type <- as.factor(rep(c("approx", "exact"), each = N))
#df <- rbind(df1, df2)
#df$type <- type

#p1 <- ggplot(df, aes(x = age, y = f1, group = type, color = type)) +
#  geom_line()
# p2 <-  ggplot(df, aes(x=age,y=f2,group=z,color=type)) +geom_line()
# p3 <-  ggplot(df, aes(x=age,y=f3,group=id,color=type)) +geom_line()

p1 <- ggplot(df1, aes(x = age, y = f1)) + geom_line(color="firebrick") + 
  geom_point(data=dat, mapping=aes(x=age, y=y))
