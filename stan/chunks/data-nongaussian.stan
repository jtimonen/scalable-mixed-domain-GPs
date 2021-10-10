  // Inputs specific to non-gaussian observation model
  int<lower=2,upper=5> obs_model; // 2-5: Poisson, NB, Bin, BB
  int<lower=0> y_int[obs_model>1, num_obs]; // response variable (int)
  int<lower=1> y_num_trials[obs_model>3, num_obs]; // for Bin or BB model
  int<lower=0> prior_phi[obs_model==3, 2];
  real hyper_phi[obs_model==3, 3];
  real hyper_gamma[obs_model==5, 2];
  vector[num_obs] c_hat; // constant vector added to the GP function f
