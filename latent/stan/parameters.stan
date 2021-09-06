  real<lower=1e-12> alpha[num_comps]; // component magnitudes
  real<lower=1e-12> ell[num_ell]; // lengthscales
  real<lower=1e-12> sigma[obs_model==1];
  real<lower=1e-12> phi[obs_model==3];
  real<lower=1e-12, upper=1-1e-12> gamma[obs_model==5];
