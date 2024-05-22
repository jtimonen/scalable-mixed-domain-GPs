  // Approximation domain and dimensions
  vector<lower=0>[num_cov_cont] X_sd;
  vector<lower=0>[num_cov_cont] X_hr;
  real<lower=0> scale_bf; // factor c to determine domain width L
  int<lower=1> num_bf;    // number of basis functions
  int<lower=0> num_xi[num_comps];
  
  // Categorical decomposition stuff
  int<lower=0> C_sizes[num_comps];
  int<lower=0> C_ranks[num_comps];
  int<lower=0> C_rsp[num_comps]; // element wise product of C_sizes and C_ranks
  vector[sum(C_rsp)] C_vecs;
  real C_vals[sum(C_ranks)];
