  // Compute total signal by summing components
    vector STAN_vectorsum(vector[] vecs, data int L){
    int num_vecs = size(vecs);
    vector[L] s = rep_vector(0, L);
    for (j in 1:num_vecs){
      s += vecs[j];
    }
    return(s);
  }
  
  // Create vector with elements 1, ..., N
  vector STAN_seq_len(int N){
    vector[N] v = rep_vector(1.0, N);
    for(n in 2:N) { 
      v[n] = n;
    }
    return(v);
  }
  
  // Same as rep(x, times=N) in R
  vector STAN_rep_vector_times(vector x, int N) {
    return to_vector(rep_matrix(x, N));
  }
  
    // Same as rep(x, each=N) in R
  vector STAN_rep_vector_each(vector x, int N) {
    return to_vector(transpose(rep_matrix(x, N)));
  }
  
  // Wrap to functions into one
  vector STAN_seq_len_rep_times(int M, int N) {
    return STAN_rep_vector_times(STAN_seq_len(M), N);
  }
  
    // Wrap two functions into one
  vector STAN_seq_len_rep_each(int M, int N) {
    return STAN_rep_vector_each(STAN_seq_len(M), N);
  }
  