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