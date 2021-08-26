// Evaluate m'th eigenfunction of the Dirichlet boundary value problem
// - x = points where to evaluate the function
// - m = index of the function
// - L = domain width
//vector STAN_phi_eq(vector x, int m, data real L){
//  return sin(pi()/(2.0*L)*m*(x+L))/sqrt(L);
//}

matrix STAN_PHI_eq(vector x, vector seq_M, real L) {
  int N = num_elements(x);
  int M = num_elements(seq_M);
  matrix[N,M] A = rep_matrix(pi()/(2*L)*(x+L), M);
  return sin(diag_post_multiply(A, seq_M))/sqrt(L);
}

// Compute diagonal of diagonal matrix Lambda
vector STAN_diag_spd_eq(real alpha, real ell, data int M, data real L,
    int num_reps){
  vector[M*num_reps] Lambda;
  for(m in 1:M) {
    int i1 = 1 + (m-1)*num_reps;
    int i2 = m*num_reps;
    real v = sqrt(2*pi()) * ell * exp(-0.5*(ell*pi()/2/L)^2 * m^2);
    Lambda[i1:i2] = rep_vector(alpha*sqrt(v+1e-9), num_reps);
  }
  return(Lambda);
}
