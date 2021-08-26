// Evaluate m'th eigenfunction of the Dirichlet boundary value problem
// - x = points where to evaluate the function
// - m = index of the function
// - L = domain width
//vector STAN_phi_eq(vector x, int m, data real L){
//  return sin(pi()/(2.0*L)*m*(x+L))/sqrt(L);
//}

matrix STAN_PHI_eq(vector x, data vector seq_M, real L) {
  int N = num_elements(x);
  int M = num_elements(seq_M);
  matrix[N,M] A = rep_matrix(pi()/(2*L)*(x+L), M);
  return sin(diag_post_multiply(A, seq_M))/sqrt(L);
}

// Compute diagonal of diagonal matrix Lambda
vector STAN_diag_spd_eq(real alpha, real ell, data vector seq_M, data real L){
  return sqrt(2*pi()) * ell * exp(-0.5*(ell*pi()/2/L)^2 * seq_M .* seq_M);
}
