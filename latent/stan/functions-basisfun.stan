// Evaluate m'th eigenfunction of the Dirichlet boundary value problem
// - x = points where to evaluate the function
// - m = index of the function
// - L = domain width
vector STAN_phi(vector x, int m, data real L){
  return sin(pi()/(2.0*L)*m*(x+L))/sqrt(L);
}

// Compute diagonal of diagonal matrix Lambda
vector STAN_diag_spd_eq(real alpha, real ell, data int M, data real L){
  vector[M] Lambda;
  for(m in 1:M) {
    real v = sqrt(2*pi()) * ell * exp(-0.5*(ell*pi()/2/L)^2 * m^2);
    Lambda[m] = alpha*sqrt(v+1e-9);
  }
  return(Lambda);
}
