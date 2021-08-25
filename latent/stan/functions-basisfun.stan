// Evaluate m'th eigenfunction of the Dirichlet boundary value problem
// - x = points where to evaluate the function
// - m = index of the function
// - L = domain width
vector STAN_phi(vector x, int m, data real L){
  real A = inv(sqrt(L));
  real B = pi()*m/(2.0*L);
  return(A*sin(B*(x+L)));
}

// Evaluate m'th eigenvalue of the Dirichlet boundary value problem
// - m = index of the function
// - L = domain width
real STAN_lambda(int m, data real L){
  real A = pi()*m/(2.0*L);
  return(square(A));
}

// Evaluate spectral density of the exponentiated quadratic kernel
// - w = frequency
// - ell = lengthscale of the kernel
real STAN_spd_eq(real w, real ell){
  real A = ell*sqrt(2.0*pi());
  real B = - 0.5*square(w*ell);
  return(A*exp(B) + 1e-9);
}

// Compute diagonal of diagonal matrix Lambda
vector STAN_lambda_matrix(real ell, data int M, data real L){
  vector[M] Lambda;
  for(m in 1:M) {
    Lambda[m] = STAN_spd_eq((m*pi())/(2.0*L), ell);
  }
  return(Lambda);
}
