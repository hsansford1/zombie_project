#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
std::vector<double> simulate_zombies_rcpp(vec& params) {
  
  double delta = params[0];
  double zeta = params[1];
  double beta = params[2];
  double alpha = params[3];
  double pi = params[4];
  double rho = params[5];

  int N = 500;             // total population
  int time_horizon = 25;
  double dt = 0.05;
  int n = time_horizon/dt;

  vec S(n+1, fill::zeros);
  vec I(n+1, fill::zeros);
  vec Z(n+1, fill::zeros);
  vec R(n+1, fill::zeros);

  // Initial conditions
  S[0] = N - 1;
  Z[0] = 1;

  // Define model dynamics (ODEs) and solve by Euler's method
  for(int i = 0; i < n; i++){

    if(S[i] < 0){ S[i] = 0; }
    if(I[i] < 0){ I[i] = 0; }
    if(Z[i] < 0){ Z[i] = 0; }
    if(R[i] < 0){ R[i] = 0; }

    S[i+1] = S[i] + dt*(pi*S[i] - beta*S[i]*Z[i] - delta*S[i]);
    I[i+1] = I[i] + dt*(beta*S[i]*Z[i] - rho*I[i] - delta*I[i]);
    Z[i+1] = Z[i] + dt*(rho*I[i] - alpha*S[i]*Z[i] + zeta*R[i]);
    R[i+1] = R[i] + dt*(alpha*S[i]*Z[i] + delta*(S[i]+I[i]) - zeta*R[i]);
  }

  uvec idx = conv_to<uvec>::from(regspace(0, 1/dt, 500));
  vec daily_S = S(idx);
  vec daily_Z = Z(idx);
  
  vec dailies = join_cols(daily_S, daily_Z);
  return as< std::vector<double> >(wrap(dailies));
}

// [[Rcpp::export]]
std::vector<double> simulate_zombies_seed_rcpp(vec& params) {
  
  // arma_rng::set_seed(params[0]);
  double delta = params[1];
  double zeta = params[2];
  double beta = params[3];
  double alpha = params[4];
  double pi = params[5];
  double rho = params[6];
  
  int N = 500;             // total population
  int time_horizon = 25;
  double dt = 0.05;
  int n = time_horizon/dt;
  
  vec S(n+1, fill::zeros);
  vec I(n+1, fill::zeros);
  vec Z(n+1, fill::zeros);
  vec R(n+1, fill::zeros);
  
  // Initial conditions
  S[0] = N - 1;
  Z[0] = 1;
  
  // Define model dynamics (ODEs) and solve by Euler's method
  for(int i = 0; i < n; i++){
    
    if(S[i] < 0){ S[i] = 0; }
    if(I[i] < 0){ I[i] = 0; }
    if(Z[i] < 0){ Z[i] = 0; }
    if(R[i] < 0){ R[i] = 0; }
    
    S[i+1] = S[i] + dt*(pi*S[i] - beta*S[i]*Z[i] - delta*S[i]);
    I[i+1] = I[i] + dt*(beta*S[i]*Z[i] - rho*I[i] - delta*I[i]);
    Z[i+1] = Z[i] + dt*(rho*I[i] - alpha*S[i]*Z[i] + zeta*R[i]);
    R[i+1] = R[i] + dt*(alpha*S[i]*Z[i] + delta*(S[i]+I[i]) - zeta*R[i]);
  }
  
  uvec idx = conv_to<uvec>::from(regspace(0, 1/dt, 500));
  vec daily_S = S(idx);
  vec daily_Z = Z(idx);
  
  vec dailies = join_cols(daily_S, daily_Z);
  return as< std::vector<double> >(wrap(dailies));
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
// 
// # /*** R
// # t1 <- simulate_zombies_rcpp(rep(0.1,6))
// # t2 <- simulate_zombies_seed_rcpp(rep(0.1,7))
// # all.equal(t1, t2)
// #  */
