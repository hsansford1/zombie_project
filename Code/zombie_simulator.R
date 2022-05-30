
## Function is in this format so that it can be used with `EasyABC` package!!

## Function to simulate zombie epidemic:
## params is a vector containing the following parameters in order:
#
  # delta: background (non-zombie related) death rate
  # zeta: zombie resurrection rate
  # beta: zombie transmission rate
  # alpha: zombie destruction rate
  # pi: birth rate
  # rho: latent infection rate

simulate_zombies <- function(params){
  
  delta <- params[1] 
  zeta <- params[2] 
  beta <- params[3]
  alpha <- params[4] 
  pi <- params[5]
  rho <- params[6]

  N <- 500             # total population
  time_horizon <- 25
  dt <- 0.05
  n <- time_horizon/dt
  t <- seq(0, time_horizon, dt)
  
  S <- rep(0, n+1)
  I <- rep(0, n+1)
  Z <- rep(0, n+1)
  R <- rep(0, n+1)
  
  # Initial conditions
  S[1] <- N -1
  Z[1] <- 1
  
  # Define model dynamics (ODEs) and solve by Euler's method
  for (i in 1:n){
    
    if (S[i] < 0){
      S[i] = 0
    }
    if (I[i] < 0){
      I[i] = 0
    }
    if (Z[i] < 0){
      Z[i] = 0
    }
    if (R[i] < 0){
      R[i] = 0
    }
    
    S[i+1] <- S[i] + dt*(pi*S[i] - beta*S[i]*Z[i] - delta*S[i])
    I[i+1] <- I[i] + dt*(beta*S[i]*Z[i] - rho*I[i] - delta*I[i])
    Z[i+1] <- Z[i] + dt*(rho*I[i] - alpha*S[i]*Z[i] + zeta*R[i])
    R[i+1] <- R[i] + dt*(alpha*S[i]*Z[i] + delta*(S[i]+I[i]) - zeta*R[i])
  }

  daily_S <- S[seq(1,501,20)]
  daily_Z <- Z[seq(1,501,20)]
  # return(list('S'=S, 'I'=I, 'Z'=Z, 'R'=R, 't'=t))
  return(c(daily_S,daily_Z))
}

