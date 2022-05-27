
simulate_zombies <- function(delta, zeta, beta, alpha, pi, rho, time_horizon, dt){
  
  # delta: background (non-zombie related) death rate
  # zeta: zombie resurrection rate
  # beta: zombie transmission rate
  # alpha: zombie destruction rate
  # pi: birth rate
  # rho: latent infection rate
  # time_horizon: length of time to run dynamics
  # dt: time-step for numerical solutions
  
  N <- 500             # total population
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

  
  return(list('S'=S, 'I'=I, 'Z'=Z, 'R'=R, 't'=t))
}

epidemic <- simulate_zombies(delta=0.0001, zeta=0.0001, beta=0.0095, alpha=0.0005, 
                             pi=0.0001, rho=0.05,  time_horizon=25, dt=0.05)

plot(epidemic$t, floor(epidemic$S), ylim=c(0,500), type='l', col='blue', lwd=2)
lines(epidemic$t, floor(epidemic$Z), col='red', lwd=2)
legend('topright', legend=c('Susceptibles', 'Zombies'), col = c('blue', 'red'), lty=1, lwd=2)


  t <- 1:time_horizon
  s <- s[(1/dt)*(1:time_horizon)]
  z <- z[(1/dt)*(1:time_horizon)]
  r <- r[(1/dt)*(1:time_horizon)]