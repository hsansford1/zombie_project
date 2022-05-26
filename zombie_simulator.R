
simulate_zombies <- function(delta, zeta, beta, alpha, time_horizon, dt){
  
  # delta: background (non-zombie related) death rate
  # zeta: zombie resurrection rate
  # beta: zombie transmission rate
  # alpha: zombie destruction rate
  # time_horizon: length of time to run dynamics
  # dt: time-step for numerical solutions
  
  N <- 1000       # total population
  n <- time_horizon/dt
  s <- rep(0, n+1)
  z <- rep(0, n+1)
  r <- rep(0, n+1)
  
  # Initial conditions
  s[1] <- N
  z[1] <- 0
  r[1] <- 0
  t <- seq(0, time_horizon, dt)
  
  # Define model dynamics (ODEs) and solve by Euler's method
  for (i in 1:n){
    s[i+1] <- s[i] + dt*(-beta*s[i]*z[i]) # assuming birth rate = background death rate
    z[i+1] <- z[i] + dt*((beta - alpha)*s[i]*z[i] + zeta*r[i])
    r[i+1] <- r[i] + dt*(alpha*s[i]*z[i] + delta*s[i] - zeta*r[i])
    if ((s[i] < 0) | (s[i] > N)){
      break
    }
    if ((z[i] > N) | (z[i] < 0)){
      break
    }
    if ((r[i] < 0) | (r[i] > N)){
      break
    }
  }
  
  plot(t, s, 'blue')
  lines(t, z, 'red')
  legend(c('Susceptibles', 'Zombies'))
}