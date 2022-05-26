
simulate_zombies <- function(delta, zeta, beta, alpha, pi, time_horizon, dt){
  
  # delta: background (non-zombie related) death rate
  # zeta: zombie resurrection rate
  # beta: zombie transmission rate
  # alpha: zombie destruction rate
  # pi: birth rate
  # time_horizon: length of time to run dynamics
  # dt: time-step for numerical solutions
  
  N <- 1000             # total population
  n <- time_horizon/dt
  c <- 0.1              # human/zombie contact probability
  
  s <- rep(0, n+1)
  z <- rep(0, n+1)
  r <- rep(0, n+1)
  
  # Initial conditions
  z[1] <- 10
  s[1] <- N - z[1]
  r[1] <- 0
  t <- seq(0, time_horizon, dt)
  
  # Define model dynamics (ODEs) and solve by Euler's method
  for (i in 1:n){
    
    s[i+1] <- s[i] + dt*(pi*s[i] - beta*s[i]*z[i]*c - delta*s[i])
    z[i+1] <- z[i] + dt*((beta - alpha)*s[i]*z[i]*c + zeta*r[i])
    r[i+1] <- r[i] + dt*(alpha*s[i]*z[i]*c + delta*s[i] - zeta*r[i])
    
    if (s[i] < 0){
      s[i:n+1] <- 0
      for (j in i:n){
        z[j+1] <- z[j] + dt*(zeta*r[j])
        r[j+1] <- r[j] + dt*(-zeta*r[j])
        if (r[j] < 0){
          r[j:n+1] <- 0
          z[j:n+1] <- z[j]
        }
      }
      break
    }
    
    if (z[i] < 0){
      z[i+1] <- dt(zeta*r[i])
      s[i+1] <- s[i] + dt((pi - delta)*s[i])
      r[i+1] <- r[i] + dt(delta*s[i] - zeta*r[i])
    }
    if (r[i] < 0){
      z[i+1] <- z[i] + dt*((beta - alpha)*s[i]*z[i]*c)
      r[i+1] <- dt*(alpha*s[i]*z[i]*c + delta*s[i])
      break
    }
  }
  
  return(list('s'=s, 'z'=z, 'r'=r, 't'=t))
}

epidemic <- simulate_zombies(0.01, 0.05, 0.3, 0.4, 0.02, 50, 0.1)
plot(epidemic$t, epidemic$s, ylim=c(0,1200), type='l', col='blue')
lines(epidemic$t, epidemic$z, col='red')
