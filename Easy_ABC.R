source('zombie_simulator.R')

epidemic <- simulate_zombies(c(0.0001, 0.0001, 0.0095, 0.0005, 0.0001, 0.05))
t <- seq(0, 25, 0.05)
plot(t, floor(epidemic[1:501]), ylim=c(0,500), type='l', col='blue', lwd=2)
lines(t, floor(epidemic[502:1002]), col='red', lwd=2)
legend('topright', legend=c('Susceptibles', 'Zombies'), col = c('blue', 'red'), lty=1, lwd=2)

library(EasyABC)

# define prior distribution
priors <- list(c('unif', 0, 0.01), c('unif', 0, 0.01), c('unif', 0, 0.1),
               c('unif', 0, 0.01), c('unif', 0, 0.01), c('unif', 0, 0.1))

n <- 5000 # number of simulations to be performed
p <- 0.01  #proportion of simulations to be retained

ABC_rej <- ABC_rejection(model = simulate_zombies, prior=priors, nb_simul=n, 
                         summary_stat_target = epidemic, tol=p, progress_bar=TRUE)

ABC_rej$computime
ABC_rej$param
par(mfrow=c(2,3))
hist(ABC_rej$param[,1], breaks=10)
hist(ABC_rej$param[,2], breaks=10)
hist(ABC_rej$param[,3], breaks=10)
hist(ABC_rej$param[,4], breaks=10)
hist(ABC_rej$param[,5], breaks=10)
hist(ABC_rej$param[,6], breaks=10)

