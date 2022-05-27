source('zombie_simulator.R')
library(jcolors)
colours <- jcolors('pal3')
library(scales)
library(EasyABC)

## Simulate 'true' epidemic data

true_params <- c(0.0001, 0.0001, 0.0095, 0.0005, 0.0001, 0.05)
epidemic <- simulate_zombies(true_params)
t <- seq(0, 25, 0.05)
plot(t, floor(epidemic[1:501]), ylim=c(0,500), type='l', col='blue', lwd=2)
lines(t, floor(epidemic[502:1002]), col='red', lwd=2)
legend('topright', legend=c('Susceptibles', 'Zombies'), col = c('blue', 'red'), lty=1, lwd=2)

## function for plotting posterior histograms
histogram_plot <- function(param, title){
  
  par(mfrow=c(2,3))
  hist(param[,1], breaks=10, col = colours[2], border=colours[2], 
     main='', xlab=expression(delta))
  abline(v=true_params[1], col=colours[3], lwd=2)
  hist(param[,2], breaks=10, col = colours[2], border=colours[2],
     main='', xlab=expression(zeta))
  abline(v=true_params[2], col=colours[3], lwd=2)
  hist(param[,3], breaks=10, col = colours[2], border=colours[2],
     main='', xlab=expression(beta))
  abline(v=true_params[3], col=colours[3], lwd=2)
  hist(param[,4], breaks=10, col = colours[2], border=colours[2],
     main='', xlab=expression(alpha))
  abline(v=true_params[4], col=colours[3], lwd=2)
  hist(param[,5], breaks=10, col = colours[2], border=colours[2],
     main='', xlab=expression(Pi))
  abline(v=true_params[5], col=colours[3], lwd=2)
  hist(param[,6], breaks=10, col = colours[2], border=colours[2],
     main='', xlab=expression(rho))
  abline(v=true_params[6], col=colours[3], lwd=2)
  mtext(title, line=-2, outer=TRUE)
}

##--------------------------------------------------------------
## REJECTION ABC
##-----------------------------------------------------------------

# define prior distribution
priors <- list(c('unif', 0, 0.001), c('unif', 0, 0.001), c('unif', 0, 0.1),
               c('unif', 0, 0.001), c('unif', 0, 0.001), c('unif', 0, 0.1))

n <- 5000 # number of simulations to be performed
p <- 0.005  #proportion of simulations to be retained

ABC_rej <- ABC_rejection(model = simulate_zombies, prior=priors, nb_simul=n, 
                         summary_stat_target = epidemic, tol=p, progress_bar=TRUE)

ABC_rej$computime
ABC_rej$param



histogram_plot(ABC_rej$param, 'Posterior distributions from Rejection ABC')

par(mfrow=c(1,1))
plot(t, floor(epidemic[1:501]), ylim=c(0,500), type='l', col=colours[2], lwd=3,
      ylab='Population Values', xlab='Time', main='Accepted simulations from rejection ABC algorithm')
lines(t, floor(epidemic[502:1002]), col=colours[3], lwd=3)
for (i in 1:20){
  rej_epidemic <- simulate_zombies(ABC_rej$param[i,])
  lines(t, floor(rej_epidemic[1:501]), col=alpha(colours[2], 0.4), lwd=1)
  lines(t, floor(rej_epidemic[502:1002]), col=alpha(colours[3], 0.4), lwd=1)
}
legend('topright', legend = c('True Susceptibles', 'True Zombies', 
                              'Simulated Susceptibles', 'Simulated Zombies'), 
       col = c(colours[2], colours[3], alpha(colours[2], 0.4), alpha(colours[3], 0.4)),
       lty=1, lwd=c(3,3,1,1))


plot(t, floor(epidemic[1:501]), ylim=c(0,500), type='l', col=colours[2], lwd=3,
     ylab='Population Values', xlab='Time', main='Mean and confidence band from rejection ABC')
lines(t, floor(epidemic[502:1002]), col=colours[3], lwd=3)

lines(t, apply(rej$ss, 2, mean), col = colours[1], lwd=3)
