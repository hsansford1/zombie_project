source('zombie_simulator.R')
source('plotting_functions.R')
library(jcolors)
colours <- jcolors('pal3')
library(scales)
library(EasyABC)


## Simulate 'true' epidemic data

true_params <- c(0.0001, 0.0001, 0.0095, 0.0005, 0.0001, 0.05)
true_epidemic <- simulate_zombies(true_params)
t <- seq(0, 25, 0.05)
plot(t, floor(true_epidemic[1:501]), ylim=c(0,500), type='l', col=colours[2], lwd=2)
lines(t, floor(true_epidemic[502:1002]), col=colours[3], lwd=2)
legend('topright', legend=c('Susceptibles', 'Zombies'), col = colours[2:3], lty=1, lwd=2)


##--------------------------------------------------------------
## REJECTION ABC
##-----------------------------------------------------------------

# define prior distribution
priors <- list(c('unif', 0, 0.001), c('unif', 0, 0.001), c('unif', 0, 0.1),
               c('unif', 0, 0.001), c('unif', 0, 0.001), c('unif', 0, 0.1))

n <- 5000 # number of simulations to be performed
p <- 0.005  #proportion of simulations to be retained

# run ABC rejection algorithm
ABC_rej <- ABC_rejection(model = simulate_zombies, prior=priors, nb_simul=n, 
                         summary_stat_target = epidemic, tol=p, progress_bar=TRUE)

# algorithm run-time
ABC_rej$computime

# plot posterior distributions
histogram_plot(ABC_rej$param, 'Posterior distributions from Rejection ABC')

# Plot accepted sampled trajectories
plot_trajectories(true_epidemic, ABC_rej, title='Accepted simulations from rejection ABC algorithm')


parnames <- c('delta', 'zeta', 'beta', 'alpha', 'pi', 'rho')
colnames(ABC_rej$param) <- parnames
abc_rej <- abc(true_epidemic, ABC_rej$param, ABC_rej$stats, tol=1, method="rejection")

# plot confidence bands
plot(t, floor(epidemic[1:501]), ylim=c(0,500), type='l', col=colours[2], lwd=3,
     ylab='Population Values', xlab='Time', main='Mean and confidence band from rejection ABC')
lines(t, floor(epidemic[502:1002]), col=colours[3], lwd=3)

lines(t, apply(abc_rej$ss[,1:501], 2, mean), col = colours[1], lwd=3, lty=2)
polygon(c(t,rev(t)), c(apply(abc_rej$ss[,1:501], 2, quantile, prob=0.05), 
                       rev(apply(abc_rej$ss[,1:501], 2, quantile, prob= 0.95))), 
                        col=alpha(colours[1], 0.25), border = NA)
lines(t, apply(abc_rej$ss[,502:1002], 2, mean), col = colours[4], lwd=3, lty=2)
polygon(c(t,rev(t)), c(apply(abc_rej$ss[,502:1002], 2, quantile, prob=0.05), 
                       rev(apply(abc_rej$ss[,502:1002], 2, quantile, prob= 0.95))), 
        col=alpha(colours[4], 0.25), border = NA)

lines(t, floor(epidemic[1:501]), col=colours[2], lwd=3)
lines(t, floor(epidemic[502:1002]), col=colours[3], lwd=3)


