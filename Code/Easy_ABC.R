source('zombie_simulator.R')
source('plotting_functions.R')
library(jcolors)
colours <- jcolors('pal3')
library(scales)
library(EasyABC)
library(abc)

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
                         summary_stat_target = true_epidemic, tol=p, progress_bar=TRUE)

# algorithm run-time
ABC_rej$computime

# plot posterior distributions
histogram_plot(ABC_rej$param, 'Posterior distributions from Rejection ABC')

# Plot accepted sampled trajectories
plot_trajectories(true_epidemic, t, ABC_rej$param, title='Accepted simulations from rejection ABC algorithm')

# plot mean and 90% confidence bands
plot_conf_bands(true_epidemic, t, ABC_rej, title='Mean and confidence band from rejection ABC')

# plot mean +/- one standard deviation
plot_sd_bands(true_epidemic, t, abc_rej, title='Mean and standard deviation from rejection ABC')

# use features of `abc` package
parnames <- c('delta', 'zeta', 'beta', 'alpha', 'pi', 'rho')
colnames(ABC_rej$param) <- parnames
abc_rej <- abc(true_epidemic, ABC_rej$param, ABC_rej$stats, tol=1, method="rejection")

# epsilon values for accepted simulations
abc_rej$dist

#------------------------------------------------------------------------------
# SEQUENTIAL ABC (ABC-SMC)
#--------------------------------------------------------------------------

#sequence of tolerance levels:
tolerance <- c(50, 20, 5, 2)
#number of simulations to obtain below the tolerance level at each iteration:
n <- 100

ABC_Beaumont <- ABC_sequential(method="Beaumont", model=simulate_zombies, prior=priors,
                              nb_simul=n, summary_stat_target=true_epidemic,
                              tolerance_tab=tolerance, progress_bar = TRUE, inside_prior = TRUE)

ABC_Beaumont$nsim # The number of model simulations performed.
ABC_Beaumont$computime # The computing time to perform the simulations
ABC_Beaumont$epsilon

# Take the 25 best simulations
dist <- sqrt(apply((ABC_Beaumont$stats- true_epidemic)^2, 1, sum))
idx <-sort(dist, index.return=TRUE)$ix
acc_params_beaumont <- ABC_Beaumont$param[idx[1:25],]

# plot posterior distributions
histogram_plot(ABC_Beaumont$param, 'Posterior distributions from PMC-ABC')

# Plot accepted sampled trajectories
plot_trajectories(true_epidemic, t, acc_params_beaumont, title='Accepted simulations from PMC-ABC algorithm')

# plot mean and 90% confidence bands
plot_conf_bands(true_epidemic, t, ABC_Beaumont, title='Mean and confidence band from PMC-ABC')

# plot mean +/- one standard deviation
plot_sd_bands(true_epidemic, t, ABC_Beaumont, title='Mean and standard deviation from PMC-ABC')
