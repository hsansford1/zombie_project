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



## Plot accepted sampled trajectories
plot_trajectories <- function(true_epidemic, t,  ABC_output, title){
  par(mfrow=c(1,1))
  plot(t, floor(true_epidemic[1:501]), ylim=c(0,500), type='l', col=colours[2], lwd=3,
       ylab='Population Values', xlab='Time', main=title)
  lines(t, floor(true_epidemic[502:1002]), col=colours[3], lwd=3)
  for (i in 1:225){
    rej_epidemic <- simulate_zombies(ABC_output$param[i,])
    lines(t, floor(rej_epidemic[1:501]), col=alpha(colours[2], 0.4), lwd=1)
    lines(t, floor(rej_epidemic[502:1002]), col=alpha(colours[3], 0.4), lwd=1)
  }
  legend('topright', legend = c('True Susceptibles', 'True Zombies', 
                                'Simulated Susceptibles', 'Simulated Zombies'), 
         col = c(colours[2], colours[3], alpha(colours[2], 0.4), alpha(colours[3], 0.4)),
         lty=1, lwd=c(3,3,1,1))
}


## Plot confidence bands
plot_conf_bands <- function(true_epidemic, t, abc_output, title){
  plot(t, floor(true_epidemic[1:501]), ylim=c(0,500), type='l', col=colours[2], lwd=3,
       ylab='Population Values', xlab='Time', main=title)
  lines(t, floor(true_epidemic[502:1002]), col=colours[3], lwd=3)
  
  lines(t, apply(abc_output$ss[,1:501], 2, mean), col = colours[1], lwd=3, lty=2)
  polygon(c(t,rev(t)), c(apply(abc_output$ss[,1:501], 2, quantile, prob=0.05), 
                         rev(apply(abc_output$ss[,1:501], 2, quantile, prob= 0.95))), 
          col=alpha(colours[1], 0.2), border = NA)
  lines(t, apply(abc_output$ss[,502:1002], 2, mean), col = colours[4], lwd=3, lty=2)
  polygon(c(t,rev(t)), c(apply(abc_output$ss[,502:1002], 2, quantile, prob=0.05), 
                         rev(apply(abc_output$ss[,502:1002], 2, quantile, prob= 0.95))), 
          col=alpha(colours[4], 0.2), border = NA)
  legend('topright', legend = c('True Susceptibles', 'True Zombies', 
                                'Mean Simulated Susceptibles', 'Mean Simulated Zombies'),
         col = colours[c(2,3,1,4)], lty=c(1,1,2,2), lwd=3)
}


## Plot mean +/- one standard deviation
plot_conf_bands2 <- function(true_epidemic, t, abc_output, title){
  plot(t, floor(true_epidemic[1:501]), ylim=c(0,500), type='l', col=colours[2], lwd=3,
       ylab='Population Values', xlab='Time', main=title)
  lines(t, floor(true_epidemic[502:1002]), col=colours[3], lwd=3)
  
  lines(t, apply(abc_output$ss[,1:501], 2, mean), col = colours[1], lwd=3, lty=2)
  polygon(c(t,rev(t)), c(apply(abc_output$ss[,1:501], 2, mean) - apply(abc_output$ss[,1:501], 2, sd), 
                         rev(apply(abc_output$ss[,1:501], 2, mean) + apply(abc_output$ss[,1:501], 2, sd))), 
          col=alpha(colours[1], 0.2), border = NA)
  lines(t, apply(abc_output$ss[,502:1002], 2, mean), col = colours[4], lwd=3, lty=2)
  polygon(c(t,rev(t)), c(apply(abc_output$ss[,502:1002], 2, mean) - apply(abc_output$ss[,502:1002], 2, sd), 
                         rev(apply(abc_output$ss[,502:1002], 2, mean) + apply(abc_output$ss[,502:1002], 2, sd))), 
          col=alpha(colours[4], 0.2), border = NA)
  legend('topright', legend = c('True Susceptibles', 'True Zombies', 
                                'Mean Simulated Susceptibles', 'Mean Simulated Zombies'),
         col = colours[c(2,3,1,4)], lty=c(1,1,2,2), lwd=3)
}

