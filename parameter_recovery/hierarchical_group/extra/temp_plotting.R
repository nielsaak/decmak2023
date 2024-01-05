dev.off()
ylim=c(0,1)

par(mfrow=c(3,2), xpd=NA)
plot(density(prior_samples$BUGSoutput$sims.list$alpha_a_rew),ylim=ylim,main="a_rew", col="blue")
lines(density(samples$BUGSoutput$sims.list$alpha_a_rew),col="red")
lines(density(rnorm(12000,0,1/sqrt(1))), col="grey")
lines(density(rnorm(12000,0,1/sqrt(.1))), col="grey", lty=2)

plot(density(prior_samples$BUGSoutput$sims.list$alpha_a_pun),ylim=ylim,main="a_pun", col="blue")
lines(density(samples$BUGSoutput$sims.list$alpha_a_pun),col="red")
lines(density(rnorm(12000,0,1/sqrt(1))), col="grey")
lines(density(rnorm(12000,0,1/sqrt(.1))), col="grey", lty=2)

plot(density(prior_samples$BUGSoutput$sims.list$alpha_omega_f),ylim=ylim,main="omega_f", col="blue")
lines(density(samples$BUGSoutput$sims.list$alpha_omega_f),col="red")
lines(density(rnorm(12000,0,1/sqrt(1))), col="grey")
lines(density(rnorm(12000,0,1/sqrt(.1))), col="grey", lty=2)

plot(density(prior_samples$BUGSoutput$sims.list$alpha_omega_p),ylim=ylim,main="omega_p", col="blue")
lines(density(samples$BUGSoutput$sims.list$alpha_omega_p),col="red")
lines(density(rnorm(12000,0,1/sqrt(1))), col="grey")
lines(density(rnorm(12000,0,1/sqrt(.1))), col="grey", lty=2)

plot(density(prior_samples$BUGSoutput$sims.list$alpha_K),ylim=ylim,main="K", col="blue")
lines(density(samples$BUGSoutput$sims.list$alpha_K),col="red")
lines(density(rnorm(12000,0,1/sqrt(1))), col="grey")
lines(density(rnorm(12000,0,1/sqrt(.1))), col="grey", lty=2)

