#install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline, tidyverse)

set.seed(1983)

setwd("/work/NielsAalundKrogsgaard#7447/Exam/decmak2023/parameter_recovery/hierarchical_group")

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#load control data
data <- read_csv("../../data/proctor_et_al_2014_combined.csv")

#----------prepare data for jags models - want trial x subject arrays for choice, gain, and loss ----
# identify and count unique subject IDs
# concat the two sessions (40 to 80 trials)
data$trial_new <- ifelse(data$session == 2, data$trial + 40, data$trial)


# filter data
human_data <- data %>% 
  filter(species == "human") %>% 
  select(id, trial_new, choice, payoff)
# identify and count unique subject IDs
subIDs_human <- unique(human_data$id)
nsubs_human <- length(subIDs_human)

# filter data
monkey_data <- data %>% 
  filter(species == "chimp") %>% 
  filter(!is.na(payoff)) %>% 
  select(id, trial_new, choice, payoff)

# identify and count unique subject IDs
subIDs_monkey <- unique(monkey_data$id)
nsubs_monkey <- length(subIDs_monkey)

ntrials_max <- 80

# SCALE THE REWARDS DOWN SO WE DON'T NEED SUCH WIDE PRIORS
# all choices (x) and outcomes (X)
x_raw_monkey <- monkey_data$choice
X_raw_monkey <- monkey_data$payoff #note the sign!!!!!

x_raw_human <- human_data$choice
X_raw_human <- human_data$payoff #note the sign!!!!!

#--- assign choices and outcomes in trial x sub matrix

#different number of trials across subjects. We'll need to fix this by padding arrays of < 100
#this is just so we can make the array
#then we'll also need to record number of valid trials for each sub, 
#then run the JAGS model on only valid trials

# empty arrays to fill
ntrials_monkey <- array(0,c(nsubs_monkey))
x_monkey <- array(0,c(nsubs_monkey,ntrials_max))
X_monkey <- array(0,c(nsubs_monkey,ntrials_max))

ntrials_human <- array(0,c(nsubs_human))
x_human <- array(0,c(nsubs_human,ntrials_max))
X_human <- array(0,c(nsubs_human,ntrials_max))

# make control data matrices
for (s in 1:ntrials_monkey) {
  
  #record n trials for subject s
  ntrials_monkey[s] <- length(x_raw_monkey[monkey_data$id==subIDs_monkey[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub_monkey <- x_raw_monkey[monkey_data$id==subIDs_monkey[s]] 
  length(x_sub_monkey) <- ntrials_max
  
  X_sub_monkey <- X_raw_monkey[monkey_data$id==subIDs_monkey[s]] 
  length(X_sub_monkey) <- ntrials_max
  
  # assign arrays
  x_monkey[s,] <- x_sub_monkey
  X_monkey[s,] <- X_sub_monkey
  
}

# make humanoid data matrices
for (s in 1:nsubs_human) {
  
  #record n trials for subject s
  ntrials_human[s] <- length(x_raw_human[human_data$id==subIDs_human[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub_human <- x_raw_human[human_data$id==subIDs_human[s]] 
  length(x_sub_human) <- ntrials_max
  
  X_sub_human <- X_raw_human[human_data$id==subIDs_human[s]] 
  length(X_sub_human) <- ntrials_max
  
  # assign arrays
  x_human[s,] <- x_sub_human
  X_human[s,] <- X_sub_human
  
}



# set up jags and run jags model on one subject
data <- list("x_monkey","X_monkey","ntrials_monkey","nsubs_monkey",
             "x_human","X_human","ntrials_human","nsubs_human") 
params<-c("alpha_a_rew","alpha_K","alpha_omega_f","alpha_omega_p",
          "mu_a_rew","mu_K","mu_omega_f","mu_omega_p")#,
#"a_rew_monkey","a_pun_monkey","K_monkey","omega_f_monkey","omega_p_monkey",
#"a_rew_human","a_pun_human","K_human","omega_f_human","omega_p_human")
start_time = Sys.time()

samples <- jags.parallel(data, inits=NULL, params,
                model.file ="ORL_compare.txt",
                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)

end_time = Sys.time()
end_time - start_time

samples

BF_effect <- NULL
BF_null <- NULL

# savage dickey plot
plot(density(rnorm(12000,0,1/sqrt(1))),ylim=c(0,1),main=" ")
lines(density(samples$BUGSoutput$sims.list$alpha_a_rew),col="red")

fit.posterior <- logspline(samples$BUGSoutput$sims.list$alpha_a_rew)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF_effect$alpha_a_rew <- null.prior/null.posterior
BF_null$alpha_a_rew <- null.posterior/null.prior

# savage dickey plot
plot(density(rnorm(12000,0,1/sqrt(1))),ylim=c(0,1),main=" ")
lines(density(samples$BUGSoutput$sims.list$alpha_omega_f),col="red")

fit.posterior <- logspline(samples$BUGSoutput$sims.list$alpha_omega_f)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF_effect$alpha_omega_f <- null.prior/null.posterior
BF_null$alpha_omega_f <- null.posterior/null.prior



par(mfrow=c(2,2))
plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="a_rew")
lines(density(samples$BUGSoutput$sims.list$alpha_a_rew),col="red")

plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="omega_f")
lines(density(samples$BUGSoutput$sims.list$alpha_omega_f),col="red")

plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="omega_p")
lines(density(samples$BUGSoutput$sims.list$alpha_omega_p),col="red")

plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="K")
lines(density(samples$BUGSoutput$sims.list$alpha_K),col="red")















############### Look at cumulative balance to justify learning ##############

# calculate cumulative value
xcum_monkey <- array(0,c(nsubs_monkey,80))
for (i in 1:nsubs_monkey) {
  xcum_monkey[i,] <- cumsum(X_monkey[i,]) 
}

xcum_human <- array(0,c(nsubs_human,80))
for (i in 1:nsubs_human) {
  xcum_human[i,] <- cumsum(X_human[i,]) 
}

plot(colMeans(xcum_monkey),ylim=c(-5,100),type='l',lwd=2,  
     xlab = "Trial", ylab = "Cumulative Balance") # should be -1500, 800 if not scaled
lines(colMeans(xcum_human),col="red",lwd=2)
legend(x='topright', legend=c('healthy controls', 'heroine', 'ampetamine'))

############### compare final means ##############

# set up jags and run jags model on one subject
data <- list("xcum_monkey","nsubs_monkey",
             "xcum_human","nsubs_human") 
params<-c("mu","alpha","Smu_monkey","Smu_human")
temp_samples <- jags.parallel(data, inits=NULL, params,
                model.file ="balance_compare.txt",
                n.chains=3, n.iter=15000, n.burnin=1000, n.thin=3, n.cluster=4)


# savage dickey plot
plot(density(rnorm(12000,0,1/sqrt(1))),main=" ")
lines(density(temp_samples$BUGSoutput$sims.list$alpha),col="red")

fit.posterior <- logspline(temp_samples$BUGSoutput$sims.list$alpha)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF <- null.prior/null.posterior
BF

# savage dickey plot
plot(density(rnorm(12000,0,1/sqrt(1))),main="mu")
lines(density(temp_samples$BUGSoutput$sims.list$mu),col="red")

fit.posterior <- logspline(temp_samples$BUGSoutput$sims.list$mu)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF <- null.prior/null.posterior
BF


