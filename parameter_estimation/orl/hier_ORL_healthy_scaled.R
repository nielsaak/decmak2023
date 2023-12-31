install.packages("pacman")
pacman::p_load(R2jags, parallel)

set.seed(1983)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#----------getting the data
# data from this paper: https://www.frontiersin.org/articles/10.3389/fpsyg.2014.00849/full
# available here: https://figshare.com/articles/dataset/IGT_raw_data_Ahn_et_al_2014_Frontiers_in_Psychology/1101324

#load control data
ctr_data <- read.table("data/IGTdata_healthy_control.txt",header=TRUE)

#----------prepare data for jags models - want trial x subject arrays for choice and gain & loss ----
# identify and count unique subject IDs
subIDs <- unique(ctr_data$subjID)
nsubs <- length(subIDs)
ntrials_max <- 100

# all choices (x) and outcomes (X)
x_raw <- ctr_data$deck
X_raw <- ctr_data$gain + ctr_data$loss #note the sign!

#--- assign choices and outcomes in trial x sub matrix

#different number of trials across subjects. We'll need to fix this by padding arrays of < 100
#this is just so we can make the array
#then we'll also need to record number of valid trials for each sub, 
#then run the JAGS model on only valid trials

# empty arrays to fill
ntrials_all <- array(0,c(nsubs))
x_all <- array(0,c(nsubs,ntrials_max))
X_all <- array(0,c(nsubs,ntrials_max))

for (s in 1:nsubs) {
  
  #record n trials for subject s
  ntrials_all[s] <- length(x_raw[ctr_data$subjID==subIDs[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub <- x_raw[ctr_data$subjID==subIDs[s]] 
  length(x_sub) <- ntrials_max
  
  X_sub <- X_raw[ctr_data$subjID==subIDs[s]] 
  length(X_sub) <- ntrials_max
  
  # assign arrays
  x_all[s,] <- x_sub
  X_all[s,] <- X_sub
  
}

# Scaling the payoffs (cuz the learning parameter becomes less relevant for very large payoffs/losses)
X_all <- X_all/100


#----------testing our data curation by running JAGS on one subject

# Now we'll fit one subject just to make sure everything works

x <- x_all[1,]
X <- X_all[1,]

ntrials <- ntrials_all[1]

# set up jags and run jags model on one subject
data <- list("x","X","ntrials") 
params<-c("a_rew","a_pun","K","theta","omega_f","omega_p","p")
samples <- jags.parallel(data, inits=NULL, params,
                model.file ="ORL.txt",
                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)

# let's look at the posteriors for the parameters
par(mfrow=c(3,2))
plot(density(samples$BUGSoutput$sims.list$a_rew))
plot(density(samples$BUGSoutput$sims.list$a_pun))
plot(density(samples$BUGSoutput$sims.list$theta))
plot(density(samples$BUGSoutput$sims.list$K))
plot(density(samples$BUGSoutput$sims.list$omega_f))
plot(density(samples$BUGSoutput$sims.list$omega_p))

# Question: how would you expect the data to look on the basis of these posteriors?


###########################################################
#---------- run the hierarchical model on controls --------
###########################################################

x <- x_all
X <- X_all

ntrials <- ntrials_all

# set up jags and run jags model
data <- list("x","X","ntrials","nsubs") 
# NB! we're not tracking theta cuz we're not modelling it in order reduce complexity a bit (hence, we're just setting it to 1 in "hier_ORL.txt")
params<-c("mu_a_rew","mu_a_pun","mu_K","mu_theta","mu_omega_f","mu_omega_p") 

start_time = Sys.time()
samples <- jags.parallel(data, inits=NULL, params,
                         model.file ="hier_ORL.txt",
                         n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)
end_time = Sys.time()
end_time - start_time

par(mfrow=c(3,2))
plot(density(samples$BUGSoutput$sims.list$mu_a_rew))
plot(density(samples$BUGSoutput$sims.list$mu_a_pun))
plot(density(samples$BUGSoutput$sims.list$mu_K))
plot(density(samples$BUGSoutput$sims.list$mu_theta))
plot(density(samples$BUGSoutput$sims.list$mu_omega_f))
plot(density(samples$BUGSoutput$sims.list$mu_omega_p))

# xlim scaled to Haines et al.
par(mfrow=c(3,2))
plot(density(samples$BUGSoutput$sims.list$mu_a_rew), xlim=c(0,0.4))
plot(density(samples$BUGSoutput$sims.list$mu_a_pun), xlim=c(0,0.08))
plot(density(samples$BUGSoutput$sims.list$mu_K), xlim=c(0,2))
plot(density(samples$BUGSoutput$sims.list$mu_theta), xlim=c(0,5))
plot(density(samples$BUGSoutput$sims.list$mu_omega_f), xlim=c(0,5))
plot(density(samples$BUGSoutput$sims.list$mu_omega_p), xlim=c(-4,2))

