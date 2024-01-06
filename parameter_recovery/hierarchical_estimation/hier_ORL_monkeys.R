#install.packages("pacman")
pacman::p_load(R2jags, parallel, tidyverse)

set.seed(1983)

setwd("/work/NielsAalundKrogsgaard#7447/Exam/decmak2023/parameter_recovery/hierarchical_estimation")

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#----------getting the data
# data from this paper: https://www.frontiersin.org/articles/10.3389/fpsyg.2014.00849/full
# available here: https://figshare.com/articles/dataset/IGT_raw_data_Ahn_et_al_2014_Frontiers_in_Psychology/1101324



#load control data
monkey_data <- read_csv("../../data/proctor_et_al_2014_combined.csv")

#----------prepare data for jags models - want trial x subject arrays for choice and gain & loss ----
# concat the two sessions (40 to 80 trials)
monkey_data$trial_new <- ifelse(monkey_data$session == 2, monkey_data$trial + 40, monkey_data$trial)

# filter data
monkey_data <- monkey_data %>% 
  filter(species == "chimp") %>% 
  filter(!is.na(payoff)) %>% 
  select(id, trial_new, choice, payoff)
# identify and count unique subject IDs
subIDs <- unique(monkey_data$id)
nsubs <- length(subIDs)
ntrials_max <- 80

# all choices (x) and outcomes (X)
x_raw <- monkey_data$choice
X_raw <- monkey_data$payoff #note the sign!

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
  ntrials_all[s] <- length(x_raw[monkey_data$id==subIDs[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub <- x_raw[monkey_data$id==subIDs[s]] 
  length(x_sub) <- ntrials_max
  
  X_sub <- X_raw[monkey_data$id==subIDs[s]] 
  length(X_sub) <- ntrials_max
  
  # assign arrays
  x_all[s,] <- x_sub
  X_all[s,] <- X_sub
  
}

# Scaling the payoffs (cuz the learning parameter becomes less relevant for very large payoffs/losses)
# X_all <- X_all/100


#----------testing our data curation by running JAGS on one subject

# # Now we'll fit one subject just to make sure everything works
# 
# x <- x_all[1,]
# X <- X_all[1,]
# 
# ntrials <- ntrials_all[1]
# 
# # set up jags and run jags model on one subject
# data <- list("x","X","ntrials") 
# params<-c("a_rew","K","omega_f","omega_p","p")
# samples <- jags.parallel(data, inits=NULL, params,
#                 model.file ="hier_ORL_wo_theta.txt",
#                 n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)
# 
# # let's look at the posteriors for the parameters
# par(mfrow=c(3,2))
# plot(density(samples$BUGSoutput$sims.list$a_rew))
# plot(density(samples$BUGSoutput$sims.list$a_pun))
# plot(density(samples$BUGSoutput$sims.list$theta))
# plot(density(samples$BUGSoutput$sims.list$K))
# plot(density(samples$BUGSoutput$sims.list$omega_f))
# plot(density(samples$BUGSoutput$sims.list$omega_p))

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
params<-c("mu_a_rew","mu_K","mu_omega_f","mu_omega_p") 

start_time = Sys.time()
samples <- jags.parallel(data, inits=NULL, params,
                         model.file ="hier_ORL_wo_theta.txt",
                         n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)
end_time = Sys.time()
end_time - start_time

par(mfrow=c(2,2))
plot(density(samples$BUGSoutput$sims.list$mu_a_rew))
#plot(density(samples$BUGSoutput$sims.list$mu_a_pun))
plot(density(samples$BUGSoutput$sims.list$mu_K))
#plot(density(samples$BUGSoutput$sims.list$mu_theta))
plot(density(samples$BUGSoutput$sims.list$mu_omega_f))
plot(density(samples$BUGSoutput$sims.list$mu_omega_p))


#----------Posterior predictive checks of descriptive accuracy

# let's see how the model goes for more than 1 subject. Let's run this on all subjects
pred_success <- array(nsubs)

start_time = Sys.time()

for (s in 1:nsubs) {
  
  x <- x_all[s, ]
  X <- X_all[s, ]
  
  ntrials <- ntrials_all[s]
  
  # set up jags and run jags model on one subject
  data <- list("x","X","ntrials") 
  params<-c("a_rew","K","omega_f","omega_p","p")
  temp_samples <- jags.parallel(data, inits=NULL, params,
                                model.file ="ORL_wo_theta.txt",
                                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=3)
  
  p_post <- temp_samples$BUGSoutput$sims.list$p
  
  x_predict <- array(ntrials)
  
  for (t in 1:ntrials) {
    p_predict <- c(
      MPD(p_post[,t,1]),
      MPD(p_post[,t,2])
    )
    
    x_predict[t] <- which.max(p_predict)
    
  }
  
  pred_success[s] <- sum(x_predict==x[1:ntrials]) # only comparing with trials for which we have choices
  print(s)
  
}

end_time = Sys.time()
end_time - start_time

pred_success_adjust <- pred_success/ntrials_all

avg_pred <- mean(pred_success_adjust)

# plotting code courtesy of Mia
pred_df <- data.frame(pred_success_adjust)
pred_df$sub <- 1:length(pred_success_adjust) # rownames(pred_df) # creating a subject index
pred_df$avg <- mean(pred_df$pred_success)
pred_df$std <- sd(pred_df$pred_success)
pred_df$chance <- .5
ggplot(pred_df, aes(sub, pred_success_adjust)) +
  geom_point() +
  geom_line(aes(y=chance), linetype="dashed", color = "black") +
  geom_ribbon(aes(xmin = -Inf, xmax = Inf, ymin = avg - std, ymax = avg + std), fill = "pink", alpha = 0.6) + 
  geom_line(aes(y=avg)) + 
  ylim(0,1)
ggsave("output/descriptive_adequacy_monkeys.png", width = 10, height = 5)

save.image(file="output/hier_orl_monkeys_data.RData") 