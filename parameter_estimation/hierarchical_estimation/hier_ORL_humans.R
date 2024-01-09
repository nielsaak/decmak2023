#install.packages("pacman")
pacman::p_load(R2jags, parallel, tidyverse, ggpubr)

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
human_data <- read_csv("../../data/proctor_et_al_2014_combined.csv")

#----------prepare data for jags models - want trial x subject arrays for choice and gain & loss ----
# concat the two sessions (40 to 80 trials)
human_data$trial_new <- ifelse(human_data$session == 2, human_data$trial + 40, human_data$trial)

# filter data
human_data <- human_data %>% 
  filter(species == "human") %>% 
  filter(!is.na(payoff)) %>% 
  select(id, trial_new, choice, payoff)
# identify and count unique subject IDs
subIDs <- unique(human_data$id)
nsubs <- length(subIDs)
ntrials_max <- 40

# all choices (x) and outcomes (X)
x_raw <- human_data$choice
X_raw <- human_data$payoff #note the sign!

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
  ntrials_all[s] <- length(x_raw[human_data$id==subIDs[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub <- x_raw[human_data$id==subIDs[s]] 
  length(x_sub) <- ntrials_max
  
  X_sub <- X_raw[human_data$id==subIDs[s]] 
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
params<-c("mu_a_rew","mu_K","mu_omega_f","mu_omega_p", "lambda_a_rew", "lambda_K", "lambda_omega_f", "lambda_omega_p") 

start_time = Sys.time()
samples <- jags.parallel(data, inits=NULL, params,
                         model.file ="hier_ORL_wo_theta.txt",
                         n.chains=4, n.iter=10000, n.burnin=3000, n.thin=1, n.cluster=4)
end_time = Sys.time()
end_time - start_time

png(file="output/posterior_estimates_humans.png",
    width=1600, height=1000, res = 200, pointsize = 10)

par(mfrow=c(2,2))
plot(density(samples$BUGSoutput$sims.list$mu_a_rew), main = expression(""*mu[a[rew]]))
abline(v=0, lty = 2)
#plot(density(samples$BUGSoutput$sims.list$mu_a_pun))
plot(density(samples$BUGSoutput$sims.list$mu_K), main = expression(""*mu[K]))
abline(v=0, lty = 2)
#plot(density(samples$BUGSoutput$sims.list$mu_theta))
plot(density(samples$BUGSoutput$sims.list$mu_omega_f), main = expression(""*mu[omega[F]]))
abline(v=0, lty = 2)
plot(density(samples$BUGSoutput$sims.list$mu_omega_p), main = expression(""*mu[omega[P]]))
abline(v=0, lty = 2)

dev.off()

pacman::p_load(bayesplot)
samples_mcmc <- as.mcmc(samples)
mcmc_trace(samples_mcmc, pars = c("mu_a_rew", "mu_K", "mu_omega_p", "mu_omega_f")) +
  ggplot2::scale_color_discrete() +
  theme_minimal()
ggsave("output/traceplot_humans_mu.png", width = 10, height = 5)

samples_mcmc <- as.mcmc(samples)
mcmc_trace(samples_mcmc, pars = c("lambda_a_rew", "lambda_K", "lambda_omega_f", "lambda_omega_p")) +
  ggplot2::scale_color_discrete() +
  theme_minimal()
ggsave("output/traceplot_humans_lambda.png", width = 10, height = 5)

mpd_a_rew <- MPD(samples$BUGSoutput$sims.list$mu_a_rew)
mpd_K <- MPD(samples$BUGSoutput$sims.list$mu_K)
mpd_omega_f <- MPD(samples$BUGSoutput$sims.list$mu_omega_f)
mpd_omega_p <- MPD(samples$BUGSoutput$sims.list$mu_omega_p)








#----------Posterior predictive checks of descriptive accuracy

# let's see how the model goes for more than 1 subject. Let's run this on all subjects
pred_success <- array(nsubs)
x_predict_saved <- array(NA, dim = c(40, 10))

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
                                n.chains=4, n.iter=10000, n.burnin=3000, n.thin=1, n.cluster=4)
  
  p_post <- temp_samples$BUGSoutput$sims.list$p
  
  x_predict <- array(ntrials)
  
  for (t in 1:ntrials) {
    p_predict <- c(
      MPD(p_post[,t,1]),
      MPD(p_post[,t,2])
    )
    
    x_predict[t] <- which.max(p_predict)
    
  }
  
  x_predict_saved[,s] <- x_predict
  pred_success[s] <- sum(x_predict==x[1:ntrials]) # only comparing with trials for which we have choices
  print(s)
  
}

end_time = Sys.time()
end_time - start_time

pred_success_adjust <- pred_success/ntrials_all

avg_pred <- mean(pred_success_adjust)
sd_pred <- sd(pred_success_adjust)

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
  ylim(0,1) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Subject") + #Setting x label
  ylab("Prediction Accuracy") + #Setting y label
  ggtitle(label = "Per Subject Predictive Accuracy: Humans")
ggsave("output/descriptive_adequacy_humans.png", width = 10, height = 5)

save.image(file="output/hier_orl_humans_data.RData") 
#load(file="output/hier_orl_humans_data.RData") 
