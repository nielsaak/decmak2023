#install.packages("pacman")
pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm)

set.seed(1983)

### NB! Don't forget to set your working directory
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/work/decmak2023/parameter_recovery/hierarchical_parameter_recovery")

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#------ create task environment -------------------
# NB! mod(ntrials, nstruct) (aka. ntrials %% nstruct) must be 0
ntrials <- 80 # total number of trials in our payoff structure
nstruct <- 10 # size of our subdivisions for pseudorandomization
# freq <- 0.5 # probability of our frequent losses (we have losses half of the time)
# infreq <- 0.1 # probability of our infrequent losses (we have losses 1/10th of the time)
# bad_r <- 100 # "bad" winnings
# bad_freq_l <- -250 # "bad" frequent loss
# bad_infreq_l <- -1250 # "bad" infrequent loss
# good_r <- 50 # "good" winnings
# good_freq_l <- -50 # "good" frequent loss
# good_infreq_l <- -250 # "good" infrequent loss

# # Bad frequent
# A_R <- rep(bad_r, nstruct) # we win on every trials
# A_L <- c(rep(bad_freq_l, nstruct*freq),rep(0,nstruct*(1-freq))) # we have losses half of the time
# 
# # Bad infrequent
# B_R <- rep(bad_r, nstruct)
# B_L <- c(rep(bad_infreq_l, nstruct*infreq),rep(0,nstruct*(1-infreq))) # we have losses 1/10th of the time
# 
# # Good frequent
# C_R <- rep(good_r, nstruct)
# C_L <- c(rep(good_freq_l, nstruct*freq),rep(0,nstruct*(1-freq)))
# 
# # Good infrequent
# D_R <- rep(good_r, nstruct)
# D_L <- c(rep(good_infreq_l, nstruct*infreq),rep(0,nstruct*(1-infreq)))
# 
# # create the pseudorandomized full payoff structure
# A <- array(NA,ntrials) # setting up and empty array to be filled
# B <- array(NA,ntrials)
# C <- array(NA,ntrials)
# D <- array(NA,ntrials)
# for (i in 1:(ntrials/nstruct)) {
#   A[(1+(i-1)*nstruct):(i*nstruct)] <- (A_R + sample(A_L)) # randomly shuffling the loss-array for every ten trials (and adding those losses to the winnings)
#   B[(1+(i-1)*nstruct):(i*nstruct)] <- (B_R + sample(B_L))
#   C[(1+(i-1)*nstruct):(i*nstruct)] <- (C_R + sample(C_L))
#   D[(1+(i-1)*nstruct):(i*nstruct)] <- (D_R + sample(D_L))
# }

# "less generic" code
# A <- c()
# B <- c()
# C <- c()
# D <- c()
# for (i in 1:10) {
#   A <- (append(A,A_R + sample(A_L)))
#   B <- (append(B,B_R + sample(B_L)))
#   C <- (append(C,C_R + sample(C_L)))
#   D <- (append(D,D_R + sample(D_L)))
# }

A <- c()

for (i in 1:8) {
  A <- c(A, sample(c(rep(2,5),rep(3,5)), 10, replace = FALSE, prob = c(rep(0.1, 10))))
}

B <- c()

for (i in 1:8) {
  B <- c(B, sample(c(rep(3,1),rep(0,4), rep(1,4), rep(6,1)), 10, replace = FALSE, prob = c(rep(0.1, 10))))
}

# payoff <- cbind(A,B,C,D)/100 # combining all four decks as columns with each 100 trials - dividing our payoffs by 100 to make the numbers a bit easier to work with
payoff <- cbind(A,B) # combining all four decks as columns with each 100 trials - dividing our payoffs by 100 to make the numbers a bit easier to work with
#payoff <- (payoff-mean(payoff))/sd(payoff) # z-scoring payoffs to assess the effect of the different kinds of scaling

# let's look at the payoff
colSums(payoff) # the two bad decks should sum to -25 (i.e. -2500), and the two good ones to 25 (i.e. 2500)

###--------------Run full parameter recovery -------------
niterations <- 100 # fewer because it takes too long
nsubs <- 9 # mimicking the data structure from Ahn et al.
ntrials_all <- rep(ntrials, nsubs) # all 48 simulated subs have 100 trials each

# mu
true_mu_a_rew <- array(NA,c(niterations))
true_mu_K <- array(NA,c(niterations))
#true_mu_theta <- array(NA,c(niterations))
true_mu_omega_f <- array(NA,c(niterations))
true_mu_omega_p <- array(NA,c(niterations))

infer_mu_a_rew <- array(NA,c(niterations))
infer_mu_K <- array(NA,c(niterations))
#infer_mu_theta <- array(NA,c(niterations))
infer_mu_omega_f <- array(NA,c(niterations))
infer_mu_omega_p <- array(NA,c(niterations))

# sigma (SD for R) / lambda (precision for JAGS)
true_lambda_a_rew <- array(NA,c(niterations))
true_lambda_K <- array(NA,c(niterations))
#true_lambda_theta <- array(NA,c(niterations))
true_lambda_omega_f <- array(NA,c(niterations))
true_lambda_omega_p <- array(NA,c(niterations))

infer_lambda_a_rew <- array(NA,c(niterations))
infer_lambda_K <- array(NA,c(niterations))
#infer_lambda_theta <- array(NA,c(niterations))
infer_lambda_omega_f <- array(NA,c(niterations))
infer_lambda_omega_p <- array(NA,c(niterations))

start_time = Sys.time()
for (i in 1:niterations) {
  ntrials <- ntrials_all
  
  # let's see how robust the model is. Does it recover all sorts of values?
  mu_a_rew <- runif(1,0,1)
  mu_K <- runif(1,0,2)
  #mu_theta <- runif(1,.2,2) # could also just be a set value (e.g. 1) to simplify the model a bit
  mu_omega_f <- runif(1,-2,2)
  mu_omega_p <- runif(1,-2,2)
  
  sigma_a_rew <- runif(1,0,0.1)
  sigma_K <- runif(1,0,0.2)
  #sigma_theta <- runif(1,0,0.2) # if theta is just a set value (e.g. 1), then this isn't relevant anymore
  sigma_omega_f <- runif(1,0,0.4)
  sigma_omega_p <- runif(1,0,0.4)
  
  # sigma_a_rew <- runif(1,0,.5)
  # sigma_a_pun <- runif(1,0,.5)
  # sigma_K <- runif(1,0,.5)
  # sigma_theta <- runif(1,0,.5)
  # sigma_omega_f <- runif(1,0,.5)
  # sigma_omega_p <- runif(1,0,.5)
  
  source('hier_ORL_sim_wo_theta.R')
  ORL_sims <- hier_ORL_sim(payoff,nsubs,ntrials,mu_a_rew,
                           mu_K,mu_omega_f,mu_omega_p,
                           sigma_a_rew,sigma_K,
                           sigma_omega_f,sigma_omega_p)
  
  x <- ORL_sims$x
  X <- ORL_sims$X
  
  # set up jags and run jags model
  data <- list("x","X","ntrials","nsubs") 
  params<-c("mu_a_rew",
            "mu_K","mu_omega_f","mu_omega_p",
            "lambda_a_rew","lambda_K",
            "lambda_omega_f","lambda_omega_p")
  samples <- jags.parallel(data, inits=NULL, params,
                           model.file ="hier_ORL_wo_theta.txt", n.chains=3, 
                           n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)
  
  # mu
  true_mu_a_rew[i] <- mu_a_rew
  true_mu_K[i] <- mu_K
  #true_mu_theta[i] <- mu_theta
  true_mu_omega_f[i] <- mu_omega_f
  true_mu_omega_p[i] <- mu_omega_p
  
  # find maximum a posteriori
  Y <- samples$BUGSoutput$sims.list
  infer_mu_a_rew[i] <- MPD(Y$mu_a_rew)
  infer_mu_K[i] <- MPD(Y$mu_K)
  #infer_mu_theta[i] <- MPD(Y$mu_theta)
  infer_mu_omega_f[i] <- MPD(Y$mu_omega_f)
  infer_mu_omega_p[i] <- MPD(Y$mu_omega_p)
  
  # lambda
  true_lambda_a_rew[i] <- sigma_a_rew
  true_lambda_K[i] <- sigma_K
  #true_lambda_theta[i] <- sigma_theta
  true_lambda_omega_f[i] <- sigma_omega_f
  true_lambda_omega_p[i] <- sigma_omega_p
  
  # find maximum a posteriori
  infer_lambda_a_rew[i] <- MPD(Y$lambda_a_rew)
  infer_lambda_K[i] <- MPD(Y$lambda_K)
  #infer_lambda_theta[i] <- MPD(Y$lambda_theta)
  infer_lambda_omega_f[i] <- MPD(Y$lambda_omega_f)
  infer_lambda_omega_p[i] <- MPD(Y$lambda_omega_p)
  
  print(i)
  end_time = Sys.time()
  print(end_time - start_time)
}

end_time = Sys.time()
end_time - start_time

#load("hier_ORL_recovery_wo_theta.RData")

# let's look at some scatter plots
# plotting code courtesy of Lasse
source('recov_plot.R')
pl1 <- recov_plot(true_mu_a_rew, infer_mu_a_rew, plot_lab_1 = expression("True "*mu[arew]), plot_lab_2 = expression("Inferred "*mu[arew]), 'smoothed linear fit', title=expression(mu[arew])) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) 
pl4 <- recov_plot(true_mu_K, infer_mu_K, plot_lab_1 = expression("True "*mu[K]), plot_lab_2 = expression("Inferred "*mu[K]), 'smoothed linear fit', title=expression(mu[K])) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0,1.7))
pl5 <- recov_plot(true_mu_omega_f, infer_mu_omega_f, plot_lab_1 = expression("True "*mu[omega][F]), plot_lab_2 = expression("Inferred "*mu[omega][F]), 'smoothed linear fit', title=expression(mu[omega][F])) +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2,3)) 
pl6 <- recov_plot(true_mu_omega_p, infer_mu_omega_p, plot_lab_1 = expression("True "*mu[omega][P]), plot_lab_2 = expression("Inferred "*mu[omega][P]), 'smoothed linear fit', title=expression(mu[omega][P])) +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2,3)) 

plot_1 <- ggarrange(pl1, pl4, pl5, pl6, ncol = 2, nrow = 2)
print(plot_1)

# Saving in a different (weird) way to ensure symbol formatting
png(filename = "output/recovery_mu_80_trials_100_iter_theta1.png", type = "cairo-png", height = 1200, width = 1200, res = 170)
print(plot_1)
dev.off()

#ggsave(plot = plot_1, "output/recovery_mu_80_trials_100_iter_theta1.png")

pl1 <- recov_plot(true_lambda_a_rew, sqrt(1/infer_lambda_a_rew), plot_lab_1 = expression("True "*lambda[arew]), plot_lab_2 = expression("Inferred "*lambda[arew]), 'smoothed linear fit', title=expression(lambda[arew])) +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(0,0.5)) 
pl4 <- recov_plot(true_lambda_K, sqrt(1/infer_lambda_K), plot_lab_1 = expression("True "*lambda[K]), plot_lab_2 = expression("Inferred "*lambda[K]), 'smoothed linear fit', title=expression(lambda[K])) +
  coord_cartesian(xlim = c(0, 0.2), ylim = c(0,3)) 
pl5 <- recov_plot(true_lambda_omega_f, sqrt(1/infer_lambda_omega_f), plot_lab_1 = expression("True "*lambda[omega][F]), plot_lab_2 = expression("Inferred "*lambda[omega][F]), 'smoothed linear fit', title=expression(lambda[omega][F])) +
  coord_cartesian(xlim = c(0, 0.4), ylim = c(0,1)) 
pl6 <- recov_plot(true_lambda_omega_p, sqrt(1/infer_lambda_omega_p), plot_lab_1 = expression("True "*lambda[omega][P]), plot_lab_2 = expression("Inferred "*lambda[omega][P]), 'smoothed linear fit', title=expression(lambda[omega][P])) +
  coord_cartesian(xlim = c(0, 0.4), ylim = c(0,3.5))

plot_2 <- ggarrange(pl1, pl4, pl5, pl6, ncol = 2, nrow = 2)
print(plot_2)

# Saving in a different (weird) way to ensure symbol formatting
png(filename = "output/recovery_lambda_80_trials_100_iter_theta1.png", type = "cairo-png", height = 1200, width = 1200, res = 170)
print(plot_2)
dev.off()

#ggsave(plot = plot_2, "output/recovery_lambda_80_trials_100_iter_theta1.png", width = 10, height = 20)

traceplot(samples)

save.image(file = "hier_ORL_recovery_wo_theta.RData")