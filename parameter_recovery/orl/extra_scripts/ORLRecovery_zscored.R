install.packages("pacman")
pacman::p_load(hesim, extraDistr, R2jags, parallel, ggpubr)

set.seed(1982)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}


#------ create task environment -------------------
# NB! mod(ntrials, nstruct) (aka. ntrials %% nstruct) must be 0
ntrials <- 100 # total number of trials in our payoff structure
nstruct <- 10 # size of our subdivisions for pseudorandomization
freq <- 0.5 # probability of our frequent losses (we have losses half of the time)
infreq <- 0.1 # probability of our infrequent losses (we have losses 1/10th of the time)
bad_r <- 100 # "bad" winnings
bad_freq_l <- -250 # "bad" frequent loss
bad_infreq_l <- -1250 # "bad" infrequent loss
good_r <- 50 # "good" winnings
good_freq_l <- -50 # "good" frequent loss
good_infreq_l <- -250 # "good" infrequent loss

# Bad frequent
A_R <- rep(bad_r, nstruct) # we win on every trials
A_L <- c(rep(bad_freq_l, nstruct*freq),rep(0,nstruct*(1-freq))) # we have losses half of the time

# Bad infrequent
B_R <- rep(bad_r, nstruct)
B_L <- c(rep(bad_infreq_l, nstruct*infreq),rep(0,nstruct*(1-infreq))) # we have losses 1/10th of the time

# Good frequent
C_R <- rep(good_r, nstruct)
C_L <- c(rep(good_freq_l, nstruct*freq),rep(0,nstruct*(1-freq)))

# Good infrequent
D_R <- rep(good_r, nstruct)
D_L <- c(rep(good_infreq_l, nstruct*infreq),rep(0,nstruct*(1-infreq)))

# create the pseudorandomized full payoff structure
A <- array(NA,ntrials) # setting up and empty array to be filled
B <- array(NA,ntrials)
C <- array(NA,ntrials)
D <- array(NA,ntrials)
for (i in 1:(ntrials/nstruct)) {
  A[(1+(i-1)*nstruct):(i*nstruct)] <- (A_R + sample(A_L)) # randomly shuffling the loss-array for every ten trials (and adding those losses to the winnings)
  B[(1+(i-1)*nstruct):(i*nstruct)] <- (B_R + sample(B_L))
  C[(1+(i-1)*nstruct):(i*nstruct)] <- (C_R + sample(C_L))
  D[(1+(i-1)*nstruct):(i*nstruct)] <- (D_R + sample(D_L))
}

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

payoff <- cbind(A,B,C,D) # combining all four decks as columns with each 100 trials 
payoff <- (payoff-mean(payoff))/sd(payoff) # z-scoring payoffs to assess the effect of the different kinds of scaling

# let's look at the payoff
colSums(payoff) # the two bad decks should sum to -25 (i.e. -2500), and the two good ones to 25 (i.e. 2500)
#----------------------------------------------------


#-------test ORL delta function and jags script ---------

#---set params

a_rew <- .3
a_pun <- .3
K <- 2
theta <- 2
omega_f <- .7
omega_p <- .7

# ntrials <- 100

source("ORL.R")
ORL_sims <- ORL(payoff,ntrials,a_rew,a_pun,K,theta,omega_f,omega_p)

par(mfrow=c(2,2))
plot(ORL_sims$Ev[,1])
plot(ORL_sims$Ev[,2])
plot(ORL_sims$Ev[,3])
plot(ORL_sims$Ev[,4])

x <- ORL_sims$x
X <- ORL_sims$X

# set up jags and run jags model
data <- list("x","X","ntrials") 
params<-c("a_rew","a_pun","K","theta","omega_f","omega_p")
samples <- jags.parallel(data, inits=NULL, params,
                model.file ="ORL.txt", n.chains=3, 
                n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=64)


###--------------Run full parameter recovery -------------
niterations <- 100 # fewer because it takes too long

true_a_rew <- array(NA,c(niterations))
true_a_pun <- array(NA,c(niterations))
true_K <- array(NA,c(niterations))
true_theta <- array(NA,c(niterations))
true_omega_f <- array(NA,c(niterations))
true_omega_p <- array(NA,c(niterations))

infer_a_rew <- array(NA,c(niterations))
infer_a_pun <- array(NA,c(niterations))
infer_K <- array(NA,c(niterations))
infer_theta <- array(NA,c(niterations))
infer_omega_f <- array(NA,c(niterations))
infer_omega_p <- array(NA,c(niterations))

start_time = Sys.time()

for (i in 1:niterations) {
  
  # let's see how robust the model is. Does it recover all sorts of values?
  a_rew <- runif(1,0,1)
  a_pun <- runif(1,0,1)
  K <- runif(1,0,2)
  theta <- runif(1,.2,3) # could also just be a set value (e.g. 1) to simplify the model a bit
  omega_f <- runif(1,-2,2)
  omega_p <- runif(1,-2,2)
  
  ORL_sims <- ORL(payoff,ntrials,a_rew,a_pun,K,theta,omega_f,omega_p)
  
  x <- ORL_sims$x
  X <- ORL_sims$X
  
  # set up jags and run jags model
  data <- list("x","X","ntrials") 
  params<-c("a_rew","a_pun","K","theta","omega_f","omega_p")
  samples <- jags.parallel(data, inits=NULL, params,
                  model.file ="ORL.txt", n.chains=3, 
                  n.iter=3000, n.burnin=1000, n.thin=1, n.cluster=64)
  
  
  true_a_rew[i] <- a_rew
  true_a_pun[i] <- a_pun
  true_K[i] <- K
  true_theta[i] <- theta
  true_omega_f[i] <- omega_f
  true_omega_p[i] <- omega_p
  
  # find maximum a posteriori
  Y <- samples$BUGSoutput$sims.list
  infer_a_rew[i] <- MPD(Y$a_rew)
  infer_a_pun[i] <- MPD(Y$a_pun)
  infer_K[i] <- MPD(Y$K)
  infer_theta[i] <- MPD(Y$theta)
  infer_omega_f[i] <- MPD(Y$omega_f)
  infer_omega_p[i] <- MPD(Y$omega_p)

  print(i)
  
}

end_time = Sys.time()
end_time - start_time

# let's look at some scatter plots

par(mfrow=c(3,2))
plot(true_a_rew,infer_a_rew)
plot(true_a_pun,infer_a_pun)
plot(true_K,infer_K)
plot(true_theta,infer_theta)
plot(true_omega_f,infer_omega_f)
plot(true_omega_p,infer_omega_p)

# plotting code courtesy of Lasse
source('216377/Module3/recov_plot.R')
pl1 <- recov_plot(true_a_rew, infer_a_rew, c("true a_rew", "infer a_rew"), 'smoothed linear fit')
pl2 <- recov_plot(true_a_pun, infer_a_pun, c("true a_pun", "infer a_pun"), 'smoothed linear fit')
pl3 <- recov_plot(true_K, infer_K, c("true K", "infer K"), 'smoothed linear fit')
pl4 <- recov_plot(true_theta, infer_theta, c("true theta", "infer theta"), 'smoothed linear fit')
pl5 <- recov_plot(true_omega_f, infer_omega_f, c("true omega_f", "infer omega_f"), 'smoothed linear fit')
pl6 <- recov_plot(true_omega_p, infer_omega_p, c("true omega_p", "infer omega_p"), 'smoothed linear fit')
ggarrange(pl1, pl2, pl3, pl4, pl5, pl6)

# for investigating multi-colinearity

# par(mfrow=c(2,2))
# plot(true_a_rew,true_a_pun)
# plot(infer_a_rew,infer_a_pun)
# plot(true_omega_f,true_omega_p)
# plot(infer_omega_f,infer_omega_p)
# 
# par(mfrow=c(2,2))
# plot(true_a_rew,true_omega_f)
# plot(infer_a_rew,infer_omega_f)
# plot(true_a_rew,true_omega_p)
# plot(infer_a_rew,infer_omega_p)









