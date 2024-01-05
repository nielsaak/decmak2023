# ORL <- function(payoff,ntrials,a_rew,a_pun,K,theta,omega_f,omega_p) {
ORL <- function(payoff,ntrials,a_rew,K,theta,omega_f,omega_p) {
  
  # arrays to populate for simulation
  x <- array(NA,c(ntrials))
  X <- array(NA,c(ntrials))

  Ev_update <- array(NA,c(ntrials,2))
  Ev <- array(NA,c(ntrials,2))
    
  signX <- array(NA,c(ntrials))
  Ef_cho <- array(NA,c(ntrials,2))
  Ef_not <- array(NA,c(ntrials,2))
  Ef <- array(NA,c(ntrials,2))
  
  PS <- array(NA,c(ntrials,2))
  
  V <- array(NA,c(ntrials,2))
  
  exp_p <- array(NA,c(ntrials,2))
  p <- array(NA,c(ntrials,2))
  
  # free parameters - turn back on when constructing
  #a_rew <- .3
  #a_pun <- .3
  #K <- 3
  #theta <- 3
  #omega_f <- .7
  #omega_p <- .7
  
  x[1] <- rcat(1,c(.5, .5))
  
  X[1] <- payoff[1, x[1]]
  
  Ev[1,] <- rep(0,2)
  
  Ef[1,] <- rep(0,2)
  
  PS[1,] <- rep(1,2)
  
  for (t in 2:ntrials) {
    
    #this is important mention this as constructing model
    #signX[t] <- ifelse(X[t-1]<3,-1,1)
    signX[t] <- ifelse(X[t-1]<0.5,-1,1) # for z-scored pay-off
    
    for (d in 1:2) {
      
      # -------- Updating expected values ------------------------
      # Ev_update[t,d] <- ifelse(X[t-1]>=0,
      #                           Ev[t-1,d] + a_rew*((X[t-1]) - Ev[t-1,d]), 
      #                           Ev[t-1,d] + a_pun*((X[t-1]) - Ev[t-1,d])
      Ev_update[t,d] <- Ev[t-1,d] + a_rew*((X[t-1]) - Ev[t-1,d])
                            
      Ev[t,d] <- ifelse(d==x[t-1],Ev_update[t,d],Ev[t-1,d])
      
      # -------- Updating expected frequencies ------------------------
      #update expected frequencies for ALL decks - AS IF THEY WERE ALL CHOSEN
      # Ef_cho[t,d] <- ifelse(X[t-1]>=0, 
      #                         Ef[t-1,d] + a_rew*(signX[t] - Ef[t-1,d]),
      #                         Ef[t-1,d] + a_pun*(signX[t] - Ef[t-1,d])
      # )
      Ef_cho[t,d] <- Ef[t-1,d] + a_rew*(signX[t] - Ef[t-1,d])
      
      #update expected frequencies for ALL decks - AS IF THEY WERE ALL UNCHOSEN. 
      # Ef_not[t,d] <- ifelse(X[t-1]>=0, 
      #                         Ef[t-1,d] + a_pun*(-(signX[t]/3) - Ef[t-1,d]),
      #                         Ef[t-1,d] + a_rew*(-(signX[t]/3) - Ef[t-1,d])
      # ) 
      Ef_not[t,d] <- Ef[t-1,d] + a_rew*(-(signX[t]) - Ef[t-1,d])
      
      #copy appropriate values to ef variable
      Ef[t,d] <- ifelse(d==x[t-1],Ef_cho[t,d],Ef_not[t,d])  
      
      #-----------Perseverance----------------------------------
      #ifelse needed to disctiminate chosen and unchosen decks
      PS[t,d] <- ifelse(x[t-1]==d,1/(1+K),PS[t-1,d]/(1+K))
      
      #-----------Valence model------------------------------
      V[t,d] <- Ev[t,d] + Ef[t,d]*omega_f + PS[t,d]*omega_p
      
      #----------softmax part 1-------------
      exp_p[t,d] <- exp(theta*V[t,d])
      
    }
    
    #----------softmax part 2-------------
    for (d in 1:2) {
      p[t,d] <- exp_p[t,d]/sum(exp_p[t,])
    }
      
    x[t] <- rcat(1,p[t,])
    
    X[t] <- payoff[t,x[t]]
    
  }
  
  result <- list(x=x,
                 X=X,
                 Ev=Ev,
                 Ef=Ef,
                 PS=PS)
  
  return(result)
  
  
  #turn back on when building
  #par(mfrow=c(2,2))
  #plot(Ev[,1])
  #plot(Ev[,2])
  #plot(Ev[,3])
  #plot(Ev[,4])
  #plot(x)
  
}
