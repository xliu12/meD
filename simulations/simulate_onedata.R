

library(mvtnorm)
library(tidyverse)
library(glue)
library(SuperLearner)
library(origami)


devtools::load_all("meD")

GenData <- function(
    sim.seed = 123, 
    n = 500, 
    Mfamily = "binomial", # "gaussian", 
    Yfamily = "binomial", # "gaussian",
    quadratic.R = FALSE,
    quadratic.tt = FALSE,
    quadratic.M = FALSE,
    quadratic.Y = FALSE,
    corr_c = TRUE,
    num_c = 20
) { 
  set.seed(seed = sim.seed)
  n <- n
  # baseline covariates ----
  # num_c <- 5
  cov_c <- diag(1, nrow = num_c)
  if (corr_c) {
    cov_c[upper.tri(cov_c)|lower.tri(cov_c)] <- .2
  }
  C_dat <- as.data.frame(rmvnorm(n, sigma = diag(1, nrow = num_c))) 
  
  # Subgroup indicator ----
  C2 <- model.matrix(~.-1 + ., C_dat)
  if (quadratic.R == TRUE) {
    C2 <- (C_dat^2 - 1) / sqrt(2)
  }
  bc <- sqrt(0.26/(1-0.26)/ncol(C2)) / rep(1, ncol(C2))
  
  R <- 1*(as.matrix(C2) %*% bc + rnorm(n) > 0)
  
  # var(as.matrix(C2) %*% bc) / var(as.matrix(C2) %*% bc + rnorm(n))
  
  # Treatment ------
  rc_dat <- model.matrix(~.-1 + ., data.frame(R=R, C_dat=C_dat))
  if (quadratic.tt == TRUE) {
    C2 <- (C_dat^2 - 1) / sqrt(2)
    rc_dat <- model.matrix(~.-1 + ., data.frame(R=R, C_dat=C2))
  }
  
  br <- sqrt(0.13/(1-0.13)/ncol(rc_dat)) / rep(1, ncol(rc_dat))
  
  tt <- 1*(rc_dat %*% br + rnorm(n) > 0) 
  # var(rc_dat[,2:6] %*% br[2:6]) / (var(rc_dat %*% br)+1)
  
  # Intermediate ----
  trc_dat <- model.matrix(~.-1 + ., data.frame(tt=tt, R=R, ttR = tt*R, C_dat=C_dat))
  if (quadratic.M == TRUE) {
    C2 <- (C_dat^2 - 1) / sqrt(2)
    trc_dat <- model.matrix(~.-1 + ., data.frame(tt=tt, R=R, ttR = tt*R, C_dat=C2))
  }
  
  btt <- c( -0.3, 0.5, 1, sqrt(0.26/(1-0.26)/ncol(C_dat))/rep(1, ncol(C_dat)) )
  names(btt) <- colnames(trc_dat)
  
  
  # var(tts %*% btt)
  if (Mfamily == "gaussian") {
    if (quadratic.M == TRUE) {
      C2 <- data.frame(C_dat, (C_dat^2 )  )
      trc_dat <- model.matrix(~.-1 + ., data.frame(tt=tt, R=R, ttR = tt*R, C_dat=C2))
      btt <- c(0,0,0.3, rep(sqrt(0.26/(1-0.26)/ncol(C2)), ncol(C2)))
      
    }
    names(btt) <- colnames(trc_dat)
    btt[c("tt","R")] <- 0; btt["ttR"] <- 0.3
    
    M <- trc_dat %*% btt + rnorm(n)
    var(trc_dat[,4:8] %*% btt[4:8])/var(M)
  }
  if (Mfamily == "binomial") {
    if (quadratic.M == TRUE) {
      C2 <- data.frame(C_dat, (C_dat^2 )  )
      trc_dat <- model.matrix(~.-1 + ., data.frame(tt=tt, R=R, ttR = tt*R, C_dat=C2))
      btt <- c(0,0,0.3, rep(sqrt(0.26/(1-0.26)/ncol(C2)), ncol(C2)))
      
    }
    names(btt) <- colnames(trc_dat)
    btt[c("tt","R")] <- 0; btt["ttR"] <- 0.3
    
    M <- 1*(trc_dat %*% btt + rnorm(n) > 0) 
    # rbinom(n, 1, prob = plogis(trc_dat %*% btt))
  }
  
  # Outcome -----
  mtrc_dat <- model.matrix(~.-1 + ., data.frame(M=M, tt=tt, R=R, ttM = M*tt, ttR = tt*R, RM = R*M, ttRM = tt*R*M, C_dat=C_dat))
  if (quadratic.Y == TRUE) {
    C2 <- (C_dat^2 - 1) / sqrt(2)
    mtrc_dat <- model.matrix(~.-1 + ., data.frame(M=M, tt=tt, R=R, ttM = M*tt, ttR = tt*R, RM = R*M, ttRM = tt*R*M, C_dat=C2))
  }
  bm <- c(c(0.6, -1.2, -1.4),  0.7, 0.8, 0.9, 1.1, sqrt(0.26/(1-0.26)/ncol(C_dat)) / rep(1, ncol(C_dat)))
  
  
  names(bm) <- colnames(mtrc_dat)
  
  if (Yfamily == "gaussian") {
    if (quadratic.Y == TRUE) {
      C2 <- data.frame(C_dat, (C_dat^2 )  )
      mtrc_dat <- model.matrix(~.-1 + ., data.frame(M=M, tt=tt, R=R, ttM = M*tt, ttR = tt*R, RM = R*M, ttRM = tt*R*M, C_dat=C2))
      
      bm <- c(c(0.4, 0, 0),  0, 0.3, 0, 0.3, sqrt(0.26/(1-0.26)/ncol(C2)) / rep(1, ncol(C2)))
    }
    names(bm) <- colnames(mtrc_dat)
    bm[c("tt","R", "ttM", "RM")] <-0 
    bm[c("M", "ttR", "ttRM")] <- c(0.3, 0.3, 0.3)
    
    Y <- mtrc_dat %*% bm + rnorm(n)
    # var(mtrc_dat[,8:12]%*% bm[8:12])/var(Y)
  }
  if (Yfamily == "binomial") {
    if (quadratic.Y == TRUE) {
      C2 <- data.frame(C_dat, (C_dat^2 )  )
      mtrc_dat <- model.matrix(~.-1 + ., data.frame(M=M, tt=tt, R=R, ttM = M*tt, ttR = tt*R, RM = R*M, ttRM = tt*R*M, C_dat=C2))
      
      bm <- c(c(0.4, 0, 0),  0, 0.3, 0, 0.3, sqrt(0.26/(1-0.26)/ncol(C2)) / rep(1, ncol(C2)))
    }
    names(bm) <- colnames(mtrc_dat)
    bm[c("tt","R", "ttM", "RM")] <-0 
    bm[c("M", "ttR", "ttRM")] <- c(0.3, 0.3, 0.3)
    
    Y <- 1*(mtrc_dat %*% bm + rnorm(n) > 0) 
    #rbinom(n, 1, prob = plogis(mtrc_dat %*% bm))
  }
  
  datobs <- data.frame(Y=Y, Mdat=M, tt=tt, R=R, Cdat=C_dat)
  
  
  # True values 
  ## continuous M & Y 
  if ( (Yfamily == "gaussian") & (Mfamily == "gaussian") ) {
    # linear Y models
    meanC_dat <- C_dat #matrix(rep(0, ncol(C_dat)), nrow = 1)
    # when tt = 1 
    # Yt1r1
    t1r1_dat <- as.data.frame(trc_dat) %>% mutate(tt=1,R=1,ttR=1)#model.matrix(~.-1 + ., data.frame(tt=1, R=1, ttR = 1, C_dat=meanC_dat))
    meanMt1r1 <- as.numeric(as.matrix(t1r1_dat) %*% btt)
    meanMt1r1_dat <- as.data.frame(mtrc_dat) %>% 
      mutate(M=meanMt1r1, tt=1,R=1) %>% mutate(ttR=tt*R, ttM=tt*M, RM=R*M,ttRM=tt*R*M) %>% as.matrix()
    #model.matrix(~.-1 + ., data.frame(M=meanMt1r1, tt=1, R=1, ttM = meanMt1r1*1, ttR = 1, RM = 1*meanMt1r1, ttRM = 1*1*meanMt1r1, C_dat=meanC_dat))
    meanYt1r1 <- as.numeric(meanMt1r1_dat %*% bm) %>% mean(.)
    
    # Yt1r1.Mt1r0
    t1r0_dat <- as.data.frame(trc_dat) %>% 
      mutate(tt=1,R=0) %>% mutate(ttR=tt*R) %>% as.matrix()
    #model.matrix(~.-1 + ., data.frame(tt=1, R=0, ttR = 0, C_dat=meanC_dat))
    meanMt1r0 <- as.numeric(t1r0_dat %*% btt)
    meanMt1r0_dat <- as.data.frame(mtrc_dat) %>% 
      mutate(M=meanMt1r0, tt=1,R=1) %>% mutate(ttR=tt*R, ttM=tt*M, RM=R*M,ttRM=tt*R*M) %>% as.matrix()
      #model.matrix(~.-1 + ., data.frame(M=meanMt1r0, tt=1, R=1, ttM = meanMt1r0*1, ttR = 1, RM = 1*meanMt1r0, ttRM = 1*1*meanMt1r0, C_dat=meanC_dat))
    meanYt1r1.Mt1r0 <- as.numeric(meanMt1r0_dat %*% bm) %>% mean(.)
    
    # Yt1r0
    # t1r0_dat <- as.data.frame(trc_dat) %>% mutate(tt=1,R=0) %>% mutate(ttR=tt*R) %>% as.matrix() #model.matrix(~.-1 + ., data.frame(tt=1, R=0, ttR = 0, C_dat=meanC_dat))
    # meanMt1r0 <- as.numeric(t1r0_dat %*% btt)
    meanMt1r0_dat <- as.data.frame(mtrc_dat) %>% 
      mutate(M=meanMt1r0, tt=1,R=0) %>% mutate(ttR=tt*R, ttM=tt*M, RM=R*M,ttRM=tt*R*M) %>% as.matrix() # model.matrix(~.-1 + ., data.frame(M=meanMt1r0, tt=1, R=0, ttM = meanMt1r0*1, ttR = 0, RM = 0*meanMt1r0, ttRM = 1*0*meanMt1r0, C_dat=meanC_dat))
    meanYt1r0 <- as.numeric(meanMt1r0_dat %*% bm) %>% mean(.)
    
    # when tt = 0 
    # Yt0r1
    t0r1_dat <- as.data.frame(trc_dat) %>% mutate(tt=0,R=1) %>% mutate(ttR=tt*R) %>% as.matrix() #model.matrix(~.-1 + ., data.frame(tt=0, R=1, ttR = 0, C_dat=meanC_dat))
    meanMt0r1 <- as.numeric(t0r1_dat %*% btt)
    meanMt0r1_dat <- as.data.frame(mtrc_dat) %>% 
      mutate(M=meanMt0r1, tt=0,R=1) %>% mutate(ttR=tt*R, ttM=tt*M, RM=R*M,ttRM=tt*R*M) %>% as.matrix() # model.matrix(~.-1 + ., data.frame(M=meanMt0r1, tt=0, R=1, ttM = meanMt0r1*0, ttR = 0, RM = 1*meanMt0r1, ttRM = 0*1*meanMt0r1, C_dat=meanC_dat))
    meanYt0r1 <- as.numeric(meanMt0r1_dat %*% bm) %>% mean(.)
    
    # Yt0r1.Mt0r0
    t0r0_dat <- as.data.frame(trc_dat) %>% mutate(tt=0,R=0) %>% mutate(ttR=tt*R) %>% as.matrix() #model.matrix(~.-1 + ., data.frame(tt=0, R=0, ttR = 0, C_dat=meanC_dat))
    meanMt0r0 <- as.numeric(t0r0_dat %*% btt)
    meanMt0r0_dat <- as.data.frame(mtrc_dat) %>% 
      mutate(M=meanMt0r0, tt=0,R=1) %>% mutate(ttR=tt*R, ttM=tt*M, RM=R*M,ttRM=tt*R*M) %>% as.matrix() # model.matrix(~.-1 + ., data.frame(M=meanMt0r0, tt=0, R=1, ttM = meanMt0r0*0, ttR = 0, RM = 1*meanMt0r0, ttRM = 0*1*meanMt0r0, C_dat=meanC_dat))
    meanYt0r1.Mt0r0 <- as.numeric(meanMt0r0_dat %*% bm) %>% mean(.)
    
    # Yt0r0
    # t0r0_dat <- model.matrix(~.-1 + ., data.frame(tt=0, R=0, ttR = 0, C_dat=meanC_dat))
    # meanMt0r0 <- as.numeric(t0r0_dat %*% btt) %>% mean(.)
    meanMt0r0_dat <- as.data.frame(mtrc_dat) %>% 
      mutate(M=meanMt0r0, tt=0,R=0) %>% mutate(ttR=tt*R, ttM=tt*M, RM=R*M,ttRM=tt*R*M) %>% as.matrix() # model.matrix(~.-1 + ., data.frame(M=meanMt0r0, tt=0, R=0, ttM = meanMt0r0*0, ttR = 0, RM = 0*meanMt0r0, ttRM = 0*0*meanMt0r0, C_dat=meanC_dat))
    meanYt0r0 <- as.numeric(meanMt0r0_dat %*% bm) %>% mean(.)
    
    
  }
  
  ## binomial M or Y 
  if ( (Mfamily == "binomial") ) {
    
    # when tt = 1 
    # Yt1r1
    t1r1_dat <- data.frame(trc_dat); t1r1_dat$tt <- 1;  t1r1_dat$R <- 1
    t1r1_dat$ttR <- t1r1_dat$tt * t1r1_dat$R
    # model.matrix(~.-1 + ., data.frame(tt=1, R=1, ttR = 1, C_dat=C_dat))
    pMt1r1 <- pnorm(as.matrix(t1r1_dat) %*% btt, lower.tail = T)
    #plogis(t1r1_dat %*% btt)
    m1t1r1_dat <- data.frame(mtrc_dat); m1t1r1_dat$M <- 1; m1t1r1_dat$tt <- 1; m1t1r1_dat$R <- 1
    m1t1r1_dat$ttM <- m1t1r1_dat$tt * m1t1r1_dat$M
    m1t1r1_dat$ttR <- m1t1r1_dat$tt * m1t1r1_dat$R
    m1t1r1_dat$RM <- m1t1r1_dat$R * m1t1r1_dat$M
    m1t1r1_dat$ttRM <- m1t1r1_dat$tt * m1t1r1_dat$R * m1t1r1_dat$M
    #model.matrix(~.-1 + ., data.frame(M=1, tt=1, R=1, ttM = 1*1, ttR = 1, RM = 1*1, ttRM = 1*1*1, C_dat=C_dat))
    m0t1r1_dat <- data.frame(mtrc_dat); m0t1r1_dat$M <- 0; m0t1r1_dat$tt <- 1; m0t1r1_dat$R <- 1
    m0t1r1_dat$ttM <- m0t1r1_dat$tt * m0t1r1_dat$M
    m0t1r1_dat$ttR <- m0t1r1_dat$tt * m0t1r1_dat$R
    m0t1r1_dat$RM <- m0t1r1_dat$R * m0t1r1_dat$M
    m0t1r1_dat$ttRM <- m0t1r1_dat$tt * m0t1r1_dat$R * m0t1r1_dat$M
    #model.matrix(~.-1 + ., data.frame(M=0, tt=1, R=1, ttM = 1*0, ttR = 1, RM = 1*0, ttRM = 1*1*0, C_dat=C_dat))
    if (Yfamily == "binomial") {
      Yt1r1 <- pnorm(as.matrix(m1t1r1_dat) %*% bm) * pMt1r1 + pnorm(as.matrix(m0t1r1_dat) %*% bm) * (1-pMt1r1)
    }
    if (Yfamily == "gaussian") {
      Yt1r1 <- (as.matrix(m1t1r1_dat) %*% bm) * pMt1r1 + (as.matrix(m0t1r1_dat) %*% bm) * (1-pMt1r1)
    }
    meanYt1r1 <- mean(Yt1r1)
    
    # Yt1r1.Mt1r0
    t1r0_dat <- data.frame(trc_dat); t1r0_dat$tt <- 1;  t1r0_dat$R <- 0
    t1r0_dat$ttR <- t1r0_dat$tt * t1r0_dat$R
    # model.matrix(~.-1 + ., data.frame(tt=1, R=0, ttR = 0, C_dat=C_dat))
    pMt1r0 <- pnorm(as.matrix(t1r0_dat) %*% btt)
    if (Yfamily == "binomial") {
      Yt1r1.Mt1r0 <- pnorm(as.matrix(m1t1r1_dat) %*% bm) * pMt1r0 + pnorm(as.matrix(m0t1r1_dat) %*% bm) * (1-pMt1r0)
    }
    if (Yfamily == "gaussian") {
      Yt1r1.Mt1r0 <- (as.matrix(m1t1r1_dat) %*% bm) * pMt1r0 + (as.matrix(m0t1r1_dat) %*% bm) * (1-pMt1r0)
    }
    meanYt1r1.Mt1r0 <- mean(Yt1r1.Mt1r0)
    
    
    # Yt1r0
    m1t1r0_dat <- data.frame(mtrc_dat); 
    m1t1r0_dat$M <- 1; m1t1r0_dat$tt <- 1; m1t1r0_dat$R <- 0
    m1t1r0_dat$ttM <- m1t1r0_dat$tt * m1t1r0_dat$M
    m1t1r0_dat$ttR <- m1t1r0_dat$tt * m1t1r0_dat$R
    m1t1r0_dat$RM <- m1t1r0_dat$R * m1t1r0_dat$M
    m1t1r0_dat$ttRM <- m1t1r0_dat$tt * m1t1r0_dat$R * m1t1r0_dat$M
    #model.matrix(~.-1 + ., data.frame(M=1, tt=1, R=0, ttM = 1*1, ttR = 0, RM = 0*1, ttRM = 1*0*1, C_dat=C_dat))
    m0t1r0_dat <- data.frame(mtrc_dat); 
    m0t1r0_dat$M <- 0; m0t1r0_dat$tt <- 1; m0t1r0_dat$R <- 0
    m0t1r0_dat$ttM <- m0t1r0_dat$tt * m0t1r0_dat$M
    m0t1r0_dat$ttR <- m0t1r0_dat$tt * m0t1r0_dat$R
    m0t1r0_dat$RM <- m0t1r0_dat$R * m0t1r0_dat$M
    m0t1r0_dat$ttRM <- m0t1r0_dat$tt * m0t1r0_dat$R * m0t1r0_dat$M
    #model.matrix(~.-1 + ., data.frame(M=0, tt=1, R=0, ttM = 1*0, ttR = 0, RM = 0*0, ttRM = 1*0*0, C_dat=C_dat))
    if (Yfamily == "binomial") {
      Yt1r0 <- pnorm(as.matrix(m1t1r0_dat) %*% bm) * pMt1r0 + pnorm(as.matrix(m0t1r0_dat) %*% bm) * (1-pMt1r0)
    }
    if (Yfamily == "gaussian") {
      Yt1r0 <- (as.matrix(m1t1r0_dat) %*% bm) * pMt1r0 + (as.matrix(m0t1r0_dat) %*% bm) * (1-pMt1r0)
    }
    
    meanYt1r0 <- mean(Yt1r0)
    
    # when tt = 0 
    # Yt0r1
    t0r1_dat <- data.frame(trc_dat); 
    t0r1_dat$tt <- 0;  t0r1_dat$R <- 1
    t0r1_dat$ttR <- t0r1_dat$tt * t0r1_dat$R
    
    pMt0r1 <- pnorm(as.matrix(t0r1_dat) %*% btt, lower.tail = T)
    
    m1t0r1_dat <- data.frame(mtrc_dat); 
    m1t0r1_dat$M <- 1; m1t0r1_dat$tt <- 0; m1t0r1_dat$R <- 1
    m1t0r1_dat$ttM <- m1t0r1_dat$tt * m1t0r1_dat$M
    m1t0r1_dat$ttR <- m1t0r1_dat$tt * m1t0r1_dat$R
    m1t0r1_dat$RM <- m1t0r1_dat$R * m1t0r1_dat$M
    m1t0r1_dat$ttRM <- m1t0r1_dat$tt * m1t0r1_dat$R * m1t0r1_dat$M
    
    m0t0r1_dat <- data.frame(mtrc_dat); 
    m0t0r1_dat$M <- 0; m0t0r1_dat$tt <- 0; m0t0r1_dat$R <- 1
    m0t0r1_dat$ttM <- m0t0r1_dat$tt * m0t0r1_dat$M
    m0t0r1_dat$ttR <- m0t0r1_dat$tt * m0t0r1_dat$R
    m0t0r1_dat$RM <- m0t0r1_dat$R * m0t0r1_dat$M
    m0t0r1_dat$ttRM <- m0t0r1_dat$tt * m0t0r1_dat$R * m0t0r1_dat$M
    if (Yfamily == "binomial") {
      Yt0r1 <- pnorm(as.matrix(m1t0r1_dat) %*% bm) * pMt0r1 + pnorm(as.matrix(m0t0r1_dat) %*% bm) * (1-pMt0r1)
    }
    if (Yfamily == "gaussian") {
      Yt0r1 <- (as.matrix(m1t0r1_dat) %*% bm) * pMt0r1 + (as.matrix(m0t0r1_dat) %*% bm) * (1-pMt0r1)
    }
    meanYt0r1 <- mean(Yt0r1)
    
    # Yt0r1.Mt0r0
    t0r0_dat <- data.frame(trc_dat); 
    t0r0_dat$tt <- 0;  t0r0_dat$R <- 0
    t0r0_dat$ttR <- t0r0_dat$tt * t0r0_dat$R
    
    pMt0r0 <- pnorm(as.matrix(t0r0_dat) %*% btt, lower.tail = T)
    if (Yfamily == "binomial") {
      Yt0r1.Mt0r0 <- pnorm(as.matrix(m1t0r1_dat) %*% bm) * pMt0r0 + pnorm(as.matrix(m0t0r1_dat) %*% bm) * (1-pMt0r0)
    }
    if (Yfamily == "gaussian") {
      Yt0r1.Mt0r0 <- (as.matrix(m1t0r1_dat) %*% bm) * pMt0r0 + (as.matrix(m0t0r1_dat) %*% bm) * (1-pMt0r0)
    }
    meanYt0r1.Mt0r0 <- mean(Yt0r1.Mt0r0)
    
    
    # Yt0r0
    m1t0r0_dat <- data.frame(mtrc_dat); 
    m1t0r0_dat$M <- 1; m1t0r0_dat$tt <- 0; m1t0r0_dat$R <- 0
    m1t0r0_dat$ttM <- m1t0r0_dat$tt * m1t0r0_dat$M
    m1t0r0_dat$ttR <- m1t0r0_dat$tt * m1t0r0_dat$R
    m1t0r0_dat$RM <- m1t0r0_dat$R * m1t0r0_dat$M
    m1t0r0_dat$ttRM <- m1t0r0_dat$tt * m1t0r0_dat$R * m1t0r0_dat$M
    
    m0t0r0_dat <- data.frame(mtrc_dat); 
    m0t0r0_dat$M <- 0; m0t0r0_dat$tt <- 0; m0t0r0_dat$R <- 0
    m0t0r0_dat$ttM <- m0t0r0_dat$tt * m0t0r0_dat$M
    m0t0r0_dat$ttR <- m0t0r0_dat$tt * m0t0r0_dat$R
    m0t0r0_dat$RM <- m0t0r0_dat$R * m0t0r0_dat$M
    m0t0r0_dat$ttRM <- m0t0r0_dat$tt * m0t0r0_dat$R * m0t0r0_dat$M
    
    if (Yfamily == "binomial") {
      Yt0r0 <- pnorm(as.matrix(m1t0r0_dat) %*% bm) * pMt0r0 + pnorm(as.matrix(m0t0r0_dat) %*% bm) * (1-pMt0r0)
    }
    if (Yfamily == "gaussian") {
      Yt0r0 <- (as.matrix(m1t0r0_dat) %*% bm) * pMt0r0 + (as.matrix(m0t0r0_dat) %*% bm) * (1-pMt0r0)
    }
    
    meanYt0r0 <- mean(Yt0r0)
    
  }
  out <- mget( ls(), envir = environment() )
  
  return(out)
}





# Simulation----------------------------

condition_all <- filter(data.frame(expand.grid(
  n = c(3000, 1000, 500, 300),
  # generating data
  quadratic.R = c(T, F),
  quadratic.tt = c(T, F),
  quadratic.M = c(T, F),
  quadratic.Y = c(T, F),
  Yfamily = c("gaussian"), 
  Mfamily = c("gaussian"),
  num_c = c(5), 
  corr_c = c(T)
)), (quadratic.R+quadratic.tt+quadratic.M+quadratic.Y == 0)|(quadratic.R+quadratic.tt+quadratic.M+quadratic.Y == 4)) %>% 
  bind_rows(filter(data.frame(expand.grid(
    n = c(1000, 500, 300),
    # generating data
    quadratic.R = c(T, F),
    quadratic.tt = c(T, F),
    quadratic.M = c(T, F),
    quadratic.Y = c(T, F),
    Yfamily = c("binomial"), 
    Mfamily = c("binomial"),
    num_c = c(20),
    corr_c = c(F)
  )), (quadratic.R+quadratic.tt+quadratic.M+quadratic.Y == 0)  |
    ((quadratic.R+quadratic.tt+quadratic.Y==0)&(quadratic.M==1))  |
    ((quadratic.R+quadratic.tt+quadratic.M==0)&(quadratic.Y==1)) |
    ((quadratic.R+quadratic.tt==2)&(quadratic.M+quadratic.Y==0))  |
    ((quadratic.R+quadratic.tt==2)&(quadratic.M+quadratic.Y==2))))



condition <- condition_all

datseeds <- c(sample(1:1e6, 1000), sample(1e6:2e6, 1000))

iseed <-1
cond <- 1



OneData <- function(iseed = 1, cond = 1){
  
  gen_data <- GenData(
    sim.seed = datseeds[iseed], 
    n = condition$n[cond],
    quadratic.R = condition$quadratic.R[cond],
    quadratic.tt = condition$quadratic.tt[cond],
    quadratic.M = condition$quadratic.M[cond],
    quadratic.Y = condition$quadratic.Y[cond],
    Mfamily = condition$Mfamily[cond], 
    Yfamily = condition$Yfamily[cond],
    num_c = condition$num_c[cond],
    corr_c = condition$corr_c[cond] 
  )
  data_in <- gen_data$datobs
  data_in$id <- 1:nrow(data_in)
  Yname <- grep("^Y", colnames(data_in), value = T)
  Rname <- grep("^R", colnames(data_in), value = T)
  ttname <- grep("^tt", colnames(data_in), value = T)
  Cnames <- grep("^Cdat", colnames(data_in), value = T)
  Mnames <- grep("^Mdat", colnames(data_in), value = T)

  
  # estimation -----
  
  # Study 1 with parametric estimation of the models
  if (condition$Mfamily[cond]!="gaussian" & condition$Yfamily[cond]!="gaussian") {
    mr_estimates <- MedMod(
      data = data_in,
      Yname, ttname, Rname, Mnames = Mnames, 
      Cnames = Cnames,
      estimator = c("onestep", "tml"), 
      nuisance_estimation = c("SL.ranger", "SL.gam", "SL.glm"),
      num_folds = 5
    )
    
    reg_estimates <- MedMod(
      data = data_in,
      Yname, ttname, Rname, Mnames = Mnames, 
      Cnames = Cnames,
      estimator = c("onestep", "tml"), 
      nuisance_estimation = c("SL.glm"),
      num_folds = 5
    )
  }
  
  # Study 2 with the regression-based estimation
  if (condition$Mfamily[cond]=="gaussian" & condition$Yfamily[cond]=="gaussian") {
    mr_estimates <- MedMod(
      data = data_in,
      Yname, ttname, Rname, Mnames = Mnames, 
      Cnames = Cnames,
      estimator = c("onestep", "tml"), 
      nuisance_estimation = c("SL.xgboost.modified", "SL.nnet"),
      num_folds = 3 # three-fold cross-fitting due to time consideration; different number of folds (e.g., 3, 4, or 5) led to similar results based on pilot simulations with 200 replications.
    )
    
    reg_estimates <- MedMod(
      data = data_in,
      Yname, ttname, Rname, Mnames = Mnames, 
      Cnames = Cnames,
      estimator = c("reg")
    )
  }
  
  
  
  # true values 
  gen_largeN <- GenData(sim.seed = datseeds[iseed], 
                        n = 1e5,
                        quadratic.R = condition$quadratic.R[cond],
                        quadratic.tt = condition$quadratic.tt[cond],
                        quadratic.M = condition$quadratic.M[cond],
                        quadratic.Y = condition$quadratic.Y[cond],
                        Mfamily = condition$Mfamily[cond], 
                        Yfamily = condition$Yfamily[cond],
                        num_c = condition$num_c[cond],
                        corr_c = condition$corr_c[cond] )
  
  true_values <- data.frame(Yt1r1 = gen_largeN$meanYt1r1, Yt1r0 = gen_largeN$meanYt1r0, Yt1r1.Mt1r0 = gen_largeN$meanYt1r1.Mt1r0, Yt0r1 = gen_largeN$meanYt0r1, Yt0r0 = gen_largeN$meanYt0r0, Yt0r1.Mt0r0 = gen_largeN$meanYt0r1.Mt0r0)
  rm(gen_largeN)
  true_effects <- effect(true_values)
  
  
  out <- data.frame(
    condition[cond, ],
    effect = names(true_effects),
    # true
    true_val = t(true_effects),
    # estimate
    mr_estimates[,-1],
    reg_estimates[,-1, drop=F],
    
    row.names = NULL
  ) 
  
  return(out)
}

