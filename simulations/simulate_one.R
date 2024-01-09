
# setwd("/work/08878/xliu19/ls6/meMO/simulation")
library(mvtnorm)
# library(cubature)
library(parallel)
library(tidyverse)


GenData <- function(
    sim.seed = 123, 
    n = 500, 
    Mfamily = "binomial", # "gaussian", 
    Yfamily = "binomial", # "gaussian",
    quadratic.R = FALSE,
    quadratic.M = FALSE,
    quadratic.Y = FALSE,
    num_c = 5
) { 
  set.seed(seed = sim.seed)
  n <- n
  # baseline covariates ----
  # num_c <- 5
  C_dat <- as.data.frame(rmvnorm(n, sigma = diag(1, nrow = num_c))) 
  
  # Subgroup indicator ----
  C2 <- model.matrix(~.-1 + ., C_dat)
  # if (quadratic.R == TRUE) {
  #   C2 <- (C_dat^2 - 1) / sqrt(2)
  # } 
  bc <- sqrt(0.26/(1-0.26)/ncol(C2)) / rep(1, ncol(C2))
  if (num_c > 5) {
    bc <- sqrt(0.26/(1-0.26)) / c(1: ncol(C2))^2
  }
  
  R <- 1*(as.matrix(C2) %*% bc + rnorm(n) > 0)
  
  # Treatment ------
  rc_dat <- model.matrix(~.-1 + ., data.frame(R=R, C_dat=C_dat))
  if (quadratic.R == TRUE) {
    C2 <- (C_dat^2 - 1) / sqrt(2)
    rc_dat <- model.matrix(~.-1 + ., data.frame(R=R, C_dat=C2))
  }
  
  br <- sqrt(0.13/(1-0.13)/ncol(rc_dat)) / rep(1, ncol(rc_dat))
  if (num_c > 5) {
    br <- sqrt(0.13/(1-0.13)) / c(1: ncol(rc_dat))^2
  }
  
  tt <- 1*(rc_dat %*% br + rnorm(n) > 0) 
  
  # Intermediate ----
  trc_dat <- model.matrix(~.-1 + ., data.frame(tt=tt, R=R, ttR = tt*R, C_dat=C_dat))
  if (quadratic.M == TRUE) {
    C2 <- (C_dat^2 - 1) / sqrt(2)
    trc_dat <- model.matrix(~.-1 + ., data.frame(tt=tt, R=R, ttR = tt*R, C_dat=C2))
  }
  btt <- c( -0.3, 0.5, 1, sqrt(0.26/(1-0.26)/ncol(C_dat))/rep(1, ncol(C_dat)) )
  if (num_c > 5) {
    btt <- c( -0.3, 0.5, 1, c(sqrt(0.26/(1-0.26))/c(1:ncol(C_dat)))^2 ) 
  }
  names(btt) <- colnames(trc_dat)
  
  # var(tts %*% btt)
  if (Mfamily == "gaussian") {
    M <- trc_dat %*% btt + rnorm(n)
  }
  if (Mfamily == "binomial") {
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
  
  if (num_c > 5) {
    bm <- c(c(0.6, -1.2, -1.4),  0.7, 0.8, 0.9, 1.1, c(sqrt(0.26/(1-0.26))/c(1:ncol(C_dat)))^2 )
  }
  
  names(bm) <- colnames(mtrc_dat)
  
  if (Yfamily == "gaussian") {
    Y <- mtrc_dat %*% bm + rnorm(n)
  }
  if (Yfamily == "binomial") {
    Y <- 1*(mtrc_dat %*% bm + rnorm(n) > 0) 
    #rbinom(n, 1, prob = plogis(mtrc_dat %*% bm))
  }
  
  datobs <- data.frame(Y=Y, Mdat=M, tt=tt, R=R, Cdat=C_dat)
  
  
  # True Values -----------------------------
  if ( (Yfamily == "gaussian") & (Mfamily == "gaussian") ) {
    # linear Y models
    meanC_dat <- t(c(0,0,0))
    # when tt = 1 ----
    # Yt1r1
    t1r1_dat <- model.matrix(~.-1 + ., data.frame(tt=1, R=1, ttR = 1, C_dat=meanC_dat))
    meanMt1r1 <- as.numeric(t1r1_dat %*% btt)
    meanMt1r1_dat <- model.matrix(~.-1 + ., data.frame(M=meanMt1r1, tt=1, R=1, ttM = meanMt1r1*1, ttR = 1, RM = 1*meanMt1r1, ttRM = 1*1*meanMt1r1, C_dat=meanC_dat))
    meanYt1r1 <- as.numeric(meanMt1r1_dat %*% bm)
    
    # Yt1r1.Mt1r0
    t1r0_dat <- model.matrix(~.-1 + ., data.frame(tt=1, R=0, ttR = 0, C_dat=meanC_dat))
    meanMt1r0 <- as.numeric(t1r0_dat %*% btt)
    meanMt1r0_dat <- model.matrix(~.-1 + ., data.frame(M=meanMt1r0, tt=1, R=1, ttM = meanMt1r0*1, ttR = 1, RM = 1*meanMt1r0, ttRM = 1*1*meanMt1r0, C_dat=meanC_dat))
    meanYt1r1.Mt1r0 <- as.numeric(meanMt1r0_dat %*% bm)
    
    # Yt1r0
    t1r0_dat <- model.matrix(~.-1 + ., data.frame(tt=1, R=0, ttR = 0, C_dat=meanC_dat))
    meanMt1r0 <- as.numeric(t1r0_dat %*% btt)
    meanMt1r0_dat <- model.matrix(~.-1 + ., data.frame(M=meanMt1r0, tt=1, R=0, ttM = meanMt1r0*1, ttR = 0, RM = 0*meanMt1r0, ttRM = 1*0*meanMt1r0, C_dat=meanC_dat))
    meanYt1r0 <- as.numeric(meanMt1r0_dat %*% bm)
    
    # when tt = 0 ----
    # Yt0r1
    t0r1_dat <- model.matrix(~.-1 + ., data.frame(tt=0, R=1, ttR = 0, C_dat=meanC_dat))
    meanMt0r1 <- as.numeric(t0r1_dat %*% btt)
    meanMt0r1_dat <- model.matrix(~.-1 + ., data.frame(M=meanMt0r1, tt=0, R=1, ttM = meanMt0r1*0, ttR = 0, RM = 1*meanMt0r1, ttRM = 0*1*meanMt0r1, C_dat=meanC_dat))
    meanYt0r1 <- as.numeric(meanMt0r1_dat %*% bm)
    
    # Yt0r1.Mt0r0
    t0r0_dat <- model.matrix(~.-1 + ., data.frame(tt=0, R=0, ttR = 0, C_dat=meanC_dat))
    meanMt0r0 <- as.numeric(t0r0_dat %*% btt)
    meanMt0r0_dat <- model.matrix(~.-1 + ., data.frame(M=meanMt0r0, tt=0, R=1, ttM = meanMt0r0*0, ttR = 0, RM = 1*meanMt0r0, ttRM = 0*1*meanMt0r0, C_dat=meanC_dat))
    meanYt0r1.Mt0r0 <- as.numeric(meanMt0r0_dat %*% bm)
    
    # Yt0r0
    t0r0_dat <- model.matrix(~.-1 + ., data.frame(tt=0, R=0, ttR = 0, C_dat=meanC_dat))
    meanMt0r0 <- as.numeric(t0r0_dat %*% btt)
    meanMt0r0_dat <- model.matrix(~.-1 + ., data.frame(M=meanMt0r0, tt=0, R=0, ttM = meanMt0r0*0, ttR = 0, RM = 0*meanMt0r0, ttRM = 0*0*meanMt0r0, C_dat=meanC_dat))
    meanYt0r0 <- as.numeric(meanMt0r0_dat %*% bm)
    
    
  }
  
  if ( (Mfamily == "binomial") ) {
    
    # when tt = 1 ----
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
    
    # when tt = 0 ----
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


# simdata_in <- dat.tt.oneX(n=500)
# aggregate(cbind(M, Y) ~ R, data = simdata_in, mean)



# Simulation----------------------------

condition1 <- data.frame(expand.grid(
  n = c(300, 500, 1000, 5000),
  # generating data
  quadratic.R = c(F, T),
  quadratic.M = c(F, T),
  quadratic.Y = c(F, T),
  Yfamily = c("binomial"), # c("gaussian"),
  num_c = c(20),
  Fit =  c("mlr")
)) 
adding1 <- data.frame(expand.grid(
  n = c(300, 500, 1000, 5000),
  # generating data
  quadratic.R = c(F, T),
  quadratic.M = c(F, T),
  quadratic.Y = c(F, T),
  Yfamily = c("binomial"),
  num_c = c(20),
  Fit =  c("linear") #c("lme")#
))
condition_all <- data.frame(do.call(rbind, list(condition1, adding1)))
condition <- condition_all[(condition_all$quadratic.R + condition_all$quadratic.M + condition_all$quadratic.Y) %in% c(1, 3, 0), ]
# condition <- condition[condition$fixed==T, ]


set.seed(12)
datseeds <- sample(1:1e6, 1000)

iseed <-1
cond <- 1



OneData <- function(iseed = 1, cond = 1){
  
  Yfamily <- as.character(condition$Yfamily[cond])
  num_c <- condition$num_c[cond]
  
  gen_data <- GenData(
    sim.seed = datseeds[iseed], 
    n = condition$n[cond],
    quadratic.R = condition$quadratic.R[cond],
    quadratic.M = condition$quadratic.M[cond],
    quadratic.Y = condition$quadratic.Y[cond],
    Mfamily = "binomial", 
    Yfamily = Yfamily, # "binomial",
    num_c = num_c
  )
  data_in <- gen_data$datobs
  data_in$id <- 1:nrow(data_in)
  data_in$obs_weights_Y <- rep(1, nrow(data_in))
  Cnames <- grep("Cdat", colnames(data_in), value = T)
  #colnames(select(data_in, starts_with("Cdat")))
  # Xnames <- colnames(select(data_in, starts_with("Xdat")))
  Mnames <- grep("Mdat", colnames(data_in), value = T)
  #colnames(select(data_in, starts_with("Mdat")))
  
  # add interactions
  data_in <- data_in %>% mutate(
    ttR = data_in$tt * data_in$R,
    ttM = data_in$tt*data_in[, Mnames],
    RM = data_in$R*data_in[, Mnames],
    ttRM = ttR*data_in[, Mnames]
  )
  data_in <- do.call(data.frame, data_in)
  ttMnames <- paste0("ttM.", Mnames)
  RMnames <- paste0("RM.", Mnames)
  ttRMnames <- paste0("ttRM.", Mnames)
  colnames(data_in)[(ncol(data_in)+1-3*length(Mnames)): ncol(data_in)] <- c(ttMnames, RMnames, ttRMnames)
  
  
  
  # estimators -----
  Fit <- condition$Fit[cond]
  
  
  crossfit_out <- crossfit.onestep(
    cv_folds = 5L,
    data_in = data_in, 
    Cnames = Cnames, 
    Mnames = Mnames,
    # Xnames = Xnames,
    Fit = Fit, # "linear",
    Yfamily =  Yfamily
  )
  
  tmle_out <- tmle.medMO(
    crossfit_onestep = crossfit_out,
    data_in = data_in,
    Cnames = Cnames,
    # Xnames = Xnames,
    Mnames = Mnames
  )
  
  
  
  # true values -----------------------------
  gen_largeN <- GenData(sim.seed = datseeds[iseed], n = 1e4,
                        quadratic.R = condition$quadratic.R[cond],
                        quadratic.M = condition$quadratic.M[cond],
                        quadratic.Y = condition$quadratic.Y[cond],
                        Mfamily = "binomial", 
                        Yfamily = Yfamily, # "binomial",
                        num_c = num_c )
  true_values <- c(Yt1r1 = gen_largeN$meanYt1r1, Yt1r0 = gen_largeN$meanYt1r0, Yt1r1.Mt1r0 = gen_largeN$meanYt1r1.Mt1r0, 
                   Yt0r1 = gen_largeN$meanYt0r1, Yt0r0 = gen_largeN$meanYt0r0, Yt0r1.Mt0r0 = gen_largeN$meanYt0r1.Mt0r0)
  true_values <- data.frame(t(true_values))
  true_effects <- within(true_values, {
    true_toD <- (Yt1r1 - Yt0r1) - (Yt1r0 - Yt0r0)
    true_meD <- (Yt1r1 - Yt1r1.Mt1r0) - (Yt0r1 - Yt0r1.Mt0r0)
    true_reD <- (Yt1r1.Mt1r0 - Yt1r0) - (Yt0r1.Mt0r0 - Yt0r0)
  })
  
  # crossfit_out$estimates / true_effects - 1
  # tmle_out$tml_estimates /true_effects-1
  
  out <- list(
    results = data.frame(
      condition[cond, ],
      effect = c("Yt1r1", "Yt1r0", "Yt1r1.Mt1r0", 
                 "Yt0r1", "Yt0r0", "Yt0r1.Mt0r0", 
                 "reD",  
                 "meD", 
                 "toD"),
      # true
      true_val = t(true_effects),
      # tml
      tml_estimate = tmle_out$tml_estimates,
      tml_interval = t(tmle_out$tml_intervals),
      # onestep
      onestep_estimate = crossfit_out$estimates,
      onestep_interval = t(crossfit_out$z_intervals),
      stderr = crossfit_out$stderrors
      # weighting-based
      , weighting_estimate = crossfit_out$weighting_estimates
    ) #, 
    # crossfit_out = crossfit_out, 
    # tmle_out = tmle_out
  )
  
  return(out)
}


# Parallel setup -----
num_reps <- 20 
mc_cores <- 20
seedseq <- seq(from=1, to=1000, by=num_reps)

jobconds <- c(1:nrow(condition))
