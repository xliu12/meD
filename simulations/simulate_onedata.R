
library(mvtnorm)
library(parallel)
library(tidyverse)


GenData <- function(
    sim.seed = 123,
    n = 500,
    Mfamily = "binomial",
    Yfamily = "binomial",
    quadratic.R = FALSE,
    quadratic.M = FALSE,
    quadratic.Y = FALSE,
    num_c = 20
) {
  set.seed(seed = sim.seed)
  n <- n
  # baseline covariates ----
  C_dat <- as.data.frame(rmvnorm(n, sigma = diag(1, nrow = num_c)))

  # Subgroup indicator ----
  C2 <- model.matrix(~.-1 + ., C_dat)
  if (quadratic.R == TRUE) {
    C2 <- (C_dat^2 - 1) / sqrt(2)
  }
  bc <- sqrt(0.35/ncol(C2)) / rep(1, ncol(C2))

  R <- 1*(as.matrix(C2) %*% bc + rnorm(n) > 0)

  # Treatment ------
  rc_dat <- model.matrix(~.-1 + ., data.frame(R=R, C_dat=C_dat))
  if (quadratic.R == TRUE) {
    C2 <- (C_dat^2 - 1) / sqrt(2)
    rc_dat <- model.matrix(~.-1 + ., data.frame(R=R, C_dat=C2))
  }

  br <- sqrt(0.15/ncol(rc_dat)) / rep(1, ncol(rc_dat))

  tt <- 1*(rc_dat %*% br + rnorm(n) > 0)

  # Mediator ----
  trc_dat <- model.matrix(~.-1 + ., data.frame(tt=tt, R=R, ttR = tt*R, C_dat=C_dat))
  if (quadratic.M == TRUE) {
    C2 <- (C_dat^2 - 1) / sqrt(2)
    trc_dat <- model.matrix(~.-1 + ., data.frame(tt=tt, R=R, ttR = tt*R, C_dat=C2))
  }
  btt <- c( -0.3, 0.5, 1, sqrt(0.35/ncol(C_dat))/rep(1, ncol(C_dat)) )

  names(btt) <- colnames(trc_dat)

  M <- 1*(trc_dat %*% btt + rnorm(n) > 0)

  # Outcome -----
  mtrc_dat <- model.matrix(~.-1 + ., data.frame(M=M, tt=tt, R=R, ttM = M*tt, ttR = tt*R, RM = R*M, ttRM = tt*R*M, C_dat=C_dat))
  if (quadratic.Y == TRUE) {
    C2 <- (C_dat^2 - 1) / sqrt(2)
    mtrc_dat <- model.matrix(~.-1 + ., data.frame(M=M, tt=tt, R=R, ttM = M*tt, ttR = tt*R, RM = R*M, ttRM = tt*R*M, C_dat=C2))
  }
  bm <- c(c(0.6, -1.2, -1.4),  0.7, 0.8, 0.9, 1.1, sqrt(0.35/ncol(C_dat)) / rep(1, ncol(C_dat)))

  names(bm) <- colnames(mtrc_dat)
  Y <- 1*(mtrc_dat %*% bm + rnorm(n) > 0)

  datobs <- data.frame(Y=Y, Mdat=M, tt=tt, R=R, Cdat=C_dat)


  # True Values -----------------------------

  if ( (Mfamily == "binomial") ) {

    # when tt = 1 ----
    # Yt1r1
    t1r1_dat <- data.frame(trc_dat); t1r1_dat$tt <- 1;  t1r1_dat$R <- 1
    t1r1_dat$ttR <- t1r1_dat$tt * t1r1_dat$R
    pMt1r1 <- pnorm(as.matrix(t1r1_dat) %*% btt, lower.tail = T)
    m1t1r1_dat <- data.frame(mtrc_dat); m1t1r1_dat$M <- 1; m1t1r1_dat$tt <- 1; m1t1r1_dat$R <- 1
    m1t1r1_dat$ttM <- m1t1r1_dat$tt * m1t1r1_dat$M
    m1t1r1_dat$ttR <- m1t1r1_dat$tt * m1t1r1_dat$R
    m1t1r1_dat$RM <- m1t1r1_dat$R * m1t1r1_dat$M
    m1t1r1_dat$ttRM <- m1t1r1_dat$tt * m1t1r1_dat$R * m1t1r1_dat$M

    m0t1r1_dat <- data.frame(mtrc_dat); m0t1r1_dat$M <- 0; m0t1r1_dat$tt <- 1; m0t1r1_dat$R <- 1
    m0t1r1_dat$ttM <- m0t1r1_dat$tt * m0t1r1_dat$M
    m0t1r1_dat$ttR <- m0t1r1_dat$tt * m0t1r1_dat$R
    m0t1r1_dat$RM <- m0t1r1_dat$R * m0t1r1_dat$M
    m0t1r1_dat$ttRM <- m0t1r1_dat$tt * m0t1r1_dat$R * m0t1r1_dat$M

    if (Yfamily == "binomial") {
      Yt1r1 <- pnorm(as.matrix(m1t1r1_dat) %*% bm) * pMt1r1 + pnorm(as.matrix(m0t1r1_dat) %*% bm) * (1-pMt1r1)
    }

    meanYt1r1 <- mean(Yt1r1)

    # Yt1r1.Mt1r0
    t1r0_dat <- data.frame(trc_dat); t1r0_dat$tt <- 1;  t1r0_dat$R <- 0
    t1r0_dat$ttR <- t1r0_dat$tt * t1r0_dat$R
    pMt1r0 <- pnorm(as.matrix(t1r0_dat) %*% btt)
    if (Yfamily == "binomial") {
      Yt1r1.Mt1r0 <- pnorm(as.matrix(m1t1r1_dat) %*% bm) * pMt1r0 + pnorm(as.matrix(m0t1r1_dat) %*% bm) * (1-pMt1r0)
    }

    meanYt1r1.Mt1r0 <- mean(Yt1r1.Mt1r0)

    # Yt1r0
    m1t1r0_dat <- data.frame(mtrc_dat);
    m1t1r0_dat$M <- 1; m1t1r0_dat$tt <- 1; m1t1r0_dat$R <- 0
    m1t1r0_dat$ttM <- m1t1r0_dat$tt * m1t1r0_dat$M
    m1t1r0_dat$ttR <- m1t1r0_dat$tt * m1t1r0_dat$R
    m1t1r0_dat$RM <- m1t1r0_dat$R * m1t1r0_dat$M
    m1t1r0_dat$ttRM <- m1t1r0_dat$tt * m1t1r0_dat$R * m1t1r0_dat$M

    m0t1r0_dat <- data.frame(mtrc_dat);
    m0t1r0_dat$M <- 0; m0t1r0_dat$tt <- 1; m0t1r0_dat$R <- 0
    m0t1r0_dat$ttM <- m0t1r0_dat$tt * m0t1r0_dat$M
    m0t1r0_dat$ttR <- m0t1r0_dat$tt * m0t1r0_dat$R
    m0t1r0_dat$RM <- m0t1r0_dat$R * m0t1r0_dat$M
    m0t1r0_dat$ttRM <- m0t1r0_dat$tt * m0t1r0_dat$R * m0t1r0_dat$M

    if (Yfamily == "binomial") {
      Yt1r0 <- pnorm(as.matrix(m1t1r0_dat) %*% bm) * pMt1r0 + pnorm(as.matrix(m0t1r0_dat) %*% bm) * (1-pMt1r0)
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

    meanYt0r1 <- mean(Yt0r1)

    # Yt0r1.Mt0r0
    t0r0_dat <- data.frame(trc_dat);
    t0r0_dat$tt <- 0;  t0r0_dat$R <- 0
    t0r0_dat$ttR <- t0r0_dat$tt * t0r0_dat$R

    pMt0r0 <- pnorm(as.matrix(t0r0_dat) %*% btt, lower.tail = T)
    if (Yfamily == "binomial") {
      Yt0r1.Mt0r0 <- pnorm(as.matrix(m1t0r1_dat) %*% bm) * pMt0r0 + pnorm(as.matrix(m0t0r1_dat) %*% bm) * (1-pMt0r0)
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

    meanYt0r0 <- mean(Yt0r0)

  }
  out <- mget( ls(), envir = environment() )

  return(out)
}


# Simulation ----------------------------

condition1 <- data.frame(expand.grid(
  n = c(300, 500, 1000),
  # generating data
  quadratic.R = c(F, T),
  quadratic.M = c(F, T),
  quadratic.Y = c(F, T),
  Yfamily = c("binomial"),
  num_c = c(20),
  Fit =  c("linear")
))
adding1 <- data.frame(expand.grid(
  n = c(300, 500, 1000),
  # generating data
  quadratic.R = c(F, T),
  quadratic.M = c(F, T),
  quadratic.Y = c(F, T),
  Yfamily = c("binomial"),
  num_c = c(20),
  Fit =  c("mlr")
))
condition_all <- data.frame(do.call(rbind, list(condition1, adding1)))
condition <- condition_all[(condition_all$quadratic.R + condition_all$quadratic.M + condition_all$quadratic.Y) %in% c(1, 3, 0), ]


set.seed(12)
datseeds <- sample(1:1e6, 1000)

iseed <-1
cond <- 1

devtools::load_all("pkg/meD")

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
    Yfamily = Yfamily,
    num_c = num_c
  )
  data_in <- gen_data$datobs
  Cnames <- grep("Cdat", colnames(data_in), value = T)
  Mnames <- grep("Mdat", colnames(data_in), value = T)

  # estimators -----

  if (condition$Fit[cond] == "mlr") {
    nuisance_estimation <- c("SL.ranger", "SL.glm", "SL.gam")
  }
  if (condition$Fit[cond] == "linear") {
    nuisance_estimation <- c("SL.glm")
  }

  estimates <- meDreD(data_in, Mnames, Cnames, nuisance_estimation, Yfamily = "binomial", estimator = c("onestep", "tml"))

  # true values -----------------------------
  gen_largeN <- GenData(sim.seed = datseeds[iseed], n = 1e5,
                        quadratic.R = condition$quadratic.R[cond],
                        quadratic.M = condition$quadratic.M[cond],
                        quadratic.Y = condition$quadratic.Y[cond],
                        Mfamily = "binomial",
                        Yfamily = Yfamily,
                        num_c = num_c )
  true_values <- c(Yt1r1 = gen_largeN$meanYt1r1, Yt1r0 = gen_largeN$meanYt1r0, Yt1r1.Mt1r0 = gen_largeN$meanYt1r1.Mt1r0,
                   Yt0r1 = gen_largeN$meanYt0r1, Yt0r0 = gen_largeN$meanYt0r0, Yt0r1.Mt0r0 = gen_largeN$meanYt0r1.Mt0r0)
  true_values <- data.frame(t(true_values))
  true_effects <- within(true_values, {
    true_toD <- (Yt1r1 - Yt0r1) - (Yt1r0 - Yt0r0)
    true_meD <- (Yt1r1 - Yt1r1.Mt1r0) - (Yt0r1 - Yt0r1.Mt0r0)
    true_reD <- (Yt1r1.Mt1r0 - Yt1r0) - (Yt0r1.Mt0r0 - Yt0r0)
  })

  out <- list(estimates = estimates, true_effects = true_effects)

  return(out)
}


# Parallel setup -----
num_reps <- 20
mc_cores <- 20
seedseq <- seq(from=1, to=1000, by=num_reps)

jobconds <- c(1:nrow(condition))
