# library(mlr3)
# library(mlr3learners)
library(SuperLearner)
library(ranger)
library(glmnet)
# library(gbm)
library(xgboost)
library(hal9001)
#library(bartMachine)
# library(origami)
library(gam)
library(earth)
#  ZINB model: Zero-inflated and Overdispersed Data https://cran.r-project.org/web/packages/mpath/vignettes/static_german.pdf
# library(pscl)
# library(zic)
# library(mpath) #  zipath(): Fit zero-inflated count data linear model with lasso  (or elastic net), snet or mnet regularization
# library(flexsurv)
# # devtools::install_github("wuziyueemory/twostageSL") 
# library(twostageSL) # super learner for zero inflated outcome
# library(earth) # Multivariate adaptive regression splines (MARS) [https://bradleyboehmke.github.io/HOML/mars.html]
# listWrappers()
# SL.gbm; For maximum accuracy one might try at least the following models: glmnet, randomForest, XGBoost, SVM, and bartMachine. 
# We specify family = binomial() because we are predicting a binary outcome, aka classification. With a continuous outcome we would specify family = gaussian().
# library(data.table)
# library(tidyverse)

SL.hal.modified <- function(...) {
  SL.hal9001(..., 
             # tuning
             max_degree = 1, smoothness_orders = 1, 
             num_knots = c(25, 15, 10), #num_knots =c(25, 10, 5),
             # reduce_basis = NULL,
             X_unpenalized = NULL
  )
}

bound_precision <- function(vals, tol = 1e-5) {
  vals[vals < tol] <- tol
  vals[vals > 1 - tol] <- 1 - tol
  return(vals)
}
bound_propensity <- function(vals, bounds = c(0.01, 0.99)) {
  # assertthat::assert_that(!(max(vals) > 1 || min(vals) < 0))
  vals[vals < bounds[1]] <- bounds[1]
  vals[vals > bounds[2]] <- bounds[2]
  return(vals)
}


fitting.R <- function(train_data, valid_data, 
                      rmodel = "R ~ tt + C",
                      Cnames, Mnames,
                      # Znames = Znames,
                      # Xnames,
                      SL_library = c("SL.glm")
) {
  RMnames <- paste0("RM.", Mnames)
  ttMnames <- paste0("ttM.", Mnames)
  ttRMnames <- paste0("ttRM.", Mnames)
  
  
  if (rmodel == "R ~ M + tt + C") {
    cov_names <- c(Mnames, "tt", Cnames)
  } 
  if (rmodel == "R ~ M * tt + C") {
    cov_names <- c(Mnames, "tt", ttMnames, Cnames)
  } 
  if (rmodel == "R ~ tt + C") {
    cov_names <- c("tt", Cnames)
  }
  
  set.seed(9999)
  
  sl_fit <- SuperLearner(
    Y = train_data$R, 
    X = train_data[, cov_names], 
    # obsWeights = train_data[, obs_weights],
    family = binomial(), 
    SL.library = SL_library
  )
  sl_pred_train <- predict(sl_fit, train_data[, cov_names], onlySL = TRUE)$pred
  sl_pred_valid <- predict(sl_fit, valid_data[, cov_names], onlySL = TRUE)$pred
  
  
  out <- list(
    r1_pred_train = sl_pred_train,
    r1_pred_valid = sl_pred_valid,
    r1_fit = sl_fit
  )
  
  return(out)
}


fitting.tt <- function(train_data, valid_data, 
                       tmodel = "tt ~ C", 
                       Cnames, Mnames,
                       # Znames,
                       Xnames,
                       SL_library = c("SL.glm")
) {
  
  RMnames <- paste0("RM.", Mnames)
  ttMnames <- paste0("ttM.", Mnames)
  ttRMnames <- paste0("ttRM.", Mnames)
  
  
  if (tmodel == "tt ~ C") {
    cov_names <- c(Cnames)
  }
  
  set.seed(9999)
  
  
  sl_fit <- SuperLearner(
    Y = train_data$tt, 
    X = train_data[, cov_names], 
    # obsWeights = train_data[, obs_weights],
    family = binomial(), 
    
    SL.library = SL_library
  )
  sl_pred_train <- predict(sl_fit, train_data[, cov_names], onlySL = TRUE)$pred
  sl_pred_valid <- predict(sl_fit, valid_data[, cov_names], onlySL = TRUE)$pred
  
  
  out <- list(
    t1_pred_train = sl_pred_train,
    t1_pred_valid = sl_pred_valid,
    t1_fit = sl_fit
  )
  
  return(out)
}



fitting.M <- function(train_data, valid_data, 
                      mmodel = "M ~ tt * R + C", 
                      Cnames, Mnames,
                      # Znames,
                      # Xnames,
                      SL_library = c("SL.glm")
                      
) {
  
  RMnames <- paste0("RM.", Mnames)
  ttMnames <- paste0("ttM.", Mnames)
  ttRMnames <- paste0("ttRM.", Mnames)
  
  if (mmodel == "M ~ tt + R + C") {
    cov_names <- c("tt", "R", Cnames)
  }
  if (mmodel == "M ~ tt * R + C"){
    cov_names <- c("tt", "R", "ttR", Cnames)
  }
  
  set.seed(9999)
  Mfamily <- ifelse(length(unique(train_data[, Mnames]))==2, "binomial", "gaussian")
  sl_fit <- SuperLearner(
    Y = train_data[, Mnames], 
    X = train_data[, cov_names], 
    # obsWeights = train_data[, obs_weights],
    family = Mfamily, 
    SL.library = SL_library
  )
  sl_pred_train <- predict(sl_fit, train_data[, cov_names], onlySL = TRUE)$pred
  sl_pred_valid <- predict(sl_fit, valid_data[, cov_names], onlySL = TRUE)$pred
  
  out <- list(
    m1_pred_train = sl_pred_train,
    m1_pred_valid = sl_pred_valid,
    m1_fit = sl_fit
  )
  
  return(out)
}


fitting.Y <- function(train_data, valid_data, 
                      ymodel = "Y ~ M * tt * R + C", 
                      Cnames, 
                      Mnames,
                      # Znames = Znames,
                      # Xnames, 
                      SL_library = c("SL.glm"), 
                      Yfamily = "gaussian"
                      # Yfamily = c(familyzero = "binomial", familypositive = "gaussian"),
                      # SL_library =list(stage1=c("SL.glm","SL.glmnet", "SL.ranger"), stage2=c("SL.glm","SL.glmnet", "SL.ranger","SL.coxph"))
                      # Yfamily = "zeroinfl", SL_library = "SL.zeroinfl" or SL_library == "SL.zinflasso"
) {
  
  RMnames <- paste0("RM.", Mnames)
  ttMnames <- paste0("ttM.", Mnames)
  ttRMnames <- paste0("ttRM.", Mnames)
  
  if (ymodel == "Y ~ M * tt * R + C"){
    cov_names <- c(Mnames, "tt", "R", "ttR", ttMnames, RMnames, ttRMnames, 
                   Cnames)
  } 
  if (ymodel == "Y ~ M + tt + R + C"){
    cov_names <- c(Mnames, "tt", "R", Cnames)
    
  } 
  
  if (ymodel == "Y ~ tt * R + C"){
    cov_names <- c("tt", "R", "ttR", Cnames)
  } 
  if (ymodel == "Y ~ tt + R + C"){
    cov_names <- c("tt", "R", Cnames)
  } 
  
  set.seed(9999)
  
  if (length(Yfamily) == 2) {
    # use the package "twostageSL" (separate fitting the zero component and the positive component): stage-1 for fitting Pr(Y>0 | X); stage-2 for fitting the conditional mean given non-zero, E[Y | Y>0, X]
    sl_fit <- twostageSL(
      Y = train_data$Y, 
      X = train_data[, cov_names],
      twostage = TRUE,
      # library.1stage = SL_library,# library.1stage	Candidate prediction algorithms in standard super learner.
      family.1 = binomial, family.2 = Yfamily[2], # gaussian,
      library.2stage = SL_library,
      library.1stage = SL_library[[2]], family.single = "gaussian"
    )
    sl_pred_train <- predict(sl_fit, train_data[, cov_names], onlySL = TRUE)$pred
    sl_pred_valid <- predict(sl_fit, valid_data[, cov_names], onlySL = TRUE)$pred
    
  }
  
  # use the package "mpath" (zero-inflated linear model with lasso) or "pscl" (Zero-inflated Count Data Regression)
  if (Yfamily == "zeroinfl") {
    dat <- cbind(
      Y = train_data$Y,
      X = train_data[, cov_names]
    )
    dat <- data.frame(dat)
    sl_fit <- zeroinfl(Y ~.|., data=dat, dist="poisson") # "negbin"
    
    sl_pred_train <- predict(sl_fit, train_data[, cov_names], type = "response")
    sl_pred_valid <- predict(sl_fit, valid_data[, cov_names], type = "response")
  }
  if (Yfamily == "zinflasso") {
    dat <- cbind(
      Y = train_data$Y,
      X = train_data[, cov_names]
    )
    dat <- data.frame(dat)
    sl_fit <- zipath(Y~.|.,data = dat, family = "poisson", nlambda=5)
    cvsl_fit <- cv.zipath(Y~.|.,data = dat, family = "poisson", nlambda=5, plot=F)
    
    sl_pred_train <- predict(sl_fit, newdata=train_data[, cov_names], type = "response", which=cvsl_fit$lambda.which)
    sl_pred_valid <- predict(sl_fit, valid_data[, cov_names], type = "response", which=cvsl_fit$lambda.which)
  }
  
  if ((Yfamily %in% c("gaussian", "binomial")) & (length(Yfamily) == 1)) {
    sl_fit <- SuperLearner(
      Y = train_data$Y, 
      X = train_data[, cov_names], 
      # newX = valid_data[, cov_names],
      family = Yfamily, 
      # obsWeights = train_data[, obs_weights],
      
      SL.library = SL_library
    )
    sl_pred_train <- predict(sl_fit, train_data[, cov_names], onlySL = TRUE)$pred
    sl_pred_valid <- predict(sl_fit, valid_data[, cov_names], onlySL = TRUE)$pred
    
  }
  
  
  out <- list(
    y_pred_train = sl_pred_train,
    y_pred_valid = sl_pred_valid,
    y_fit = sl_fit
  )
  
  return(out)
}



fitting.u_pseudo <- function(train_data, valid_data, 
                             umodel = "u_pseudo ~ X + t1 + r1 + C", 
                             Cnames, Mnames, Xnames,
                             SL_library = c("SL.glm"),
                             y_out.mtrxc, 
                             x_out.mtrc, 
                             x_out.trc, 
                             t_out.rc,
                             t_out.mrc,
                             r_out.c, 
                             r_out.mc
) {
  
  if (umodel == "u_pseudo ~ X + t1 + r1 + C"){
    cov_names <- c(Xnames, Cnames)
  } 
  
  set.seed(9999)
  
  # calculating psedo outcome
  train_data_r1 <- train_data; train_data_r1$R <- 1; train_data_r1$tt <- 1
  valid_data_r1 <- valid_data; valid_data_r1$R <- 1; valid_data_r1$tt <- 1
  
  y_pred_train_r1 <- predict(y_out.mtrxc$y_fit, train_data_r1[, y_out.mtrxc$y_fit$varNames], onlySL = TRUE)$pred
  y_pred_valid_r1 <- predict(y_out.mtrxc$y_fit, valid_data_r1[, y_out.mtrxc$y_fit$varNames], onlySL = TRUE)$pred
  
  train_data_r0 <- train_data; train_data_r0$R <- 0; train_data_r0$tt <- 1
  valid_data_r0 <- valid_data; valid_data_r0$R <- 0; valid_data_r0$tt <- 1
  
  # w(r,t,x,m,c) -----
  # , for training ----
  # for x
  x1.mtrc_pred_train_r1 <- predict(x_out.mtrc$x1_fit, train_data_r1[, x_out.mtrc$x1_fit$varNames])$pred
  
  x1.trc_pred_train_r1 <- predict(x_out.trc$x1_fit, train_data_r1[, x_out.trc$x1_fit$varNames])$pred
  
  # for r, given t=1
  r1.c_pred_train <- (r_out.c$r1_pred_train)
  r1.mc_pred_train <- (predict(r_out.mc$r1_fit, train_data_r1[, r_out.mc$r1_fit$varNames])$pred )
  
  # for tt
  t1.rc_pred_train_r0 <- predict(t_out.rc$t1_fit, train_data_r0[, t_out.rc$t1_fit$varNames])$pred
  t1.rc_pred_train_r1 <- predict(t_out.rc$t1_fit, train_data_r1[, t_out.rc$t1_fit$varNames])$pred
  t1.mrc_pred_train_r0 <- predict(t_out.mrc$t1_fit, train_data_r0[, t_out.mrc$t1_fit$varNames])$pred
  t1.mrc_pred_train_r1 <- predict(t_out.mrc$t1_fit, train_data_r1[, t_out.mrc$t1_fit$varNames])$pred
  
  # odds_rt_train <- with(train_data_r1, {
  #   (r1.c_pred_train/(1-r1.c_pred_train)) * 
  #     ((1-r1.mc_pred_train)/(r1.mc_pred_train)) *
  #     (t1.rc_pred_train_r1/t1.rc_pred_train_r0) * 
  #     (t1.mrc_pred_train_r0/t1.mrc_pred_train_r1)
  # })
  # calculate w
  
  w_mtrxc_train <- with(train_data_r1, {
    
    den <- (1-r1.c_pred_train)*(r1.mc_pred_train)*t1.rc_pred_train_r0*t1.mrc_pred_train_r1 * ( x1.mtrc_pred_train_r1*as.numeric(train_data_r1[, Xnames]==1) + (1 - x1.mtrc_pred_train_r1)*as.numeric(train_data_r1[, Xnames]==0) )
    den <- bound_propensity(den)
    
    nume <- r1.c_pred_train * (1-r1.mc_pred_train) * t1.rc_pred_train_r1 * t1.mrc_pred_train_r0 * ( x1.trc_pred_train_r1*as.numeric(train_data_r1[, Xnames]==1) + (1 - x1.trc_pred_train_r1)*as.numeric(train_data_r1[, Xnames]==0) )
    
    nume / den
    # odds_rt_train * I(train_data_r1[, Xnames]==1)* 
    #   (x1.trc_pred_train_r1 / x1.mtrc_pred_train_r1) + 
    #   odds_rt_train * I(train_data_r1[, Xnames]==0)* 
    #   ((1 - x1.trc_pred_train_r1) / (1 - x1.mtrc_pred_train_r1)) 
  })
  
  
  
  u_pseudo_train <- y_pred_train_r1 * w_mtrxc_train
  
  # , for valid -----
  # for x
  x1.mtrc_pred_valid_r1 <- predict(x_out.mtrc$x1_fit, valid_data_r1[, x_out.mtrc$x1_fit$varNames])$pred
  
  x1.trc_pred_valid_r1 <- predict(x_out.trc$x1_fit, valid_data_r1[, x_out.trc$x1_fit$varNames])$pred
  
  # for r, given t=1
  r1.c_pred_valid <- r_out.c$r1_pred_valid
  r1.mc_pred_valid <- predict(r_out.mc$r1_fit, valid_data_r1[, r_out.mc$r1_fit$varNames])$pred 
  
  # for tt
  t1.rc_pred_valid_r0 <- predict(t_out.rc$t1_fit, valid_data_r0[, t_out.rc$t1_fit$varNames])$pred
  t1.rc_pred_valid_r1 <- predict(t_out.rc$t1_fit, valid_data_r1[, t_out.rc$t1_fit$varNames])$pred
  t1.mrc_pred_valid_r0 <- predict(t_out.mrc$t1_fit, valid_data_r0[, t_out.mrc$t1_fit$varNames])$pred
  t1.mrc_pred_valid_r1 <- predict(t_out.mrc$t1_fit, valid_data_r1[, t_out.mrc$t1_fit$varNames])$pred
  
  # odds_rt_valid <- with(valid_data_r1, {
  #   (r1.c_pred_valid/(1-r1.c_pred_valid)) * 
  #     ((1-r1.mc_pred_valid)/(r1.mc_pred_valid)) *
  #     (t1.rc_pred_valid_r1/t1.rc_pred_valid_r0) * 
  #     (t1.mrc_pred_valid_r0/t1.mrc_pred_valid_r1)
  # })
  # calculate w
  w_mtrxc_valid <- with(valid_data_r1, {
    
    den <- (1-r1.c_pred_valid)*(r1.mc_pred_valid)*t1.rc_pred_valid_r0*t1.mrc_pred_valid_r1 * ( x1.mtrc_pred_valid_r1*as.numeric(valid_data_r1[, Xnames]==1) + (1 - x1.mtrc_pred_valid_r1)*as.numeric(valid_data_r1[, Xnames]==0) )
    den <- bound_propensity(den)
    
    nume <- r1.c_pred_valid * (1-r1.mc_pred_valid) * t1.rc_pred_valid_r1 * t1.mrc_pred_valid_r0 * ( x1.trc_pred_valid_r1*as.numeric(valid_data_r1[, Xnames]==1) + (1 - x1.trc_pred_valid_r1)*as.numeric(valid_data_r1[, Xnames]==0) )
    
    nume / den
    # odds_rt_valid * I(valid_data_r1[, Xnames]==1)* 
    #   (x1.trc_pred_valid_r1 / x1.mtrc_pred_valid_r1) + 
    #   odds_rt_valid * I(valid_data_r1[, Xnames]==0)* 
    #   ((1 - x1.trc_pred_valid_r1) / (1 - x1.mtrc_pred_valid_r1)) 
  })
  
  u_pseudo_valid <- y_pred_valid_r1 * w_mtrxc_valid
  
  
  # fit u_pseudo ----------------------------
  
  u_train_data <- train_data # should here only be the subsample with R==1?
  u_train_data$u_pseudo_train <- u_pseudo_train 
  
  u_train_data <- u_train_data[u_train_data$R==1 & u_train_data$tt==1, ]
  
  
  sl_fit <- SuperLearner(
    Y = u_train_data$u_pseudo_train, 
    X = u_train_data[, cov_names], 
    family = "gaussian", 
    # obsWeights = u_train_data[, obs_weights],
    SL.library = SL_library
  )
  sl_pred_train <- predict(sl_fit, u_train_data[, cov_names], onlySL = TRUE)$pred
  
  u_valid_data <- valid_data
  u_valid_data$u_pseudo_valid <- u_pseudo_valid
  
  sl_pred_valid <- predict(sl_fit, u_valid_data[, cov_names], onlySL = TRUE)$pred
  
  out <- list(
    y_pred_valid_r1 = y_pred_valid_r1,
    u_pseudo_valid = u_pseudo_valid,
    w_mtrxc_valid = w_mtrxc_valid,
    
    u_pred_train = sl_pred_train,
    u_pred_valid = sl_pred_valid,
    u_fit = sl_fit
  )
  
  return(out)
}


fitting.v_pseudo <- function(train_data, valid_data, 
                             vmodel = "v_pseudo ~ t1 + r0 + C", 
                             alltrain = FALSE,
                             Cnames, Mnames, # Xnames,
                             SL_library = c("SL.glm"),
                             y_out.mtrc, 
                             # x_out.mtrc, 
                             # x_out.trc, 
                             t_out.c,
                             r_out.tc, 
                             r_out.mtc,
                             Yfamily 
) {
  RMnames <- paste0("RM.", Mnames)
  ttMnames <- paste0("ttM.", Mnames)
  ttRMnames <- paste0("ttRM.", Mnames)
  
  
  cov_names <- c(Cnames)
  
  
  set.seed(9999)
  
  RMnames <- paste0("RM.", Mnames)
  ttMnames <- paste0("ttM.", Mnames)
  ttRMnames <- paste0("ttRM.", Mnames)
  
  # when tt = 1 -----
  if (vmodel %in% c("v_pseudo ~ t1 + r0 + C", "y_M ~ t1 + r0 + C", "y_M ~ t1 + r1 + C") ) { 
    
    train_data_r0 <- train_data; train_data_r0$R <- 0; train_data_r0$tt <- 1
    train_data_r0[, "ttR"] <- with(train_data_r0, {tt*R})
    train_data_r0[, RMnames] <- train_data_r0$R * train_data_r0[, Mnames]
    train_data_r0[, ttMnames] <- train_data_r0$tt * train_data_r0[, Mnames]
    train_data_r0[, ttRMnames] <- train_data_r0$tt * train_data_r0$R * train_data_r0[, Mnames]
    
    train_data_r1 <- train_data; train_data_r1$tt <- 1; train_data_r1$R <- 1
    train_data_r1[, "ttR"] <- with(train_data_r1, {tt*R})
    train_data_r1[, RMnames] <- train_data_r1$R * train_data_r1[, Mnames]
    train_data_r1[, ttMnames] <- train_data_r1$tt * train_data_r1[, Mnames]
    train_data_r1[, ttRMnames] <- train_data_r1$tt * train_data_r1$R * train_data_r1[, Mnames]
    
    valid_data_r0 <- valid_data; valid_data_r0$tt <- 1; valid_data_r0$R <- 0
    valid_data_r0[, "ttR"] <- with(valid_data_r0, {tt*R})
    valid_data_r0[, RMnames] <- valid_data_r0$R * valid_data_r0[, Mnames]
    valid_data_r0[, ttMnames] <- valid_data_r0$tt * valid_data_r0[, Mnames]
    valid_data_r0[, ttRMnames] <- valid_data_r0$tt * valid_data_r0$R * valid_data_r0[, Mnames]
    
    valid_data_r1 <- valid_data; valid_data_r1$tt <- 1; valid_data_r1$R <- 1
    valid_data_r1[, "ttR"] <- with(valid_data_r1, {tt*R})
    valid_data_r1[, RMnames] <- valid_data_r1$R * valid_data_r1[, Mnames]
    valid_data_r1[, ttMnames] <- valid_data_r1$tt * valid_data_r1[, Mnames]
    valid_data_r1[, ttRMnames] <- valid_data_r1$tt * valid_data_r1$R * valid_data_r1[, Mnames]
    
    
    # creating outcome and training data 
    
    # for v_pseudo -------
    if (vmodel == "v_pseudo ~ t1 + r0 + C") {
      # for training 
      
      y_pred_train_r1 <- predict(y_out.mtrc$y_fit, train_data_r1[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_train_r1 <- bound_precision(y_pred_train_r1)
      }
      v_pseudo_train <- y_pred_train_r1
      
      #for valid 
      
      y_pred_valid_r1 <- predict(y_out.mtrc$y_fit, valid_data_r1[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_valid_r1 <- bound_precision(y_pred_valid_r1)
      }
      
      v_pseudo_valid <- y_pred_valid_r1
      
      # fit v_pseudo
      v_train_data <- train_data 
      v_train_data$v_pseudo_train <- v_pseudo_train
      if (alltrain == FALSE) {
        # here only use the subsample with R==0, differing from the rvalue for pseudo outcome 
        v_train_data <- v_train_data[v_train_data$R==0 & v_train_data$tt==1, ]
      }
    }
    
    
    # for E[Y | tt=1, r=0] -----
    if (vmodel == "y_M ~ t1 + r0 + C") {
      # for training 
      
      y_pred_train_r0 <- predict(y_out.mtrc$y_fit, train_data_r0[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_train_r0 <- bound_precision(y_pred_train_r0)
      }
      v_pseudo_train <- y_pred_train_r0
      
      #for valid 
      
      y_pred_valid_r0 <- predict(y_out.mtrc$y_fit, valid_data_r0[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_valid_r0 <- bound_precision(y_pred_valid_r0)
      }
      
      v_pseudo_valid <- y_pred_valid_r0
      
      # fit v_pseudo
      v_train_data <- train_data 
      v_train_data$v_pseudo_train <- v_pseudo_train
      # here only use the subsample with R==0
      if (alltrain == FALSE) {
        v_train_data <- v_train_data[v_train_data$R==0 & v_train_data$tt==1, ]
      }
      
    }
    
    
    # for E[Y | tt=1, r=1] -------
    if (vmodel == "y_M ~ t1 + r1 + C") {
      # for training 
      
      y_pred_train_r1 <- predict(y_out.mtrc$y_fit, train_data_r1[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_train_r1 <- bound_precision(y_pred_train_r1)
      }
      v_pseudo_train <- y_pred_train_r1
      
      #for valid 
      
      y_pred_valid_r1 <- predict(y_out.mtrc$y_fit, valid_data_r1[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_valid_r1 <- bound_precision(y_pred_valid_r1)
      }
      
      v_pseudo_valid <- y_pred_valid_r1
      
      # fit v_pseudo
      v_train_data <- train_data 
      v_train_data$v_pseudo_train <- v_pseudo_train
      # here only use the subsample with R==1, the same as the rvalue for pseudo outcome 
      if (alltrain == FALSE) {
        v_train_data <- v_train_data[v_train_data$R==1 & v_train_data$tt==1, ]
      }
    }
    
  }
  
  #  when tt = 0 -----
  if (vmodel %in% c("v_pseudo ~ t0 + r0 + C", "y_M ~ t0 + r0 + C", "y_M ~ t0 + r1 + C") ) { 
    
    train_data_r0 <- train_data; train_data_r0$R <- 0; train_data_r0$tt <- 0
    train_data_r0[, "ttR"] <- with(train_data_r0, {tt*R})
    train_data_r0[, RMnames] <- train_data_r0$R * train_data_r0[, Mnames]
    train_data_r0[, ttMnames] <- train_data_r0$tt * train_data_r0[, Mnames]
    train_data_r0[, ttRMnames] <- train_data_r0$tt * train_data_r0$R * train_data_r0[, Mnames]
    
    train_data_r1 <- train_data; train_data_r1$R <- 1; train_data_r1$tt <- 0
    train_data_r1[, "ttR"] <- with(train_data_r1, {tt*R})
    train_data_r1[, RMnames] <- train_data_r1$R * train_data_r1[, Mnames]
    train_data_r1[, ttMnames] <- train_data_r1$tt * train_data_r1[, Mnames]
    train_data_r1[, ttRMnames] <- train_data_r1$tt * train_data_r1$R * train_data_r1[, Mnames]
    
    valid_data_r0 <- valid_data; valid_data_r0$R <- 0; valid_data_r0$tt <- 0
    valid_data_r0[, "ttR"] <- with(valid_data_r0, {tt*R})
    valid_data_r0[, RMnames] <- valid_data_r0$R * valid_data_r0[, Mnames]
    valid_data_r0[, ttMnames] <- valid_data_r0$tt * valid_data_r0[, Mnames]
    valid_data_r0[, ttRMnames] <- valid_data_r0$tt * valid_data_r0$R * valid_data_r0[, Mnames]
    
    valid_data_r1 <- valid_data; valid_data_r1$R <- 1; valid_data_r1$tt <- 0
    valid_data_r1[, "ttR"] <- with(valid_data_r1, {tt*R})
    valid_data_r1[, RMnames] <- valid_data_r1$R * valid_data_r1[, Mnames]
    valid_data_r1[, ttMnames] <- valid_data_r1$tt * valid_data_r1[, Mnames]
    valid_data_r1[, ttRMnames] <- valid_data_r1$tt * valid_data_r1$R * valid_data_r1[, Mnames]
    
    
    
    # creating outcome and training data 
    
    # for v_pseudo -------
    if (vmodel == "v_pseudo ~ t0 + r0 + C") {
      # for training 
      
      y_pred_train_r1 <- predict(y_out.mtrc$y_fit, train_data_r1[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_train_r1 <- bound_precision(y_pred_train_r1)
      }
      v_pseudo_train <- y_pred_train_r1
      
      #for valid 
      
      y_pred_valid_r1 <- predict(y_out.mtrc$y_fit, valid_data_r1[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_valid_r1 <- bound_precision(y_pred_valid_r1)
      }
      
      v_pseudo_valid <- y_pred_valid_r1
      
      # fit v_pseudo
      v_train_data <- train_data 
      v_train_data$v_pseudo_train <- v_pseudo_train
      # here only use the subsample with R==0, differing from the rvalue for pseudo outcome 
      if (alltrain == FALSE) {
        v_train_data <- v_train_data[v_train_data$R==0 & v_train_data$tt==0, ]
      }
      
    }
    
    
    # for E[Y | tt=0, r=0] -----
    if (vmodel == "y_M ~ t0 + r0 + C") {
      # for training 
      
      y_pred_train_r0 <- predict(y_out.mtrc$y_fit, train_data_r0[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_train_r0 <- bound_precision(y_pred_train_r0)
      }
      v_pseudo_train <- y_pred_train_r0
      
      #for valid 
      
      y_pred_valid_r0 <- predict(y_out.mtrc$y_fit, valid_data_r0[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_valid_r0 <- bound_precision(y_pred_valid_r0)
      }
      
      v_pseudo_valid <- y_pred_valid_r0
      
      # fit v_pseudo
      v_train_data <- train_data 
      v_train_data$v_pseudo_train <- v_pseudo_train
      # here only use the subsample with R==0
      if (alltrain == FALSE) {
        v_train_data <- v_train_data[v_train_data$R==0 & v_train_data$tt==0, ]
      }
      
    }
    
    
    # for E[Y | tt=0, r=1] -------
    if (vmodel == "y_M ~ t0 + r1 + C") {
      # for training 
      
      y_pred_train_r1 <- predict(y_out.mtrc$y_fit, train_data_r1[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_train_r1 <- bound_precision(y_pred_train_r1)
      }
      v_pseudo_train <- y_pred_train_r1
      
      #for valid 
      
      y_pred_valid_r1 <- predict(y_out.mtrc$y_fit, valid_data_r1[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
      if (Yfamily == "binomial") {
        y_pred_valid_r1 <- bound_precision(y_pred_valid_r1)
      }
      
      v_pseudo_valid <- y_pred_valid_r1
      
      # fit v_pseudo
      v_train_data <- train_data 
      v_train_data$v_pseudo_train <- v_pseudo_train
      # here only use the subsample with R==1, the same as the rvalue for pseudo outcome 
      if (alltrain==FALSE) {
        v_train_data <- v_train_data[v_train_data$R==1 & v_train_data$tt==0, ]
      }
      
    }
  }
  
  # fitting ------
  # cov_names <- c(Cnames)
  if (alltrain==TRUE) {
    cov_names <- c("tt", "R", "ttR", Cnames)
  }
  if (alltrain==FALSE) {
    cov_names <- c(Cnames)
  }
  suppressWarnings(
    sl_fit <- SuperLearner(
      Y = v_train_data$v_pseudo_train, 
      X = v_train_data[, cov_names], 
      # family = Yfamily, #
      # family = "binomial",
      # obsWeights = v_train_data[, obs_weights],
      SL.library = SL_library #if Yfamily=="binomial", SL_library=c("SL.glm.scaledY", "SL.gam.scaledY", "SL.glmnet.scaledY", "SL.ranger.scaledY")  
    )
  )
  # only covariates C are used in prediction, which are the same across the datasets, _data_r0, _data_r1, etc.
  sl_pred_train <- predict(sl_fit, train_data[, cov_names], onlySL = TRUE)$pred
  sl_pred_valid <- predict(sl_fit, valid_data[, cov_names], onlySL = TRUE)$pred
  
  
  out <- list(
    v_pseudo_valid = v_pseudo_valid, 
    # y_pred_valid_r1m0 = y_pred_valid_r1m0,
    # y_pred_valid_r1m1 = y_pred_valid_r1m1,
    v_pred_train = sl_pred_train,
    v_pred_valid = sl_pred_valid,
    v_fit = sl_fit
  )
  
  return(out)
}


SL.gam.scaledY <- function(..., family =  binomial() ) { #"quasibinomial"
  SL.gam(..., family =  binomial() )
}

SL.glm.scaledY <- function(..., family =  binomial() ) {
  SL.glm(..., family =  binomial() )
}

SL.earth.scaledY <- function(..., family = binomial() ) { #"quasibinomial"
  SL.earth(..., family = binomial())
}

SL.glmnet.scaledY <- function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, 
                               nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  # .SL.require("glmnet")
  propY <- as.matrix(data.frame(p0 = 1-Y, p1=Y))
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = propY, 
                             family = "binomial" 
                             # weights = obsWeights, 
                             # lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             # alpha = alpha, nlambda = nlambda, ...
  )
  pred1 <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                     "lambda.min", "lambda.1se"))
  
  # fitCVgau <- glmnet::cv.glmnet(x = (X), y = Y, 
  #                            family = "gaussian" 
  #                            # weights = obsWeights, 
  #                            # lambda = NULL, type.measure = loss, nfolds = nfolds, 
  #                            # alpha = alpha, nlambda = nlambda, ...
  # )
  
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- c("SL.glmnet")
  out <- list(pred = pred1, fit = fit)
  return(out)
}


SL.ranger.scaledY <- function (Y, X, newX, family= "gaussian", obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), 
                               write.forest = TRUE, probability = FALSE, 
                               min.node.size = 5, #ifelse(family$family == "gaussian", 5, 1), 
                               replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), 
                               num.threads = 1, verbose = T, ...) 
{
  # .SL.require("ranger")
  # if (family$family == "binomial") {
  #   Y = as.factor(Y)
  # }
  # if ( length(unique(Y))==2 ) {
  #   Y = as.factor(Y)
  # }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X), 
                        splitrule = NULL, # "beta",#
                        num.trees = num.trees, mtry = mtry, 
                        min.node.size = 5, probability = FALSE, 
                        replace = replace, sample.fraction = sample.fraction, 
                        case.weights = obsWeights, write.forest = write.forest, 
                        num.threads = num.threads, 
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  # if (family$family == "binomial") {
  #   pred = pred[, "1"]
  # }
  
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}


SL.nnet.modified <- function(..., family = binomial()) {
  SL.nnet(..., family = binomial(),
          entropy = TRUE   )
}


SL.xgboost.scaledY <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
          max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
          nthread = 1, verbose = 0, save_period = NULL, ...) 
{
  # .SL.require("xgboost")
  # if (packageVersion("xgboost") < "0.6") 
  #   stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  
  #if (family$family == "gaussian") {
    # if (packageVersion("xgboost") >= "1.1.1.1") {
    #   objective <- "reg:squarederror"
    # }
    # else {
    #   objective <- "reg:linear"
    # }
    model = xgboost::xgboost(
      data = xgmat, 
      objective = "reg:logistic",
      nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
      eta = shrinkage, verbose = verbose, nthread = nthread, 
      params = params, save_period = save_period
      )
  #}
  
  
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}
