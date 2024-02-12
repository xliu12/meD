library(SuperLearner)
library(origami)

source("fitmodels.R")

eif.onefold <- function(fold,
                        data_in, 
                        Cnames, 
                        Mnames,
                        # Fit = "linear", 
                        Fit = "mlr",
                        # Znames,
                        # Xnames,
                        Yfamily = "binomial",
                        cv = TRUE
) {
  if(cv==TRUE) {
    # # make training and validation data
    train_data <- origami::training(data_in)
    valid_data <- origami::validation(data_in)
    
    fold <- origami::fold_index()
    
  }
  if(cv==FALSE) {
    train_data <- valid_data <- data_in
    fold <- 0
  }
  
  
  # Fit models for components of the eif -----------------
  set.seed(12222)
  
  if (Fit == "linear") { 
    SL_user.vfit <- SL_user <- c("SL.glm") 
    }
  if (Fit != "linear") { 
    SL_user.vfit <- SL_user <- c("SL.glm", "SL.xgboost", "SL.ranger") 
  }
  
  if (Yfamily == "binomial") {
    if (Fit == "linear") { 
      SL_user.vfit <- c("SL.glm.scaledY") 
      }
    if (Fit != "linear") { 
      SL_user.vfit <-  c("SL.glm.scaledY", "SL.xgboost", "SL.ranger") 
    }
  }
  
  r_out.tc <- fitting.R(
    train_data = train_data, valid_data = valid_data,
    rmodel = "R ~ tt + C",
    Cnames = Cnames,
    Mnames = Mnames,
    # Znames = Znames,
    # Xnames = Xnames,
    SL_library = SL_user # c("SL.glm", "SL.ranger", "SL.glmnet")
  )
  r_out.mtc <- fitting.R(
    train_data = train_data, valid_data = valid_data,
    rmodel = "R ~ M * tt + C",
    Cnames = Cnames,
    Mnames = Mnames,
    # Znames = Znames,
    # Xnames = Xnames, 
    SL_library = SL_user #c("SL.glm", "SL.ranger", "SL.glmnet")
  )
  
  
  # fit tt models
  t_out.c <- fitting.tt(
    train_data = train_data, valid_data = valid_data,
    tmodel = "tt ~ C", 
    Cnames = Cnames,
    Mnames = Mnames,
    # Znames = Znames,
    # Xnames = Xnames, 
    SL_library = SL_user #c("SL.glm", "SL.ranger", "SL.glmnet")
  )
  
  
  # # fit m models 
  # m_out.trc <- fitting.M(
  #   train_data = train_data, valid_data = valid_data,
  #   mmodel = "M ~ tt * R + C", 
  #   Cnames = Cnames,
  #   Mnames = Mnames,
  #   SL_library = c("SL.glm", "SL.ranger", "SL.glmnet")
  # )
  
  
  # fit y models
  y_out.mtrc <- fitting.Y(
    train_data = train_data, valid_data = valid_data,
    # ymodel = "Y ~ M + tt + R + X + C", 
    ymodel = "Y ~ M * tt * R + C",
    Cnames = Cnames, 
    Mnames = Mnames,
    # Znames = Znames,
    # Xnames = Xnames, 
    Yfamily = Yfamily,
    SL_library = SL_user #c("SL.glm", "SL.ranger", "SL.glmnet")
  )
  
  # directly estimating E(Y|t,r,C)
  y_out.trc <- fitting.Y(
    train_data = train_data, valid_data = valid_data,
    # ymodel = "Y ~ tt + R + C",
    ymodel = "Y ~ tt * R + C",
    Cnames = Cnames,
    Mnames = Mnames,
    # Znames = Znames,
    # Xnames = Xnames,
    Yfamily = Yfamily,
    SL_library = SL_user #c("SL.glm", "SL.ranger", "SL.glmnet") #, "SL.ranger", "SL.glmnet"
  )
  
  
  # fit pseudo outcome models
  
  v_out.t1r0c <- fitting.v_pseudo(
    train_data = train_data, valid_data = valid_data,
    vmodel = "v_pseudo ~ t1 + r0 + C", 
    alltrain = FALSE, #TRUE,
    Cnames = Cnames, 
    Mnames = Mnames,
    # Znames = Znames,
    # Xnames = Xnames, 
    y_out.mtrc = y_out.mtrc, 
    t_out.c = t_out.c,
    r_out.tc = r_out.tc, 
    r_out.mtc = r_out.mtc,
    Yfamily = Yfamily,
    SL_library = SL_user.vfit #c("SL.glm", "SL.ranger", "SL.glmnet")
  )
  
  v_out.t0r0c <- fitting.v_pseudo(
    train_data = train_data, valid_data = valid_data,
    vmodel = "v_pseudo ~ t0 + r0 + C", 
    alltrain = FALSE, #TRUE,
    Cnames = Cnames, 
    Mnames = Mnames,
    # Znames = Znames,
    # Xnames = Xnames, 
    y_out.mtrc = y_out.mtrc, 
    t_out.c = t_out.c,
    r_out.tc = r_out.tc, 
    r_out.mtc = r_out.mtc,
    Yfamily = Yfamily,
    SL_library = SL_user.vfit #c("SL.glm", "SL.ranger", "SL.glmnet")
  )
  
  # for validation data
  # !!(updated) valid_data_r1,r0 when tt=1 -------
  
  RMnames <- paste0("RM.", Mnames)
  ttMnames <- paste0("ttM.", Mnames)
  ttRMnames <- paste0("ttRM.", Mnames)
  
  # weights ----------------
  # for r
  r1.tc_pred_valid <- r_out.tc$r1_pred_valid
  r1.mtc_pred_valid <- r_out.mtc$r1_pred_valid
  # for tt
  t1.c_pred_valid <- t_out.c$t1_pred_valid
  
  components.tval <- function(tval = 1) {
    valid_data_r0 <- valid_data; valid_data_r0$tt <- tval; valid_data_r0$R <- 0
    valid_data_r0[, "ttR"] <- with(valid_data_r0, {tt*R})
    valid_data_r0[, RMnames] <- valid_data_r0$R * valid_data_r0[, Mnames]
    valid_data_r0[, ttMnames] <- valid_data_r0$tt * valid_data_r0[, Mnames]
    valid_data_r0[, ttRMnames] <- valid_data_r0$tt * valid_data_r0$R * valid_data_r0[, Mnames]
    
    valid_data_r1 <- valid_data; valid_data_r1$tt <- tval; valid_data_r1$R <- 1
    valid_data_r1[, "ttR"] <- with(valid_data_r1, {tt*R})
    valid_data_r1[, RMnames] <- valid_data_r1$R * valid_data_r1[, Mnames]
    valid_data_r1[, ttMnames] <- valid_data_r1$tt * valid_data_r1[, Mnames]
    valid_data_r1[, ttRMnames] <- valid_data_r1$tt * valid_data_r1$R * valid_data_r1[, Mnames]
    
    # calculate tilting covariate
    tval.c_pred_valid <- 1*(tval==1)*t1.c_pred_valid + 1*(tval==0)*(1 - t1.c_pred_valid)
    H_y.mtrc <- with(valid_data, {
      1*(tt==tval)*(R==1) * (1-r1.mtc_pred_valid) / bound_propensity(tval.c_pred_valid * r1.mtc_pred_valid * (1-r1.tc_pred_valid))
    })
    
    y.mtrc_valid_r1 <- predict(y_out.mtrc$y_fit, valid_data_r1[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
    
    
    # EIF scores calculating ----
    west_y <- with(valid_data, {
      1*(tt==tval)*(R==1) * H_y.mtrc * (Y) / 
        mean( 1*(tt==tval)*(R==1) * H_y.mtrc )
    })
    
    eif_y <- with(valid_data, {
      1*(tt==tval)*(R==1) * H_y.mtrc * (Y - y.mtrc_valid_r1) / 
        mean( 1*(tt==tval)*(R==1) * H_y.mtrc )
    })
    
    if (tval == 1) {
      v_pseudo_valid <- v_out.t1r0c$v_pseudo_valid
      v_pred_valid <- est_reg <- v_out.t1r0c$v_pred_valid
    } 
    if (tval == 0) {
      v_pseudo_valid <- v_out.t0r0c$v_pseudo_valid
      v_pred_valid <- est_reg <- v_out.t0r0c$v_pred_valid
    }
    
    # v_pseudo_valid <- y_pred_valid_r1
    
    H_v.trc <- with(valid_data, {
      1*(tt==tval)*(R==0) / bound_propensity(tval.c_pred_valid * (1 - r1.tc_pred_valid))
    })
    eif_v <- with(valid_data, {
      1*(tt==tval)*(R==0) * H_v.trc * (v_pseudo_valid - v_pred_valid) / 
        mean( 1*(tt==tval)*(R==0)*H_v.trc )
    })
    
    eif <- eif_y + eif_v + v_pred_valid
    
    # Y( tval ) | r=1 -----------
    y.trc_valid_r1 <- predict(y_out.trc$y_fit, valid_data_r1[, y_out.trc$y_fit$varNames], onlySL = TRUE)$pred
    
    # y.trc_valid_r1 <- 1*(tval==1)*(y_out.t1r1c$v_pred_valid) + 1*(tval==0)*(y_out.t0r1c$v_pred_valid)
    
    H_y1.trc <- with(valid_data, {
      1*(tt==tval)*(R==1) / bound_propensity(tval.c_pred_valid * (r1.tc_pred_valid))
    })
    
    # doubly robust estimator
    eif_y1_c <- with(valid_data, {
      1*(tt==tval)*(R==1) * H_y1.trc * (Y - y.trc_valid_r1) / 
        mean( 1*(tt==tval)*(R==1) * H_y1.trc )
    })
    
    eif_y1 <- eif_y1_c + y.trc_valid_r1
    
    # regression estimator
    est_reg1 <- y.trc_valid_r1
    # weighted estimator
    west_y1 <- with(valid_data, {
      1*(tt==tval)*(R==1) * H_y1.trc * (Y) / 
        mean( 1*(tt==tval)*(R==1) * H_y1.trc )
    })
    
    
    
    # Y( tval ) | r=0 -----------
    y.trc_valid_r0 <- predict(y_out.trc$y_fit, valid_data_r0[, y_out.trc$y_fit$varNames], onlySL = TRUE)$pred
    
    # y.trc_valid_r0 <- 1*(tval==1)*(y_out.t1r0c$v_pred_valid) + 1*(tval==0)*(y_out.t0r0c$v_pred_valid)
    
    H_y0.trc <- with(valid_data, {
      1*(tt==tval)*(R==0) / bound_propensity(tval.c_pred_valid * (1-r1.tc_pred_valid))
    })
    
    # doubly robust estimator
    eif_y0_c <- with(valid_data, {
      1*(tt==tval)*(R==0) * H_y0.trc * (Y - y.trc_valid_r0) / 
        mean( 1*(tt==tval)*(R==0) * H_y0.trc )
    })
    eif_y0 <- eif_y0_c + y.trc_valid_r0
    
    # regression estimator 
    est_reg0 <- y.trc_valid_r0
    # weighted estimator
    west_y0 <- with(valid_data, {
      1*(tt==tval)*(R==0) * H_y0.trc * (Y) / 
        mean( 1*(tt==tval)*(R==0) * H_y0.trc )
    })
    
    
    components <- data.frame(
      valid_id = valid_data$id,
      r1.tc_pred_valid = r1.tc_pred_valid,
      r1.mtc_pred_valid = r1.mtc_pred_valid,
      t1.c_pred_valid = t1.c_pred_valid,
      
      y.mtrc_valid_r1 = y.mtrc_valid_r1,
      y.trc_valid_r0 = y.trc_valid_r0,
      y.trc_valid_r1 = y.trc_valid_r1,
      v_pred_valid = v_pred_valid,
      # eif & estimator # grep("eif", ls(), value = TRUE) # grep("est", ls(), value = TRUE)
      # Y(1, Gm(tt=1)|0 ), eif & est
      eif = eif, est_reg = est_reg, west_y = west_y,
      # eif for Y(tt=1) | r1
      eif_y1 = eif_y1, est_reg1 = est_reg1, west_y1 = west_y1,
      # eif for Y(tt=1) | r0
      eif_y0 = eif_y0, est_reg0 = est_reg0, west_y0 = west_y0
    )
    
    return(components)
  }
  
  components_tt1 <- components.tval(tval = 1)
  components_tt0 <- components.tval(tval = 0)
  # which(components_tt0$valid_id != components_tt1$valid_id)
  # which(components_tt0$r1.mtc_pred_valid != components_tt1$r1.mtc_pred_valid)
  
  # for tml ----------------
  
  
  # output ----
  eif_out <- list(
    components_tt1 = components_tt1,
    components_tt0 = components_tt0,
    # fold IDs
    fold = fold# origami::fold_index()
  )
  
  return(eif_out)
}


# Crossfit ------------------------------------
crossfit.onestep <- function(
    cv_folds = 5L,
    data_in, 
    Cnames, 
    Mnames,
    # Xnames,
    Fit = "mlr", #Fit = "linear", 
    Yfamily = "binomial"
) {
  
  set.seed(12345)
  
  if(cv_folds < 1) {
    cv_eif_out <- eif.onefold(fold = 0,
                              data_in = data_in, 
                              Cnames = Cnames, 
                              Mnames = Mnames,
                              # Xnames = Xnames,
                              Fit = Fit, # Fit = "mlr"
                              Yfamily = Yfamily,
                              cv = FALSE
    )
    
    cv_components_tt1 <- cv_eif_out[["components_tt1"]]
    # do.call(rbind, cv_eif_out[["components_tt1"]])
    cv_components_tt1 <- cv_components_tt1[order(cv_components_tt1$valid_id), ]
    
    cv_components_tt0 <- cv_eif_out[["components_tt0"]]
    cv_components_tt0 <- cv_components_tt0[order(cv_components_tt0$valid_id), ]
  }
  
  if (cv_folds > 1) {
    # create cross-validation folds
    folds <- origami::make_folds(data_in,
                                 fold_fun = origami::folds_vfold,
                                 V = cv_folds
    )
    
    cv_eif_out <- origami::cross_validate(
      cv_fun = eif.onefold, 
      folds = folds,
      data_in = data_in, 
      Cnames = Cnames, 
      Mnames = Mnames,
      # Xnames = Xnames,
      Fit = Fit, # Fit = "mlr"
      Yfamily = Yfamily,
      cv = TRUE,
      use_future = FALSE, .combine = FALSE
    )
    
    obs_valid_idx <- do.call(c, lapply(folds, `[[`, "validation_set"))
    # obs_valid_idx[1:4]
    # cv_components$valid_id[1:4]
    # cv_components <- cv_components[order(obs_valid_idx), ]
    
    cv_components_tt1 <- do.call(rbind, cv_eif_out[["components_tt1"]])
    cv_components_tt1 <- cv_components_tt1[order(cv_components_tt1$valid_id), ]
    
    cv_components_tt0 <- do.call(rbind, cv_eif_out[["components_tt0"]])
    cv_components_tt0 <- cv_components_tt0[order(cv_components_tt0$valid_id), ]
    
  }
  
  eifscores <- data.frame(
    Yt1r1 = cv_components_tt1$eif_y1,
    Yt1r0 = cv_components_tt1$eif_y0,
    Yt1r1.Mt1r0 = cv_components_tt1$eif,
    Yt0r1 = cv_components_tt0$eif_y1,
    Yt0r0 = cv_components_tt0$eif_y0,
    Yt0r1.Mt0r0 = cv_components_tt0$eif
  )
  
  # confidence intervals based on asymptotic normality
  effects_scores <- within(eifscores, {
    toD_scores <- (Yt1r1 - Yt0r1) - (Yt1r0 - Yt0r0)
    meD_scores <- (Yt1r1 - Yt1r1.Mt1r0) - (Yt0r1 - Yt0r1.Mt0r0)
    reD_scores <- (Yt1r1.Mt1r0 - Yt1r0) - (Yt0r1.Mt0r0 - Yt0r0)
  })
  
  estimates <- colMeans(effects_scores)
  stderrors <- apply(effects_scores, 2, FUN = function(score) {
    as.numeric(sqrt(var(score) / nrow(data_in)))
  } )
  z_intervals <- apply(effects_scores, 2, FUN = function(score) {
    mean(score) + c(qnorm(.025), qnorm(.975))*as.numeric(sqrt(var(score) / nrow(data_in)))
  } )
  
  
  set.seed(1212)
  odds <- function(x) { mean(x)/(1-mean(x)) } # marginal effect, following Nguyen et al., 2016 [ https://doi.org/10.1080/10705511.2015.1062730]
  
  boot.one <- function(b=0) {
    if (b == 0) { # original sample
      bootid <- 1:nrow(data_in)
    } else {
      bootid <- sample(1:nrow(data_in), nrow(data_in), replace = TRUE)
    }
    
    booteif <- eifscores[bootid, ]
    toD <- with(booteif, {
      mean( (Yt1r1 - Yt0r1) - (Yt1r0 - Yt0r0) )
    })
    
    meD <- with(booteif, {
      mean((Yt1r1 - Yt1r1.Mt1r0) - (Yt0r1 - Yt0r1.Mt0r0))
    })
    
    reD <- with(booteif, {
      mean((Yt1r1.Mt1r0 - Yt1r0) - (Yt0r1.Mt0r0 - Yt0r0))
    })
    
    bootone <- c(toD, meD, reD)
    names(bootone) <- c("toD", "meD", "reD")
    
    if (Yfamily == "binomial") {
      toD_or <- with(booteif, { 
        (odds(Yt1r1) / odds(Yt0r1)) / (odds(Yt1r0) / odds(Yt0r0))
      })
      
      meD_or <- with(booteif, { 
        (odds(Yt1r1) / odds(Yt1r1.Mt1r0)) / (odds(Yt0r1) / odds(Yt0r1.Mt0r0))
      })
      
      reD_or <- with(booteif, { 
        (odds(Yt1r1.Mt1r0) / odds(Yt1r0)) / (odds(Yt0r1.Mt0r0) / odds(Yt0r0))
      })
      
      bootone <- c(toD, meD, reD, toD_or, meD_or, reD_or)
      names(bootone) <- c("toD", "meD", "reD", "toD_or", "meD_or", "reD_or")
      
    }
    
    return(bootone)
  }
  # nboot_out <- lapply(0:1000, boot.one)
  # nboot_out <- do.call(rbind, nboot_out)
  # original_est <- nboot_out[1, ]
  # boot_intervals <- apply(nboot_out, 2, quantile, prob = c(0.05, 0.95))
  
  
  
  # weighting-based estimators ----------------
  westscores <- data.frame(
    Yt1r1 = cv_components_tt1$west_y1,
    Yt1r0 = cv_components_tt1$west_y0,
    Yt1r1.Mt1r0 = cv_components_tt1$west_y,
    Yt0r1 = cv_components_tt0$west_y1,
    Yt0r0 = cv_components_tt0$west_y0,
    Yt0r1.Mt0r0 = cv_components_tt0$west_y
  )
  
  effects_west <- within(westscores, {
    toD_scores <- (Yt1r1 - Yt0r1) - (Yt1r0 - Yt0r0)
    meD_scores <- (Yt1r1 - Yt1r1.Mt1r0) - (Yt0r1 - Yt0r1.Mt0r0)
    reD_scores <- (Yt1r1.Mt1r0 - Yt1r0) - (Yt0r1.Mt0r0 - Yt0r0)
  })
  
  weighting_estimates <- colMeans(effects_west)
  
  # output
  out <- list(
    cv_components_tt1 = cv_components_tt1,
    cv_components_tt0 = cv_components_tt0,
    effects_scores = effects_scores,
    estimates = estimates,
    stderrors = stderrors,
    z_intervals = z_intervals
    , weighting_estimates = weighting_estimates
    #, boot_intervals = boot_intervals
  )
  
  return(out)
}














# TMLE ----------------------------------------
tmle.medMO <- function(
    crossfit_onestep = crossfit_out,
    data_in, 
    Cnames, Mnames #, Xnames
){
  
  tml_data <- data_in[, c("id", "R", "tt", Mnames, "Y", "obs_weights_Y")]
  
  # list2env(as.list(cv_components), envir = environment())
  
  tilt.tval <- function(tval = 1) {
    if (tval == 1) {
      cv_components <- crossfit_onestep$cv_components_tt1
    }
    if (tval == 0) {
      cv_components <- crossfit_onestep$cv_components_tt0
    }
    
    # make sure rows of the components and data for the same individuals
    cv_components <- cv_components[order(cv_components$valid_id), ]
    # sum(cv_components$valid_id != data_in$id)
    
    tml_data <- merge(tml_data, cv_components, by.x = "id", by.y = "valid_id")
    colnames(tml_data)
    range_y <- with(tml_data, {
      range(Y, y.mtrc_valid_r1, y.trc_valid_r0, y.trc_valid_r1, v_pred_valid)
    })
    min_y <- range_y[1] 
    max_y <- range_y[2] 
    
    
    # scaled outcome, pseudo outcome, predicted outcome
    tml_data$Yscaled <- bound_precision((tml_data$Y - min_y) / (max_y - min_y))
    tml_data$yscaled.mtrc_valid_r1 <- bound_precision((tml_data$y.mtrc_valid_r1 - min_y) / (max_y - min_y)) 
    
    tml_data$yscaled.trc_valid_r0 <- bound_precision((tml_data$y.trc_valid_r0 - min_y) / (max_y - min_y))  
    tml_data$yscaled.trc_valid_r1 <- bound_precision((tml_data$y.trc_valid_r1 - min_y) / (max_y - min_y))  
    
    tml_data$vscaled_pred_valid <- bound_precision((tml_data$v_pred_valid - min_y) / (max_y - min_y)) 
    
    # calculate tilting covariate
    tml_data$tval.c_pred_valid <- with(tml_data, {
      tval.c_pred_valid <- 1*(tval==1)*t1.c_pred_valid + 1*(tval==0)*(1 - t1.c_pred_valid)
    })
    
    tml_data$H_y.mtrc <- with(tml_data, {
      1*(tt==tval)*(R==1) * (1-r1.mtc_pred_valid) / bound_propensity(tval.c_pred_valid * r1.mtc_pred_valid * (1-r1.tc_pred_valid))
    })
    
    
    max_iter <- 3
    tiltmod_tol <- 10
    
    # targeting y.mtrc----
    se_eif <- as.numeric(sqrt(var(tml_data$eif) / length(tml_data$eif)))
    tilt_stop_crit <- se_eif / log( length(tml_data$eif) )
    
    # Iterative targeting 
    yfit_score <- Inf
    eif_stop_crit <- FALSE
    n_iter <- 0
    while ((!eif_stop_crit) && n_iter <= max_iter) {
      if (abs(mean(yfit_score)) > tilt_stop_crit) {
        
        
        y_up_fit <- glm(
          "Yscaled ~ 1 + offset(qlogis(yscaled.mtrc_valid_r1))", 
          weights = H_y.mtrc,  
          family="quasibinomial", start = 0,
          data = tml_data
        ) 
        y_up_fit$df.null
        if (is.na(coef(y_up_fit))) {
          y_up_fit$coefficients <- 0
        } else if (!y_up_fit$converged || abs(max(stats::coef(y_up_fit))) >
                   tiltmod_tol) {
          y_up_fit$coefficients <- 0
        }
        y_up_coef <- unname(stats::coef(y_up_fit))
        
        tml_data$yscaled.mtrc_valid_r1 <- with(tml_data, {
          plogis(qlogis(yscaled.mtrc_valid_r1) + y_up_coef)
        })
        
        # same as unstab_eif_y
        yfit_score <- with(tml_data, {
          H_y.mtrc * (Yscaled - yscaled.mtrc_valid_r1)
        })
      } else { # (abs(mean(yfit_score)) <= tilt_stop_crit)
        yfit_score <- 0
      }
      
      eif_stop_crit <- (abs(mean(yfit_score)) < tilt_stop_crit)
      
      n_iter <- n_iter + 1
    }
    
    # targeting v----
    tml_data$vscaled_pseudo <- bound_precision(tml_data$yscaled.mtrc_valid_r1)
    
    tml_data$H_v.trc <- with(tml_data, {
      1*(tt==tval)*(R==0) / bound_propensity(tval.c_pred_valid * (1 - r1.tc_pred_valid))
    })
    
    # Iterative targeting 
    vfit_score <- Inf
    eif_stop_crit <- FALSE
    n_iter <- 0
    while ((!eif_stop_crit) && n_iter <= max_iter) {
      if (abs(mean(vfit_score)) > tilt_stop_crit) {
        v_up_fit <- glm(
          "vscaled_pseudo ~ 1 + offset(qlogis(vscaled_pred_valid))",
          weights = H_v.trc,
          family="quasibinomial", start = 0, data = tml_data
        )
        
        if (is.na(coef(v_up_fit))) {
          v_up_fit$coefficients <- 0
        } else if (!v_up_fit$converged || abs(max(stats::coef(v_up_fit))) >
                   tiltmod_tol) {
          v_up_fit$coefficients <- 0
        }
        v_up_coef <- unname(stats::coef(v_up_fit))
        
        tml_data$vscaled_pred_valid <- with(tml_data, {
          plogis(qlogis(vscaled_pred_valid) + v_up_coef)
        })
        
        # check fit score
        vfit_score <- with(tml_data, {
          H_v.trc * (vscaled_pseudo - vscaled_pred_valid)
        })
      } else {
        vfit_score <- 0
      }
      
      eif_stop_crit <- (abs(mean(vfit_score)) < tilt_stop_crit)
      
      n_iter <- n_iter + 1
    }
    tml_data$vscaled_tmle <- unname(predict(v_up_fit, type = "response"))
    # tml_data$vscaled_pred_valid 
    # v_tmle <- plogis(qlogis(v_pred_valid) + coef(v_trc_up_fit))
    
    
    # ~~~~~~~~~~~~~~~~~~~~-------------------
    
    
    # tmle y1_c -----------
    se_eif <- as.numeric(sqrt(var(tml_data$eif_y1) / length(tml_data$eif_y1) ))
    tilt_stop_crit <- se_eif / log(length(tml_data$eif_y1))
    
    # tilting covariate
    tml_data$H_y1.trc <- with(tml_data, {
      1*(tt==tval)*(R==1) / bound_propensity(tval.c_pred_valid * (r1.tc_pred_valid))
    })
    
    # iterative targeting 
    eif_stop_crit <- FALSE
    n_iter <- 0
    y_c_fit_score <- Inf
    
    while ( (!eif_stop_crit) && n_iter <= max_iter ) {
      if (abs(mean(y_c_fit_score)) > tilt_stop_crit) {
        
        y_c_up_fit <- glm(
          "Yscaled ~ 1 + offset(qlogis(yscaled.trc_valid_r1))", 
          weights = H_y1.trc, # only the overlap part
          family="quasibinomial", start = 0, data=tml_data
        ) 
        y_c_up_fit$df.null
        
        if (is.na(stats::coef(y_c_up_fit))) {
          y_c_up_fit$coefficients <- 0
        } else if (!y_c_up_fit$converged || abs(max(stats::coef(y_c_up_fit))) > tiltmod_tol) {
          y_c_up_fit$coefficients <- 0
        }
        y_c_up_coef <- unname(stats::coef(y_c_up_fit))
        
        tml_data$yscaled.trc_valid_r1 <- plogis(qlogis(tml_data$yscaled.trc_valid_r1) + y_c_up_coef)
        
        # same as unstab_eif_y
        y_c_fit_score <- with(tml_data, {
          H_y1.trc * (Yscaled - yscaled.trc_valid_r1)
        })
        
      } else {
        y_c_fit_score <- 0
      }
      
      eif_stop_crit <- all( abs(c(mean(y_c_fit_score) )) < tilt_stop_crit)
      
      n_iter <- n_iter + 1
    }
    
    # y1_tmle <- predict(y_zc_up_fit, type = "response") # alternatively, 
    tml_data$y1scaled_tmle <- tml_data$yscaled.trc_valid_r1 # y1mean_c_up_r1
    
    
    # tmle y0_c -----------
    
    se_eif <- as.numeric(sqrt(var(tml_data$eif_y0) / length(tml_data$eif_y0) ))
    tilt_stop_crit <- se_eif / log(length(tml_data$eif_y0))
    
    # tilting covariate
    tml_data$H_y0.trc <- with(tml_data, {
      1*(tt==tval)*(R==0) / bound_propensity(tval.c_pred_valid * (1-r1.tc_pred_valid))
    })
    
    eif_stop_crit <- FALSE
    n_iter <- 0
    y_c_fit_score <- Inf
    
    # iterative targeting
    while ( (!eif_stop_crit) && n_iter <= max_iter ) {
      if (abs(mean(y_c_fit_score)) > tilt_stop_crit) {
        
        y_c_up_fit <- glm(
          "Yscaled ~ 1 + offset(qlogis(yscaled.trc_valid_r0))", 
          weights = H_y0.trc,  
          family="quasibinomial", start = 0, data=tml_data
        ) 
        y_c_up_fit$df.null
        
        if (is.na(stats::coef(y_c_up_fit))) {
          y_c_up_fit$coefficients <- 0
        } else if (!y_c_up_fit$converged || abs(max(stats::coef(y_c_up_fit))) > tiltmod_tol) {
          y_c_up_fit$coefficients <- 0
        }
        y_c_up_coef <- unname(stats::coef(y_c_up_fit))
        
        tml_data$yscaled.trc_valid_r0 <- plogis(qlogis(tml_data$yscaled.trc_valid_r0) + y_c_up_coef)
        
        # same as unstab_eif_y
        y_c_fit_score <- with(tml_data, {
          H_y0.trc * (Yscaled - yscaled.trc_valid_r0)
        })
        
      } else {
        y_c_fit_score <- 0
      }
      
      eif_stop_crit <- all( abs(c(mean(y_c_fit_score) )) < tilt_stop_crit)
      
      n_iter <- n_iter + 1
    }
    
    # y0_tmle <- predict(y_c_up_fit, type = "response") # alternatively, 
    tml_data$y0scaled_tmle <- tml_data$yscaled.trc_valid_r0 
    
    # scale back the outcome parameters -----
    
    tml_data$v_tmle <- scale_from_unit(tml_data$vscaled_tmle, max_orig = max_y, min_orig = min_y)
    tml_data$y1_tmle <- scale_from_unit(tml_data$y1scaled_tmle, max_orig = max_y, min_orig = min_y)
    tml_data$y0_tmle <- scale_from_unit(tml_data$y0scaled_tmle, max_orig = max_y, min_orig = min_y)
    
    return(tml_data)
  }
  
  tml_data_tt1 <- tilt.tval(tval = 1)
  tml_data_tt0 <- tilt.tval(tval = 0)
  
  
  # tml effect estimates -------
  tml_scores <- data.frame(
    Yt1r1 = tml_data_tt1$y1_tmle,
    Yt1r0 = tml_data_tt1$y0_tmle,
    Yt1r1.Mt1r0 = tml_data_tt1$v_tmle,
    Yt0r1 = tml_data_tt0$y1_tmle,
    Yt0r0 = tml_data_tt0$y0_tmle,
    Yt0r1.Mt0r0 = tml_data_tt0$v_tmle
  )
  tml_effects <- within(tml_scores, {
    toD_scores <- (Yt1r1 - Yt0r1) - (Yt1r0 - Yt0r0)
    meD_scores <- (Yt1r1 - Yt1r1.Mt1r0) - (Yt0r1 - Yt0r1.Mt0r0)
    reD_scores <- (Yt1r1.Mt1r0 - Yt1r0) - (Yt0r1.Mt0r0 - Yt0r0)
  })
  tml_estimates <- colMeans(tml_effects)
  
  
  # same asymptotic behavior as one-step estimator
  cv_components_tt1 <- crossfit_onestep$cv_components_tt1
  cv_components_tt0 <- crossfit_onestep$cv_components_tt0
  eifscores <- data.frame(
    Yt1r1 = cv_components_tt1$eif_y1,
    Yt1r0 = cv_components_tt1$eif_y0,
    Yt1r1.Mt1r0 = cv_components_tt1$eif,
    Yt0r1 = cv_components_tt0$eif_y1,
    Yt0r0 = cv_components_tt0$eif_y0,
    Yt0r1.Mt0r0 = cv_components_tt0$eif
  )
  effects_scores <- within(eifscores, {
    toD_scores <- (Yt1r1 - Yt0r1) - (Yt1r0 - Yt0r0)
    meD_scores <- (Yt1r1 - Yt1r1.Mt1r0) - (Yt0r1 - Yt0r1.Mt0r0)
    reD_scores <- (Yt1r1.Mt1r0 - Yt1r0) - (Yt0r1.Mt0r0 - Yt0r0)
  })
  
  stderrors <- apply(effects_scores, 2, FUN = function(score) {
    as.numeric(sqrt(var(score) / nrow(data_in)))
  } )
  
  # tml estimates with standard errors of the eif scores
  tml_intervals <- t(sapply(c(qnorm(.025), qnorm(.975)), FUN = function(qz) {
    tml_estimates + qz*stderrors
  }))
  
  
  # output
  tmle_out <- list(
    tml_effects = tml_effects,
    tml_estimates = tml_estimates,
    tml_intervals = tml_intervals,
    stderrors = stderrors
  )
  return(tmle_out)
}

