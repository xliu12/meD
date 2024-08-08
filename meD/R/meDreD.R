


# Estimating the effects --------------------------------------
meDreD <- function(
    data_in, 
    Yname, ttname, Rname,
    Mnames,
    Cnames,
    Yfamily = "binomial",
    estimator,
    Yt1only=FALSE,
    nuisance_estimation, #"linear"
    num_folds = 5
) {
  
  if ("reg" %in% estimator) {
    reg_estimate <- reg.medmod(data_in = data_in, 
                                Yname = Yname, ttname = ttname, Rname = Rname,
                                Mnames = Mnames,
                                Cnames = Cnames, Yt1only=Yt1only)
    
    results <- reg_estimate
  }
  
  
  
  data_in$id <- 1:nrow(data_in)
  data_in$obs_weights_Y <- rep(1, nrow(data_in))

  data_in <- data_in %>% 
    rename(Y=Yname, R=Rname, tt=ttname)
  # add interactions among the key variables
  data_in <- data_in %>% 
  mutate(
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

  # estimators

  # if ((!"tml"%in%estimator) & (!"onestep"%in%estimator)) {
  
  if (("tml"%in%estimator) | ("onestep"%in%estimator)) {
    crossfit_out <- crossfit.onestep(
      cv_folds = num_folds,
      data_in = data_in,
      Cnames = Cnames,
      Mnames = Mnames,
      Fit = nuisance_estimation,
      Yfamily =  Yfamily,
      Yt1only = Yt1only
    )
    
    tmle_out <- tmle.medMO(
      crossfit_onestep = crossfit_out,
      data_in = data_in,
      Cnames = Cnames,
      Mnames = Mnames,
      Yt1only = Yt1only
    )
  }
  

  if ((!"tml"%in%estimator) & ("onestep"%in%estimator)) {
    results = data.frame(
      effect = c("Yt1r1", "Yt1r0", "Yt1r1.Mt1r0",
                 "Yt0r1", "Yt0r0", "Yt0r1.Mt0r0",
                 "reD",
                 "meD",
                 "toD"),

      # onestep
      onestep_estimate = crossfit_out$estimates,
      onestep_interval = t(crossfit_out$z_intervals)
    )
  }
  if (("tml"%in%estimator) & (!"onestep"%in%estimator)) {
    results = data.frame(
      effect = c("Yt1r1", "Yt1r0", "Yt1r1.Mt1r0",
                 "Yt0r1", "Yt0r0", "Yt0r1.Mt0r0",
                 "reD",
                 "meD",
                 "toD"),
      # tml
      tml_estimate = tmle_out$tml_estimates,
      tml_interval = t(tmle_out$tml_intervals)
    )
  }

  if ( ("tml"%in%estimator) & ("onestep"%in%estimator) ) {
    results = data.frame(
      effect = c("Yt1r1", "Yt1r0", "Yt1r1.Mt1r0",
                 "Yt0r1", "Yt0r0", "Yt0r1.Mt0r0",
                 "reD",
                 "meD",
                 "toD"),

      # onestep
      onestep_estimate = crossfit_out$estimates,
      onestep_interval = t(crossfit_out$z_intervals),
      # tml
      tml_estimate = tmle_out$tml_estimates,
      tml_interval = t(tmle_out$tml_intervals)
    )
  }
  
  
  row.names(results) <- NULL


  return(results)
}









crossfit.onestep <- function(
    cv_folds = 5L,
    data_in,
    Cnames,
    Mnames,
    Fit, # e.g., Fit = c("SL.glm")
    Yfamily = "binomial", Yt1only=FALSE
) {

  set.seed(12345)

  if(cv_folds <= 1) {
    cv_eif_out <- eif.onefold(fold = 0,
                              data_in = data_in,
                              Cnames = Cnames,
                              Mnames = Mnames,
                              Fit = Fit,
                              Yfamily = Yfamily,
                              cv = FALSE
    )

    cv_components_tt1 <- cv_eif_out[["components_tt1"]]
    # do.call(rbind, cv_eif_out[["components_tt1"]])
    cv_components_tt1 <- cv_components_tt1[order(cv_components_tt1$valid_set), ]

    cv_components_tt0 <- cv_eif_out[["components_tt0"]]
    cv_components_tt0 <- cv_components_tt0[order(cv_components_tt0$valid_set), ]
  }

  # if (cv_folds > 1) {
  #   # create cross-validation folds
  #   folds <- origami::make_folds(data_in,
  #                                fold_fun = origami::folds_vfold,
  #                                V = cv_folds
  #   )
  #   
  #   
  # 
  #   cv_eif_out <- origami::cross_validate(
  #     cv_fun = eif.onefold,
  #     folds = folds,
  #     data_in = data_in,
  #     Cnames = Cnames,
  #     Mnames = Mnames,
  #     Fit = Fit,
  #     Yfamily = Yfamily,
  #     cv = TRUE,
  #     use_future = FALSE, .combine = FALSE
  #   )
  # 
  #   obs_valid_idx <- do.call(c, lapply(folds, `[[`, "validation_set"))
  # 
  #   cv_components_tt1 <- do.call(rbind, cv_eif_out[["components_tt1"]])
  #   cv_components_tt1 <- cv_components_tt1[order(cv_components_tt1$valid_id), ]
  # 
  #   cv_components_tt0 <- do.call(rbind, cv_eif_out[["components_tt0"]])
  #   cv_components_tt0 <- cv_components_tt0[order(cv_components_tt0$valid_id), ]
  # 
  # }
  
  if(cv_folds > 1) {
    data_in <- data_in %>%
      group_by(R, tt) %>% mutate(K = cur_group_id())
    fold_K <- lapply(unique(data_in$K), FUN = function(k=1) {

      if (nrow(data_in[data_in$K==k, ]) >= 1) {
        fk <- origami::make_folds(data_in[data_in$K==k, ],
                                  fold_fun = origami::folds_vfold,
                                  V = cv_folds)
        fold_k <- fk
        for(v in 1:cv_folds) {
          fold_k[[v]]$validation_set <- data_in$id[data_in$K==k][fk[[v]]$validation_set]
          fold_k[[v]]$training_set <- data_in$id[data_in$K==k][fk[[v]]$training_set]
        }
      }

     return(fold_k)
    } )

    folds <- origami::make_folds(data_in,
                                 fold_fun = origami::folds_vfold,
                                 V = cv_folds)
    v <- 1
    for(v in 1:cv_folds) {
      folds[[v]]$validation_set <- unlist(lapply(1:length(fold_K), FUN = function(k=1) {
        fold_K[[k]][[v]]$validation_set
      }))
      folds[[v]]$training_set <- unlist(lapply(1:length(fold_K), FUN = function(k=1) {
        fold_K[[k]][[v]]$training_set
      }))
    }
    
    eif_tt1 <- eif_tt0 <- NULL
    for(v in 1:cv_folds) {
      eif_out <- eif.onefold(
        v = v,
        folds = folds,
        data_in = data_in,
        Cnames = Cnames,
        Mnames = Mnames,
        Fit = Fit,
        Yfamily = Yfamily,
        cv = TRUE
      )
      
      eif_tt1 <- rbind(eif_tt1, eif_out$components_tt1)
      eif_tt0 <- rbind(eif_tt0, eif_out$components_tt0)
    }
    cv_components_tt1 <- eif_tt1[order(eif_tt1$valid_set), ] # make sure it is in the same order of data_in
    cv_components_tt0 <- eif_tt0[order(eif_tt0$valid_set), ]
    # sum(cv_components$valid_set - data_in$id)
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
  # effects_scores <- within(eifscores, {
  #   toD_scores <- (Yt1r1 - Yt0r1) - (Yt1r0 - Yt0r0)
  #   meD_scores <- (Yt1r1 - Yt1r1.Mt1r0) - (Yt0r1 - Yt0r1.Mt0r0)
  #   reD_scores <- (Yt1r1.Mt1r0 - Yt1r0) - (Yt0r1.Mt0r0 - Yt0r0)
  # })
  effects_scores <- effect(eifscores, Yt1only = Yt1only)

  estimates <- colMeans(effects_scores)
  stderrors <- apply(effects_scores, 2, FUN = function(score) {
    as.numeric(sqrt(var(score) / nrow(data_in)))
  } )
  z_intervals <- apply(effects_scores, 2, FUN = function(score) {
    mean(score) + c(qnorm(.025), qnorm(.975))*as.numeric(sqrt(var(score) / nrow(data_in)))
  } )


  # output
  out <- list(
    cv_components_tt1 = cv_components_tt1,
    cv_components_tt0 = cv_components_tt0,
    effects_scores = effects_scores,
    estimates = estimates,
    stderrors = stderrors,
    z_intervals = z_intervals
  )

  return(out)
}


tmle.medMO <- function(
    crossfit_onestep = crossfit_out,
    data_in,
    Cnames, Mnames, Yt1only=FALSE
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
    # cv_components <- cv_components[order(cv_components$valid_id), ]

    # tml_data <- merge(tml_data, cv_components, by.x = "id", by.y = "valid_id")
    tml_data <- merge(tml_data, cv_components, by.x = "id", by.y = "valid_set")
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
    }) %>% bound_ipw(.)

    max_iter <- 3
    tiltmod_tol <- 10

    # targeting y.mtrc
    se_eif <- as.numeric(sqrt(var(tml_data$eif) / length(tml_data$eif)))
    tilt_stop_crit <- se_eif / log( length(tml_data$eif) )

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

    # targeting v
    tml_data$vscaled_pseudo <- bound_precision(tml_data$yscaled.mtrc_valid_r1)

    tml_data$H_v.trc <- with(tml_data, {
      1*(tt==tval)*(R==0) / bound_propensity(tval.c_pred_valid * (1 - r1.tc_pred_valid))
    }) %>% bound_ipw(.)

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

    # tmle y1_c
    se_eif <- as.numeric(sqrt(var(tml_data$eif_y1) / length(tml_data$eif_y1) ))
    tilt_stop_crit <- se_eif / log(length(tml_data$eif_y1))

    # tilting covariate
    tml_data$H_y1.trc <- with(tml_data, {
      1*(tt==tval)*(R==1) / bound_propensity(tval.c_pred_valid * (r1.tc_pred_valid))
    }) %>% bound_ipw(.)

    eif_stop_crit <- FALSE
    n_iter <- 0
    y_c_fit_score <- Inf

    while ( (!eif_stop_crit) && n_iter <= max_iter ) {
      if (abs(mean(y_c_fit_score)) > tilt_stop_crit) {

        y_c_up_fit <- glm(
          "Yscaled ~ 1 + offset(qlogis(yscaled.trc_valid_r1))",
          weights = H_y1.trc,
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
    tml_data$y1scaled_tmle <- tml_data$yscaled.trc_valid_r1


    # tmle y0_c

    se_eif <- as.numeric(sqrt(var(tml_data$eif_y0) / length(tml_data$eif_y0) ))
    tilt_stop_crit <- se_eif / log(length(tml_data$eif_y0))

    # tilting covariate
    tml_data$H_y0.trc <- with(tml_data, {
      1*(tt==tval)*(R==0) / bound_propensity(tval.c_pred_valid * (1-r1.tc_pred_valid))
    }) %>% bound_ipw(.)

    eif_stop_crit <- FALSE
    n_iter <- 0
    y_c_fit_score <- Inf

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

    # scale back the outcome parameters

    tml_data$v_tmle <- scale_from_unit(tml_data$vscaled_tmle, max_orig = max_y, min_orig = min_y)
    tml_data$y1_tmle <- scale_from_unit(tml_data$y1scaled_tmle, max_orig = max_y, min_orig = min_y)
    tml_data$y0_tmle <- scale_from_unit(tml_data$y0scaled_tmle, max_orig = max_y, min_orig = min_y)

    return(tml_data)
  }

  tml_data_tt1 <- tilt.tval(tval = 1)
  tml_data_tt0 <- tilt.tval(tval = 0)

  # tml effect estimates
  tml_scores <- data.frame(
    Yt1r1 = tml_data_tt1$y1_tmle,
    Yt1r0 = tml_data_tt1$y0_tmle,
    Yt1r1.Mt1r0 = tml_data_tt1$v_tmle,
    Yt0r1 = tml_data_tt0$y1_tmle,
    Yt0r0 = tml_data_tt0$y0_tmle,
    Yt0r1.Mt0r0 = tml_data_tt0$v_tmle
  )
  # tml_effects <- within(tml_scores, {
  #   toD_scores <- (Yt1r1 - Yt0r1) - (Yt1r0 - Yt0r0)
  #   meD_scores <- (Yt1r1 - Yt1r1.Mt1r0) - (Yt0r1 - Yt0r1.Mt0r0)
  #   reD_scores <- (Yt1r1.Mt1r0 - Yt1r0) - (Yt0r1.Mt0r0 - Yt0r0)
  # })
  tml_effects <- effect(tml_scores, Yt1only = Yt1only)
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








# Fitting nuisance functions --------


fitting.R <- function(train_data, valid_data,
                      rmodel = "R ~ tt + C",
                      Cnames, Mnames,
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
                      SL_library = c("SL.glm"),
                      Yfamily = "binomial"
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

  out <- list(
    y_pred_train = sl_pred_train,
    y_pred_valid = sl_pred_valid,
    y_fit = sl_fit
  )

  return(out)
}



fitting.v_pseudo <- function(train_data, valid_data,
                             vmodel = "v_pseudo ~ t1 + r0 + C",
                             alltrain = FALSE,
                             Cnames, Mnames,
                             SL_library = c("SL.glm"),
                             y_out.mtrc,
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
        v_train_data <- v_train_data[v_train_data$R==0 & v_train_data$tt==1, ]
      }
    }


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
      if (alltrain == FALSE) {
        v_train_data <- v_train_data[v_train_data$R==0 & v_train_data$tt==1, ]
      }

    }


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

      if (alltrain == FALSE) {
        v_train_data <- v_train_data[v_train_data$R==1 & v_train_data$tt==1, ]
      }
    }

  }

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

      if (alltrain == FALSE) {
        v_train_data <- v_train_data[v_train_data$R==0 & v_train_data$tt==0, ]
      }

    }


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
      if (alltrain == FALSE) {
        v_train_data <- v_train_data[v_train_data$R==0 & v_train_data$tt==0, ]
      }

    }


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
      if (alltrain==FALSE) {
        v_train_data <- v_train_data[v_train_data$R==1 & v_train_data$tt==0, ]
      }

    }
  }

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
      SL.library = SL_library
    )
  )

  sl_pred_train <- predict(sl_fit, train_data[, cov_names], onlySL = TRUE)$pred
  sl_pred_valid <- predict(sl_fit, valid_data[, cov_names], onlySL = TRUE)$pred


  out <- list(
    v_pseudo_valid = v_pseudo_valid,
    v_pred_train = sl_pred_train,
    v_pred_valid = sl_pred_valid,
    v_fit = sl_fit
  )

  return(out)
}


