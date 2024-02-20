
#' Estimation of the mediated moderation (meD) and remaining moderation (reD)
#' @param data_in a \code{data.frame} containing the observed data.
#'    In "data_in", column "tt" is the treatment assignment ("tt" is coded as 0 for individuals in the control condition and as 1 for individuals in the treatment condition);
#'    column "R" is a dummy indicator of the subgroup status.
#'    column "Y" is the outcome.
#' @param Mnames a character vector of the names of the columns in "data_in" that correspond to mediators (M).
#' @param Cnames a character vector of the names of the columns in "data_in" that correspond to baseline covariates (C).
#' @param Yfamily a character specifying the link function for a generalized linear model for the outcome. Currently supported options include "gaussian" and "binomial".
#' @param estimator a character specifying the estimator for the mediated moderation (meD) and remaining moderation (reD). Currently supported options include "onestep" and "tml", denoting the one-step estimator and targeted minimum loss estimator, respectively.
#' @param nuisance_estimation a character vector specifying the method for estimating the nuisance models. Currently supported options include the list of functions included in the "SuperLearner" package (<https://cran.r-project.org/package=SuperLearner>; can be found with the function listWrappers() of the "SuperLearner" package).
#' @param num_folds the number of folds used for the cross-fitting procedure.
#'
#' @return A data frame containing the estimates of the mediated moderation (meD) and remaining moderation (reD), as well as the total moderation (toD) and the average potential outcomes in the effect definitions. These include the average potential outcomes of each subgroup under each treatment assignment (Yt1r1, Yt1r0, Yt0r1, Yt0r0), and the average potential outcomes of the comparison subgroup if the potential mediator distribution were shifted to be the same as that of the reference subgroup (Yt1r1.Mt1r0, Yt0r1.Mt0r0). The output also includes the 0.95 confidence intervals constructed based on the asymptotic variance estimates (column names ending with "_interval.1" and "_interval.2").

#'
#'
#' @export
#'
#' @examples
#'  # data(data_in)
#'  # data_in <- read.csv("data/data_in.csv", header = TRUE)
#'  # Mnames <- grep("Mdat", colnames(data_in), value = T)
#'  # Cnames <- grep("Cdat", colnames(data_in), value = T)
#'  # out <- meDreD(
#'  # data_in = data_in,
#'  # Mnames = Mnames,
#'  # Cnames = Cnames,
#'  # nuisance_estimation = "SL.glm",
#'  # Yfamily = "binomial",
#'  # estimator = "onestep",
#'  # num_folds = 5)

#'
#'



# Estimating the effects --------------------------------------
meDreD <- function(
    data_in,
    Mnames,
    Cnames,
    nuisance_estimation = c("SL.glm"), #"linear"
    Yfamily = "binomial",
    estimator = "onestep",
    num_folds = 5
) {

  data_in$id <- 1:nrow(data_in)
  data_in$obs_weights_Y <- rep(1, nrow(data_in))
  Cnames <- grep("Cdat", colnames(data_in), value = T)
  Mnames <- grep("Mdat", colnames(data_in), value = T)

  # add interactions among the key variables
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

  # estimators

  crossfit_out <- crossfit.onestep(
    cv_folds = num_folds,
    data_in = data_in,
    Cnames = Cnames,
    Mnames = Mnames,
    Fit = nuisance_estimation,
    Yfamily =  Yfamily
  )

  tmle_out <- tmle.medMO(
    crossfit_onestep = crossfit_out,
    data_in = data_in,
    Cnames = Cnames,
    Mnames = Mnames
  )

  if (!"tml"%in%estimator & "onestep"%in%estimator) {
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
  if ("tml"%in%estimator & !"onestep"%in%estimator) {
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

  if ( "tml"%in%estimator & "onestep"%in%estimator ) {
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


  return(results)
}









crossfit.onestep <- function(
    cv_folds = 5L,
    data_in,
    Cnames,
    Mnames,
    Fit, # e.g., Fit = c("SL.glm")
    Yfamily = "binomial"
) {

  set.seed(12345)

  if(cv_folds < 1) {
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
      Fit = Fit,
      Yfamily = Yfamily,
      cv = TRUE,
      use_future = FALSE, .combine = FALSE
    )

    obs_valid_idx <- do.call(c, lapply(folds, `[[`, "validation_set"))

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
    Cnames, Mnames
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
    })

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
    })

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
    })

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




# Calculating the scores for each individual observation ------


eif.onefold <- function(fold,
                        data_in,
                        Cnames,
                        Mnames,
                        Fit,
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

  # Fit models for components of the eif
  set.seed(12222)

  SL_user <- Fit

  r_out.tc <- fitting.R(
    train_data = train_data, valid_data = valid_data,
    rmodel = "R ~ tt + C",
    Cnames = Cnames,
    Mnames = Mnames,
    SL_library = SL_user
  )
  r_out.mtc <- fitting.R(
    train_data = train_data, valid_data = valid_data,
    rmodel = "R ~ M * tt + C",
    Cnames = Cnames,
    Mnames = Mnames,
    SL_library = SL_user
  )

  # fit tt models
  t_out.c <- fitting.tt(
    train_data = train_data, valid_data = valid_data,
    tmodel = "tt ~ C",
    Cnames = Cnames,
    Mnames = Mnames,
    SL_library = SL_user
  )

  # fit y models
  y_out.mtrc <- fitting.Y(
    train_data = train_data, valid_data = valid_data,
    ymodel = "Y ~ M * tt * R + C",
    Cnames = Cnames,
    Mnames = Mnames,
    Yfamily = Yfamily,
    SL_library = SL_user
  )

  # directly estimating E(Y|t,r,C)
  y_out.trc <- fitting.Y(
    train_data = train_data, valid_data = valid_data,
    ymodel = "Y ~ tt * R + C",
    Cnames = Cnames,
    Mnames = Mnames,
    Yfamily = Yfamily,
    SL_library = SL_user
  )

  # fit pseudo outcome models

  v_out.t1r0c <- fitting.v_pseudo(
    train_data = train_data, valid_data = valid_data,
    vmodel = "v_pseudo ~ t1 + r0 + C",
    alltrain = FALSE,
    Cnames = Cnames,
    Mnames = Mnames,
    y_out.mtrc = y_out.mtrc,
    t_out.c = t_out.c,
    r_out.tc = r_out.tc,
    r_out.mtc = r_out.mtc,
    Yfamily = Yfamily,
    SL_library = SL_user
  )

  v_out.t0r0c <- fitting.v_pseudo(
    train_data = train_data, valid_data = valid_data,
    vmodel = "v_pseudo ~ t0 + r0 + C",
    alltrain = FALSE,
    Cnames = Cnames,
    Mnames = Mnames,
    y_out.mtrc = y_out.mtrc,
    t_out.c = t_out.c,
    r_out.tc = r_out.tc,
    r_out.mtc = r_out.mtc,
    Yfamily = Yfamily,
    SL_library = SL_user
  )

  # for validation data

  RMnames <- paste0("RM.", Mnames)
  ttMnames <- paste0("ttM.", Mnames)
  ttRMnames <- paste0("ttRM.", Mnames)

  # weights
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


    # scores calculating
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

    H_v.trc <- with(valid_data, {
      1*(tt==tval)*(R==0) / bound_propensity(tval.c_pred_valid * (1 - r1.tc_pred_valid))
    })
    eif_v <- with(valid_data, {
      1*(tt==tval)*(R==0) * H_v.trc * (v_pseudo_valid - v_pred_valid) /
        mean( 1*(tt==tval)*(R==0)*H_v.trc )
    })

    eif <- eif_y + eif_v + v_pred_valid

    # Y( tval ) | r=1
    y.trc_valid_r1 <- predict(y_out.trc$y_fit, valid_data_r1[, y_out.trc$y_fit$varNames], onlySL = TRUE)$pred

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


    # Y( tval ) | r=0
    y.trc_valid_r0 <- predict(y_out.trc$y_fit, valid_data_r0[, y_out.trc$y_fit$varNames], onlySL = TRUE)$pred

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


  # output ----
  eif_out <- list(
    components_tt1 = components_tt1,
    components_tt0 = components_tt0,
    # fold IDs
    fold = fold# origami::fold_index()
  )

  return(eif_out)
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
scale_to_unit <- function(vals) {
  vals_scaled <- (vals - min(vals)) / (max(vals) - min(vals))
  vals_scaled <- bound_precision(vals_scaled)
  return(vals_scaled)
}
scale_from_unit <- function(scaled_vals, max_orig, min_orig) {
  vals_orig <- (scaled_vals * (max_orig - min_orig)) + min_orig
  return(vals_orig)
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
