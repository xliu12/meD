
# Calculating the scores for each individual observation 



eif.onefold <- function(v = 1,  #fold,
                        folds,
                        data_in,
                        Cnames,
                        Mnames,
                        Fit,
                        Yfamily = "binomial",
                        cv = TRUE
) {
  if(cv==TRUE) {
    # # make training and validation data
    # train_data <- origami::training(data_in)
    # valid_data <- origami::validation(data_in)
    # fold <- origami::fold_index()
    train_data <- data_in[folds[[v]]$training_set, ]
    valid_data <- data_in[folds[[v]]$validation_set, ]
    # fold.ind <- origami::fold_index()
    fold.ind <- v
    valid_set <- folds[[v]]$validation_set
    valid_set_id <- data_in[valid_set, "id"]
  }
  if(cv==FALSE) {
    train_data <- valid_data <- data_in
    fold.ind <- 0
    valid_set <- 1:nrow(valid_data)
    valid_set_id <- data_in[valid_set, "id"]
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
    H_y.mtrc <- bound_ipw(H_y.mtrc)
    
    y.mtrc_valid_r1 <- predict(y_out.mtrc$y_fit, valid_data_r1[, y_out.mtrc$y_fit$varNames], onlySL = TRUE)$pred
    
    
    # scores calculating
    # west_y <- with(valid_data, {
    #   1*(tt==tval)*(R==1) * H_y.mtrc * (Y) /
    #     mean( 1*(tt==tval)*(R==1) * H_y.mtrc )
    # })
    
    eif_y <- with(valid_data, {
      1*(tt==tval)*(R==1) * H_y.mtrc * (Y - y.mtrc_valid_r1) /
        mean( 1*(tt==tval)*(R==1) * H_y.mtrc )
    })
    
    if (tval == 1) {
      v_pseudo_valid <- v_out.t1r0c$v_pseudo_valid
      v_pred_valid <- v_out.t1r0c$v_pred_valid
    }
    if (tval == 0) {
      v_pseudo_valid <- v_out.t0r0c$v_pseudo_valid
      v_pred_valid <- v_out.t0r0c$v_pred_valid
    }
    
    H_v.trc <- with(valid_data, {
      1*(tt==tval)*(R==0) / bound_propensity(tval.c_pred_valid * (1 - r1.tc_pred_valid))
    })
    H_v.trc <- bound_ipw(H_v.trc)
    
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
    H_y1.trc <- bound_ipw(H_y1.trc)
    
    # doubly robust estimator
    eif_y1_c <- with(valid_data, {
      1*(tt==tval)*(R==1) * H_y1.trc * (Y - y.trc_valid_r1) /
        mean( 1*(tt==tval)*(R==1) * H_y1.trc )
    })
    
    eif_y1 <- eif_y1_c + y.trc_valid_r1
    
    # Y( tval ) | r=0
    y.trc_valid_r0 <- predict(y_out.trc$y_fit, valid_data_r0[, y_out.trc$y_fit$varNames], onlySL = TRUE)$pred
    
    H_y0.trc <- with(valid_data, {
      1*(tt==tval)*(R==0) / bound_propensity(tval.c_pred_valid * (1-r1.tc_pred_valid))
    })
    H_y0.trc <- bound_ipw(H_y0.trc)
    
    # doubly robust estimator
    eif_y0_c <- with(valid_data, {
      1*(tt==tval)*(R==0) * H_y0.trc * (Y - y.trc_valid_r0) /
        mean( 1*(tt==tval)*(R==0) * H_y0.trc )
    })
    eif_y0 <- eif_y0_c + y.trc_valid_r0
    
    
    
    components <- data.frame(
      # valid_id = valid_data$id,
      fold.ind = rep(v, nrow(valid_data)),
      valid_set = valid_set,
      valid_set_id = valid_set_id,
      # nuisance
      r1.tc_pred_valid = r1.tc_pred_valid,
      r1.mtc_pred_valid = r1.mtc_pred_valid,
      t1.c_pred_valid = t1.c_pred_valid,
      
      y.mtrc_valid_r1 = y.mtrc_valid_r1,
      y.trc_valid_r0 = y.trc_valid_r0,
      y.trc_valid_r1 = y.trc_valid_r1,
      v_pred_valid = v_pred_valid,
      
      # Y(1, Gm(tt=1)|0 ), eif & est
      eif = eif, 
      # eif for Y(tt=1) | r1
      eif_y1 = eif_y1, 
      # eif for Y(tt=1) | r0
      eif_y0 = eif_y0 
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
    # fold = fold# origami::fold_index()
    fold = v
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
bound_ipw <- function(vals, bound=1) {
  qupper <- quantile(vals, prob = bound)
  vals[vals>qupper] <- qupper
  vals
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
