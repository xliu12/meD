

reg.medmod <- function(data_in, 
                       Yname, ttname, Rname,
                       Mnames,
                       Cnames, 
                       Yt1only = FALSE) {
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
  
  Mreg <- lm(M ~ ., data = data.frame(M=data_in[[Mnames]], 
                                 data_in[, c("tt", "R", "ttR", Cnames)]))
  hat.btt <- coef(Mreg)
  
  Yreg <- lm(Y ~ ., data = data.frame(Y=data_in[[Yname]], 
                                      data_in[, c(Mnames, "tt", "R", "ttR", ttMnames, RMnames, ttRMnames, Cnames)]))
  hat.bm <- coef(Yreg)
  
  
  # linear Y models
  meanC_dat <- data_in[, Cnames] # matrix(rep(0, length(Cnames)), nrow = 1)
  
  
  # parametric bootstrap interval -----
  nboot <- 500
  set.seed(12345)
  each.boot <- function(boot=1) {
    if (boot==0) {hat.btt <- coef(Mreg)
    hat.bm <- coef(Yreg) }
    if (boot > 0) {
      hat.btt <- drop(rmvnorm(1, hat.btt, vcov(Mreg)) )
      hat.bm <- drop(rmvnorm(1, hat.bm, vcov(Yreg)) )
    }
    
    # when tt = 1 
    # Yt1r1
    t1r1_dat <- model.matrix(~., data.frame(tt=1, R=1, ttR = 1, C_dat=meanC_dat))
    meanMt1r1 <- as.numeric(t1r1_dat %*% hat.btt)
    meanMt1r1_dat <- model.matrix(~., data.frame(M=meanMt1r1, tt=1, R=1, ttR = 1, ttM = meanMt1r1*1, RM = 1*meanMt1r1, ttRM = 1*1*meanMt1r1, C_dat=meanC_dat))
    meanYt1r1 <- as.numeric(meanMt1r1_dat %*% hat.bm) %>% mean(.)
    
    # Yt1r1.Mt1r0
    t1r0_dat <- model.matrix(~., data.frame(tt=1, R=0, ttR = 0, C_dat=meanC_dat))
    meanMt1r0 <- as.numeric(t1r0_dat %*% hat.btt)
    meanMt1r0_dat <- model.matrix(~., data.frame(M=meanMt1r0, tt=1, R=1, ttR = 1, ttM = meanMt1r0*1, RM = 1*meanMt1r0, ttRM = 1*1*meanMt1r0, C_dat=meanC_dat))
    meanYt1r1.Mt1r0 <- as.numeric(meanMt1r0_dat %*% hat.bm) %>% mean(.)
    
    # Yt1r0
    t1r0_dat <- model.matrix(~., data.frame(tt=1, R=0, ttR = 0, C_dat=meanC_dat))
    meanMt1r0 <- as.numeric(t1r0_dat %*% hat.btt)
    meanMt1r0_dat <- model.matrix(~., data.frame(M=meanMt1r0, tt=1, R=0, ttR = 0, ttM = meanMt1r0*1, RM = 0*meanMt1r0, ttRM = 1*0*meanMt1r0, C_dat=meanC_dat))
    meanYt1r0 <- as.numeric(meanMt1r0_dat %*% hat.bm) %>% mean(.)
    
    # when tt = 0 
    # Yt0r1
    t0r1_dat <- model.matrix(~., data.frame(tt=0, R=1, ttR = 0, C_dat=meanC_dat))
    meanMt0r1 <- as.numeric(t0r1_dat %*% hat.btt)
    meanMt0r1_dat <- model.matrix(~., data.frame(M=meanMt0r1, tt=0, R=1, ttR = 0, ttM = meanMt0r1*0, RM = 1*meanMt0r1, ttRM = 0*1*meanMt0r1, C_dat=meanC_dat))
    meanYt0r1 <- as.numeric(meanMt0r1_dat %*% hat.bm) %>% mean(.)
    
    # Yt0r1.Mt0r0
    t0r0_dat <- model.matrix(~., data.frame(tt=0, R=0, ttR = 0, C_dat=meanC_dat))
    meanMt0r0 <- as.numeric(t0r0_dat %*% hat.btt)
    meanMt0r0_dat <- model.matrix(~., data.frame(M=meanMt0r0, tt=0, R=1, ttR = 0, ttM = meanMt0r0*0, RM = 1*meanMt0r0, ttRM = 0*1*meanMt0r0, C_dat=meanC_dat))
    meanYt0r1.Mt0r0 <- as.numeric(meanMt0r0_dat %*% hat.bm) %>% mean(.)
    ## when balanced_M_nott == TRUE, we have meanMt0r0=meanMt0r1, and therefore  meanYt0r1=meanYt0r1.Mt0r0
    
    # Yt0r0
    t0r0_dat <- model.matrix(~., data.frame(tt=0, R=0, ttR = 0, C_dat=meanC_dat))
    meanMt0r0 <- as.numeric(t0r0_dat %*% hat.btt)
    meanMt0r0_dat <- model.matrix(~., data.frame(M=meanMt0r0, tt=0, R=0, ttR = 0, ttM = meanMt0r0*0, RM = 0*meanMt0r0, ttRM = 0*0*meanMt0r0, C_dat=meanC_dat))
    meanYt0r0 <- as.numeric(meanMt0r0_dat %*% hat.bm) %>% mean(.)
    
    estimates <- data.frame(Yt1r1 = meanYt1r1, Yt1r0 = meanYt1r0, Yt1r1.Mt1r0 = meanYt1r1.Mt1r0, 
                            Yt0r1 = meanYt0r1, Yt0r0 = meanYt0r0, Yt0r1.Mt0r0 = meanYt0r1.Mt0r0)
    estimates
  }
  set.seed(12345)
  boots <- purrr::map_df(1:4, each.boot) 

  
  reg_estimate <- effect(each.boot(boot = 0), Yt1only = Yt1only)
  
  boots_eff <- effect(boots, Yt1only = Yt1only)
  
  reg_result <- data.frame(effect = colnames(reg_estimate),
    reg_estimate=t(reg_estimate), 
  reg_interval.1 = apply(boots_eff, 2, quantile, 0.025, na.rm=TRUE),
  reg_interval.2 = apply(boots_eff, 2, quantile, 0.975, na.rm=TRUE))
  
  reg_result
}