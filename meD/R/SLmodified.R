
SL.gam.modified <- function (Y, X, newX, family, obsWeights, deg.gam = 2, cts.num = 4, 
                             ...) 
{
  if ("mgcv" %in% loadedNamespaces()) 
    warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
  # cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  cts.x <- grepl("C", colnames(X))
  if (sum(!cts.x) > 0) {
    gam.model <- as.formula(paste("Y~", paste(paste("s(", 
                                                    colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, 
                                                    ")", sep = ""), collapse = "+"), "+", paste(colnames(X[, 
                                                                                                           !cts.x, drop = FALSE]), collapse = "+")))
  }
  else {
    gam.model <- as.formula(paste("Y~", paste(paste("s(", 
                                                    colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, 
                                                    ")", sep = ""), collapse = "+")))
  }
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(paste("Y~", paste(colnames(X), 
                                              collapse = "+"), sep = ""))
  }
  fit.gam <- gam::gam(gam.model, data = X, family = family, 
                      control = gam::gam.control(maxit = 50, bf.maxit = 50), 
                      weights = obsWeights)
  if (packageVersion("gam") >= "1.15") {
    pred <- gam::predict.Gam(fit.gam, newdata = newX, type = "response")
  }
  else {
    stop("This SL.gam wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
  }
  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam")
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


SL.ranger.modified <- function(...) {
  SL.ranger(..., num.trees = 100)
}


SL.xgb <- function(...) {
  SL.xgboost(..., ntrees = 100)
}
SL.xgboost.modified <- function(...) {
  SL.xgboost(..., ntrees = 100)
}
