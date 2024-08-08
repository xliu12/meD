
#' Estimation of the mediated moderation and remaining moderation
#' @param data a \code{data.frame} containing the observed data.
#'    In "data", column "tt" is the treatment assignment ("tt" is coded as 0 for individuals in the control condition and as 1 for individuals in the treatment condition);
#'    column "R" is a dummy indicator of the subgroup status.
#'    column "Y" is the outcome.
#' @param Yname a character string of the name of the column in "data" that corresponds to outcome variable (Y).
#' @param ttname a character string of the name of the column in "data" that corresponds to a dummy coded treatment indicator (tt). 
#' @param Rname a character string of the name of the column in "data" that corresponds to a dummy coded subgroup indicator (R). 
#' @param Mnames a character vector of the names of the columns in "data" that correspond to mediators (M).
#' @param Cnames a character vector of the names of the columns in "data" that correspond to baseline covariates (C).
#' @param estimator a character specifying the estimator to use. Currently supported options include "onestep" and "tml", implementing the one-step estimator and targeted minimum loss estimator, respectively.
#' @param nuisance_estimation a character vector specifying estimation methods for the nuisance models. Currently supported options include methods in the "SuperLearner" package, available via the "listWrappers()" function (<https://cran.r-project.org/package=SuperLearner>). 
#' @param num_folds the number of folds to use for the cross-fitting procedure.
#'
#' @return A data frame containing the results (estimates and 95% confidence intervals) for the mediated moderation (MedMod) and remaining moderation (RemainMod), which decompose the total moderation (TotMod). The output also includes results for  the average potential outcomes in their definitions. 
#' These include the average potential outcomes of each subgroup under each treatment assignment (Yt1r1, Yt1r0, Yt0r1, Yt0r0), and the average potential outcomes of the comparison subgroup if the potential mediator distribution were shifted to be the same as that of the reference subgroup (Yt1r1.Mt1r0, Yt0r1.Mt0r0). In the output data frame, the confidence intervals are the columns with names ending with "_interval.1" and "_interval.2".

#'
#'
#' @export
#'
#' @examples
#'  # data(data_in)
#'  # data_in <- read.csv("data/data_in.csv", header = TRUE)
#'  # Mnames <- grep("Mdat", colnames(data_in), value = T)
#'  # Cnames <- grep("Cdat", colnames(data_in), value = T)
#'  # out <- MedMod(
#'  # data = data_in, Yname = "Y", ttname = "tt", Rname = "R", 
#'  # Mnames = Mnames,
#'  # Cnames = Cnames)

#'
#'


MedMod <- function(
    data, 
    Yname, ttname, Rname,
    Mnames,
    Cnames,
    estimator = c("onestep"), 
    nuisance_estimation = c("SL.glm", "SL.nnet"), 
    num_folds = 5
) {
  
  data_in <- data
  Yfamily <- ifelse(length(unique(data_in[[Yname]]))>2, "gaussian", "binomial")
  
  estimates <- meDreD(
    data_in = data_in, 
    Yname = Yname, 
    ttname = ttname, 
    Rname = Rname,
    Mnames = Mnames,
    Cnames = Cnames,
    Yfamily = Yfamily,
    estimator = estimator,
    Yt1only=FALSE,
    nuisance_estimation = nuisance_estimation, 
    num_folds = num_folds
  )
  
  estimates
}