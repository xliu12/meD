# The estimation in the illustrative example was implemented using the R package.
devtools::install_github(repo = "xliu12/meD", subdir = "pkg/meD", build_manual = TRUE)
?meD::meDreD
# The dataset used comes from the Linking the Interests of Families and Teachers (LIFT) Study, which is not available.
# data_in

estimates <- meDreD(
  data_in, # data frame containing the variables;
             # column "tt" is the treatment assignment ("tt" is 0 for the control condition and 1 for the treatment condition);
             # column "R" is a dummy indicator of the subgroup status;
             # column "Y" is the outcome.
  Mnames, # column names of mediators
  Cnames, # column names of covariates
  nuisance_estimation = c("SL.glm", "SL.xgboost", "SL.ranger"), # the function names included in the "SuperLearner" package; <https://cran.r-project.org/package=SuperLearner>
  Yfamily = "binomial", # binary outcome
  estimator = c("onestep") # name of the estimator
)
