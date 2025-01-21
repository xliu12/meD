# meD
The R package `meD` is used to assess causal mediation mechanisms underlying the moderation in the effect of treatment between two different subgroups using the multiply robust estimators incorporating machine learning methods.

Using the R package, a user can decompose a moderated effect of treatment into (i) the mediated moderation through a mediator and (ii) the remaining moderation not through the mediator. For details, please see the accompanying methods paper:

Liu, X., Eddy, J. M., & Martinez, C. R. (2025). Causal estimands and multiply robust estimation of mediated-moderation. Multivariate Behavioral Research. https://doi.org/10.1080/00273171.2024.2444949


## Installation

To install `meD` directly from Github:

```
remotes::install_github(repo = "xliu12/meD", subdir = "meD")
```


## Example

```
# Data ----------------------

# import the data
data(data_example)

head(data_example)

# tt: Treatment. In the example, the treatment is the treatment assignment (1 = intervention; 0 = control).
# R: Subgroup. In the example, the subgroup is gender (1 = boy; 0 = girl).
# Mdat.M_mediator1, Mdat.M_mediator2: Mediator. In the example, the mediator is a vector of intermediate risk factors. 
# Y: Outcome. In the example, the outcome is tobacco use (1 = any use; 0 = no use).
# Cdat.V1,...,Cdat.V5: Baseline (i.e., pre-treatment) covariate. 

Mnames <- grep("^Mdat", colnames(data_example), value = TRUE)
Cnames <- grep("^Cdat", colnames(data_example), value = TRUE)

library(tidyverse)
library(mvtnorm)

library(SuperLearner)
# see available methods to estimate the models (e.g., SL.glm runs generalized linear model)
SuperLearner::listWrappers()

# run ---------------------------------

library(meD)

estimates <- MedMod(
  data = data_example,
  Yname = "Y", ttname = "tt", Rname = "R", Mnames = Mnames, 
  Cnames = Cnames,
  estimator = c("onestep"), 
  nuisance_estimation = c("SL.glm", "SL.xgboost", "SL.ranger"), 
  num_folds = 5
)
estimates

```
