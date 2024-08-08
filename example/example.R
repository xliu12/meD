# Data ----------------------

# import the data
load("data_example.RData")

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

# run ---------------------------------

devtools::install_github(repo = "xliu12/meD", subdir = "meD")
library(meD)
# devtools::load_all("meD") # download the package and then we can install it in this way as well
library(SuperLearner)
# see available methods to estimate the models (e.g., SL.glm runs generalized linear model)
SuperLearner::listWrappers()

estimates <- MedMod(
  data = data_example,
  Yname = "Y", ttname = "tt", Rname = "R", Mnames = Mnames, 
  Cnames = Cnames,
  estimator = c("onestep"), 
  nuisance_estimation = c("SL.glm", "SL.nnet"),
  num_folds = 5
)
estimates
# focal/reference subgroup labels
focal_subgroup_r1 <- "focal"
ref_subgroup_r0 <- "reference"
# treatment assignment labels
intervention_t1 <- "under the intervention condition"
control_t0 <- "under the control condition"

plot_results <- estimates %>%
  filter(effect %in% c("Yt1r1","Yt1r0","Yt0r1","Yt0r0", "Yt1r1.Mt1r0", "Yt0r1.Mt0r0")) %>% 
  mutate(
    type = case_when(
      effect %in% c("Yt1r1","Yt1r0","Yt0r1","Yt0r0") ~ "As natural",
      effect %in% c("Yt1r1.Mt1r0", "Yt0r1.Mt0r0") ~ "If shifted to that of the other subgroup (with the same baseline covariate values)"),
    subgroup = case_when(
      effect %in% c("Yt1r1.Mt1r0", "Yt1r1", "Yt0r1.Mt0r0","Yt0r1") ~ focal_subgroup_r1,
      effect %in% c("Yt1r0", "Yt0r0") ~ ref_subgroup_r0
    ),
    `Treatment assignment` = factor(ifelse(str_detect(effect, "Yt1"), intervention_t1, control_t0), levels = c(intervention_t1, control_t0))
  )


library(ggpattern)
library(glue)
# outcome and mediator labels
Yname <- "risk of tobacco use outcome"
Mname <- "Intermediate risk factors (Mediator)"

fig_results <- plot_results %>% 
  ggplot(aes(y = onestep_estimate, x = subgroup, pattern = type, fill = type) ) +
  geom_bar_pattern(
    aes(linetype = type),
    fill = "white", 
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.05,
    pattern_spacing = 0.1,
    pattern_key_scale_factor = 0.3,
    alpha = 0.7, stat = "identity", color = "black", linewidth =0.4,
    position = position_dodge()) +
  geom_errorbar( aes(ymax = onestep_interval.2, ymin = onestep_interval.1
                     ,linetype = type),
                 position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.9), 
                 width = 0.4, linewidth =0.6 ) +
  scale_y_continuous(glue("Estimate of {Yname}")) +
  scale_x_discrete("Subgroup") + 
  scale_fill_discrete(Mname) +
  scale_color_discrete(Mname) +
  scale_linetype_discrete(Mname) +
  scale_pattern_discrete(Mname)+
  facet_grid(. ~ `Treatment assignment`, labeller = label_value) +
  # theme
  theme_bw() +
  theme(panel.grid.minor = element_line(linewidth = 0),
        panel.grid.major.x = element_line(linewidth = 0),
        panel.grid.major.y = element_line(linewidth = 0.5, lineend = "round", color = "grey", linetype = "longdash"),
        # axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0),
        plot.title = element_text(size = 12, face = "plain"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.position = "bottom",
        legend.spacing.x = unit(0.2, "mm"),
        legend.key.height = unit(10, "mm"),
        legend.key.width = unit(10, "mm"),
        legend.key.size = unit(10, "mm"))

fig_results

