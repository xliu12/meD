library(haven)
library(ggpattern)

library(tidyverse)
library(mice)
library(mitml)
library(miceadds)
library(fastDummies)
library(missForest)
options(dplyr.print_max = 100, pillar.print_max = 50, dplyr.width = Inf)

# load("~/Library/CloudStorage/Box-Box/Labs/Mediated/R_mediated/lift_tmp.RData")
# load("~/Library/CloudStorage/Box-Box/Labs/Mediated/R_mediated/lift_others2.RData")


load("~/Library/CloudStorage/Box-Box/Labs/Mediated/R_mediated/lift1228.RData") 

# Tables, Figures ----------------------------
# / theme figure -----------
fig.theme <- 
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


# tables ------------




# for r1 male ----
tmle_Rmale_out <- mlr.tmle_Rmale
crossfit_Rmale_out <- mlr.crossfit_Rmale

tab_Rmale <- data.frame(
  subgroupFM_r1 = "male", 
  subgroup_r1 = "boy",
  treatment_t1 = "treatment",
  outcome = "tobacco use in grade 12",
  mediated_moderator = "risk factors and deviant peer associations (grade 10)",
  effect = c("Yt1r1", "Yt1r0", "Yt1r1.Mt1r0", 
             "Yt0r1", "Yt0r0", "Yt0r1.Mt0r0", 
             "remaining moderation",  
             "mediated moderation", 
             "total moderation"), 
  tml_estimate = tmle_Rmale_out$tml_estimates,
  tml_95CI = paste0("(", round(tmle_Rmale_out$tml_intervals, 2)[1, ], ", ",
                    round(tmle_Rmale_out$tml_intervals, 2)[2, ], ")"), 
  tml_interval = t(tmle_Rmale_out$tml_intervals),
  onestep_estimate = crossfit_Rmale_out$estimates,
  onestep_95CI = paste0("(", round(crossfit_Rmale_out$z_intervals, 2)[1, ], ", ",
                        round(crossfit_Rmale_out$z_intervals, 2)[2, ], ")"), 
  onestep_interval = t(crossfit_Rmale_out$z_intervals),
  stderr = crossfit_Rmale_out$stderrors
)

tab_Rmale <- tab_Rmale %>%
  mutate(effect = factor(effect, levels = c("Yt1r1",  "Yt0r1", "Yt1r0", "Yt0r0", "Yt1r1.Mt1r0", "Yt0r1.Mt0r0", "total moderation", "mediated moderation",  "remaining moderation"))) %>%
  arrange(effect)

# write_csv(tab_Rmale, "Tables_Figures/tab_Rmale.csv")


# write_csv(tab_Rmale, "Tables_Figures/tab_Rmale_Mrisk8dvp8.csv")
write_csv(tab_Rmale, "Tables_Figures/tab_Rmale_Mrisk8dvp8_1228.csv")




# plot ----


# plotdf_Rfe <- tab_Rfe %>%
plotdf_Rmale <- tab_Rmale %>%
  filter(
    effect %in% c("Yt1r1","Yt1r0","Yt0r1","Yt0r0", "Yt1r1.Mt1r0", "Yt0r1.Mt0r0")
  ) %>% mutate(
    type = case_when(
      effect %in% c("Yt1r1","Yt1r0","Yt0r1","Yt0r0") ~ "As natural",
      effect %in% c("Yt1r1.Mt1r0", "Yt0r1.Mt0r0") ~ "If shifted to that of the other subgroup (with the same baseline covariate values)"
    ),
    subgroup = case_when(
      effect %in% c("Yt1r1.Mt1r0", "Yt1r1", "Yt0r1.Mt0r0","Yt0r1") ~ subgroup_r1,
      effect %in% c("Yt1r0", "Yt0r0") ~ ifelse(subgroup_r1=="boy", "girl", "boy")
      # effect %in% c("Yt1r0", "Yt0r0") ~ ifelse(subgroup_r1=="female", "male", "female")
    ),
    `Treatment assignment` = factor(ifelse(str_detect(effect, "Yt1"), "under the intervention condition", "under the control condition"), levels = c("under the intervention condition", "under the control condition"))
  )


# fig_female <- ggplot(data = plotdf_Rfe,
fig_male <- ggplot(data = plotdf_Rmale,
                   aes(y = onestep_estimate, x = subgroup, 
                 pattern = type #, color = type 
                 , fill = type
                 ) ) +
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
  # geom_point( aes(shape = type) ) +
  geom_errorbar( aes(ymax = onestep_interval.2, ymin = onestep_interval.1
                     ,linetype = type
                     ),
                 position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.9), #"dodge", 
                 width = 0.4, linewidth =0.6 ) +
  scale_y_continuous("Estimated risk of later tobacco use (Y)", limits = c(0, 1)) +
  # scale_x_discrete("Treatment assignment in grade 5") + 
  scale_x_discrete("Subgroup") + 
  # ggtitle("Tobacco use (Y)") +
  # annotate("text", x = 0.1, y = 0.95, label = "The results were adjusted for observed baseline covariates.") +
  scale_fill_discrete("Intermediate risk factors (Mediator)") +
  scale_color_discrete("Intermediate risk factors (Mediator)") +
  scale_linetype_discrete("Intermediate risk factors (Mediator)") +
  scale_pattern_discrete("Intermediate risk factors (Mediator)")+
  facet_grid(. ~ `Treatment assignment`, labeller = label_value) +
  fig.theme


# for r1 male ----
fig_male
# ggsave("fig_Rmale.pdf", path = "./Tables_Figures", width = 12, height = 7)
# ggsave("fig_Rmale_Mdvp456agg456.pdf", path = "./Tables_Figures", width = 12, height = 7)
# ggsave("fig_Rmale_Mdvp456.pdf", path = "./Tables_Figures", width = 12, height = 7)
# ggsave("fig_Rmale_Mdvp8.pdf", path = "./Tables_Figures", width = 7, height = 8)

# ggsave("fig_Rmale_Mrisk8dvp8.pdf", path = "./Tables_Figures", width =8, height = 7.5)


ggsave("fig_Rboy_Mrisk8dvp8.pdf", path = "./Tables_Figures", width =9, height = 7.5)


save.image("example_tabsfigs.RData")

# NOT Used ----------------
# for r1 female ----
tmle_Rfemale_out <- tmle_Rfe_Mrisk8dvp8 # tmle_Rfe_Mdvp8 # tmle_Rfe_Mdvp456 
crossfit_Rfemale_out <- crossfit_Rfe_Mrisk8dvp8 # crossfit_Rfe_Mdvp8 # crossfit_Rfe_Mdvp456agg456

tab_Rfe <- data.frame(
  subgroup_r1 = "female", treatment_t1 = "treatment",
  outcome = "tobacco use in grade 12",
  mediated_moderator = "risk factors and deviant peer association (grade 10)",
  effect = c("Yt1r1", "Yt1r0", "Yt1r1.Mt1r0", 
             "Yt0r1", "Yt0r0", "Yt0r1.Mt0r0", 
             "remaining moderation",  
             "mediated moderation", 
             "total moderation"), 
  tml_estimate = tmle_Rfemale_out$tml_estimates,
  tml_95CI = paste0("(", round(tmle_Rfemale_out$tml_intervals, 2)[1, ], ", ",
                    round(tmle_Rfemale_out$tml_intervals, 2)[2, ], ")"), 
  tml_interval = t(tmle_Rfemale_out$tml_intervals),
  onestep_estimate = crossfit_Rfemale_out$estimates,
  onestep_95CI = paste0("(", round(crossfit_Rfemale_out$z_intervals, 2)[1, ], ", ",
                        round(crossfit_Rfemale_out$z_intervals, 2)[2, ], ")"), 
  onestep_interval = t(crossfit_Rfemale_out$z_intervals),
  stderr = crossfit_Rfemale_out$stderrors
) 

tab_Rfe <- tab_Rfe %>%
  mutate(effect = factor(effect, levels = c("Yt1r1",  "Yt0r1", "Yt1r0", "Yt0r0", "Yt1r1.Mt1r0", "Yt0r1.Mt0r0", "total moderation",  "mediated moderation", "remaining moderation"))) %>%
  arrange(effect)

# write_csv(tab_Rfe, "Tables_Figures/tab_Rfemale.csv")
# write_csv(tab_Rfe, "Tables_Figures/tab_Rfemale_Mdvp456agg456.csv")
# write_csv(tab_Rfe, "Tables_Figures/tab_Rfemale_Mdvp456.csv")
# write_csv(tab_Rfe, "Tables_Figures/tab_Rfemale_Mdvp8.csv")
write_csv(tab_Rfe, "Tables_Figures/tab_Rfemale_Mrisk8dvp8.csv")

# for r1 female ----
fig_female

# ggsave("fig_Rfemale.pdf", path = "./Tables_Figures", width = 12, height = 7)
# ggsave("fig_Rfemale_Mdvp456agg456.pdf", path = "./Tables_Figures", width = 12, height = 7)
# ggsave("fig_Rfemale_Mdvp456.pdf", path = "./Tables_Figures", width = 12, height = 7)
ggsave("fig_Rfemale_Mdvp8.pdf", path = "./Tables_Figures", width = 7, height = 8)
