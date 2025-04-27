################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Generate Fig 4.2 (Sec 4.4)
# @Author - Katherine Jia
################################################################################
# <I> Calculate relative changes in incidence, ANC prev, and new DX rate -------

for(d in 1:10) {
  
  scenario <- all_scenarios[d]
  
  # (1.1) Get model outputs for each scenario
  
  data <- get(paste0("true_dynamics_", scenario))
  
  # (1.2) Calculate relative changes
  ## Create an empty df
  
  relative_change <- data.frame(
    step = 0:40,
    rel_change_prev = rep(0, 41),
    rel_change_rate_new_diag = rep(0, 41),
    rel_change_lambda = rep(0, 41)
  )
  
  ## Fill in the results
  
  for (j in 1:40) {
    relative_change$rel_change_prev[1 + j] <- (data$ANC_prevalence[1 + j] - data$ANC_prevalence[1]) / data$ANC_prevalence[1] * 100
    relative_change$rel_change_rate_new_diag[1 + j] <- (data$rate_self_report_new_diag[1 + j] - data$rate_self_report_new_diag[1]) / data$rate_self_report_new_diag[1] * 100
    relative_change$rel_change_lambda[1 + j] <- (data$lambda[1 + j] - data$lambda[1]) / data$lambda[1] * 100
  }
  
  # (1.3) Combine the relative changes with model outputs
  
  data_relative_change <- data %>%
    left_join(relative_change, by="step")
  
  assign(paste0("data_relative_change_", scenario), data_relative_change)
  
}

# <II> Plot Fig 4.2 ------------------------------------------------------------
# NOTE: Only Scenario 0 (with 0% vs 50% non-disclosure) are plotted
# (2.1) Make a data frame

plot_relative_change_scenario0_0_vs_50_df <- data_relative_change_scenario0_00 %>%
  select(step, rel_change_prev, rel_change_rate_new_diag, rel_change_lambda) %>%
  left_join(data_relative_change_scenario0_05 %>% select(step, rel_change_rate_new_diag_05 = rel_change_rate_new_diag), by="step") %>%
  pivot_longer(-step) %>%
  mutate(label = case_when(name=="rel_change_prev"~"3. ANC HIV prevalence",
                           name=="rel_change_rate_new_diag"~"2a. Reported new diagnoses \n rate (full disclosure)",
                           name=="rel_change_rate_new_diag_05"~"2b. Reported new diagnoses \n rate (50% non-disclosure)",
                           name=="rel_change_lambda"~"1. Incidence"))

# (2.2) Plot

plot_relative_change_scenario0_0_vs_50 <- ggplot(plot_relative_change_scenario0_0_vs_50_df) + 
  geom_line(aes(x = step/4, y = value, group = label, lty = label)) + 
  lims(x = c(0, 14.9),
       y = c(-100, 0.1)) +
  scale_linetype_manual(guide="none", values=c("solid", "twodash", "longdash", "dotted")) +
  geom_text_repel(data = subset(plot_relative_change_scenario0_0_vs_50_df, step==40), 
            aes(label = label, x = step/4, y = value), 
            segment.color = 'transparent',
            hjust=0,
            direction = "y",
            size=4,
            na.rm = TRUE) +
  scale_colour_discrete(guide = 'none')  +    
  theme_minimal() +
  labs(x = "year",
       y = "% change compared to baseline",
       col = " ")  +
theme(
  legend.position = "bottom",
  legend.direction = "vertical",
  text = element_text(size = 12),
  legend.text = element_text(size = 12)) 

# (2.3) Save plot

  ggsave("3_figures/Fig4_2.png", plot_relative_change_scenario0_0_vs_50, width = 8, height = 7.5, dpi=300)

# (2.4) Report first time to observe a 25% drop 
# ANC HIV prevalence 
  
time_to_25_drop_ANC_prev <- data_relative_change_scenario0_00 %>% 
  filter(rel_change_prev < - 25 ) %>%
  first() %>%
  pull(step) 

(time_to_25_drop_ANC_prev - 4) / 4

# New DX rate (full disclosure)

time_to_25_drop_new_diag <- data_relative_change_scenario0_00 %>% 
  filter(rel_change_rate_new_diag < - 25 ) %>%
  first() %>%
  pull(step) 

(time_to_25_drop_new_diag - 4) / 4

# New DX rate (50% non-disclosure)

time_to_25_drop_new_diag_05 <- data_relative_change_scenario0_05 %>% 
  filter(rel_change_rate_new_diag < - 25 ) %>%
  first() %>%
  pull(step) 

(time_to_25_drop_new_diag_05 - 4 + 1)/4

