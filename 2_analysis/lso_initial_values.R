################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Specify pre-baseline incidence and initial prevalence for DREAMS and non-DREAMS districts in Lesotho
# See Appendix C.5 for details
# @Author - Katherine Jia
################################################################################

# <I> Import national AGYW incidence (per 100 person-year) estimates from UNAIDS ---------------------
# Note: publicly available on https://aidsinfo.unaids.org/

lso_15to24_national <- read_csv("1_data/lso_15to24_national.csv")               # incidence per 100 person-year

New_HIV_Infections_National <- lso_15to24_national %>%
  mutate(time = as.numeric(Year) + 0.25) %>%
  filter(time >= 2006, time < 2021) %>%
  dplyr::select(time, lesotho_2024_02_20)

# Approx the national AGYW incidence by fitting a linear time trend

New_HIV_Infections_National_lm <- lm(lesotho_2024_02_20~time, data=New_HIV_Infections_National)

New_HIV_Infections_National_lm_pred <- data.frame(time=rep(2006:2020, each=4)+1:4/4,
                                                  linear_female_15_24_incidence_estimates = predict(New_HIV_Infections_National_lm, newdata = list(time=rep(2006:2020, each=4)+1:4/4)))

# <II> Approx incidence in DREAMS and non-DREAMS district from 2006 to 2020 ----
# (2.1) Import PHIA estimates on prevalence and % viraemia among AGYW ----------

lso_survey_hiv_indicators <- read_csv("1_data/lso_survey_hiv_indicators.csv")   # Note: dataset not publicly available

# (2.2) Calculate the prevalence ratio of viraemia among AGYW ------------------
# Extract indicator

district_prev <- lso_survey_hiv_indicators %>%
  filter(sex=="female",
         area_id %in% c("LSO_1_1",
                        "LSO_1_3",
                        "LSO_1_4",
                        "LSO_1_5",
                        "LSO_1_6"),
         indicator %in% c("prevalence", "viral_suppression_plhiv"),
         age_group=="Y015_024")

# Calculate `n_unsuppressed` (# AGYW with viraemia by districts) in 2017 

lso_prop_viraemia <- district_prev %>%
  dplyr::select(survey_id, survey_mid_calendar_quarter, area_id, area_name, sex, age_group, indicator, estimate) %>%
  filter(survey_id=="LSO2017PHIA",
         indicator %in% c("prevalence", "viral_suppression_plhiv")) %>% 
  pivot_wider(names_from = indicator, values_from = estimate) %>%
  mutate(unsuppressed = prevalence * (1-viral_suppression_plhiv),
         female_15_24 = case_when(area_name=="Maseru"~56872,
                                  area_name=="Leribe"~33720,
                                  area_name=="Berea"~26463,
                                  area_name=="Mafeteng"~16674,
                                  area_name=="Mohale's Hoek"~15251),
         group = case_when(area_name=="Maseru"~"DREAMS",
                           area_name=="Leribe"~"non-DREAMS",
                           area_name=="Berea"~"DREAMS",
                           area_name=="Mafeteng"~"non-DREAMS",
                           area_name=="Mohale's Hoek"~"non-DREAMS")) %>%
  mutate(n_unsuppressed = prevalence * (1-viral_suppression_plhiv) * female_15_24) 

# Calculate the overall prevalence of viraemia among AGYW in (non-)DREAMS districts

prop_DREAMS <- sum(lso_prop_viraemia[lso_prop_viraemia$group=="DREAMS", "n_unsuppressed"])/sum(lso_prop_viraemia[lso_prop_viraemia$group=="DREAMS", "female_15_24"])

prop_nonDREAMS <- sum(lso_prop_viraemia[lso_prop_viraemia$group=="non-DREAMS", "n_unsuppressed"])/sum(lso_prop_viraemia[lso_prop_viraemia$group=="non-DREAMS", "female_15_24"])

# Calculate the national prevalence of viraemia among AGYW in Lesotho

prop_overall <- lso_survey_hiv_indicators %>%
  filter(survey_id=="LSO2017PHIA",
         sex=="female",
         area_id=="LSO",
         indicator %in% c("prevalence", "viral_suppression_plhiv"),
         age_group=="Y015_024") %>%
  dplyr::select(survey_id, survey_mid_calendar_quarter, area_id, area_name, sex, age_group, indicator, estimate) %>%
  pivot_wider(names_from = indicator, values_from = estimate) %>%
  mutate(unsuppressed = prevalence * (1-viral_suppression_plhiv),
         female_15_24 = 200243) %>%
  mutate(n_unsuppressed = prevalence * (1-viral_suppression_plhiv) * female_15_24,
         prop_lso = n_unsuppressed / 200243) %>%
  pull(prop_lso)

# (2.3) scaled national AGYW incidence by the ratio of viraemia prevalence in AGYW within DREAMS (or non-DREAMS) districts to the national viraemia prevalence in AGYW in 2017

assumed_incidence <-
  New_HIV_Infections_National_lm_pred %>%
  mutate(incidence_AGYW_DREAMS = linear_female_15_24_incidence_estimates * prop_DREAMS / prop_overall,
         incidence_AGYW_nonDREAMS = linear_female_15_24_incidence_estimates * prop_nonDREAMS / prop_overall) %>%
  filter(time < 2020.75)

# (2.4) Plot Fig. C3 -----------------------------------------------------------

plot_lso_inc <- ggplot() +
  geom_point(data = New_HIV_Infections_National, aes(x = time,
                                           y = lesotho_2024_02_20*10)) +
  geom_line(data = assumed_incidence,
            aes(x = time,
                y = linear_female_15_24_incidence_estimates*10,
                lty = "1. AGYW (national)")) +
  geom_line(data = assumed_incidence,
            aes(
              x = time,
              y = incidence_AGYW_DREAMS*10,
              lty = "2. AGYW (DREAMS)"
            )) +
  geom_line(data = assumed_incidence,
            aes(
              x = time,
              y = incidence_AGYW_nonDREAMS*10,
              lty = "3. AGYW (non-DREAMS)"
            )) +
  theme_minimal() +
  labs(y = "incidence per 1,000 person-year",
       x = "",
       col="",
       lty="") +
  lims(x = c(2004.75, 2020.5),
       y = c(0, 50)) 

ggsave("3_figures/FigC_3.png", plot_lso_inc, width = 8, height = 5, units = "in", dpi=300)

# <III> Run the model forward to get prevalence at 2016 ------------------------
## Get number of timesteps from 2006Q1 to 2020Q2

N_t <- assumed_incidence %>%
  filter(time < 2020.75) %>%
  pull(time) %>%
  length()

arms <- c("intervention", "control")

# (4.2) Simulate forward from 2006Q1 to 2020Q2 ---------------------------------

for (n in 1:length(arms)) {
    
    arm <- arms[n]
    
    if (arm == "intervention") {
      lambda <- c(assumed_incidence$incidence_AGYW_DREAMS/400, NA)              # convert the unit from per 100 person-year to per person-quarter
      pop_size <- 56872 + 26463
    }
    
    if (arm == "control") {
      lambda <- c(assumed_incidence$incidence_AGYW_nonDREAMS/400, NA)           # convert the unit from per 100 person-year to per person-quarter
      pop_size <- 33720 + 16674 + 15251
    }
    
    non_disc_DX <- approx(c(0, N_t),
                          c(0.1, 0.1),
                          n = N_t + 1)$y
    
    delta <- approx(c(0, N_t),
                    c(0.1914059, 0.1914059),
                    n = N_t + 1)$y
    
    tau_pos <- rep(0.0143811 * 1.55, N_t + 1)
    tau <- rep(0.0143811, N_t + 1)
    mu_pos <- 0.025
    
    # create empty model compartments
    S <-
      UnDX <-
      DX <-
      ANC_TOTAL <-
      ANC_POS <-
      ANC_NEW_POS_TRUE <-
      ANC_NEW_POS_OBSERVED <- ANC_KN_POS_OBSERVED <- rep(0, N_t + 1)
    
    # (4.2.2) initial values (conditional probabilities)
    
    if (arm == "intervention") {
      # prevalence in 2006 is approximated by estimates from 2004 DHS 
      
      weighted_prev_2004 <- (district_prev[district_prev$survey_id=="LSO2004DHS"&district_prev$area_id=="LSO_1_1","estimate"]*56872 + 
                             district_prev[district_prev$survey_id=="LSO2004DHS"&district_prev$area_id=="LSO_1_4","estimate"]*26463)/(56872+26463)
      initial_pos <- as.numeric(weighted_prev_2004)
      initial_DX_given_pos <- 0 
    }
    if (arm == "control") {
      # prevalence in 2006 is approximated by DHS 2004
      
      weighted_prev_2004 <- (district_prev[district_prev$survey_id=="LSO2004DHS"&district_prev$area_id=="LSO_1_3","estimate"]*33720 +
                             district_prev[district_prev$survey_id=="LSO2004DHS"&district_prev$area_id=="LSO_1_5","estimate"]*16674 +
                             district_prev[district_prev$survey_id=="LSO2004DHS"&district_prev$area_id=="LSO_1_6","estimate"]*15251)/(33720+16674+15251)
      initial_pos <- as.numeric(weighted_prev_2004)
      initial_DX_given_pos <- 0
    }
    
    # (4.2.3) initial values (joint probabilities)
    
    initial_UnDX <- initial_pos * (1 - initial_DX_given_pos)
    initial_DX <- initial_pos * initial_DX_given_pos
    
    S[1] <- pop_size * (1 - initial_pos)
    UnDX[1] <- pop_size * initial_UnDX
    DX[1] <- pop_size * initial_DX
    ANC_TOTAL[1] <- tau[1] * S[1] + tau_pos[1] * (UnDX[1] + DX[1])
    ANC_POS[1] <- tau_pos[1] * (UnDX[1] + DX[1])
    ANC_NEW_POS_OBSERVED[1] <- tau_pos[1] * (UnDX[1] + non_disc_DX[1] * DX[1])
    ANC_KN_POS_OBSERVED[1] <- tau_pos[1] * (1 - non_disc_DX[1]) * DX[1]
    ANC_NEW_POS_TRUE[1] <- tau_pos[1] * UnDX[1]
    
    # (4.3) solve the model 
    
    for (i in 1:N_t) {
      
      S[i + 1] <- S[i] - lambda[i] * S[i] + mu_pos * (UnDX[i] + DX[i])
      
      UnDX[i + 1] <- UnDX[i] + lambda[i] * S[i] - (delta[i] + tau_pos[i]) * UnDX[i] - mu_pos * UnDX[i]
      
      DX[i + 1] <- DX[i] + (delta[i] + tau_pos[i]) * UnDX[i] - mu_pos * DX[i]
      
      ANC_TOTAL[i + 1] <- tau[i + 1] * S[i + 1] + tau_pos[i + 1] * (UnDX[i + 1] + DX[i + 1])
      
      ANC_POS[i + 1] <- tau_pos[i + 1] * (UnDX[i + 1] + DX[i + 1])
      
      ANC_NEW_POS_OBSERVED[i + 1] <- tau_pos[i + 1] * (UnDX[i + 1] + non_disc_DX[i + 1] * DX[i + 1])
      
      ANC_KN_POS_OBSERVED[i + 1] <- tau_pos[i + 1] * (1 - non_disc_DX[i + 1]) * DX[i + 1]
      
      ANC_NEW_POS_TRUE[i + 1] <- tau_pos[i + 1] * UnDX[i + 1]
      
    }
    
    true_dynamics <- data.frame(
      step = 0:N_t,
      lambda = lambda,
      S = S,
      UnDX = UnDX,
      DX = DX,
      ANC_TOTAL = ANC_TOTAL,
      ANC_POS = ANC_POS,
      ANC_NEW_POS_OBSERVED = ANC_NEW_POS_OBSERVED,
      ANC_NEW_POS_TRUE = ANC_NEW_POS_TRUE,
      ANC_prevalence = ANC_POS / ANC_TOTAL,
      rate_self_report_new_diag = ANC_NEW_POS_OBSERVED / (ANC_TOTAL - ANC_KN_POS_OBSERVED)
    )
    
    # (4.5) SAVE 
    
    assign(paste0("lso_2006_", arm),
           true_dynamics)
}
