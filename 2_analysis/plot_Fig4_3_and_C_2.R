################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Generate Fig 4.3 (Sec 4.5 & 4.6)
# @Author - Katherine Jia
################################################################################
# <I> INDICES ------------------------------------------------------------------
# method1: inferring incidence from ANC HIV prevalence
# method2: inferring incidence from new DX rate

all_methods = c("method1", "method2")

# <II> SET PARAMETERS ----------------------------------------------------------
# Table 4.2: These parameters were correctly specified by the model

tau_S_CORRECT <- tau

mu_POS_CONSTANT <- 0.0225

# Table 4.2: In Scenarios 2 to 4, one of these ancillary factors changed over time, but the model incorrectly assumes that all ancillary factors were constant over time.

tau_pos_ASSUMED <- rep(0.04 * 0.9, N_t + 1)

non_disc_DX_ASSUMED <- rep(INITIAL_NON_DISC_DX, N_t + 1)

delta_ASSUMED <- rep(0.04, N_t + 1)

# <III> BACK-CACULATION TO INFER INCIDENCE -------------------------------------

for(j in 7:10){ # NOTE: In Sec 4.5 & 4.6, we considered Scenarios 1 to 4 (Table 4.2) only, which corresponds to `j in 7:10`
  
  scenario <- all_scenarios[j]
  
  # Get simulated data under each scenario
  
  data <- get(paste0("true_dynamics_", scenario))
  
  prev_t <- data$ANC_prevalence[1:41] 
  
  rate_self_report_new_diag_t <- data$rate_self_report_new_diag[1:41] 
  
  # Back-calculate incidence
  
  for(k in 1:2){
    
    method <- all_methods[k]
    
    # (3.1) Set initial values (NOTE: the model has the correct initial values)
    
    S_MODELED <- UnDX_MODELED <- DX_MODELED <- ANC_TOTAL_MODELED <- ANC_POS_MODELED <- ANC_NEW_POS_OBSERVED_MODELED <- ANC_KN_POS_OBSERVED_MODELED <- rep(0, 41)
    lambda_MODELED <- rep(0, 40)
    
    S_MODELED[1] <- true_dynamics$S[1] 
    UnDX_MODELED[1] <- true_dynamics$UnDX[1]
    DX_MODELED[1] <- true_dynamics$DX[1]
    
    ANC_TOTAL_MODELED[1] <- tau_S_CORRECT[1] * S_MODELED[1] + tau_pos_ASSUMED[1] * (UnDX_MODELED[1] + DX_MODELED[1])
    ANC_POS_MODELED[1] <- tau_pos_ASSUMED[1] * (UnDX_MODELED[1] + DX_MODELED[1])
    ANC_NEW_POS_OBSERVED_MODELED[1] <- tau_pos_ASSUMED[1] * (UnDX_MODELED[1] + non_disc_DX_ASSUMED[1] * DX_MODELED[1])
    ANC_KN_POS_OBSERVED_MODELED[1] <- tau_pos_ASSUMED[1] * (1 - non_disc_DX_ASSUMED[1]) * DX_MODELED[1]
    
    # (3.2) Back-calculate incidence using ANC HIV prevalence
    
    if(method=="method1"){
      
      for(i in 1:40){
        
        # See equation C1 in Appendix C.1
        
        C <- tau_S_CORRECT[i+1] * (S_MODELED[i] + mu_POS_CONSTANT * (UnDX_MODELED[i] + DX_MODELED[i]))
        y <- tau_S_CORRECT[i+1] * S_MODELED[i]
        A <- tau_pos_ASSUMED[i+1] * ((1 - mu_POS_CONSTANT) * DX_MODELED[i] + (tau_pos_ASSUMED[i] + delta_ASSUMED[i]) * UnDX_MODELED[i])
        B <- tau_pos_ASSUMED[i+1] * (1 - mu_POS_CONSTANT - tau_pos_ASSUMED[i] - delta_ASSUMED[i]) * UnDX_MODELED[i]
        x <- (tau_pos_ASSUMED[i+1] - tau_S_CORRECT[i+1]) * S_MODELED[i]
        one_minus_prev <- (1 - prev_t[i+1])
        
        # Inferred incidence at Day i
        
        lambda_MODELED[i] <- (C - one_minus_prev*(A+B+C))/(one_minus_prev*x+y)
        
        # Predicted states at the beginning of Day i+1 
        
        S_MODELED[i+1] <- S_MODELED[i] - lambda_MODELED[i] * S_MODELED[i] + mu_POS_CONSTANT * (UnDX_MODELED[i] + DX_MODELED[i])
        
        UnDX_MODELED[i+1] <- UnDX_MODELED[i] + lambda_MODELED[i] * S_MODELED[i] - (tau_pos_ASSUMED[i] + delta_ASSUMED[i]) * UnDX_MODELED[i] - mu_POS_CONSTANT * UnDX_MODELED[i]
        
        DX_MODELED[i+1] <- DX_MODELED[i] + (tau_pos_ASSUMED[i] + delta_ASSUMED[i]) * UnDX_MODELED[i] - mu_POS_CONSTANT * DX_MODELED[i]
        
        ANC_TOTAL_MODELED[i+1] <- tau_S_CORRECT[i+1] * S_MODELED[i+1] + tau_pos_ASSUMED[i+1] * (UnDX_MODELED[i+1] + DX_MODELED[i+1])
        
        ANC_POS_MODELED[i+1] <- tau_pos_ASSUMED[i+1] * (UnDX_MODELED[i+1] + DX_MODELED[i+1])
        
        ANC_NEW_POS_OBSERVED_MODELED[i+1] <- tau_pos_ASSUMED[i+1] * UnDX_MODELED[i+1] + tau_pos_ASSUMED[i+1] * non_disc_DX_ASSUMED[i+1] * DX_MODELED[i+1]
        
        ANC_KN_POS_OBSERVED_MODELED[i+1] <- tau_pos_ASSUMED[i+1] * (1 - non_disc_DX_ASSUMED[i+1]) * DX_MODELED[i+1]
        
      }
      
      df <- data.frame(step=0:40,
                       lambda = c(lambda_MODELED, NA),
                       S=S_MODELED,
                       UnDX=UnDX_MODELED,
                       DX=DX_MODELED,
                       ANC_TOTAL = ANC_TOTAL_MODELED,
                       ANC_POS = ANC_POS_MODELED,
                       ANC_NEW_POS_OBSERVED = ANC_NEW_POS_OBSERVED_MODELED,
                       ANC_prevalence = ANC_POS_MODELED/ANC_TOTAL_MODELED,
                       rate_self_report_new_diag = ANC_NEW_POS_OBSERVED_MODELED/(ANC_TOTAL_MODELED - ANC_KN_POS_OBSERVED_MODELED))
      }
    
    # (3.3) Back-calculate incidence using new DX rate
    
    if(method=="method2"){
      
      for(i in 1:40){
        
        # See equation C2 in Appendix C.1
        
        Q <- tau_pos_ASSUMED[i + 1] * (1 - mu_POS_CONSTANT - tau_pos_ASSUMED[i] - delta_ASSUMED[i]) * UnDX_MODELED[i]
        r <- tau_pos_ASSUMED[i + 1] * S_MODELED[i]
        A <- tau_pos_ASSUMED[i+1] * ((1 - mu_POS_CONSTANT) * DX_MODELED[i] + (tau_pos_ASSUMED[i] + delta_ASSUMED[i]) * UnDX_MODELED[i])
        B <- tau_pos_ASSUMED[i+1] * (1 - mu_POS_CONSTANT - tau_pos_ASSUMED[i] - delta_ASSUMED[i]) * UnDX_MODELED[i]
        C <- tau_S_CORRECT[i+1] * (S_MODELED[i] + mu_POS_CONSTANT * (UnDX_MODELED[i] + DX_MODELED[i]))
        x <- (tau_pos_ASSUMED[i+1] - tau_S_CORRECT[i+1]) * S_MODELED[i]
        d <- rate_self_report_new_diag_t[i+1]
        
        # Inferred incidence at Day i
        
        lambda_MODELED[i] <- (Q + non_disc_DX_ASSUMED[i+1] * A - d * (non_disc_DX_ASSUMED[i+1] * A + B + C)) / (d * x - r)
        
        # Predicted states at Day i+1 
        
        S_MODELED[i+1] <- S_MODELED[i] - lambda_MODELED[i] * S_MODELED[i] + mu_POS_CONSTANT * (UnDX_MODELED[i] + DX_MODELED[i])
        
        UnDX_MODELED[i+1] <- UnDX_MODELED[i] + lambda_MODELED[i] * S_MODELED[i] - (tau_pos_ASSUMED[i] + delta_ASSUMED[i]) * UnDX_MODELED[i] - mu_POS_CONSTANT * UnDX_MODELED[i]
        
        DX_MODELED[i+1] <- DX_MODELED[i] + (tau_pos_ASSUMED[i] + delta_ASSUMED[i]) * UnDX_MODELED[i] - mu_POS_CONSTANT * DX_MODELED[i]
        
        ANC_TOTAL_MODELED[i+1] <- tau_S_CORRECT[i+1] * S_MODELED[i+1] + tau_pos_ASSUMED[i+1] * (UnDX_MODELED[i+1] + DX_MODELED[i+1])
        
        ANC_POS_MODELED[i+1] <- tau_pos_ASSUMED[i+1] * (UnDX_MODELED[i+1] + DX_MODELED[i+1])
        
        ANC_NEW_POS_OBSERVED_MODELED[i+1] <- tau_pos_ASSUMED[i+1] * (UnDX_MODELED[i+1] + non_disc_DX_ASSUMED[i+1] * DX_MODELED[i+1])
        
        ANC_KN_POS_OBSERVED_MODELED[i+1] <- tau_pos_ASSUMED[i+1] * (1 - non_disc_DX_ASSUMED[i+1]) * DX_MODELED[i+1]

      }
      
      df <- data.frame(step=0:40,
                       lambda = c(lambda_MODELED, NA),
                       S=S_MODELED,
                       UnDX=UnDX_MODELED,
                       DX=DX_MODELED,
                       ANC_TOTAL = ANC_TOTAL_MODELED,
                       ANC_POS = ANC_POS_MODELED,
                       ANC_NEW_POS_OBSERVED = ANC_NEW_POS_OBSERVED_MODELED,
                       ANC_prevalence = ANC_POS_MODELED/ANC_TOTAL_MODELED,
                       rate_self_report_new_diag = ANC_NEW_POS_OBSERVED_MODELED/(ANC_TOTAL_MODELED - ANC_KN_POS_OBSERVED_MODELED))
      
    }
    assign(paste0("inference_", scenario, "_", method), df)
  }
}

# <VI> Fig 4.3A & 4.3B PLOT CHANGES IN ANC MEASURES ----------------------------

for(j in 7:10) {

  scenario <- all_scenarios[j]
  
  data <- get(paste0("true_dynamics_", scenario))
  
  plot_ANC_prevalence <- ggplot() +
      geom_line(data = data,
                aes(y = ANC_prevalence * 100, x = step/4)) +
      labs(y = "ANC HIV prevalence (%)",
           x = "year") +
      theme_minimal() +
      lims(y = c(0, 5)) +
      theme(text = element_text(size = 12))
    
    assign(paste0("plot_ANC_prevalence_", scenario),
           plot_ANC_prevalence)
  
  plot_rate_new_DX <- ggplot() +
      geom_line(data = data,
                aes(y = rate_self_report_new_diag * 100, x = step/4)) +
      labs(y = "reported new \n diagnoses rate (%)",
           x = "year") +
      theme_minimal() +
      lims(y = c(0, 5)) +
      theme(text = element_text(size = 12))
    
    assign(paste0("plot_rate_new_DX_", scenario),
           plot_rate_new_DX)
}

# <V> Fig 4.3C PLOT INFERRED INCIDENCE -----------------------------------------

for(j in 7:10) {
  
  # (5.1) Get incidence (true, method1, method2)
  scenario <- all_scenarios[j]
  
  true_dynamics <- get(paste0("true_dynamics_", scenario))
  true_dynamics <- true_dynamics %>%
    mutate(method="1. Truth") %>%
    as.data.frame() %>% select(step, lambda, method) %>%
    filter(step<40)
  
  inference_df_method1 <- get(paste0("inference_", scenario, "_method1"))
  inference_df_method1 <- inference_df_method1 %>%
    mutate(method="2. Inferred from ANC HIV prevalence") %>%
    as.data.frame() %>% 
    select(step, lambda, method) %>%
    filter(step<40)
  
  inference_df_method2 <- get(paste0("inference_", scenario, "_method2"))
  inference_df_method2 <- inference_df_method2 %>%
    mutate(method="3. Inferred from \n reported new diagnoses rate") %>%
    as.data.frame() %>% select(step, lambda, method) %>%
    filter(step<40)
  
  # (5.2) Combine df
  lambda_df <- rbind(true_dynamics, inference_df_method1, inference_df_method2) 
  
  # (5.3) Plot
  
  plot_inferred_inc <- ggplot(lambda_df) +
    geom_line(aes(x = step/4, y = lambda * 1000, group = method, lty = method)) + 
    labs(y = "incidence per \n 1,000 person-quarter",
         x = "year",
         lty="Incidence") +
    scale_linetype_manual("Incidence", values=c("solid", "dotdash", "longdash")) +
    theme_minimal() +
    lims(y = c(0, 1.3),
         x = c(0, 9.75)) +
    theme(text = element_text(size = 12))
  
  assign(paste0("plot_inferred_inc_", scenario),
         plot_inferred_inc)

}

# <VI> ARRANGE PLOTS ----------------------------------------------------------------

row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 1: Time-varying incidence                     ", size = 5, hjust=1) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 2: Time-varying FRR                            ", size = 5, hjust=1) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 3: Time-varying non-disclosure                ", size = 5, hjust=1) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 4: Time-varying testing outside ANC            ", size = 5, hjust=1) + theme_void() 

layoutplot_fig2 <- "
aaaaaaaaaaaaaaaaa##
ggggggeeeeeellllll#
ggggggeeeeeellllll#
ggggggeeeeeellllll#
bbbbbbbbbbbbbbbbb##
iiiiiiffffffhhhhhh#
iiiiiiffffffhhhhhh#
iiiiiiffffffhhhhhh#
cccccccccccccccccc#
mmmmmmjjjjjjkkkkkk#
mmmmmmjjjjjjkkkkkk#
mmmmmmjjjjjjkkkkkk#
ddddddddddddddddddd
nnnnnnxxxxxxyyyyyy#
nnnnnnxxxxxxyyyyyy#
nnnnnnxxxxxxyyyyyy#
############zzzzz##
############zzzzz##
"

plotlist_fig2 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    g = plot_ANC_prevalence_scenario1,
    e = plot_rate_new_DX_scenario1,
    l = plot_inferred_inc_scenario1,
    i = plot_ANC_prevalence_scenario2,
    f = plot_rate_new_DX_scenario2,
    h = plot_inferred_inc_scenario2,
    m = plot_ANC_prevalence_scenario3,
    j = plot_rate_new_DX_scenario3,
    k = plot_inferred_inc_scenario3,
    n = plot_ANC_prevalence_scenario4,
    x = plot_rate_new_DX_scenario4,
    y = plot_inferred_inc_scenario4,
    z = guide_area())

fig2 <- wrap_plots(plotlist_fig2, guides = 'collect', design = layoutplot_fig2) 

# For the main results in Fig 4.3, Scenario 3 considers the case when non-disc decreases from 50% to 30%.

if(INITIAL_NON_DISC_DX==0.5){
  ggsave("3_figures/Fig4_3.png", fig2, width = 10, height = 11, dpi=300, units="in")
}

# For the suppl results in Fig C2, Scenario 3 considers the case when non-disc is 0%.

if(INITIAL_NON_DISC_DX==0){
  ggsave("3_figures/FigC_2.png", fig2, width = 10, height = 11, dpi=300, units="in")
}

# <VII> REPORT RELATIVE AND ABSOLUTE CHANGES -----------------------------------
# Report the % change in ANC measures under each scenario

change_by_endpoints_scenarios1to4 <- NULL

for(j in 7:10) {
  
  scenario <- all_scenarios[j]
  
  data <- get(paste0("data_relative_change_", scenario))
  
  change_by_endpoints <- data %>%
  filter(step %in% c(0, 40)) %>%
  mutate(scenario = scenario,
         ANC_prevalence_percent = round(ANC_prevalence * 100,1),
         rel_change_prev = round(rel_change_prev),
         rate_self_report_new_diag_percent = round(rate_self_report_new_diag * 100,1),
         rel_change_rate_new_diag = round(rel_change_rate_new_diag)) %>% 
  dplyr::select(scenario, step, ANC_prevalence_percent, rate_self_report_new_diag_percent, rel_change_prev, rel_change_rate_new_diag)

  change_by_endpoints_scenarios1to4 <- rbind(change_by_endpoints_scenarios1to4, 
                                             change_by_endpoints)
}

change_by_endpoints_scenarios1to4 

# <VIII> REPORT ANY BIASES IN INCIDENCE ----------------------------------------

inf_inc_all_scenarios <- NULL

for(j in 7:10) {
  scenario <- all_scenarios[j]

  for (k in 1:2) {
    method <- all_methods[k]

    inference_df <- get(paste0("inference_", scenario, "_", method))

    inf_inc_df <- inference_df %>%
      filter(step %in% c(0, 39)) %>%
      mutate(inc_per_1000 = round(lambda * 1000, 2),
             bias = round((1 - inc_per_1000 / 1.25)*100),
             scenario = scenario,
             method = method) %>%
      dplyr::select(scenario,
             method,
             inc_per_1000,
             bias)

    inf_inc_all_scenarios <- rbind(inf_inc_all_scenarios, inf_inc_df)
  }
}

inf_inc_all_scenarios
