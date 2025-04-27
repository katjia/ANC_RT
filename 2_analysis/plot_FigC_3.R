################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Generate Fig C.3 (Back-calculate HIV incidence when tau_UnDX = tau_S, which are different from tau_DX)
# @Author - Katherine Jia
################################################################################
# <I> INDICES ------------------------------------------------------------------
# method1: inferring incidence from ANC HIV prevalence
# method2: inferring incidence from new DX rate

all_methods = c("method1", "method2")

# <II> Calculate relative changes in incidence, ANC prev, new DX rate ----------

for(d in 1:4) {
  
  scenario <- all_scenarios[d]
  
  # (2.1) Get model outputs for each scenario
  
  data <- get(paste0("true_dynamics_", scenario, "_suppl"))
  
  # (2.2) Calculate relative changes
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
  
  # (2.3) Combine the relative changes with model outputs

  data_relative_change <- data %>%
    left_join(relative_change, by="step")
  
  assign(paste0("data_relative_change_", scenario, "_suppl"), data_relative_change)
  
}

# <III> SET PARAMETERS ----------------------------------------------------------
# These parameters were correctly specified by the model

tau_S_CORRECT <- tau

mu_POS_CONSTANT <- 0.0225

# In Scenarios 2 to 4, one of the ancillary factors changed over time, but the model incorrectly assumes that all ancillary factors were constant over time.

tau_pos_ASSUMED <- rep(0.04 * 0.9, N_t + 1)

non_disc_DX_ASSUMED <- rep(INITIAL_NON_DISC_DX, N_t + 1)

delta_ASSUMED <- rep(0.04, N_t + 1)

# <IV> BACK-CACULATION TO INTER INCIDENCE -------------------------------------

for(j in 1:4){
  
  scenario <- all_scenarios[j]
  
  data <- get(paste0("true_dynamics_", scenario, "_suppl"))
  
  prev_t <- data$ANC_prevalence[1:41] 
  
  rate_self_report_new_diag_t <- data$rate_self_report_new_diag[1:41] 
  
  for(k in 1:2){
    
    method <- all_methods[k]
    
    # (4.1) Set initial values (NOTE: the model has the correct initial values)
    
    S_MODELED <- UnDX_MODELED <- DX_MODELED <- ANC_TOTAL_MODELED <- ANC_POS_MODELED <- ANC_NEW_POS_OBSERVED_MODELED <- ANC_KN_POS_OBSERVED_MODELED <- rep(0, 41)
    lambda_MODELED <- rep(0, 40)
    
    S_MODELED[1] <- data$S[1] 
    UnDX_MODELED[1] <- data$UnDX[1]
    DX_MODELED[1] <- data$DX[1]
    
    ANC_TOTAL_MODELED[1] <- tau_S_CORRECT[1] * S_MODELED[1] + tau_S_CORRECT[1] * UnDX_MODELED[1] + tau_pos_ASSUMED[1] * DX_MODELED[1]
    ANC_POS_MODELED[1] <- tau_S_CORRECT[1] * UnDX_MODELED[1] + tau_pos_ASSUMED[1] * DX_MODELED[1]
    ANC_NEW_POS_OBSERVED_MODELED[1] <-  tau_S_CORRECT[1] * UnDX_MODELED[1] + tau_pos_ASSUMED[1] * non_disc_DX_ASSUMED[1] * DX_MODELED[1]
    ANC_KN_POS_OBSERVED_MODELED[1] <- tau_pos_ASSUMED[1] * (1 - non_disc_DX_ASSUMED[1]) * DX_MODELED[1]
    
    # (4.2) Back-calculate incidence using ANC HIV prevalence
    
    if(method=="method1"){
      
      for(i in 1:40){

        # Similar to equation C1 in Appendix C.1, but now UnDX has fertility rate tau_S
        
        C <- tau_S_CORRECT[i+1] * (S_MODELED[i] + mu_POS_CONSTANT * (UnDX_MODELED[i] + DX_MODELED[i]))
        y <- tau_S_CORRECT[i+1] * S_MODELED[i]
        A <- tau_pos_ASSUMED[i+1] * ((1 - mu_POS_CONSTANT) * DX_MODELED[i] + (tau_S_CORRECT[i] + delta_ASSUMED[i]) * UnDX_MODELED[i])
        B <- tau_S_CORRECT[i+1] * (1 - mu_POS_CONSTANT - tau_S_CORRECT[i] - delta_ASSUMED[i]) * UnDX_MODELED[i]
        one_minus_prev <- (1 - prev_t[i+1])
        
        # Inferred incidence at Day i

        lambda_MODELED[i] <- (C - one_minus_prev*(A+B+C))/y
        
        # Predicted states at the beginning of Day i+1 

        S_MODELED[i+1] <- S_MODELED[i] - lambda_MODELED[i] * S_MODELED[i] + mu_POS_CONSTANT * (UnDX_MODELED[i] + DX_MODELED[i])
        
        UnDX_MODELED[i+1] <- UnDX_MODELED[i] + lambda_MODELED[i] * S_MODELED[i] - (delta_ASSUMED[i] + tau_S_CORRECT[i]) * UnDX_MODELED[i] - mu_POS_CONSTANT * UnDX_MODELED[i]
        
        DX_MODELED[i+1] <- DX_MODELED[i] + (delta_ASSUMED[i] + tau_S_CORRECT[i]) * UnDX_MODELED[i] - mu_POS_CONSTANT * DX_MODELED[i]
        
        ANC_TOTAL_MODELED[i+1] <- tau_S_CORRECT[i+1] * S_MODELED[i+1] + tau_S_CORRECT[i+1] * UnDX_MODELED[i+1] + tau_pos_ASSUMED[i+1] * DX_MODELED[i+1]
        
        ANC_POS_MODELED[i+1] <- tau_S_CORRECT[i+1] * UnDX_MODELED[i+1] + tau_pos_ASSUMED[i+1] * DX_MODELED[i+1]
        
        ANC_NEW_POS_OBSERVED_MODELED[i+1] <- tau_S_CORRECT[i+1] * UnDX_MODELED[i+1] + tau_pos_ASSUMED[i+1] * non_disc_DX_ASSUMED[i+1] * DX_MODELED[i+1]
        
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
    
    # (4.3) Back-calculate incidence using new DX rate
    
    if(method=="method2"){
      
      for(i in 1:40){
        
        # Similar to equation C2 in Appendix C.1, but now UnDX has fertility rate tau_S
        
        Q <- tau_S_CORRECT[i + 1] * (1 - mu_POS_CONSTANT - tau_S_CORRECT[i] - delta_ASSUMED[i]) * UnDX_MODELED[i]
        r <- tau_S_CORRECT[i + 1] * S_MODELED[i]
        A <- tau_pos_ASSUMED[i+1] * ((1 - mu_POS_CONSTANT) * DX_MODELED[i] + (tau_S_CORRECT[i] + delta_ASSUMED[i]) * UnDX_MODELED[i])
        B <- tau_S_CORRECT[i+1] * (1 - mu_POS_CONSTANT - tau_S_CORRECT[i] - delta_ASSUMED[i]) * UnDX_MODELED[i]
        C <- tau_S_CORRECT[i+1] * (S_MODELED[i] + mu_POS_CONSTANT * (UnDX_MODELED[i] + DX_MODELED[i]))
        d <- rate_self_report_new_diag_t[i+1]
        
        # Inferred incidence at Day t
        
        lambda_MODELED[i] <- (Q - d * (non_disc_DX_ASSUMED[i+1] * A + B + C) + non_disc_DX_ASSUMED[i+1] * A) / (- r)
        
        # Predicted states at Day t+1
        
        S_MODELED[i+1] <- S_MODELED[i] - lambda_MODELED[i] * S_MODELED[i] + mu_POS_CONSTANT * (UnDX_MODELED[i] + DX_MODELED[i])
        
        UnDX_MODELED[i+1] <- UnDX_MODELED[i] + lambda_MODELED[i] * S_MODELED[i] - (tau_S_CORRECT[i] + delta_ASSUMED[i]) * UnDX_MODELED[i] - mu_POS_CONSTANT * UnDX_MODELED[i]
        
        DX_MODELED[i+1] <- DX_MODELED[i] + (tau_S_CORRECT[i] + delta_ASSUMED[i]) * UnDX_MODELED[i] - mu_POS_CONSTANT * DX_MODELED[i]
        
        ANC_TOTAL_MODELED[i+1] <- tau_S_CORRECT[i+1] * S_MODELED[i+1] + tau_S_CORRECT[i+1] * UnDX_MODELED[i+1] + tau_pos_ASSUMED[i+1] * DX_MODELED[i+1]
        
        ANC_POS_MODELED[i+1] <- tau_S_CORRECT[i+1] * UnDX_MODELED[i+1] + tau_pos_ASSUMED[i+1] * DX_MODELED[i+1]
        
        ANC_NEW_POS_OBSERVED_MODELED[i+1] <- tau_S_CORRECT[i+1] * UnDX_MODELED[i+1] + tau_pos_ASSUMED[i+1] * non_disc_DX_ASSUMED[i+1] * DX_MODELED[i+1]
        
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
    assign(paste0("inference_", scenario, "_", method, "_suppl"), df)
  }
}

# <V> Fig C.3A & C.3B PLOT CHANGES IN ANC MEASURES -----------------------------

for(j in 1:4) {

  scenario <- all_scenarios[j]
  
  data <- get(paste0("true_dynamics_", scenario, "_suppl"))
  
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

# <VI> Fig C.3C PLOT INFERRED INCIDENCE ----------------------------------------

for(j in 1:4) {
  
  # (6.1) Get incidence (true, method1, method2)
  scenario <- all_scenarios[j]
  
  true_dynamics <- get(paste0("true_dynamics_", scenario, "_suppl")) %>%
    mutate(method="1. Truth") %>%
    as.data.frame() %>% select(step, lambda, method) %>%
    filter(step<40)
  
  inference_df_method1 <- get(paste0("inference_", scenario, "_method1_suppl"))
  inference_df_method1 <- inference_df_method1 %>%
    mutate(method="2. Inferred from ANC prevalence") %>%
    as.data.frame() %>% select(step, lambda, method) %>%
    filter(step<40)
  
  inference_df_method2 <- get(paste0("inference_", scenario, "_method2_suppl"))
  inference_df_method2 <- inference_df_method2 %>%
    mutate(method="3. Inferred from \n reported new diagnoses rate") %>%
    as.data.frame() %>% select(step, lambda, method) %>%
    filter(step<40)
  
  # (6.2) Combine df
  lambda_df <- rbind(true_dynamics, inference_df_method1, inference_df_method2) 
  
  # (6.3) Plot
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

# <VII> ARRANGE PLOTS ----------------------------------------------------------------

row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 1: Time-varying incidence                     ", size = 5, hjust=1) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 2: Time-varying FRR                            ", size = 5, hjust=1) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 3: Time-varying non-disclosure                ", size = 5, hjust=1) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 4: Time-varying testing outside ANC            ", size = 5, hjust=1) + theme_void() 

layoutplot_figC3 <- "
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

plotlist_figC3 <-
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

figC3 <- wrap_plots(plotlist_figC3, guides = 'collect', design = layoutplot_figC3) 

ggsave("3_figures/FigC_3.png", figC3, width = 10, height = 11, dpi=300, units="in")
