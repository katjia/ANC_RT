################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Predict new DX rate and ANC prevalence using fitted values
# @Author - Katherine Jia
################################################################################

# <I> INDICES ------------------------------------------------------------------
all_hypotheses = c("hypothesis1", 
                   "hypothesis2", 
                   "hypothesis3",
                   "hypothesis4")

arms <- c("intervention", "control")

N_t <- 17                                                                       ## Number of timesteps

pop_size_DREAMS <- 56872 + 26463

pop_size_nonDREAMS <- 33720 + 16674 + 15251

# <II> Run model using fitted parameter values ---------------------------------

for (j in 1:length(all_hypotheses)) {
  
  hypothesis <- all_hypotheses[j]
  
  for (n in 1:length(arms)) {
    
    arm <- arms[n]
    
    # (2.1) Set parameters (Table C2 & C3)
    
    if (arm == "intervention") {
      lambda <- assumed_incidence_2016_to_2020$incidence_AGYW_DREAMS/400
      pop_size <- pop_size_DREAMS
    }
    
    if (arm == "control") {
      lambda <- assumed_incidence_2016_to_2020$incidence_AGYW_nonDREAMS/400
      pop_size <- pop_size_nonDREAMS
    }
    
    non_disc_DX <- rep(expit(non_intervention$par[2]), N_t + 1)
    delta <- rep(0.1914059, N_t + 1)
    tau_pos <- rep(0.0143811 * 1.55, N_t + 1)
    tau <- rep(0.0143811, N_t + 1)
    mu_pos <- 0.025
    
    # (2.2) Set parameters fitted under H1: differently changing incidence
    
    if (hypothesis == "hypothesis1") {
      if (arm == "intervention") {
        lambda <- approx(c(0, 4, N_t),
                         c(assumed_incidence_2016_to_2020$incidence_AGYW_DREAMS[1]/400, 
                           assumed_incidence_2016_to_2020$incidence_AGYW_DREAMS[5]/400, 
                           expit(hypothesis1_intervention$par)),
                         n = N_t + 1)$y
      }
    }
    
    # (2.3) Set parameters fitted under H2: differently changing FRR
    
    if (hypothesis == "hypothesis2") {
      if (arm == "intervention") {
        tau_pos <- approx(c(0, 4, N_t),
                          c(0.0143811 * 1.55, 
                            0.0143811 * 1.55, 
                            0.0143811 * hypothesis2_intervention$par), 
                          n = N_t + 1)$y
        tau <- rep(0.0143811, N_t + 1)
      }
    }
    
    # (2.4) Set parameters fitted under H3: differently changing non-disclosure
    
    if (hypothesis == "hypothesis3") {
      if (arm == "intervention") {
        non_disc_DX <- approx(c(0, 4, N_t),
                              c(expit(non_intervention$par[2]), 
                                expit(non_intervention$par[2]), 
                                expit(hypothesis3_intervention$par)),
                              n = N_t + 1)$y
      }
    }
    
    # (2.5) Set parameters fitted under H4: different;y changing testing outside ANC
    
    if (hypothesis == "hypothesis4") {
      if (arm == "intervention") {
        delta <- approx(c(0, 3, 4, N_t),
                        c(0.1914059, 
                          0.1914059, 
                          expit(hypothesis4_intervention$par), 
                          expit(hypothesis4_intervention$par)),
                        n = N_t + 1)$y
      }
    }
    
    # (2.6) Empty compartments for model compartments
    
    S <-
      UnDX <-
      DX <-
      ANC_TOTAL <-
      ANC_POS <-
      ANC_NEW_POS_TRUE <-
      ANC_NEW_POS_OBSERVED <- ANC_KN_POS_OBSERVED <- rep(0, N_t + 1)
    
    # (2.7) Initial conditions 
    # (2.7.1) Initial conditions (conditional probabilities)
    
    if (arm == "intervention") {
      initial_pos <- 1 - lso_2006_intervention[41,"S"] / (lso_2006_intervention[41,"S"] + lso_2006_intervention[41,"UnDX"] + lso_2006_intervention[41,"DX"])
      initial_DX_given_pos <- int_ini_prop_DX$par[1] 
    }
    if (arm == "control") {
      initial_pos <- 1 - lso_2006_control[41,"S"] / (lso_2006_control[41,"S"] + lso_2006_control[41,"UnDX"] + lso_2006_control[41,"DX"])
      initial_DX_given_pos <- expit(non_intervention$par[1])
    }
    
    # (2.7.2) Initial conditions (joint probabilities)
    
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
    
    # (2.8) Solve the model
    
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
    
    dynamics <- data.frame(
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
    
    # (2.9) Save
    assign(paste0("dynamics_", hypothesis, "_", arm),
           dynamics)
  }
}
