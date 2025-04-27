################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Fit to Pelletier et al data in intervention districts under four hypotheses
# @Author - Katherine Jia
################################################################################

# <I> Create a function to run the model ---------------------------------------

  fn_hypothesis_intervention <- function(hypothesis, par){
    
    # For H2 to H4, incidence in intervention districts is assumed to follow the same linear trend as in Fig C3
    
    lambda <- assumed_incidence_2016_to_2020$incidence_AGYW_DREAMS/400                    # convert from incidence per 100 person-year to per person-quarter
    
    # (1.1) Set parameters (See Tables C2 & C3) 
    
    non_disc_DX <- rep(expit(non_intervention$par[2]), N_t + 1)
    delta <- rep(0.1914059, N_t + 1)
    tau_pos <- rep(0.0143811 * 1.55, N_t + 1)
    tau <- rep(0.0143811, N_t + 1)
    mu_pos <- 0.025
    
    # (1.2) H1: differently changing incidence 
    
    if (hypothesis == "hypothesis1") {
      
      lambda_end <- par
      
      lambda <- approx(c(0, 4, N_t),
                       c(assumed_incidence_2016_to_2020$incidence_AGYW_DREAMS[1]/400, 
                         assumed_incidence_2016_to_2020$incidence_AGYW_DREAMS[5]/400, 
                         lambda_end),
                       n = N_t + 1)$y
      
    }
    
    # (1.3) H2: differently changing FRR 
    
    if (hypothesis == "hypothesis2") {
      
      FRR_end <- par
      
      tau_pos <- approx(c(0, 4, N_t),
                        c(0.0143811 * 1.55, 
                          0.0143811 * 1.55, 
                          0.0143811 * FRR_end), 
                        n = N_t + 1)$y
      tau <- rep(0.0143811, N_t + 1)
    }
    
    # (1.4) H3: differently changing non-disclosure 
    
    if (hypothesis == "hypothesis3") {
      
      non_disc_end <- par
      
      non_disc_DX <- approx(c(0, 4, N_t),
                            c(expit(non_intervention$par[2]), 
                              expit(non_intervention$par[2]), 
                              non_disc_end),
                            n = N_t + 1)$y
    }
    
    # (1.5) H4: differently changing testing outside ANC 
    
    if (hypothesis == "hypothesis4") {
      
      delta_end <- par
      
      delta <- approx(c(0, 3, 4, N_t),
                      c(0.1914059, 0.1914059, delta_end, delta_end),
                      n = N_t + 1)$y
    }
    
    # (1.6) Create empty rows for model compartments 
    
    S <-
      UnDX <-
      DX <-
      ANC_TOTAL <-
      ANC_POS <-
      ANC_NEW_POS_TRUE <-
      ANC_NEW_POS_OBSERVED <- ANC_KN_POS_OBSERVED <- rep(0, N_t + 1)
    
    # (1.7) Initial conditions 
    # (1.7.1) Initial conditions (conditional probabilities)
    
    initial_pos <- 1 - lso_2006_intervention[41,"S"] / (lso_2006_intervention[41,"S"] + lso_2006_intervention[41,"UnDX"] + lso_2006_intervention[41,"DX"])
    initial_DX_given_pos <- int_ini_prop_DX$par[1]                              # estimated from `lso_estimate_initial_DX.R`
    
    # (1.7.2) Initial conditions (joint probabilities)
    
    initial_UnDX <- initial_pos * (1 - initial_DX_given_pos)
    initial_DX <- initial_pos * initial_DX_given_pos
    
    S[1] <- pop_size_DREAMS * (1 - initial_pos)
    UnDX[1] <- pop_size_DREAMS * initial_UnDX
    DX[1] <- pop_size_DREAMS * initial_DX
    ANC_TOTAL[1] <- tau[1] * S[1] + tau_pos[1] * (UnDX[1] + DX[1])
    ANC_POS[1] <- tau_pos[1] * (UnDX[1] + DX[1])
    ANC_NEW_POS_OBSERVED[1] <- tau_pos[1] * (UnDX[1] + non_disc_DX[1] * DX[1])
    ANC_KN_POS_OBSERVED[1] <- tau_pos[1] * (1 - non_disc_DX[1]) * DX[1]
    ANC_NEW_POS_TRUE[1] <- tau_pos[1] * UnDX[1]
    
    # (1.8) Solve the model 
    
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
      ANC_KN_POS_OBSERVED = ANC_KN_POS_OBSERVED,
      ANC_NEW_POS_TRUE = ANC_NEW_POS_TRUE,
      ANC_prevalence = ANC_POS / ANC_TOTAL,
      rate_self_report_new_diag = ANC_NEW_POS_OBSERVED / (ANC_TOTAL - ANC_KN_POS_OBSERVED)
    )
    
    # (1.9) negative log likelihood 
    
    nll = - sum(dbinom(as.numeric(unlist(pelletier_long[pelletier_long$group=="DREAMS", "NewDX"])),
                       as.numeric(unlist(pelletier_long[pelletier_long$group=="DREAMS", "Denom"])),
                       dynamics$rate_self_report_new_diag,
                       log = T))
    
    return(nll)
  }
  
# <II> Optim -------------------------------------------------------------------
  # (2.1) Fit under hypothesis 1
  
  fn_hypothesis1_intervention <- function(logit_par){
    par <- expit(logit_par)
    fn_hypothesis_intervention(hypothesis = "hypothesis1", par)
  }
  
  hypothesis1_intervention <- optim(
    c(logit(0.0055/5)),
    fn_hypothesis1_intervention)
  
  # par: incidence at 2020Q2
  expit(hypothesis1_intervention$par)
  
  # (2.2) Fit under hypothesis 1
  
  fn_hypothesis2_intervention <- function(par){
    fn_hypothesis_intervention(hypothesis = "hypothesis2", par)
  }
  
  hypothesis2_intervention <- optim(
    c(1),
    fn_hypothesis2_intervention
  )
  
  hypothesis2_intervention$par
  
  # (2.3) Optim for Hypothesis 3 
  
  fn_hypothesis3_intervention <- function(logit_par){
    par <- expit(logit_par)
    fn_hypothesis_intervention(hypothesis = "hypothesis3", par)
  }
  
  hypothesis3_intervention <- optim(
    c(logit(0.05)),
    fn_hypothesis3_intervention
  )
  
  expit(hypothesis3_intervention$par)
  
  # (2.4) Optim for Hypothesis 4 
  
  fn_hypothesis4_intervention <- function(logit_par){
    par <- expit(logit_par)
    fn_hypothesis_intervention(hypothesis = "hypothesis4", par)
  }
  
  hypothesis4_intervention <- optim(
    c(logit(0.05)),
    fn_hypothesis4_intervention
  )
  
  expit(hypothesis4_intervention$par)
  
# <III> NLL --------------------------------------------------------------------
  
  hypothesis1_intervention$value
  hypothesis2_intervention$value
  hypothesis3_intervention$value
  hypothesis4_intervention$value
