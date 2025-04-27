################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Estimate % diagnosed among PLHIV in non-DREAMS districts in 2016Q1 and % non-disclosure
# @Author - Katherine Jia
################################################################################
# GOAL: Fit to Pelletier et al data in non-intervention districts to estimate % diagnosed among PLHIV in non-DREAMS districts in 2016Q1 and  % non-disclosure

# <I> INDICES ------------------------------------------------------------------
# NOTE: See Table 4.4 for descriptions on hypotheses

all_hypotheses = c("hypothesis1", 
                   "hypothesis2", 
                   "hypothesis3",
                   "hypothesis4")

arms <- c("intervention", "control")

N_t <- 17                                                                       # Number of timesteps from 2016Q1 to 2020Q2

pop_size_DREAMS <- 56872 + 26463

pop_size_nonDREAMS <- 33720 + 16674 + 15251

# Subset the incidence from 2016Q1 to 2020Q2

assumed_incidence_2016_to_2020 <- assumed_incidence %>%
  filter(time >= 2016.25, time<2020.75)

# <II> FIT THE MODEL FOR NON-INTERVENTION GROUP --------------------------------
# (2.1) Import data from Pelletier et al 

pelletier <- read_csv("1_data/pelletier.csv", col_names = FALSE)

# (2.2) Add column names and tidy the df 

colnames(pelletier) <- c("time", "NewDX_DREAMS", "Denom_DREAMS", "NewDX_nonDREAMS", "Denom_nonDREAMS")

pelletier_long <- pelletier %>%
  pivot_longer(cols = -time,
               names_to = c(".value", "group"), 
               names_sep = "\\_") %>%
  mutate(NewDX_rate=NewDX/Denom,
         arm = case_when(group=="DREAMS"~"1. Intervention districts",
                         group=="nonDREAMS"~"2. Non-intervention districts"))

# (2.3) Create a function to run model (% non-disclosure and % initial DX in 2016Q1 are variables)

run_model_non_intervention <- function(prop_non_disc_DX, 
                                       prop_initial_DX){
  
    for (n in 1:length(arms)) {
      
      lambda <- assumed_incidence_2016_to_2020$incidence_AGYW_nonDREAMS / 400
      
      # Set parameters for non-intervention districts (Table C.2 & C.3)
      
      non_disc_DX <- rep(prop_non_disc_DX, N_t + 1)
      delta <- rep(0.1914059, N_t + 1)
      tau_pos <- rep(0.0143811 * 1.55, N_t + 1)
      tau <- rep(0.0143811, N_t + 1)
      mu_pos <- 0.025
      
      # Create empty model compartments
      
      S <-
        UnDX <-
        DX <-
        ANC_TOTAL <-
        ANC_POS <-
        ANC_NEW_POS_TRUE <-
        ANC_NEW_POS_OBSERVED <- ANC_KN_POS_OBSERVED <- rep(0, N_t + 1)
      
      # Initial values (conditional probabilities)
      # Get HIV prevalence in 2006 
      
      initial_pos <- 1 - lso_2006_control[41,"S"] / (lso_2006_control[41,"S"] + lso_2006_control[41,"UnDX"] + lso_2006_control[41,"DX"])
      initial_DX_given_pos <- prop_initial_DX
      
      # Initial values (joint probabilities) 
      
      initial_UnDX <- initial_pos * (1 - initial_DX_given_pos)
      initial_DX <- initial_pos * initial_DX_given_pos
      
      S[1] <- pop_size_nonDREAMS * (1 - initial_pos)
      UnDX[1] <- pop_size_nonDREAMS * initial_UnDX
      DX[1] <- pop_size_nonDREAMS * initial_DX
      ANC_TOTAL[1] <- tau[1] * S[1] + tau_pos[1] * (UnDX[1] + DX[1])
      ANC_POS[1] <- tau_pos[1] * (UnDX[1] + DX[1])
      ANC_NEW_POS_OBSERVED[1] <- tau_pos[1] * (UnDX[1] + non_disc_DX[1] * DX[1])
      ANC_KN_POS_OBSERVED[1] <- tau_pos[1] * (1 - non_disc_DX[1]) * DX[1]
      ANC_NEW_POS_TRUE[1] <- tau_pos[1] * UnDX[1]
      
      # Solve the model
      
      for (i in 1:N_t) {
        
        S[i + 1] <- S[i] - lambda[i] * S[i] + mu_pos * (UnDX[i] + DX[i])
        
        UnDX[i + 1] <- UnDX[i] + lambda[i] * S[i] - (tau_pos[i] + delta[i]) * UnDX[i] - mu_pos * UnDX[i]
        
        DX[i + 1] <- DX[i] + (tau_pos[i] + delta[i]) * UnDX[i] - mu_pos * DX[i]
        
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
    }
    
    return(dynamics)
    
  }
  
  # (2.4) Create a function for MLE 

  fn_non_intervention <- function(par){
    
    logit_prop_initial_DX <- par[1]
    
    prop_initial_DX <- expit(logit_prop_initial_DX)
    
    logit_prop_non_disc_DX <- par[2]
    
    prop_non_disc_DX <- expit(logit_prop_non_disc_DX)
    
    dynamics <- run_model_non_intervention(prop_initial_DX=prop_initial_DX,
                                                prop_non_disc_DX=prop_non_disc_DX)
    
    nll = - sum(dbinom(as.numeric(unlist(pelletier_long[pelletier_long$group=="nonDREAMS", "NewDX"])),
                       as.numeric(unlist(pelletier_long[pelletier_long$group=="nonDREAMS", "Denom"])),
                       dynamics$rate_self_report_new_diag,
                       log = T))
    
  }
  
# (2.5) Optim 
  
  non_intervention <- optim(
    c(0, 0),
    fn_non_intervention
  )
  
# (2.6) Results 
  
  # proportion diagnosed among HIV positive AGYW in 2016Q1 in non-intervention districts
  
  expit(non_intervention$par[1])
  
  # proportion of non-disclosure of status in non-intervention districts
  
  expit(non_intervention$par[2])
