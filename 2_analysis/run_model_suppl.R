################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Fertility for HIV-undiagnosed is different from HIV-diagnosed women
# @Author - Katherine Jia
################################################################################
# <I> INDICES ------------------------------------------------------------------
# In the main analysis, fertility rate is assumed the same for HIV-undiagnosed and diagnosed women.
# Here, fertility rate is the same for HIV-undiagnosed and susceptible women, but different from HIV-diagnosed women.

all_scenarios = c("scenario1", 
                  "scenario2", 
                  "scenario3", 
                  "scenario4")

N_t <- 40                                          # number of quarters

pop_size <- 90000                                  # population size

# <II> PARAMETERS --------------------------------------------------------------
for(j in 1:length(all_scenarios)){

  # (2.1) Set parameter values -------------------------------------------------

  
  lambda <- rep(0.00125, N_t + 1)
  tau <- rep(0.04, N_t + 1)                        # fertility for susceptible and HIV-undiagnosed women
  tau_pos <- rep(0.04*0.9, N_t + 1)                # fertility for the HIV-diagnosed
  non_disc_DX <- rep(INITIAL_NON_DISC_DX, N_t + 1)
  delta <- rep(0.04, N_t + 1) 
  mu_pos <- 0.0225
  
  # (2.2) Vary parameter values by scenarios -----------------------------------
  
  scenario <- all_scenarios[j]
 
    if(scenario=="scenario1"){
      lambda <- approx(c(0, N_t-1, N_t),
                       c(0.00125, 0.00025, 0.00025),
                       n = N_t + 1)$y
    }
    
    if(scenario=="scenario2"){
      tau_pos <- approx(c(0, N_t),
                        c(0.04 * 0.9,
                          0.04 * 0.8),
                        n = N_t + 1)$y
    }
    
    if(scenario=="scenario3"){
      non_disc_DX <- approx(c(0, N_t),
                            c(INITIAL_NON_DISC_DX, FINAL_NON_DISC_DX),
                            n = N_t + 1)$y
    }
    
    if(scenario=="scenario4"){
      delta <- approx(c(0, N_t),
                      c(0.04, 0.08),
                      n = N_t + 1)$y
    }

# <III> Initial conditions -----------------------------------------------------
# (3.1) Create empty model compartments

S <- UnDX <- DX <- ANC_TOTAL <- ANC_POS <- ANC_NEW_POS_TRUE <- ANC_NEW_POS_OBSERVED <- ANC_KN_POS_OBSERVED <- rep(0, N_t + 1)

# (3.2) Initial values (Conditional probabilities)

initial_pos <- 1 - 85263.16/9e4
initial_DX_given_pos <- 3697.047/(1039.795 + 3697.047)

# (3.3) Initial conditions (Joint probabilities)

initial_UnDX <- initial_pos * (1-initial_DX_given_pos)
initial_DX <- initial_pos * initial_DX_given_pos 

S[1] <- pop_size * (1-initial_pos)
UnDX[1] <- pop_size * initial_UnDX
DX[1] <- pop_size * initial_DX
ANC_TOTAL[1] <- tau[1] * S[1] + tau[1] * UnDX[1] + tau_pos[1] * DX[1]
ANC_POS[1] <- tau[1] * UnDX[1] + tau_pos[1] * DX[1]
ANC_NEW_POS_OBSERVED[1] <- tau[1] * UnDX[1] + tau_pos[1] * non_disc_DX[1] * DX[1]
ANC_KN_POS_OBSERVED[1] <- tau_pos[1] * (1 - non_disc_DX[1]) * DX[1]
ANC_NEW_POS_TRUE[1] <- tau[1] * UnDX[1]

# <IV> Solve the model ---------------------------------------------------------

for(i in 1:N_t){

# (4.1) Process model (NOTE: Transition occurs during interval i and a new state was updated at the beginning of interval i+1.)
  
S[i+1] <- S[i] - lambda[i] * S[i] + mu_pos * (UnDX[i] + DX[i])

UnDX[i+1] <- UnDX[i] + lambda[i] * S[i] - (tau[i] + delta[i]) * UnDX[i] - mu_pos * UnDX[i]

DX[i+1] <- DX[i] + (tau[i] + delta[i]) * UnDX[i] - mu_pos * DX[i]

# (4.2) Observation model (NOTE: ANC data are observed at the beginning of interval i+1, before the next round of transitions occurs.)

ANC_TOTAL[i+1] <- tau[i+1] * S[i+1] + tau[i+1] * UnDX[i+1] + tau_pos[i+1] * DX[i+1]

ANC_POS[i+1] <- tau[i+1] * UnDX[i+1] + tau_pos[i+1] * DX[i+1]

ANC_NEW_POS_OBSERVED[i+1] <- tau[i+1] * UnDX[i+1] + tau_pos[i+1] * non_disc_DX[i+1] * DX[i+1]

ANC_KN_POS_OBSERVED[i+1] <- tau_pos[i+1] * (1 - non_disc_DX[i+1]) * DX[i+1]

ANC_NEW_POS_TRUE[i+1] <- tau[i+1] * UnDX[i+1]

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
  ANC_prevalence = ANC_POS/ANC_TOTAL,
  rate_self_report_new_diag = ANC_NEW_POS_OBSERVED/(ANC_TOTAL - ANC_KN_POS_OBSERVED)
)

# <V> SAVE ---------------------------------------------------------------------

  assign(paste0("true_dynamics_", scenario, "_suppl"), true_dynamics)

}

