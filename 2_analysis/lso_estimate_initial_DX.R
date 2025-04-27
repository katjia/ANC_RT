################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Estimate % diagnosed among PLHIV in DREAMS districts in 2016
# @Author - Katherine Jia
################################################################################
# GOAL: Fit to the Pelletier et al data in the DREAMS districts to estimate `prop_DX` (% diagnosed among PLHIV) in 2016
# <I> Set parameters -----------------------------------------------------------

pos <- pop_size_DREAMS - lso_2006_intervention[41,"S"]                          # HIV+ AGYW at 2016 Note: row #41 pertains to 2016Q1
prop_pos <- pos/pop_size_DREAMS                                                 # Prevalence among AGYW at 2016
tau <- 0.0143811                                                                # Fertility rate of negative (Table C.3)
FRR <- 1.55                                                                     # Fertility rate ratio of pos to neg (Table C.2)
tau_pos <- tau * FRR
non_disc_fitted <- expit(non_intervention$par[2])                               # Non-disclosure proportion fitted from new DX rate in non-DREAMS districts

# <II> Create a function for MLE -----------------------------------------------
cali_rate_new_DX <- function(prop_DX){
  predicted_new_DX_rate_num <- tau_pos * pos * (1-prop_DX) + non_disc_fitted * tau_pos * pos * prop_DX
  predicted_new_DX_rate_denom <- non_disc_fitted * tau_pos * pos * prop_DX + tau_pos * pos * (1-prop_DX) + tau * (pop_size_DREAMS * (1-prop_pos))
  diff <- 0.11350575 - predicted_new_DX_rate_num/predicted_new_DX_rate_denom    # 0.11350575 is reported new DX rate in DREAMS districts in 2016
  abs(diff)
}

# <III> Optim ------------------------------------------------------------------
int_ini_prop_DX <- optim(c(0.2),
                         cali_rate_new_DX)

int_ini_prop_DX$par[1]
