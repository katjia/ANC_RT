################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Source other scripts
# @Author - Katherine Jia
################################################################################

# <I> LOAD LIBRARIES -----------------------------------------------------------
library(odin)
library(TMB)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(scales)
library(readr)
library(tidyr)
library(dplyr)
library(readxl)
library(grid)
library(ggrepel)

# <II> HELPER FUNCTIONS --------------------------------------------------------

expit <- function(x){1/(1+exp(-x))}
logit <- function(x){log(x/(1-x))}
`%nin%` <- Negate(`%in%`)

# <III> SIMULATION ANALYSIS (Sec 4.4 to 4.6) -----------------------------------
# (3.1) Main results 
# In the main results, Scenario 3 considers the case when non-disc decreases from 50% to 30%.

INITIAL_NON_DISC_DX <- 0.5
FINAL_NON_DISC_DX <- 0.3

source("2_analysis/run_model.R")
source("2_analysis/plot_Fig4_2.R")
source("2_analysis/plot_Fig4_3_and_C_2.R")
source("2_analysis/plot_FigC_1.R")

# (3.2) Suppl results 
# In Appendix C3, Scenario 3 considers the case when non-disc is 0%.

INITIAL_NON_DISC_DX <- 0
FINAL_NON_DISC_DX <- 0
source("2_analysis/run_model.R")
source("2_analysis/plot_Fig4_3_and_C_2.R")

# In the main analysis, fertility rate is assumed the same for HIV-undiagnosed and diagnosed women.
# In Appendix C6, fertility rate is the same for HIV-undiagnosed and susceptible women, but different from HIV-diagnosed women.

INITIAL_NON_DISC_DX <- 0.5
FINAL_NON_DISC_DX <- 0.3
source("2_analysis/run_model_suppl.R")
source("2_analysis/plot_FigC_3.R")

# <IV> REAL DATA EXAMPLE (Sec 4.7) ---------------------------------------------
source("2_analysis/lso_initial_values.R")                                       # Specify initial prevalence for DREAMS and non-DREAMS districts
source("2_analysis/lso_fit_model_hypotheses_non_intervention.R")                # Estimate % diagnosed in 2016Q1 in non-DREAMS districts and % non-disclose 
source("2_analysis/lso_estimate_initial_DX.R")                                  # Estimate % diagnosed in 2016Q1 in DREAMS districts
source("2_analysis/lso_fit_model_hypotheses_intervention.R")                    # Estimate parameters for DREAMS districts under each of the hypotheses
source("2_analysis/lso_run_model_hypotheses_fitted.R")                          # 
source("2_analysis/lso_Fig4_5.R")
source("2_analysis/lso_Fig4_4.R")

