################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Plot fitted parameter values under each hypothesis (Fig 4.4)
# @Author - Katherine Jia
################################################################################

# <I> INDICES ------------------------------------------------------------------

all_hypotheses = c("hypothesis1", 
                   "hypothesis2", 
                   "hypothesis3",
                   "hypothesis4")

arms <- c("intervention", "control")

hypothesis_labels = c("Differently changing incidence (H1)",
                      "Differently changing FRR (H2)",
                      "Differently changing non-disclosure (H3)",
                      "Differently changing testing outside ANC (H4)")

N_t <- 17                                                                       ## Number of timesteps

# <II> Plot parameters under each hypothesis -----------------------------------

for (j in 1:length(all_hypotheses)) {
  
  hypothesis <- all_hypotheses[j]
  
  for (n in 1:length(arms)) {
    
    arm <- arms[n]
    
    # (2.1) Set parameters (Table C2 & C3)
    
    if (arm == "intervention") {
      lambda <- assumed_incidence_2016_to_2020$incidence_AGYW_DREAMS/400
    }
    
    if (arm == "control") {
      lambda <- assumed_incidence_2016_to_2020$incidence_AGYW_nonDREAMS/400
    }
    
    non_disc_DX <- rep(expit(non_intervention$par[2]), N_t + 1)
    delta <- rep(0.1914059, N_t + 1)
    tau_pos <- rep(0.0143811 * 1.55, N_t + 1)
    tau <- rep(0.0143811, N_t + 1)
    mu_pos <- 0.025

    # (2.2) H1: differently changing incidence
    if (hypothesis == "hypothesis1") {
      if (arm == "intervention") {
        
        lambda <- approx(c(0, 4, N_t),
                         c(assumed_incidence_2016_to_2020$incidence_AGYW_DREAMS[1]/400, 
                           assumed_incidence_2016_to_2020$incidence_AGYW_DREAMS[5]/400, 
                           expit(hypothesis1_intervention$par)),
                         n = N_t + 1)$y
        
      }
    }
    
    # (2.3) H2: differently changing FRR
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
    
    # (2.4) H3: differently changing non-disclosure
    if (hypothesis == "hypothesis3") {
      if (arm == "intervention") {

        non_disc_DX <- approx(c(0, 4, N_t),
                              c(expit(non_intervention$par[2]), 
                                expit(non_intervention$par[2]), 
                                expit(hypothesis3_intervention$par)),
                              n = N_t + 1)$y
      }
    }
    
    # (2.5) H4: differently changing test frequency outside ANC
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
    
    # (2.6) Compile together the parameters
    parameters <- data.frame(
      arm = arm,
      step = 0:N_t,
      Incidence = lambda,
      FRR = tau_pos/tau,
      Non_disc = non_disc_DX,
      Testing = delta
    )
    
    # (2.7) Tidy the df
    parameters_long <- parameters %>%
      pivot_longer(-c(step,arm), values_to = "Freq", names_to = "Parameters") %>%
      mutate(par_labels = case_when(Parameters == "Incidence"~"1. Incidence",
                                    Parameters == "FRR"~"2. FRR",
                                    Parameters == "Non_disc"~"3. Non-disc",
                                    Parameters == "Testing"~"4. Testing \n outside ANC"),
             district_labels = case_when(arm=="intervention" ~ "1. Intervention",
                                         arm=="control" ~ "2. Non-intervention"))
    
    assign(paste0("parameters_", hypothesis, "_", arm),
           parameters_long)
  }
  
  # (2.8) combine intervention and control dfs
  hypotheses_par <-
    rbind(get(paste0("parameters_", hypothesis, "_intervention")),
          get(paste0("parameters_", hypothesis, "_control")))
  
  assign(paste0("parameters_", hypothesis),
         hypotheses_par)
}

# <III> Plot -------------------------------------------------------------------

all_parameters <- c("1. Incidence", "2. FRR", "3. Non-disc", "4. Testing \n outside ANC")

for(i in 1:length(all_hypotheses)){
  hypothesis <- all_hypotheses[i]
  for(j in 1:4){
    
    # Set limits
    parameter_j <- all_parameters[j]
    if(j==1){
      ymax=0.006
      ymin=0
    }
    if(j==2){
      ymax=1.6
      ymin=1
    }
    if(j==3){
      ymax=0.1
      ymin=0
    }
    if(j==4){
      ymax=0.3
      ymin=0.17
    }
    
    # Plot
    plot_hypothesis_par <- ggplot() +
      geom_line(data = get(paste0("parameters_", hypothesis)) %>% filter(par_labels == parameter_j), 
                aes(x=2016+step/4, y=Freq, lty=district_labels)) +
      lims(y=c(ymin, ymax),
           x=c(2016,2020)) +
      scale_x_continuous(breaks=c(2016, 2020)) +
      theme_minimal() +
      labs(title = parameter_j) +
      theme(axis.text.x = element_text(angle = 90),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title = element_text(size=12, hjust = 0.5),
            legend.title=element_blank(),
            legend.text=element_text(size=12)) 
    
    assign(paste0("plot_parameters_", hypothesis, "_", parameter_j),
           plot_hypothesis_par)
  }
}  

# <VI> Arrange the plot ----------------------------------------------

row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="H1: Differing declines in incidence        ", size=5) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="H2: Differing declines in FRR                ", size=5) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="H3: Differing declines in non-disclosure", size=5) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="H4: Differing declines in testing outside ANC", size=5) + theme_void() 

layoutplot_fig4_4 <- "
AAAAAAAAA###############
aaaaaaggggggeeeeeehhhhhh
aaaaaaggggggeeeeeehhhhhh
aaaaaaggggggeeeeeehhhhhh
BBBBBBBBB###############
bbbbbbiiiiiiffffffpppppp
bbbbbbiiiiiiffffffpppppp
bbbbbbiiiiiiffffffpppppp
CCCCCCCCC###############
ccccccnnnnnnjjjjjjqqqqqq
ccccccnnnnnnjjjjjjqqqqqq
ccccccnnnnnnjjjjjjqqqqqq
DDDDDDDDDDD#############
ddddddmmmmmmkkkkkkxxxxxx
ddddddmmmmmmkkkkkkxxxxxx
ddddddmmmmmmkkkkkkxxxxxx
######zzzzzzzzzzzz######
"

plotlist_fig4_4 <-
  list(
    A = row1,
    B = row2,
    C = row3,
    D = row4,
    a = `plot_parameters_hypothesis1_1. Incidence`,
    g = `plot_parameters_hypothesis1_2. FRR`,
    e = `plot_parameters_hypothesis1_3. Non-disc`,
    h = `plot_parameters_hypothesis1_4. Testing \n outside ANC`,
    b = `plot_parameters_hypothesis2_1. Incidence`,
    i = `plot_parameters_hypothesis2_2. FRR`,
    f = `plot_parameters_hypothesis1_3. Non-disc`,
    p = `plot_parameters_hypothesis2_4. Testing \n outside ANC`,
    c = `plot_parameters_hypothesis3_1. Incidence`,
    n = `plot_parameters_hypothesis3_2. FRR`,
    j = `plot_parameters_hypothesis3_3. Non-disc`,
    q = `plot_parameters_hypothesis3_4. Testing \n outside ANC`,
    d = `plot_parameters_hypothesis4_1. Incidence`,
    m = `plot_parameters_hypothesis4_2. FRR`,
    k = `plot_parameters_hypothesis4_3. Non-disc`,
    x = `plot_parameters_hypothesis4_4. Testing \n outside ANC`,
    z = guide_area())

fig4_4 <- wrap_plots(plotlist_fig4_4, guides = 'collect', design = layoutplot_fig4_4) 

ggsave("3_figures/Fig4_4.png", fig4_4, width = 10, height = 12, dpi=300, units="in")
  
