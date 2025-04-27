################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Generate Fig C.1 (Sec 4.6)
# @Author - Katherine Jia
################################################################################

# <I> Plot the true versus inferred population dynamics ------------------------

for(j in 7:10) {
  
  scenario <- all_scenarios[j]
  
  true_dynamics <- get(paste0("true_dynamics_", scenario))
  
  # (4.1) Plot the true distribution
  
  true_df_long <- true_dynamics %>%
    dplyr::select(step,
           `2. Total HIV positive women (undiagnosed)` = UnDX, 
           `1. Total HIV positive women (diagnosed)` = DX) %>%
    pivot_longer(-step, names_to = "Comp", values_to = "Freq") 
  
  plot_true_pop <- ggplot(true_df_long, 
                              aes(x = step, y = Freq, fill=Comp)) +
    geom_bar(position = "stack", stat = "identity") + 
    labs(y = "number of women",
         x = " ",
         fill = " ") +
    theme_minimal() +
    scale_fill_grey(start=0.8, end=0.2) +
    theme(text = element_text(size = 12))
  
  assign(paste0("plot_true_pop_", scenario),
         plot_true_pop)
  
  # (4.2) Plot the inferred distribution 
  
  for (k in 1:2) {
    method <- all_methods[k]

    inference_df <- get(paste0("inference_", scenario, "_", method))
    
    inference_df_long <- inference_df %>%
      dplyr::select(step,
             `2. Total HIV positive women (undiagnosed)` = UnDX, 
             `1. Total HIV positive women (diagnosed)` = DX) %>%
      pivot_longer(-step, names_to = "Comp", values_to = "Freq") 
      
    plot_inferred_pop <- ggplot(inference_df_long, 
                                aes(x = step, y = Freq, fill=Comp)) +
      geom_bar(position = "stack", stat = "identity") + 
      labs(y = "number of women",
           x = " ",
           fill = " ") +
      theme_minimal() +
      scale_fill_grey(start=0.8, end=0.2) +
      theme(text = element_text(size = 12))
    
    assign(paste0("plot_inferred_pop_", scenario, "_", method),
           plot_inferred_pop)
  }
}

# <II> ARRANGE PLOTS -----------------------------------------------------------

col1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="True population dynamics") + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Population dynamics inferred from \n ANC HIV prevalence") + theme_void() 

col3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Population dynamics inferred from \n observed new diagnoses rate") + theme_void() 

row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 1", angle = 90) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 2", angle = 90) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 3", angle = 90) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="Scenario 4", angle = 90) + theme_void() 

layoutplot_Fig_C1 <- "
#AAAAAABBBBBBCCCCCC
#AAAAAABBBBBBCCCCCC
aggggggeeeeeellllll
aggggggeeeeeellllll
aggggggeeeeeellllll
biiiiiiffffffhhhhhh
biiiiiiffffffhhhhhh
biiiiiiffffffhhhhhh
cnnnnnnjjjjjjkkkkkk
cnnnnnnjjjjjjkkkkkk
cnnnnnnjjjjjjkkkkkk
dmmmmmmxxxxxxyyyyyy
dmmmmmmxxxxxxyyyyyy
dmmmmmmxxxxxxyyyyyy
######zzzzzzz######
######zzzzzzz######
"

plotlist_Fig_C1 <-
  list(
    A = col1,
    B = col2,
    C = col3,
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    g = plot_true_pop_scenario1,
    i = plot_true_pop_scenario2,
    n = plot_true_pop_scenario3,
    m = plot_true_pop_scenario4,
    e = plot_inferred_pop_scenario1_method1,
    l = plot_inferred_pop_scenario1_method2,
    f = plot_inferred_pop_scenario2_method1,
    h = plot_inferred_pop_scenario2_method2,
    j = plot_inferred_pop_scenario3_method1,
    k = plot_inferred_pop_scenario3_method2,
    x = plot_inferred_pop_scenario4_method1,
    y = plot_inferred_pop_scenario4_method2,
    z = guide_area())

Fig_C1 <- wrap_plots(plotlist_Fig_C1, guides = 'collect', design = layoutplot_Fig_C1) 

ggsave("3_figures/FigC_1.png", Fig_C1, width = 10, height = 8, dpi=300, units="in")

# <III> REPORT QUANTITATIVELY CHANGES IN POPULATION DYNAMICS -------------------

# Scenario 3: actual
true_dynamics_scenario3 %>%
  filter(step %in% c(0,40)) %>%
  mutate(prop_DX_among_pos = DX/(UnDX+DX)*100,
         prop_UnDX_among_pos = UnDX/(UnDX+DX)*100)
# Scenario 3: inferred using Method 2
inference_scenario3_method2 %>%
  filter(step %in% c(0,40)) %>%
  mutate(prop_DX_among_pos = DX/(UnDX+DX)*100,
         prop_UnDX_among_pos = UnDX/(UnDX+DX)*100)
# Scenario 4
true_dynamics_scenario4 %>%
  filter(step %in% c(0,40)) %>%
  mutate(pct_UnDX = UnDX/(S+UnDX+DX)*100)
