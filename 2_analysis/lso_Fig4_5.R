################################################################################
# @Project - Considerations for Using Reported New Diagnoses Rate at Antenatal Care Clinics to Infer Trends in HIV Incidence
# @Description - Generate Fig 4.5
# @Author - Katherine Jia
################################################################################

# <I> Plot the fit -------------------------------------------------------------

for (j in 1:length(all_hypotheses)) {
  
  hypothesis <- all_hypotheses[j]
  
# (1.1) new diagnoses rate
  
plot_fit_new_diag_rate <- ggplot() +
  geom_point(data = pelletier_long,
             aes(y = NewDX_rate * 100, x = time - 0.25, col = arm)) +
  geom_line(
    data = get(paste0("dynamics_",hypothesis,"_intervention")),
    aes(y = rate_self_report_new_diag * 100, 
        x = step / 4 + 2016,
        lty = "1. Intervention districts")) +
  geom_line(data = get(paste0("dynamics_",hypothesis,"_control")),
            aes(y = rate_self_report_new_diag * 100,
                x = step / 4 + 2016,
                lty = "2. Non-intervention districts")) +
  labs(y = "reported new \n diagnoses rate (%)",
       x = " ",
       lty = " ",
       col = " ") +
  scale_color_manual(
    values = c("black", "grey60"),
    labels = c('1. Intervention districts',
               '2. Non-intervention districts')) +
  lims(y = c(0, 15)) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

# (1.2) ANC HIV prevalence

plot_pred_ANC_prevalence <- ggplot() +
  geom_line(data = get(paste0("dynamics_",hypothesis,"_intervention")),
            aes(y = ANC_prevalence * 100, x = step / 4 + 2016, lty = "1. Intervention districts")) +
  geom_line(data = get(paste0("dynamics_",hypothesis,"_control")),
            aes(y = ANC_prevalence * 100, x = step / 4 + 2016, lty = "2. Non-intervention districts")) +
  labs(y = "ANC HIV prevalence (%)",
       x = " ",
       lty = " ",
       col = " ") +
  scale_color_manual(
    values = c("black", "grey60"),
    labels = c('1. Intervention districts',
               '2. Non-intervention districts')) +
  lims(y = c(0, 40)) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

assign(paste0("plot_fit_new_diag_rate_", hypothesis), plot_fit_new_diag_rate)
assign(paste0("plot_pred_ANC_prevalence_", hypothesis), plot_pred_ANC_prevalence)

}

# <II> ARRANGE PLOTS -----------------------------------------------------------

row1 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="H1: Differing declines in incidence                    ", size=5) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="H2: Differing declines in FRR                          ", size=5) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="H3: Differing declines in non-disclosure           ", size=5) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=0.1, y=0.1, label="H4: Differing declines in testing outside ANC  ", size=5) + theme_void() 

layoutplot_fig4_5 <- "
aaaaaa######
ggggggeeeeee
ggggggeeeeee
ggggggeeeeee
bbbbbb######
iiiiiiffffff
iiiiiiffffff
iiiiiiffffff
cccccc######
nnnnnnjjjjjj
nnnnnnjjjjjj
nnnnnnjjjjjj
dddddd######
mmmmmmkkkkkk
mmmmmmkkkkkk
mmmmmmkkkkkk
###zzzzzz###
"

plotlist_fig4_5 <-
  list(
    a = row1,
    b = row2,
    c = row3,
    d = row4,
    g = plot_fit_new_diag_rate_hypothesis1,
    i = plot_fit_new_diag_rate_hypothesis2,
    n = plot_fit_new_diag_rate_hypothesis3,
    m = plot_fit_new_diag_rate_hypothesis4,
    e = plot_pred_ANC_prevalence_hypothesis1,
    f = plot_pred_ANC_prevalence_hypothesis2,
    j = plot_pred_ANC_prevalence_hypothesis3,
    k = plot_pred_ANC_prevalence_hypothesis4,
    z = guide_area())

fig4_5 <- wrap_plots(plotlist_fig4_5, guides = 'collect', design = layoutplot_fig4_5) 

ggsave("3_figures/Fig4_5.png", fig4_5, width = 10, height = 15, dpi=300, units="in")
