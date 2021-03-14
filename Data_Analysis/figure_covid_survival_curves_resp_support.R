################################################################################
## Article Evolving changes in mortality of 13,301 critically ill adult 
##   patients with COVID-19 over eight months
##
## Survival Curves (Main Analysis)
## Leonardo S.L. Bastos (lslbastos), Pedro Kurtz 
##
################################################################################


# Library -----------------------------------------------------------------
library(tidyverse)
library(tidylog)
library(gtsummary)
library(survival)
library(survminer)
library(ggfortify)


# Obtaining main data frame of COVID-19 patients (with preparations)
source("code/Data_Analysis/descriptive_stats_covid.R")



# Defining data -----------------------------------------------------------
df_admissions_adj_resp_support <- 
    df_covid_admissions %>%
    filter(VentSupport != "none") %>% 
    droplevels() %>%  
    mutate(
        first_resp_support = case_when(
            VentSupport %in% c("only_niv", "niv_to_mv") ~ "nirs_first",
            VentSupport %in% c("only_mv") ~ "imv_first"
        ),
        first_resp_support = factor(first_resp_support, 
                                    levels = c("nirs_first", "imv_first"))
        )






# Survival Curves - Patients under adv. resp. support (Figure 2) ------------

# Period of admission (estimated using breakpoints)
period_surv_vent <- survfit(
    Surv(hosp_los_follow_up, outcome_follow_up) ~ period,
    data = df_admissions_adj_resp_support
    )


period_surv_plot_vent <- (ggsurvplot(period_surv_vent, conf.int = TRUE,
                                ggtheme = theme_bw(), 
                                censor = FALSE,
                                legend.title = "",
                                risk.table = TRUE, pval.method = TRUE,
                                legend.labs = c("Period 1", "Period 2", 
                                                "Period 3", "Period 4"),
                                xlab = "Time (days)",
                                title = "(A)", pval = TRUE, pval.size = 3,
                                surv.scale = "percent"
                                )
                          )




# Age group
age_surv_vent <- survfit(
    Surv(hosp_los_follow_up, outcome_follow_up) ~ idade_grupo ,
    data = df_admissions_adj_resp_support
    )


age_surv_plot_vent <- (ggsurvplot(age_surv_vent, conf.int = TRUE,
                             ggtheme = theme_bw(), censor = FALSE, pval = TRUE,
                             risk.table = TRUE, pval.method = TRUE, pval.size = 3,
                             legend.title = "Age", 
                             legend.labs = c("< 40",
                                             "40-49",
                                             "50-59",
                                             "60-69", 
                                             "70-79",
                                             ">= 80"),
                             xlab = "Time (days)",
                             surv.scale = "percent",
                             title = "(B)") +
                      guides(colour = guide_legend(nrow = 1)
                             )
                      ) 


# MFI
MFI_surv_vent <- survfit(
    Surv(hosp_los_follow_up, outcome_follow_up) ~ MFI_level ,
    data = df_admissions_adj_resp_support
    )


MFI_surv_plot_vent <- (ggsurvplot(MFI_surv_vent, conf.int = TRUE,
                             ggtheme = theme_bw(), censor = FALSE,
                             legend.title = "MFI", pval = TRUE, pval.method = TRUE,
                             risk.table = TRUE, pval.size = 3,
                             # legend.position = "right",
                             legend.labs = c("Non-frail",
                                             "Pre-frail",
                                             "Frail"),
                             xlab = "Time (days)",
                             surv.scale = "percent",
                             title = "(C)"
                             )
                       )


# Ventilation groups
first_support_surv_vent <- survfit(
    Surv(hosp_los_follow_up, outcome_follow_up) ~ first_resp_support,
    data = df_admissions_adj_resp_support
    )


first_support_surv_plot_vent <- (ggsurvplot(first_support_surv_vent, conf.int = TRUE,
                              ggtheme = theme_bw(), censor = FALSE, 
                              pval = TRUE, pval.method = TRUE, pval.size = 3,
                              legend.title = "First Respiratory support",
                              legend.labs = c("NIRS first",
                                              # "NIV failure",
                                              "IMV first"),
                              xlab = "Time (days)",
                              surv.scale = "percent",
                              title = "(D)",
                              risk.table = TRUE
                              )
                        )




# Combining survival plots (Figure 2)
arrange_plots_fig_vent <- 
    arrange_ggsurvplots(list(period_surv_plot_vent,  
                             MFI_surv_plot_vent,
                             age_surv_plot_vent,
                             first_support_surv_plot_vent), 
                        print = FALSE, nrow = 2, ncol = 2)


ggsave("output/main/figures/figure_2_survival_curves_ventilated_60day.pdf", 
       arrange_plots_fig_vent, 
       width = 11, height = 13.5, dpi = 800)







# finished












