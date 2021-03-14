################################################################################
## Article Evolving changes in mortality of 13,301 critically ill adult 
##   patients with COVID-19 over eight months
##
## Survival Curves (Supplementary material)
## Leonardo S.L. Bastos (lslbastos), Pedro Kurtz 
##
################################################################################

# Obtaining main data frame of COVID-19 patients (with preparations)
source("code/Data_Analysis/descriptive_stats_covid.R")


# Library -----------------------------------------------------------------
library(tidyverse)
library(tidylog)
library(gtsummary)
library(survival)
library(survminer)
library(ggfortify)





# Defining data -----------------------------------------------------------
df_admissions_adj <- 
    df_covid_admissions %>% 
    mutate(
        first_resp_support = case_when(
            VentSupport %in% c("only_niv", "niv_to_mv") ~ "niv_first",
            VentSupport %in% c("only_mv") ~ "imv_first",
            TRUE ~ "oxygen"
        ),
        first_resp_support = factor(first_resp_support, 
                                    levels = c("oxygen", "niv_first", "imv_first"))
    )



# Survival Curves - All admissions ----------------------------------------
# General survival
period_surv <- survfit(
    Surv(hosp_los_follow_up, outcome_follow_up) ~ period,
    data = df_admissions_adj 
    )


period_surv_plot <- (ggsurvplot(period_surv, conf.int = TRUE,
                                 ggtheme = theme_bw(), 
                                 censor = FALSE,
                                 legend.title = "",
                                 risk.table = TRUE, pval.method = TRUE,
                                 legend.labs = c("Period 1", "Period 2", 
                                                 "Period 3", "Period 4"),
                                 xlab = "Time (days)",
                                 surv.scale = "percent",
                                 title = "(A)", pval = TRUE, pval.size = 3
                                 )
                      )

# Age
age_surv <- survfit(
    Surv(hosp_los_follow_up, outcome_follow_up) ~ idade_grupo ,
    data = df_admissions_adj)


age_surv_plot <- (ggsurvplot(age_surv, conf.int = TRUE,
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
                      guides(colour = guide_legend(nrow = 1))
                  ) 


# MFI
MFI_surv <- survfit(
    Surv(hosp_los_follow_up, outcome_follow_up) ~ MFI_level ,
    data = df_admissions_adj)


MFI_surv_plot <- (ggsurvplot(MFI_surv, conf.int = TRUE,
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
vent_surv <- survfit(
    Surv(hosp_los_follow_up, outcome_follow_up) ~ first_resp_support ,
    data = df_admissions_adj)


vent_surv_plot <- (ggsurvplot(vent_surv, conf.int = TRUE,
                            ggtheme = theme_bw(), censor = FALSE, 
                            pval = TRUE, pval.method = TRUE, pval.size = 3,
                            legend.title = "First Respiratory Support",
                            legend.labs = c("Oxygen",
                                            "NIRS first",
                                            # "NIV failure",
                                            "IMV first"),
                            xlab = "Time (days)",
                            surv.scale = "percent",
                            title = "(D)",
                            risk.table = TRUE
                            )
                 )


# Figure - all patients - surv curves

arrange_plots_fig1 <- 
arrange_ggsurvplots(list(period_surv_plot,  
                         MFI_surv_plot,
                         age_surv_plot,
                         vent_surv_plot), print = FALSE, nrow = 2, ncol = 2)


ggsave("output/supplementary/figures_supplementary/figre_KM_all_admissions_60day.png", 
       arrange_plots_fig1, 
       width = 13, height = 14.5, dpi = 800)






# Survival curve - strafitied by type of resp support ---------------------
df_admissions_adj_resp_support <- 
    df_admissions_adj %>% 
    filter(VentSupport != "none") %>% 
    droplevels() 



# Resp support groups - failure
resp_support_group_surv <- survfit(
    Surv(hosp_los_follow_up, outcome_follow_up) ~ VentSupport ,
    data = df_admissions_adj_resp_support)


resp_support_surv_plot <- (ggsurvplot(resp_support_group_surv, conf.int = TRUE,
                              ggtheme = theme_bw(), censor = FALSE, 
                              pval = TRUE, pval.method = TRUE, pval.size = 3,
                              legend.title = "Advanced Respiratory Support",
                              legend.labs = c("NIRS only",
                                              "NIV failure",
                                              "IMV only"),
                              xlab = "Time (days)",
                              surv.scale = "percent",
                              title = "",
                              risk.table = TRUE
                              )
                         )




ggsave("output/supplementary/figures_supplementary/figure_KM_resp_support_group_60day.png", 
       arrange_ggsurvplots(list(resp_support_surv_plot), 
                           print = FALSE, nrow = 1, ncol = 1), 
       width = 6, height = 8, dpi = 800)
