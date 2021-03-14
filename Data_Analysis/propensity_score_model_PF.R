################################################################################
## Article Evolving changes in mortality of 13,301 critically ill adult 
##   patients with COVID-19 over eight months
## 
##  Patients with P/F ratio values - Propensity score estimation
## Leonardo S.L. Bastos (lslbastos), Pedro Kurtz 
##
################################################################################


# Library -----------------------------------------------------------------
library(tidyverse)
library(tidylog)
library(broom)


# Obtaining main data frame of COVID-19 patients (with preparations)
source("code/Data_Analysis/descriptive_stats_covid.R")



# User functions ----------------------------------------------------------
## Calculating goodness-of-fit metrics
perf_metrics <- function(model) {
    auc <- ModelMetrics::auc(model)
    brier <- ModelMetrics::brier(model)

    return(list(auc = auc, brier = brier))
    }


## Obtaining goodness-of-fit metrics's confidence interval with bootstrap
model_perf_boot <- function(formula, list) {
    df_results <- 
        list %>% 
        map(~glm(formula, # RRT in the first 24h of admission,
                 family = "binomial",
                 data = .)
        ) %>% 
        map(~perf_metrics(.)) %>% 
        map_dfr(
            ~tibble(
                auc = .$auc,
                brier = .$brier
            )
        ) %>% 
        summarise(
            auc_low    = round(quantile(auc, probs = 1 - 0.975), 3),
            auc_high   = round(quantile(auc, probs = 0.975), 3),
            
            brier_low  = round(quantile(brier, probs = 1 - 0.975), 3),
            brier_high = round(quantile(brier, probs = 0.975), 3)
        ) %>% 
        mutate(
            auc_ci   = paste0("(", auc_low, " - ", auc_high, ")"),
            brier_ci = paste0("(", brier_low, " - ", brier_high, ")"),
        ) %>% 
        select(auc_ci, brier_ci)
    
    return(df_results)
    
    }
    







################################################################################
# Data preparation --------------------------------------------------------
## Defining data for modeling (PS and final model)
df_model_PS_PF <-
    df_covid_admissions %>%
    filter(
        VentSupport != "none",
        !is.na(PaO2FiO2)
    ) %>% 
    droplevels() %>% 
    mutate(
        first_resp_support = case_when(
            VentSupport %in% c("only_niv", "niv_to_mv") ~ "nirs_first",
            VentSupport %in% c("only_mv") ~ "imv_first",
        ),
        first_resp_support = factor(first_resp_support, 
                                    levels = c("imv_first", "nirs_first")),
        emergency = if_else(AdmissionSource == "Emergency", 1, 0),
        period = factor(period, levels = c(2, 1, 3, 4)),
        Saps3Q = cut_number(Saps3Points, n = 4) # SAPS3 quartiles
    ) %>% 
    select(
        outcome_follow_up, 
        hosp_los_follow_up, 
        first_resp_support,
        Age,
        Gender,
        Saps3Points,
        Saps3Q,
        SofaScore, 
        idade_grupo, 
        MFI_level, 
        MFIpoints, 
        VentSupport,
        emergency,
        hypertension,
        diabetes,
        cardio_disease,
        imunossupression,
        cerebro_disease,
        obesity,
        copd_asthma,
        malignancy,
        chronic_kidney,
        tobacco,
        liver_cirrhosis,
        PaO2FiO2,
        PaO2FiO21h_level,
        IsVasopressors,
        IsRenalReplacementTherapy,
        ResourceIsVasopressors,
        ResourceIsRenalReplacementTherapy,
        period,
        HospitalCode
    ) %>% 
    mutate_at(
        c("emergency",
          "hypertension",
          "diabetes",
          "cardio_disease",
          "imunossupression",
          "hypertension",
          "diabetes",
          "cardio_disease",
          "imunossupression",
          "cerebro_disease",
          "obesity",
          "copd_asthma",
          "malignancy",
          "chronic_kidney",
          "tobacco",
          "liver_cirrhosis",
          "IsVasopressors",
          "IsRenalReplacementTherapy",
          "ResourceIsVasopressors",
          "ResourceIsRenalReplacementTherapy"),
        function(x) { return(as.factor(if_else(x == 1, "yes", "no"))) }
    ) 






################################################################################
# Propensity score --------------------------------------------------------


## Model with all variables
model_full_PS_PF <- 
    first_resp_support ~
    idade_grupo +
    Gender +
    MFI_level +
    SofaScore +
    period +
    emergency +
    hypertension +
    diabetes +
    obesity +
    copd_asthma +
    cardio_disease +
    IsVasopressors + # Vasopressor in the first 24h of admission
    IsRenalReplacementTherapy +
    PaO2FiO21h_level

ps_model_full_PS_PF <- 
    glm(model_full_PS_PF,
        family = "binomial",
        data = df_model_PS_PF
        ) 


## No significant comorbidities
model_no_comorb_PS_PF <- 
    first_resp_support ~
    idade_grupo +
    Gender +
    MFI_level +
    SofaScore +
    period +
    emergency +
    # hypertension +
    # diabetes +
    obesity +
    # copd_asthma +
    # cardio_disease +
    IsVasopressors + 
    IsRenalReplacementTherapy +
    PaO2FiO21h_level
    


ps_model_no_comorb_PS_PF <- 
    glm(model_no_comorb_PS_PF, # RRT in the first 24h of admission,
        family = "binomial",
        data = df_model_PS_PF
        ) 



# Model without any significant variables
model_no_var_PS_PF <- 
    first_resp_support ~
    idade_grupo +
    Gender +
    MFI_level +
    SofaScore +
    period +
    emergency +
    # hypertension +
    # diabetes +
    obesity +
    # copd_asthma +
    # cardio_disease +
    IsVasopressors # Vasopressor in the first 24h of admission
    # IsRenalReplacementTherapy +
    # PaO2FiO21h_level


ps_model_no_var_PS_PF <- 
    glm(model_no_var_PS_PF,
        family = "binomial",
        data = df_model_PS_PF
    ) 


## Table with model comparison (sTable 6)
df_ps_models_comparison_PS_PF <- 
    left_join(
        ps_model_full_PS_PF %>% 
            tidy(exponentiate = TRUE,  conf.int = TRUE) %>%
            mutate(
                p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
                HR_ci = paste0(round(estimate, 2), 
                               " (", round(conf.low, 2), " - ", round(conf.high, 2), ")")
            ) %>%
            select(term, HR_ci, p_adj),
        ps_model_no_comorb_PS_PF %>% 
            tidy(exponentiate = TRUE,  conf.int = TRUE) %>%
            mutate(
                p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
                HR_ci = paste0(round(estimate, 2), 
                               " (", round(conf.low, 2), " - ", round(conf.high, 2), ")")
            ) %>%
            select(term, HR_ci_no_comorb = HR_ci, p_adj_no_comorb = p_adj)
    ) %>% 
    left_join(
        ps_model_no_var_PS_PF %>% 
            tidy(exponentiate = TRUE,  conf.int = TRUE) %>%
            mutate(
                p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
                HR_ci = paste0(round(estimate, 2), 
                               " (", round(conf.low, 2), " - ", round(conf.high, 2), ")")
            ) %>%
            select(term, HR_ci_no_var = HR_ci, p_adj_no_var = p_adj)
    )

writexl::write_xlsx(df_ps_models_comparison_PS_PF, 
                    "output/supplementary/propensity_scores_supplementary/table_PS_models_PF_ratio_comparison.xlsx")



## Estimating AUC and Brier with bootstrap

set.seed(2^31 - 1)
ls_model_boot_PS_PF <- 
    df_model_PS_PF %>% 
    infer::generate(reps = 1000, type = "bootstrap") %>% 
    split(.$replicate)
    

df_ps_model_full_PS_PF <- 
    model_perf_boot(formula = model_full_PS_PF, list = ls_model_boot_PS_PF)

df_ps_model_no_comorb_PS_PF <- 
    model_perf_boot(formula = model_no_comorb_PS_PF, list = ls_model_boot_PS_PF)

df_ps_model_no_var_PS_PF <- 
    model_perf_boot(formula = model_no_var_PS_PF, list = ls_model_boot_PS_PF)



## Comparison of PS models and metrics
df_ps_models_metric_PS_PF <- 
    tibble(
        model = c(
            "ps_model_full_PS_PF",
            "ps_model_no_comorb_PS_PF",
            "ps_model_no_var_PS_PF"
            ),
        AIC = c(
            AIC(ps_model_full_PS_PF),
            AIC(ps_model_no_comorb_PS_PF),
            AIC(ps_model_no_var_PS_PF)
        ),
        BIC = c(
            BIC(ps_model_full_PS_PF),
            BIC(ps_model_no_comorb_PS_PF),
            BIC(ps_model_no_var_PS_PF)
        ),
        AUROC = c(
            pROC::auc(pROC::roc(predictor = ps_model_full_PS_PF$fitted.values, 
                                response  = ps_model_full_PS_PF$y)),
            pROC::auc(pROC::roc(predictor = ps_model_no_comorb_PS_PF$fitted.values, 
                                response  = ps_model_no_comorb_PS_PF$y)),
            pROC::auc(pROC::roc(predictor = ps_model_no_var_PS_PF$fitted.values, 
                                response  = ps_model_no_var_PS_PF$y))
        ),
        AUROC_CI = c(
            df_ps_model_full_PS_PF$auc_ci,
            df_ps_model_no_comorb_PS_PF$auc_ci,
            df_ps_model_no_var_PS_PF$auc_ci
        ),
        Brier = c(
            ModelMetrics::brier(ps_model_full_PS_PF),
            ModelMetrics::brier(ps_model_no_comorb_PS_PF),
            ModelMetrics::brier(ps_model_no_var_PS_PF)
        ),
        Brier_CI = c(
            df_ps_model_full_PS_PF$brier_ci,
            df_ps_model_no_comorb_PS_PF$brier_ci,
            df_ps_model_no_var_PS_PF$brier_ci
        )
    ) %>% 
    mutate(
        AUC_CI   = paste(round(AUROC, 3), AUROC_CI),
        Brier_CI = paste(round(Brier, 3), Brier_CI)
    ) %>% 
    select(
        model, AIC, BIC, AUC_CI, Brier_CI
    )


## Table with comparison of metrics from models (sTable 6, additional)
writexl::write_xlsx(df_ps_models_metric_PS_PF, 
                    "output/supplementary/propensity_scores_supplementary/table_PS_models_PF_ratio_metric.xlsx")



## Patients with P/F ratio variables (N= 1,963)
## The best model is the 'reduced' model without any significant variables (lowest AIC and BIC)
# df_ps_models_comparison
# model_no_var



# Finished






