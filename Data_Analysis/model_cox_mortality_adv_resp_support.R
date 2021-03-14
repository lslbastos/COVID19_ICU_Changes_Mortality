################################################################################
## Article Evolving changes in mortality of 13,301 critically ill adult 
##   patients with COVID-19 over eight months
##
## Random-effects Cox model (mortality model - main analysis)
## Leonardo S.L. Bastos (lslbastos), Pedro Kurtz 
##
################################################################################



# Library -----------------------------------------------------------------
library(tidyverse)
library(tidylog)
library(WeightIt)
library(coxme)
library(broom)
library(ehahelper)




# Obtaining main data frame of COVID-19 patients (with preparations)
source("code/Data_Analysis/descriptive_stats_covid.R")





################################################################################
# Data preparation --------------------------------------------------------
## Defining data for modeling (Cox Mortality)
## Sample for model development: patients that required advanced respiratory support
# N = 4,188

df_model <-
    df_covid_admissions %>%
    filter(
        VentSupport != "none"
    ) %>% 
    droplevels() %>% 
    mutate(
        first_resp_support = case_when(
            VentSupport %in% c("only_niv", "niv_to_mv") ~ "niv_first",
            VentSupport %in% c("only_mv") ~ "imv_first",
            ),
        first_resp_support = factor(first_resp_support, 
                                    levels = c("imv_first", "niv_first")),
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
# Propensity score estimation (final model in 'propensity_score_mode.R')
model_no_var <- 
    first_resp_support ~
    Gender +
    idade_grupo +
    MFI_level +
    SofaScore +
    period +
    emergency +
    cardio_disease +
    IsVasopressors


## Estimating Propensity scores and weights - "ATE" estimand
ps_values <- 
    weightit(
        model_no_var,
        family = "binomial",
        estimand = "ATE",
        data = df_model,
        method = "ps"
        )


## Obtaining propensity scores and IPT weights for each patients
df_cox_model_ps_ate <-
    df_model %>%
    bind_cols(
        ps_ate = ps_values$ps,
        ps_w_ate = ps_values$weights
    ) 


### Propensity Score model evaluation
## Results from the propensity score model fitting and choice are in 
## "code/2_Data_Analysis/data_model_propensity_main.R"










################################################################################
# Random-effects cox model (Frailty model) --------------------------------
## "Full model": demographics, selected comorbidities and complications
cox_model_full <- 
    Surv(hosp_los_follow_up, outcome_follow_up) ~
    Gender +
    idade_grupo +
    MFI_level +
    SofaScore +
    Saps3Q +
    emergency +
    hypertension +
    diabetes +
    obesity +
    copd_asthma +
    cardio_disease +
    ResourceIsVasopressors +
    ResourceIsRenalReplacementTherapy +
    first_resp_support +
    period +
    (1 | HospitalCode)


df_cox_model_full <-
    coxme(
        cox_model_full,
        weights = ps_w_ate,
        data = df_cox_model_ps_ate
    )


df_cox_model_full_results <-
    df_cox_model_full %>%  
    tidy(exponentiate = TRUE) %>%
    mutate(
        p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
        HR_ci = paste0(round(estimate, 2), 
                       " (", round(conf.low, 2), " - ", round(conf.high, 2), ")")
    ) %>%
    select(term, HR_ci, p_adj)







## Model w/ non-significant  comorbidities: demographics and complications  
cox_model_no_comorb <- 
    Surv(hosp_los_follow_up, outcome_follow_up) ~
    Gender +
    idade_grupo +
    MFI_level +
    SofaScore +
    Saps3Q +
    emergency +
    hypertension +
    # diabetes +
    # obesity +
    # copd_asthma +
    cardio_disease +
    ResourceIsVasopressors +
    ResourceIsRenalReplacementTherapy +
    first_resp_support +
    period +
    (1 | HospitalCode)


df_cox_model_no_comorb <-
    coxme(
        cox_model_no_comorb,
        weights = ps_w_ate,
        data = df_cox_model_ps_ate
    )


df_cox_model_no_comorb_results <- 
    df_cox_model_no_comorb %>% 
    tidy(exponentiate = TRUE) %>%
    mutate(
        p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
        HR_ci = paste0(round(estimate, 3), " (", 
                       round(conf.low, 2), " - ", round(conf.high, 2), ")")
    ) %>%
    select(term, HR_ci, p_adj)











## Comparison of estiamted models (all admissions, no sign comorb, no sign var)
writexl::write_xlsx(
    left_join(
        df_cox_model_full_results,
        df_cox_model_no_comorb_results %>%
            rename(HR_ci_comorb = HR_ci, p_adj_comorb = p_adj)
        , by = c("term" = "term")
    )
    , "output/main/model/table_model_backward_comparison.xlsx")





## Comparison of models using estimated AIC and BIC
df_models_aic_bic <- 
    tibble(
        model = c("cox_model_full", 
                  "cox_model_no_comorb"),
        AIC = c(
            2 * (df_cox_model_full$loglik[2] - df_cox_model_full$loglik[1]) - 2 * df_cox_model_full$df[1],
            2 * (df_cox_model_no_comorb$loglik[2] - df_cox_model_no_comorb$loglik[1]) - 2 * df_cox_model_no_comorb$df[1]
            ),
        BIC = c(
            2 * (df_cox_model_full$loglik[2] - df_cox_model_full$loglik[1]) - log(df_cox_model_full$n[1]) * df_cox_model_full$df[1],
            2 * (df_cox_model_no_comorb$loglik[2] - df_cox_model_no_comorb$loglik[1]) - log(df_cox_model_no_comorb$n[1]) * df_cox_model_no_comorb$df[1]
            )
    )



writexl::write_xlsx(
    df_models_aic_bic
    , "output/main/model/table_model_backward_aic_bic.xlsx"
    )



## The model with lowest AIC/BIC was 'cox_model_full' (all variables of interest)






# Finished








################################################################################
# Sensitivity Analysis: Propensity Scores - SMR-weighted -------



## Estimating Propensity scores and weights - "ATT" estimand or SMR-weighted
ps_values_att <- weightit(model_no_var,
                          family = "binomial",
                          estimand = "ATT",
                          data = df_model,
                          method = "ps")


df_cox_model_ps_att <-
    df_model %>%
    bind_cols(
        ps_att = ps_values_att$ps,
        ps_w_att = ps_values_att$weights
    ) 




## "Full model": demographics, selected comorbidities and complications 
df_cox_model_full_att <-
    coxme(
        cox_model_full,
        weights = ps_w_att,
        data = df_cox_model_ps_att
    )



df_cox_model_full_att_results <- 
    df_cox_model_full_att %>% 
    tidy(exponentiate = TRUE) %>%
    mutate(
        p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
        HR_ci = paste0(round(estimate, 3), " (", round(conf.low, 3), " - ", round(conf.high, 3), ")")
    ) %>%
    select(term, HR_ci, p_adj)


writexl::write_xlsx(
    df_cox_model_full_att_results
    , "output/main/model/table_cox_model_full_ATT.xlsx"
    )






# Sensitivity Analysis: Propensity Scores - Trim upper 5% --------

## Estimating Propensity scores - "ATT" estimand or SMR-weighted
df_cox_model_ps_ate_trim <-
    df_model %>%
    bind_cols(
        ps_ate = ps_values$ps,
        ps_w_ate = ps_values$weights
    ) %>% 
    filter(
        ps_ate <= quantile(ps_ate, probs = 0.95)
    )




## "Full model": demographics, selected comorbidities and complications - All admissions N = 13,301
df_cox_model_full_trim <-
    coxme(
        cox_model_full,
        weights = ps_w_ate,
        data = df_cox_model_ps_ate_trim
    )



df_cox_model_full_trim_results <- 
    df_cox_model_full_trim %>% 
    tidy(exponentiate = TRUE) %>%
    mutate(
        p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
        HR_ci = paste0(round(estimate, 3), " (", round(conf.low, 3), " - ", round(conf.high, 3), ")")
    ) %>%
    select(term, HR_ci, p_adj)




writexl::write_xlsx(
    df_cox_model_full_trim_results
    , "output/main/model/table_cox_model_full_ATE_trim.xlsx"
    )














# Analysis of quantitative variables for Cox model ------------------------

## Table Age groups x 60-day in-hospital mortality
df_age_mortality <- 
    df_model %>% 
    group_by(idade_grupo) %>% 
    summarise(
        deaths       = sum(outcome_follow_up),
        total        = n(),
        mortal       = sum(outcome_follow_up) / n(),
        mortal_ratio = paste0(sum(outcome_follow_up), "/", n())
    ) %>% 
    ungroup()


## Age x 60-day in-hospital mortality
plot_age_mortality <- 
    df_age_mortality %>% 
    ggplot() +
    geom_col(aes(x = idade_grupo, y = mortal)) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Age (years)", y = "60-day in-hospital mortality") +
    theme_bw()



## Table SOFA x 60-day in-hospital mortality
df_sofa_mortality <- 
    df_model %>% 
    group_by(SofaScore) %>% 
    summarise(
        deaths       = sum(outcome_follow_up),
        total        = n(),
        mortal       = sum(outcome_follow_up) / n(),
        mortal_ratio = paste0(sum(outcome_follow_up), "/", n())
    ) %>% 
    ungroup()



## SOFA x Mortality
plot_sofa_mortality <- 
    df_sofa_mortality %>% 
    ggplot() +
    geom_col(aes(x = SofaScore, y = mortal)) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "SOFA score", y = "60-day in-hospital mortality") +
    theme_bw()



## Table SAPS3 x in-hospital mortality
df_saps_mortality <- 
    df_model %>%
    mutate(Saps3decil = cut_number(Saps3Points, n = 4)) %>% 
    group_by(Saps3decil) %>% 
    summarise(
        deaths       = sum(outcome_follow_up),
        total        = n(),
        mortal       = sum(outcome_follow_up) / n(),
        mortal_ratio = paste0(sum(outcome_follow_up), "/", n())
    ) %>% 
    ungroup()


## SAPS3 x Mortality
plot_saps_mortality <- 
    df_saps_mortality %>% 
    ggplot() +
    geom_col(aes(x = Saps3decil, y = mortal)) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_x_discrete(labels = c("<= 42", "43 to 50", "51 to 61", "> 61")) +
    labs(x = "SAPS-3 (quartiles)", y = "60-day in-hospital mortality") +
    theme_bw()



ggsave("output/supplementary/figures_supplementary/figure_Age_SOFA_SAPS_60day_mortality.png", 
       ggpubr::ggarrange(
           plot_age_mortality,
           plot_saps_mortality,
           plot_sofa_mortality,
           nrow = 3, ncol = 1
       ),
       width = 4, height = 8, dpi = 800)









## Martingale's residuals for continuous variables

survminer::ggcoxfunctional(
    coxph(Surv(hosp_los_follow_up, outcome_follow_up) ~
              Age +
              SofaScore +
              Saps3Points,
          data = df_cox_model_ps_ate)
    )



# Finished