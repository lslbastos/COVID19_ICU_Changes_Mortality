################################################################################
## Article Evolving changes in mortality of 13,301 critically ill adult 
##   patients with COVID-19 over eight months
##
## Random-effects Cox model - subset of patients with P/F and imputation
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
    






# Subset of Patients With PF ratio values ---------------------------------


df_model_PF <- 
    df_model %>% 
    filter(
        !is.na(PaO2FiO2)
    )





################################################################################
# Propensity score estimation 
# (final model in 'propensity_score_mode_PF.R')
model_no_var_PF <- 
    first_resp_support ~
    idade_grupo +
    Gender +
    MFI_level +
    SofaScore +
    period +
    emergency +
    obesity +
    IsVasopressors 

## Estimating Propensity scores and weights - "ATE" estimand
ps_values_PF <- 
    weightit(
        model_no_var_PF,
        family = "binomial",
        estimand = "ATE",
        data = df_model_PF,
        method = "ps"
        )


## Obtaining propensity scores and IPT weights for each patients
df_cox_model_ps_PF <-
    df_model_PF %>%
    bind_cols(
        ps_ate = ps_values_PF$ps,
        ps_w_ate = ps_values_PF$weights
    ) 










################################################################################
# "Full model" for patients that presented PaO2/FiO2 values - N = 1,963
df_cox_model_ps_PF <- 
    df_cox_model_ps_ate %>% 
    filter(!is.na(PaO2FiO2)) %>% 
    droplevels()



cox_model_full_PF <- 
    Surv(hosp_los_follow_up, outcome_follow_up) ~
    Gender +
    idade_grupo +
    MFI_level +
    SofaScore +
    Saps3Q +
    period +
    emergency +
    hypertension +
    diabetes +
    obesity +
    copd_asthma +
    cardio_disease +
    ResourceIsVasopressors +
    ResourceIsRenalReplacementTherapy +
    first_resp_support +
    PaO2FiO21h_level +
    (1 | HospitalCode)


coxme_full_PF <-
    coxme(
        cox_model_full_PF,
        weights = ps_w_ate,
        data = df_cox_model_ps_PF
    )



df_cox_model_full_PF_results <- 
    coxme_full_PF %>% 
    tidy(exponentiate = TRUE) %>%
    mutate(
        p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
        HR_ci = paste0(round(estimate, 3), " (", round(conf.low, 3), " - ", round(conf.high, 3), ")")
    ) %>%
    select(term, HR_ci, p_adj)


# writexl::write_xlsx(
#     df_cox_model_full_PF_results
#     , "output/supplementary/model_sensitivity/table_cox_model_full_PF.xlsx"
#     )



## Model w/ non-significant  comorbidities: demographics and complications  
cox_model_no_comorb_PF <- 
    Surv(hosp_los_follow_up, outcome_follow_up) ~
    Gender +
    idade_grupo +
    MFI_level +
    SofaScore +
    Saps3Q +
    period +
    emergency +
    hypertension +
    # diabetes +
    # obesity +
    copd_asthma +
    # cardio_disease +
    ResourceIsVasopressors +
    ResourceIsRenalReplacementTherapy +
    first_resp_support +
    PaO2FiO21h_level +
    (1 | HospitalCode)


coxme_no_comorb_PF <-
    coxme(
        cox_model_no_comorb_PF,
        weights = ps_w_ate,
        data = df_cox_model_ps_PF
    )


df_cox_model_no_comorb_PF_results <- 
    coxme_no_comorb_PF %>% 
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
        df_cox_model_full_PF_results,
        df_cox_model_no_comorb_PF_results %>%
            rename(HR_ci_comorb = HR_ci, p_adj_comorb = p_adj)
        , by = c("term" = "term")
    )
    , "output/supplementary/model_sensitivity/table_model_PF_backward_comparison.xlsx")





## Comparison of models using estimated AIC and BIC
df_models_aic_bic_PF <- 
    tibble(
        model = c("cox_model_full_PF", 
                  "cox_model_no_comorb_PF"),
        AIC = c(
            2 * (coxme_full_PF$loglik[2] - coxme_full_PF$loglik[1]) - 2 * coxme_full_PF$df[1],
            2 * (coxme_no_comorb_PF$loglik[2] - coxme_no_comorb_PF$loglik[1]) - 2 * coxme_no_comorb_PF$df[1]
        ),
        BIC = c(
            2 * (coxme_full_PF$loglik[2] - coxme_full_PF$loglik[1]) - log(coxme_full_PF$n[1]) * coxme_full_PF$df[1],
            2 * (coxme_no_comorb_PF$loglik[2] - coxme_no_comorb_PF$loglik[1]) - log(coxme_no_comorb_PF$n[1]) * coxme_no_comorb_PF$df[1]
        )
    )



writexl::write_xlsx(
    df_models_aic_bic_PF
    , "output/supplementary/model_sensitivity/table_model_PF_backward_aic_bic.xlsx"
)












































###############################################################################
###############################################################################
## Imputing values of PF ratio
## Using chosen Cox mortality model patients with PF missing



# Propensity score - Imputed model ----------------------------------------

# Propensity score estimation 
# (final model in 'propensity_score_mode_PF.R')
model_no_var_PF_imp <- 
    first_resp_support ~
    idade_grupo +
    Gender +
    MFI_level +
    SofaScore +
    period +
    emergency +
    obesity +
    IsVasopressors 

## Estimating Propensity scores and weights - All patients (missing and non missing P/F Ratio)
ps_values_PF_imp <- 
    weightit(
        model_no_var_PF_imp,
        family = "binomial",
        estimand = "ATE",
        data = df_model,
        method = "ps"
    )


## Obtaining propensity scores and IPT weights for each patients
df_cox_model_ps_PF_imp <-
    df_model %>%
    bind_cols(
        ps_ate = ps_values_PF_imp$ps,
        ps_w_ate = ps_values_PF_imp$weights
    ) 
    


# Sensitivity Analysis: Imputed PaO2/FiO2  --------
library(mice)

# Missing pattern
pilot_impute <- mice(df_cox_model_ps_PF_imp, maxit = 0, m = 1)

# first predictor matrix
pred_matrix <- pilot_impute$predictorMatrix

# Defining variables for imputing PaO2/FiO2 values
pred_matrix[, colnames(pred_matrix)] <- 0 

pred_matrix["PaO2FiO2", c("Age", 
                          "idade_grupo",
                          "Gender",
                          "emergency",
                          "MFI_level",
                          "MFIpoints", 
                          "Saps3Points",
                          "SofaScore",
                          "hypertension",
                          "diabetes",
                          "cardio_disease",
                          "obesity",
                          "copd_asthma",
                          "VentSupport",
                          "first_resp_support",
                          "IsVasopressors",
                          "IsRenalReplacementTherapy",
                          "period",
                          "HospitalCode",
                          "outcome_follow_up",
                          "hosp_los_follow_up")] <- 1


# Final imputation (stantard settings)
df_imputed_mice <- mice(df_cox_model_ps_PF_imp, 
                        predictorMatrix = pred_matrix, 
                        seed = 20210107, m = 100)

# Retrieving imputed datasets and producing new "Lactate" variable
complete_datasets <- 
    complete(df_imputed_mice,
             action = "long",
             include = TRUE) %>% 
    mutate(
        PaO2FiO21h_level = case_when(
            PaO2FiO2 <= 100 ~ "severe",
            between(PaO2FiO2, 101, 200) ~ "moderate",
            between(PaO2FiO2, 201, 300) ~ "mild",
            PaO2FiO2 > 300 ~ "normal"
        )    
    ) %>% 
    mutate(
        PaO2FiO21h_level = factor(PaO2FiO21h_level,
                                  levels = c("normal", 
                                             "mild",
                                             "moderate",
                                             "severe")
        )
    ) 


# Turning datasets into MICE object
complete_datasets_mice <- as.mids(complete_datasets)



# Full model All Admissions + imputed PaO2/FiO2
df_cox_model_full_PF_imp <- 
    pool(
        with(
            complete_datasets_mice,
            expr =
                coxme(
                    Surv(hosp_los_follow_up, outcome_follow_up) ~
                        Gender +
                        idade_grupo +
                        MFI_level +
                        SofaScore +
                        Saps3Q +
                        period +
                        emergency +
                        hypertension +
                        diabetes +
                        obesity +
                        copd_asthma +
                        cardio_disease +
                        ResourceIsVasopressors +
                        ResourceIsRenalReplacementTherapy +
                        first_resp_support +
                        PaO2FiO21h_level +
                        (1 | HospitalCode),
                    weights = ps_w_ate
                )
        ) 
        , dfcom = df_cox_model_full_PF$df[1]
    ) %>%
    summary() %>%
    # tidy()
    mutate(
        HR = exp(estimate),
        HR_low = exp(estimate - abs(qt(0.025, df)) * std.error),
        HR_high = exp(estimate + abs(qt(0.025, df)) * std.error),
        HR_ci = paste0(round(HR, 2), " (", round(HR_low, 2), " - ", round(HR_high, 2), ")"),
        p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3))
    ) %>%
    select(term, HR_ci, p_adj)



writexl::write_xlsx(
    df_cox_model_full_PF_imp
    , "output/supplementary/model_sensitivity/cox_model_full_PF_imp.xlsx"
    )
    





