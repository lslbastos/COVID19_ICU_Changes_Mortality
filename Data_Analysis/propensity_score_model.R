################################################################################
## Article Evolving changes in mortality of 13,301 critically ill adult 
##   patients with COVID-19 over eight months
## 
## Propensity score estimation - (Main Analysis)
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
df_model <-
    df_covid_admissions %>%
    filter(
        VentSupport != "none"
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
model_full <- 
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
    IsRenalReplacementTherapy

ps_model_full <- 
    glm(model_full, # RRT in the first 24h of admission,
        family = "binomial",
        data = df_model
        ) 


## No significant comorbidities
model_no_comorb <- 
    first_resp_support ~
    idade_grupo +
    Gender +
    MFI_level +
    SofaScore +
    period +
    emergency +
    # hypertension +
    # diabetes +
    # obesity +
    # copd_asthma +
    cardio_disease +
    IsVasopressors + # Vasopressor in the first 24h of admission
    IsRenalReplacementTherapy


ps_model_no_comorb <- 
    glm(model_no_comorb, # RRT in the first 24h of admission,
        family = "binomial",
        data = df_model
        ) 



# Model without any significant variables
model_no_var <- 
    first_resp_support ~
    idade_grupo +
    Gender +
    MFI_level +
    SofaScore +
    period +
    emergency +
    # hypertension +
    # diabetes +
    # obesity +
    # copd_asthma +
    cardio_disease +
    IsVasopressors # Vasopressor in the first 24h of admission
# IsRenalReplacementTherapy, # RRT in the first 24h of admission,


ps_model_no_var <- 
    glm(model_no_var,
        family = "binomial",
        data = df_model
    ) 


## Table with model comparison (sTable 6)
df_ps_models_comparison <- 
    left_join(
        ps_model_full %>% 
            tidy(exponentiate = TRUE,  conf.int = TRUE) %>%
            mutate(
                p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
                HR_ci = paste0(round(estimate, 2), 
                               " (", round(conf.low, 2), " - ", round(conf.high, 2), ")")
            ) %>%
            select(term, HR_ci, p_adj),
        ps_model_no_comorb %>% 
            tidy(exponentiate = TRUE,  conf.int = TRUE) %>%
            mutate(
                p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
                HR_ci = paste0(round(estimate, 2), 
                               " (", round(conf.low, 2), " - ", round(conf.high, 2), ")")
            ) %>%
            select(term, HR_ci_no_comorb = HR_ci, p_adj_no_comorb = p_adj)
    ) %>% 
    left_join(
        ps_model_no_var %>% 
            tidy(exponentiate = TRUE,  conf.int = TRUE) %>%
            mutate(
                p_adj = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
                HR_ci = paste0(round(estimate, 2), 
                               " (", round(conf.low, 2), " - ", round(conf.high, 2), ")")
            ) %>%
            select(term, HR_ci_no_var = HR_ci, p_adj_no_var = p_adj)
    )

writexl::write_xlsx(df_ps_models_comparison, 
                    "output/supplementary/propensity_score_supplementary/table_PS_models_comparison.xlsx")



## Estimating AUC and Brier with bootstrap

set.seed(2^31 - 1)
ls_model_boot <- 
    df_model %>% 
    infer::generate(reps = 1000, type = "bootstrap") %>% 
    split(.$replicate)
    

df_ps_model_full <- 
    model_perf_boot(formula = model_full, list = ls_model_boot)

df_ps_model_no_comorb <- 
    model_perf_boot(formula = model_no_comorb, list = ls_model_boot)

df_ps_model_no_var <- 
    model_perf_boot(formula = model_no_var, list = ls_model_boot)



## Comparison of PS models and metrics
df_ps_models_metric <- 
    tibble(
        model = c(
            "ps_model_full",
            "ps_model_no_comorb",
            "ps_model_no_var"
            ),
        AIC = c(
            AIC(ps_model_full),
            AIC(ps_model_no_comorb),
            AIC(ps_model_no_var)
        ),
        BIC = c(
            BIC(ps_model_full),
            BIC(ps_model_no_comorb),
            BIC(ps_model_no_var)
        ),
        AUROC = c(
            pROC::auc(pROC::roc(predictor = ps_model_full$fitted.values, 
                                response  = ps_model_full$y)),
            pROC::auc(pROC::roc(predictor = ps_model_no_comorb$fitted.values, 
                                response  = ps_model_no_comorb$y)),
            pROC::auc(pROC::roc(predictor = ps_model_no_var$fitted.values, 
                                response  = ps_model_no_var$y))
        ),
        AUROC_CI = c(
            df_ps_model_full$auc_ci,
            df_ps_model_no_comorb$auc_ci,
            df_ps_model_no_var$auc_ci
        ),
        Brier = c(
            ModelMetrics::brier(ps_model_full),
            ModelMetrics::brier(ps_model_no_comorb),
            ModelMetrics::brier(ps_model_no_var)
        ),
        Brier_CI = c(
            df_ps_model_full$brier_ci,
            df_ps_model_no_comorb$brier_ci,
            df_ps_model_no_var$brier_ci
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
writexl::write_xlsx(df_ps_models_metric, 
                    "output/supplementary/propensity_score_supplementary/table_PS_models_metric.xlsx")



# Finished

## The best model is the 'reduced' model without any significant variables (lowest AIC and BIC)
# df_ps_models_comparison
# model_no_var





# Evaluation of propensity scores -----------------------------------------


## Obtaining PS scores for each patient from the final PS model (model_no_var)
df_model_ps_ate <-
    df_model %>%
    bind_cols(
        ps = as.numeric(predict(ps_model_no_var, type = "response")) 
    ) 




# Distribution of PS values x outcomes (NIRS first x IMV first) - sFigure 7 -----------
# Our binary variable is: 1 - NIRS first, 0 - IMV first
# Hence, propensity scores are probability in respect to receiving NIRS first

plot_ps_dist <- 
    df_model_ps_ate %>%
    ggplot() +
    geom_density(aes(x = ps, color = first_resp_support)) +
    scale_fill_discrete(labels = c("IMV first", "NIRS first"),
                        name = "First respiratory support") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = "Propensity Score (PS)", y = "Density") +
    theme_bw() +
    theme(legend.position = "top")

plot_ps_hist <- 
    df_model_ps_ate %>%
    ggplot() +
    geom_histogram(aes(x = ps, fill = first_resp_support)) +
    scale_fill_discrete(labels = c("IMV first", "NIRS first"),
                        name = "First respiratory support") +
    scale_y_continuous(labels = scales::comma_format()) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = "Propensity Score (PS)", y = "Frequency") +
    theme_bw() +
    theme(legend.position = "top")


# Exporting plots
ggsave("output/supplementary/propensity_score_supplementary/figure_ps_distribution.png",
       ggpubr::ggarrange(plot_ps_hist, 
                         plot_ps_dist, common.legend = TRUE), 
       width =  10, height = 5, dpi = 800)


# Finished






