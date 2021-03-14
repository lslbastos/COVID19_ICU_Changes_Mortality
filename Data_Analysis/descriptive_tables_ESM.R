# Library -----------------------------------------------------------------
library(tidyverse)
library(tidylog)
library(gtsummary)

# Obtaining estimated breakpoints
source("code/Data_Analysis/temporal_breakpoints_covid.R")


# Input data --------------------------------------------------------------
df_covid_admissions <-
    read_csv("input/df_all_unique_60day_28_10_2020.csv") %>% 
    filter(StatusCOVID19 == "confirmado") %>%
    mutate(
        AdmissionSource = case_when(
            AdmissionSourceName == "Emergência" ~ "Emergency",
            AdmissionSourceName == "Enfermaria / quarto" ~ "Ward/Room",
            TRUE ~ "Other locations"
        )
    ) %>% 
    mutate(
        MFI_level = factor(MFI_level,
                           levels = c("non_frail", "pre_frail", "frail"),
        ),
        idade_grupo = factor(idade_grupo,
                             levels = c("<40",
                                        "40-49",
                                        "50-59",
                                        "60-69",
                                        "70-79",
                                        ">=80")
        ),
        AdmissionSource = factor(AdmissionSource,
                                 levels = c("Emergency",
                                            "Ward/Room",
                                            "Other locations"
                                 )
        ),
        PaO2FiO21h_level = factor(PaO2FiO21h_level,
                                  levels = c("normal",
                                             "mild",
                                             "moderate",
                                             "severe")
        ),
        VentSupport = factor(VentSupport,
                             levels = c("none", "only_niv", 
                                        "niv_to_mv", "only_mv")
        )
    ) %>%
    arrange(UnitAdmissionDate) %>% 
    mutate(
        is_admission     = 1,
        Saps3Q           = cut_number(Saps3Points, n = 4), # SAPS3 quartiles
        outcome_icu_new  = if_else(outcome_icu  != "Ongoing", outcome_icu_new,  NA_real_),
        outcome_hosp_new = if_else(outcome_hosp != "Ongoing", outcome_hosp_new, NA_real_),

        icu_los_adj      = if_else(outcome_icu  != "Ongoing", icu_los_adj,      NA_real_),
        hosp_los_adj     = if_else(outcome_hosp != "Ongoing", hosp_los_adj,     NA_real_),
        icu_los_adj      = if_else(outcome_icu  != "Ongoing", icu_los_adj,      NA_real_),
        hosp_los_adj     = if_else(outcome_hosp != "Ongoing", hosp_los_adj,     NA_real_),
        icu_stay_7       = if_else(outcome_icu  != "Ongoing", icu_stay_7,       NA_real_),
        hosp_stay_7      = if_else(outcome_hosp != "Ongoing", hosp_stay_7,      NA_real_),
        
        is_ongoing_hosp  = if_else(outcome_hosp == "Ongoing", 1, 0),
        is_ongoing_icu   = if_else(outcome_icu  == "Ongoing", 1, 0),
        
        is_ventilated    = if_else(VentSupport != "none",      1, 0),
        resp_oxygen      = if_else(VentSupport == "none",      1, 0),
        resp_only_niv    = if_else(VentSupport == "only_niv",  1, 0),
        resp_niv_fail    = if_else(VentSupport == "niv_to_mv", 1, 0),
        resp_only_mv     = if_else(VentSupport == "only_mv",  1, 0),
        niv_mode         = case_when(
                                (ResourceIsNonInvasiveVentilation == 1 | IsNonInvasiveVentilation == 1)
                                    & ResourceIsHFNC == 1 ~ "both",
                                (ResourceIsNonInvasiveVentilation == 1 | IsNonInvasiveVentilation == 1) ~ "only_nppv",
                                ResourceIsHFNC == 1 ~ "only_hfnc",
                                 TRUE ~ NA_character_
                             ),
        outcome_follow_up_niv = if_else(VentSupport    == "only_niv",  outcome_follow_up, NA_real_),
        outcome_follow_up_niv_mv = if_else(VentSupport == "niv_to_mv",  outcome_follow_up, NA_real_),
        outcome_follow_up_imv = if_else(VentSupport    == "only_mv",  outcome_follow_up, NA_real_),
        niv_first = if_else(VentSupport %in% c("only_niv", "niv_to_mv"), 1, 0),
        imv_first = if_else(VentSupport %in% c("only_mv"), 1, 0)
    ) %>% 
    mutate(
        niv_mode         = factor(niv_mode, levels = c("only_nppv", "only_hfnc", "both"))
        ) %>% 
    left_join(
        bp_dates,
        by = c("UnitAdmissionDate" = "date")
        ) %>% 
    arrange(UnitAdmissionDate) %>% 
    fill(period, .direction = "up") %>% 
    mutate(
        period = ifelse(is.na(period), max(period, na.rm = TRUE) + 1, period)
    )










# Descriptive - Patients under advanced respiratory support per period (sTable 3) -------
tb_admission_mv <- 
    df_covid_admissions %>%
    filter(VentSupport != "none") %>% 
    mutate(
        Saps3Q = cut_number(Saps3Points, n = 4), # SAPS3 quartiles
        MFIpoints_2 = MFIpoints
    ) %>% 
    select(
        is_admission,
        Age,
        idade_grupo,
        Gender,
        AdmissionSource,
        MFIpoints,
        MFIpoints_2,
        MFI_level,
        Saps3Points,
        Saps3Q,
        SofaScore,
        PaO2FiO2,
        PaO2FiO21h_level,
        resp_oxygen,
        is_ventilated,
        had_niv,
        niv_mode,
        resp_only_niv,
        resp_niv_fail,
        resp_only_mv,
        ResourceIsVasopressors,
        ResourceIsRenalReplacementTherapy,
        icu_los_adj,
        hosp_los_adj,
        icu_stay_7,
        hosp_stay_7,
        outcome_follow_up,
        outcome_follow_up_niv,
        outcome_follow_up_niv_mv,
        outcome_follow_up_imv,
        outcome_icu_new,
        outcome_hosp_new,
        is_ongoing_hosp,
        period
    ) %>%
    tbl_summary(missing = "no",
                type = list(MFIpoints ~ "continuous",
                            MFIpoints_2 ~ "continuous"),
                statistic = list(MFIpoints ~ "{mean} ({sd})"),
                label = ls_labels,
                by = period
    ) %>% 
    add_overall() %>% 
    add_n()


tb_admission_mv_desc <- 
    tb_admission_mv$table_body %>% 
    filter(label != "") %>% 
    mutate(label = ifelse(row_type == "label" & variable %in% c("icu_los_adj", "hosp_los_adj",
                                                                "icu_stay_7", "hosp_stay_7",
                                                                "outcome_icu_new", "outcome_hosp_new",
                                                                "outcome_follow_up_niv", "outcome_follow_up_niv_mv",
                                                                "outcome_follow_up_imv", "PaO2FiO2"), 
                          paste0(label, " [n = ", n, "]"), label)) %>%  
    rename("Characteristics" = "label",
           # "N" = "n",
           "Overall" = "stat_0",
           "Period 1" = "stat_1",
           "Period 2" = "stat_2",
           "Period 3" = "stat_3",
           "Period 4" = "stat_4") %>% 
    # "Northeast" = "stat_2") %>%
    select(-c(row_type, n, variable))



writexl::write_xlsx(
    tb_admission_mv_desc
    , "output/supplementary/descriptive_supplementary/table_descriptive_mv_patients_28_10.xlsx"
    )









# Descriptive - Comorbidities - All admissions (sTable 4) -----
ls_labels_comorbidities <- 
    list(
        is_admission     ~ "Total, No (%)",
        n_comorb         ~ "Any comorbidities, No. (%)",
        hypertension     ~ "Hypertension",
        diabetes         ~ "Diabetes",
        cardio_disease   ~ "Cardiovascular disease",
        imunossupression ~ "Immunosupression",
        cerebro_disease  ~ "Cerebrovascular disease ",
        obesity          ~ "Obesity ",
        copd_asthma      ~ "COPD or Asthma ",
        malignancy       ~ "Malignancy",
        chronic_kidney   ~ "Chronic kidney disease",
        tobacco          ~ "Tobacco",
        liver_cirrhosis  ~ "Liver cirrhosis",
        other_comorb     ~ "Other comorbidities",
        CharlsonComorbidityIndex   ~ "Charlson Comorbidity Index (CCI), mean (SD)",
        CharlsonComorbidityIndex_2 ~ "Charlson Comorbidity Index (CCI), median (IQR)"
    )

## Comorbidities - All Admissions
tb_admission_comorb <- 
    df_covid_admissions %>%
    mutate(
        CharlsonComorbidityIndex_2 = CharlsonComorbidityIndex
    ) %>% 
    select(
        is_admission,
        n_comorb,
        hypertension,
        diabetes,
        imunossupression,
        cardio_disease,
        obesity,
        copd_asthma,
        malignancy,
        cerebro_disease,
        chronic_kidney,
        tobacco,
        liver_cirrhosis,
        other_comorb,
        CharlsonComorbidityIndex,
        CharlsonComorbidityIndex_2,
        period
    ) %>%
    tbl_summary(missing = "no",
                type = list(CharlsonComorbidityIndex ~ "continuous",
                            CharlsonComorbidityIndex_2 ~ "continuous"
                ),
                statistic = list(CharlsonComorbidityIndex ~ "{mean} ({sd})"),
                label = ls_labels_comorbidities,
                by = period
    ) %>% 
    add_overall() %>% 
    add_n()



writexl::write_xlsx(
    tb_admission_comorb$table_body %>%  
        rename("Comorbidities" = "label",
               # "N" = "n",
               "Overall" = "stat_0",
               "Period 1" = "stat_1",
               "Period 2" = "stat_2",
               "Period 3" = "stat_3",
               "Period 4" = "stat_4") %>% 
        select(-c(row_type, n, variable))
    , "output/supplementary/descriptive_supplementary/table_comorbidities_all_admissions_28_10.xlsx"
)





# Descriptive - Comorbidities - Patients Advanced Resp. Support (sTable 4) -----
tb_admission_mv_comorb <- 
    df_covid_admissions %>%
    filter(VentSupport != "none") %>% 
    mutate(
        CharlsonComorbidityIndex_2 = CharlsonComorbidityIndex
    ) %>% 
    select(
        is_admission,
        n_comorb,
        hypertension,
        diabetes,
        imunossupression,
        cardio_disease,
        obesity,
        copd_asthma,
        malignancy,
        cerebro_disease,
        chronic_kidney,
        tobacco,
        liver_cirrhosis,
        other_comorb,
        CharlsonComorbidityIndex,
        CharlsonComorbidityIndex_2,
        period
    ) %>%
    tbl_summary(missing = "no",
                type = list(CharlsonComorbidityIndex ~ "continuous",
                            CharlsonComorbidityIndex_2 ~ "continuous"
                ),
                statistic = list(CharlsonComorbidityIndex ~ "{mean} ({sd})"),
                label = ls_labels_comorbidities,
                by = period
    ) %>% 
    add_overall() %>% 
    add_n()



writexl::write_xlsx(
    tb_admission_mv_comorb$table_body %>%  
        rename("Comorbidities" = "label",
               # "N" = "n",
               "Overall" = "stat_0",
               "Period 1" = "stat_1",
               "Period 2" = "stat_2",
               "Period 3" = "stat_3",
               "Period 4" = "stat_4") %>% 
        select(-c(row_type, n, variable))
    , "output/supplementary/descriptive_supplementary/table_comorbidities_mv_patients_28_10.xlsx"
)













# Descriptive: NIRS x No NIRS - Patients Under Adv. Resp. Support (sTable 5) ------------
ls_labels_NIV <-
    list(
        is_admission      ~ "Total, No (%)",
        Age               ~ "Age, mean (SD)", 
        idade_grupo       ~ "Age group , No. (%)", 
        Gender            ~ "Sex, No. (%)", 
        AdmissionSource   ~ "Admission Source, No. (%)",
        MFI_level         ~ "Modified Frailty Index (MFI), No. (%)", 
        Saps3Points       ~ "SAPS-3, median (IQR)",
        Saps3Q            ~ "",
        SofaScore         ~ "SOFA score, median (IQR)", 
        
        
        n_comorb          ~ "Any comorbidities, No. (%)",
        hypertension      ~ "Hypertension ",
        diabetes          ~ "Diabetes",
        cardio_disease    ~ "Cardiovascular disease",
        obesity           ~ "Obesity ",
        copd_asthma       ~ "COPD or Asthma ",

        PaO2FiO2          ~ "PaO2/FiO2, median (IQR)",
        PaO2FiO21h_level  ~ "PaO2/FiO2, No. (%)",
        
        resp_oxygen       ~ "Oxygen support, No. (%)",
        resp_only_niv     ~ "Only NIRS",
        resp_niv_fail     ~ "NIRS failure",
        resp_only_mv      ~ "Invasive mechanical ventilation",
        
        IsVasopressors    ~ "Vasopressors (first 24h), No. (%)",
        IsRenalReplacementTherapy         ~ "Renal Replacement Therapy (first 24h), No. (%)", 
        ResourceIsVasopressors            ~ "Vasopressor, No. (%)", 
        ResourceIsRenalReplacementTherapy ~ "Renal Replacement Therapy, No. (%)",
        
        icu_los_adj       ~ "ICU Length-of-stay (LOS), median (IQR)", 
        hosp_los_adj      ~ "Hospital Length-of-stay (LOS), median (IQR)",
        icu_stay_7        ~ "ICU Length-of-stay (LOS) > 7 days, No. (%)", 
        hosp_stay_7       ~ "Hospital Length-of-stay (LOS) > 7 days, No. (%)",
        outcome_follow_up ~ "60-day In-hospital deaths, No. (%)",
        period            ~ "Periods, No. (%)"
    )



tb_admission_NIV <- 
    df_covid_admissions %>%
    filter(VentSupport != "none") %>% 
    select(
        is_admission,
        Age,
        idade_grupo,
        Gender,
        AdmissionSource,
        MFI_level,
        Saps3Points,
        Saps3Q,
        SofaScore,
        
        n_comorb,
        hypertension,
        diabetes,
        cardio_disease,
        obesity,
        copd_asthma,
        
        PaO2FiO2,
        PaO2FiO21h_level,
        
        resp_oxygen,
        resp_only_niv,
        resp_niv_fail,
        resp_only_mv,
        
        IsVasopressors,
        IsRenalReplacementTherapy,
        ResourceIsVasopressors,
        ResourceIsRenalReplacementTherapy,
        
        icu_los_adj,
        hosp_los_adj,
        icu_stay_7,
        hosp_stay_7,
        outcome_follow_up,
        period,
        had_niv
    ) %>%
    tbl_summary(missing = "no",
                label = ls_labels_NIV,
                by = had_niv
    ) %>% 
    add_overall() %>% 
    add_n()


tb_admission_NIV_desc <- 
    tb_admission_NIV$table_body %>% 
    filter(label != "") %>% 
    mutate(label = ifelse(row_type == "label" & variable %in% c("icu_los_adj", "hosp_los_adj",
                                                                "icu_stay_7", "hosp_stay_7",
                                                                "outcome_icu_new", "outcome_hosp_new",
                                                                "PaO2FiO2"), 
                          paste0(label, " [n = ", n, "]"), label)) %>%  
    rename("Characteristics" = "label",
           "All admissions" = "stat_0",
           "No NIRS" = "stat_1",
           "NIRS"    = "stat_2") %>% 
    select(-c(row_type, n, variable))



writexl::write_xlsx(
    tb_admission_NIV_desc
    , "output/supplementary/descriptive_supplementary/table_NIRS_use_28_10.xlsx"
)













# Descriptive: Missing pattern of PaO2/FiO2 (sTable 11) --------
ls_labels_missing <-
    list(
        is_admission      ~ "Total, No (%)",
        Age               ~ "Age, mean (SD)", 
        idade_grupo       ~ "Age group , No. (%)", 
        Gender            ~ "Sex, No. (%)", 
        AdmissionSource   ~ "Admission Source, No. (%)",
        MFI_level         ~ "Modified Frailty Index (MFI), No. (%)", 
        Saps3Points       ~ "SAPS-3, median (IQR)",
        Saps3Q            ~ "",
        SofaScore         ~ "SOFA score, median (IQR)", 

        hypertension      ~ "Hypertension ",
        diabetes          ~ "Diabetes",
        cardio_disease    ~ "Cardiovascular disease",
        obesity           ~ "Obesity ",
        copd_asthma       ~ "COPD or Asthma ",
        
        resp_oxygen       ~ "Oxygen support, No. (%)",
        is_ventilated     ~ "Advanced respiratory support, No. (%)",
        resp_only_niv     ~ "Only NIRS",
        resp_niv_fail     ~ "NIRS failure",
        resp_only_mv      ~ "Invasive mechanical ventilation",
        
        IsVasopressors            ~ "Vasopressor (first 24h), No. (%)", 
        IsRenalReplacementTherapy ~ "Renal Replacement Therapy (first 24h), No. (%)",
        
        icu_los_adj              ~ "ICU Length-of-stay (LOS), median (IQR)", 
        hosp_los_adj             ~ "Hospital Length-of-stay (LOS), median (IQR)",
        outcome_follow_up        ~ "60-day In-hospital deaths, No. (%)"
    )



## Missing pattern of PaO2/FiO2 - All admissions
tb_admission_missing <- 
    df_covid_admissions %>%
    mutate(
        missingPF = if_else(is.na(PaO2FiO21h_level), "missing", "not_missing")
    ) %>% 
    # filter(VentSupport != "none") %>% 
    droplevels() %>% 
    select(
        is_admission,
        Age,
        idade_grupo,
        Gender,
        AdmissionSource,
        MFI_level,
        Saps3Points,
        Saps3Q,
        SofaScore,
        hypertension,
        diabetes,
        cardio_disease,
        obesity,
        copd_asthma,
        resp_oxygen,
        is_ventilated,
        resp_only_niv,
        resp_niv_fail,
        resp_only_mv,
        IsVasopressors,
        IsRenalReplacementTherapy,
        period,
        icu_los_adj,
        hosp_los_adj,
        outcome_follow_up,
        missingPF
    ) %>%
    tbl_summary(missing = "no",
                label = ls_labels_missing,
                by = missingPF
    ) %>% 
    add_overall() %>% 
    add_n()



tb_admission_missing_desc <- 
    tb_admission_missing$table_body %>% 
    filter(label != "") %>% 
    mutate(label = ifelse(row_type == "label" & variable %in% c("icu_los_adj", "hosp_los_adj",
                                                                "icu_stay_7", "hosp_stay_7",
                                                                "outcome_icu_new", "outcome_hosp_new"), 
                          paste0(label, " [n = ", n, "]"), label)) %>%  
    rename("Characteristics" = "label",
           "Overall" = "stat_0",
           "Missing P/F ratio" = "stat_1",
           "Not Missing P/F ratio" = "stat_2"
    ) %>% 
    select(-c(row_type, n, variable))



writexl::write_xlsx(
    tb_admission_missing_desc
    , "output/supplementary/descriptive_supplementary/table_PF_missing_pattern_all_admissions_28_10.xlsx"
)



## Missing pattern of PaO2/FiO2 - MV patients
tb_admission_missing_MV <- 
    df_covid_admissions %>%
    mutate(
        missingPF = if_else(is.na(PaO2FiO21h_level), "missing", "not_missing")
    ) %>% 
    filter(VentSupport != "none") %>%
    droplevels() %>% 
    select(
        is_admission,
        Age,
        idade_grupo,
        Gender,
        AdmissionSource,
        MFI_level,
        Saps3Points,
        Saps3Q,
        SofaScore,
        hypertension,
        diabetes,
        cardio_disease,
        obesity,
        copd_asthma,
        resp_oxygen,
        is_ventilated,
        resp_only_niv,
        resp_niv_fail,
        resp_only_mv,
        IsVasopressors,
        IsRenalReplacementTherapy,
        period,
        icu_los_adj,
        hosp_los_adj,
        outcome_follow_up,
        missingPF
    ) %>%
    tbl_summary(missing = "no",
                label = ls_labels_missing,
                by = missingPF
    ) %>% 
    add_overall() %>% 
    add_n()



tb_admission_missing_MV_desc <- 
    tb_admission_missing_MV$table_body %>% 
    filter(label != "") %>% 
    mutate(label = ifelse(row_type == "label" & variable %in% c("icu_los_adj", "hosp_los_adj",
                                                                "icu_stay_7", "hosp_stay_7",
                                                                "outcome_icu_new", "outcome_hosp_new"), 
                          paste0(label, " [n = ", n, "]"), label)) %>%  
    rename("Characteristics" = "label",
           "Overall" = "stat_0",
           "Missing P/F ratio" = "stat_1",
           "Not Missing P/F ratio" = "stat_2"
    ) %>% 
    select(-c(row_type, n, variable))



writexl::write_xlsx(
    tb_admission_missing_MV_desc
    , "output/supplementary/descriptive_supplementary/table_PF_missing_pattern_mv_patients_28_10.xlsx"
)



# Finished