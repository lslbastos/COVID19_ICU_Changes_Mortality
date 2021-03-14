################################################################################
## Article Evolving changes in mortality of 13,301 critically ill adult 
##   patients with COVID-19 over eight months
## 
## Descriptive tables - Main Analysis
## Leonardo S.L. Bastos (lslbastos), Pedro Kurtz 
##
################################################################################



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











# Descriptive - All Patients  (Table 1) -------
ls_labels <-
    list(
        is_admission      ~ "Total, No (%)",
        Age               ~ "Age, mean (SD)", 
        idade_grupo       ~ "Age group , No. (%)", 
        Gender            ~ "Sex, No. (%)", 
        AdmissionSource   ~ "Admission Source, No. (%)",
        MFIpoints         ~ "Modified Frailty Index (MFI), mean (SD)", 
        MFIpoints_2       ~ "Modified Frailty Index (MFI), median (IQR)", 
        MFI_level         ~ "Modified Frailty Index (MFI), No. (%)", 
        Saps3Points       ~ "SAPS-3, median (IQR)",
        Saps3Q            ~ "",
        SofaScore         ~ "SOFA score, median (IQR)", 
        PaO2FiO2          ~ "PaO2/FiO2, median (IQR)",
        PaO2FiO21h_level  ~ "PaO2/FiO2, No. (%)",
        
        resp_oxygen       ~ "Oxygen support, No. (%)",
        is_ventilated     ~ "Advanced respiratory support, No. (%)",
        had_niv           ~ "Noninvasive respiratory support (NIRS), No. (%)",
        resp_only_niv     ~ "Only NIRS",
        resp_niv_fail     ~ "NIRS failure",
        niv_mode          ~ "",
        resp_only_mv      ~ "Invasive mechanical ventilation",
        
        ResourceIsVasopressors            ~ "Vasopressor, No. (%)", 
        ResourceIsRenalReplacementTherapy ~ "Renal Replacement Therapy, No. (%)",
        
        icu_los_adj              ~ "ICU Length-of-stay (LOS), median (IQR)", 
        hosp_los_adj             ~ "Hospital Length-of-stay (LOS), median (IQR)",
        icu_stay_7               ~ "ICU Length-of-stay (LOS) > 7 days, No. (%)", 
        hosp_stay_7              ~ "Hospital Length-of-stay (LOS) > 7 days, No. (%)",
        outcome_follow_up        ~ "60-day In-hospital deaths, No. (%)",
        outcome_follow_up_niv    ~ "60-day In-hospital deaths (NIV only), No. (%)",
        outcome_follow_up_niv_mv ~ "60-day In-hospital deaths (NIV failure), No. (%)",
        outcome_follow_up_imv    ~ "60-day In-hospital deaths (IMV), No. (%)",
        outcome_icu_new          ~ "ICU Deaths, No, (%)", 
        outcome_hosp_new         ~ "In-hospital deaths, No. (%)",
        is_ongoing_hosp          ~ "Ongoing patients, No. (%)"
    )




## Descriptive table object
tb_admission <- 
    df_covid_admissions %>%
    mutate(
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


tb_admission_desc <- 
    tb_admission$table_body %>% 
    filter(label != "") %>% 
    mutate(label = ifelse(row_type == "label" & variable %in% c("icu_los_adj", "hosp_los_adj",
                                                                "icu_stay_7", "hosp_stay_7",
                                                                "outcome_icu_new", "outcome_hosp_new",
                                                                "outcome_follow_up_niv", "outcome_follow_up_niv_mv",
                                                                "outcome_follow_up_imv", "PaO2FiO2"), 
                          paste0(label, " [n = ", n, "]"), label)) %>%  
    rename("Characteristics" = "label",
           "All admissions" = "stat_0",
           "Period 1" = "stat_1",
           "Period 2" = "stat_2",
           "Period 3" = "stat_3",
           "Period 4" = "stat_4") %>% 
    select(-c(row_type, n, variable))



writexl::write_xlsx(
    tb_admission_desc
    , "output/main/descriptive/table1_descriptive_all_admissions_28_10.xlsx"
    )







# Descriptive - Adv. Respiratory Support modes (Table 2)  --------
ls_labels_RS <-
    list(
        is_admission      ~ "Total, No (%)",
        Age               ~ "Age, mean (SD)", 
        idade_grupo       ~ "Age group , No. (%)", 
        Gender            ~ "Sex, No. (%)", 
        MFI_level         ~ "Modified Frailty Index (MFI), No. (%)",
        Saps3Points       ~ "SAPS-3, median (IQR)",
        Saps3Q            ~ "",
        SofaScore         ~ "SOFA score, median (IQR)",
        n_comorb          ~ "Any comorbidities, No. (%)",
        PaO2FiO2          ~ "PaO2/FiO2, median (IQR)",
        PaO2FiO21h_level  ~ "PaO2/FiO2, No. (%)",
        niv_mode          ~ "NIRS mode, No. (%)",
        ResourceIsVasopressors            ~ "Vasopressor, No. (%)", 
        ResourceIsRenalReplacementTherapy ~ "Renal Replacement Therapy, No. (%)",
        icu_los_adj       ~ "ICU Length-of-stay (LOS), median (IQR)", 
        hosp_los_adj      ~ "Hospital Length-of-stay (LOS), median (IQR)",
        icu_stay_7        ~ "ICU Length-of-stay (LOS) > 7 days, No. (%)", 
        hosp_stay_7       ~ "Hospital Length-of-stay (LOS) > 7 days, No. (%)",
        outcome_follow_up ~ "60-day In-hospital deaths, No. (%)",
        outcome_icu_new   ~ "ICU Deaths, No, (%)", 
        outcome_hosp_new  ~ "In-hospital deaths, No. (%)",
        period            ~ "Periods, No. (%)"
    )


tb_admission_RS <- 
    df_covid_admissions %>%
    filter(VentSupport != "none") %>% 
    mutate(
        Saps3Q = cut_number(Saps3Points, n = 4), # SAPS3 quartiles (reestimation for MV patients)
    ) %>% 
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
        n_comorb,
        PaO2FiO2,
        PaO2FiO21h_level,
        niv_mode,
        ResourceIsVasopressors,
        ResourceIsRenalReplacementTherapy,
        icu_los_adj,
        hosp_los_adj,
        icu_stay_7,
        hosp_stay_7,
        period,
        outcome_follow_up,
        outcome_icu_new,
        outcome_hosp_new,
        VentSupport
    ) %>%
    tbl_summary(missing = "no",
                label = ls_labels_RS,
                by = VentSupport
    ) %>% 
    add_overall() %>% 
    add_n()



tb_admission_RS_desc <- 
    tb_admission_RS$table_body %>% 
    filter(label != "") %>% 
    mutate(label = ifelse(row_type == "label" & variable %in% c("icu_los_adj", "hosp_los_adj",
                                                                "icu_stay_7", "hosp_stay_7",
                                                                "outcome_icu_new", "outcome_hosp_new",
                                                                "PaO2FiO2"), 
                          paste0(label, " [n = ", n, "]"), label)) %>%  
    rename("Characteristics" = "label",
           "Overall" = "stat_0",
           "NIRS only" = "stat_1",
           "NIRS failure" = "stat_2",
           "IMV" = "stat_3"
           ) %>% 
    select(-c(row_type, n, variable))



writexl::write_xlsx(
    tb_admission_RS_desc
    , "output/main/descriptive/table2_admission_respiratory_supports_28_10.xlsx"
)



# Finished