################################################################################
## Article Evolving changes in mortality of 13,301 critically ill adult 
##   patients with COVID-19 over eight months
##
## Breakpoint estimation for time series and Temporal Plots
## Leonardo S.L. Bastos (lslbastos), Pedro Kurtz 
##
################################################################################

# Defining default date locale
Sys.setlocale(category = "LC_ALL", locale = "english")


# Library -----------------------------------------------------------------
library(tidyverse)
library(tidylog)
library(tsibble)


# Input data --------------------------------------------------------------
df_daily <- 
    read_csv("input/df_covid_volume_time_28_10_2020.csv") %>% 
    mutate(prop_deaths = hosp_deaths_admi / new_icu_admissions)






# Breakpoints estimation - Daily Mortality rates --------------------------

# Obtaining date of first deaths
first_death_adm <- 
    df_daily %>% 
    filter(hosp_deaths > 0) %>% 
    slice(1) %>% 
    pull(date)


# Data filtered by the date of the first death
df_daily_death <- 
    df_daily %>%
    filter(date >= first_death_adm) 

# Creating TS object
ts_mortality <- 
    tsibble::tsibble(
        day = df_daily_death$date,
        y = df_daily_death$icu_deaths,
        ) 


## Estimating Breakpoints using 'strucchange'
breakpoints <- strucchange::breakpoints(ts_mortality$y ~ 1, sort = TRUE)$breakpoints


# Obtaining breakpoint dates
bp_dates <- tibble(
    date = ts_mortality$day[breakpoints],
    bp = ts_mortality$day[breakpoints],
    period = 1:length(breakpoints)
    )  



## Merging breakpoint dates into DF with daily admissions and deaths
df_daily_death_bp <- 
    df_daily %>% 
    left_join(bp_dates, by = c("date" = "date")) %>% 
    mutate(
        period_fill = as.numeric(period)
    ) %>% 
    fill(bp, .direction = "downup") %>% 
    fill(period_fill, .direction = "up") %>% 
    mutate(
        period_fill = if_else(is.na(period_fill), max(period_fill, na.rm = TRUE) + 1, period_fill),
        date_label = if_else(!is.na(period), as.character(bp), NA_character_),
        month = lubridate::month(date),
        bp_fill = if_else(!is.na(bp), as.character(bp), NA_character_)
    ) %>% 
    group_by(period_fill) %>% 
    mutate(
        date_period = paste0(min(date), "-\n", max(date))
    ) %>% 
    ungroup()


## Creating a DF with automatic labels for periods based on breakpoints
df_dates_labels <- 
    df_daily_death_bp %>% 
    distinct(bp_fill, period_fill, date_period) %>% 
    arrange(bp_fill) %>% 
    group_by(period_fill) %>% 
    slice(1) %>% 
    mutate(
        bp_fill = if_else(period_fill == 1, "2020-02-27", bp_fill), 
        date_min = str_split(date_period, "-\n", simplify = TRUE)[1],
        date_max = str_split(date_period, "-\n", simplify = TRUE)[2],
        date_position = as.Date(bp_fill) + (as.Date(date_max) - as.Date(date_min)) / 2 + 4
    ) %>% 
    mutate(
        date_period = paste0(format(as.Date(date_min), "%d/%b"), "-\n", format(as.Date(date_max), "%d/%b"))
    ) %>% 
    ungroup()







################################################################################
# Plots - Daily admissions, deaths, and mortality rates -------------------
## Plots - Cases per day

# Plot New ICU Admissions
plot_daily_new_adm <- 
    df_daily_death_bp %>% 
    ggplot() +
    geom_line(aes(x = date, y = new_icu_admissions)) +
    geom_vline(aes(xintercept = bp), linetype = "dashed") +
    geom_smooth(aes(x = date, y = new_icu_admissions), size = 0.5, se = FALSE) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    labs(title = "New ICU admissions", x = "", y = "") +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))

# Plot ICU deaths
plot_daily_new_deaths <- 
    df_daily_death_bp %>% 
    ggplot() +
    geom_line(aes(x = date, y = icu_deaths)) +
    geom_vline(aes(xintercept = bp), linetype = "dashed") +
    geom_smooth(aes(x = date, y = icu_deaths), size = 0.5, se = FALSE) +
    ylim(c(0, NA)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    labs(title = "ICU deaths", x = "", y = "") +
    theme_bw() +
    theme(legend.position = "none") +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))

# Plot Patients in ICU
plot_daily_hosp <- 
    df_daily_death_bp %>% 
    ggplot() +
    geom_line(aes(x = date, y = total_icu)) +
    geom_vline(aes(xintercept = bp), linetype = "dashed") +
    geom_smooth(aes(x = date, y = total_icu), size = 0.5, se = FALSE) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    geom_text(data = df_dates_labels,
              aes(x = date_position, y = 10, label = date_period), size = 2.7) +
    labs(title = "Patients in the ICU", x = "", y = "") +
    theme_bw() +
    theme(legend.position = "none") +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))


## Daily in-hospital mortality
plot_daily_mortal <-
    df_daily_death_bp %>% 
    ggplot() +
    geom_line(aes(x = date, y = prop_deaths)) +
    geom_vline(aes(xintercept = bp), linetype = "dashed") +
    geom_smooth(aes(x = date, y = prop_deaths), size = 0.5, se = FALSE) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    labs(title = "In-hospital mortality rate", x = "", y = "") +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10))





## Combining Plots (Figure 1)
plot_comb_volume <- 
    ggpubr::annotate_figure(
        ggpubr::ggarrange(plot_daily_hosp,
                          plot_daily_new_adm,
                          plot_daily_new_deaths,
                          plot_daily_mortal,
                          nrow = 2, ncol = 2,
                          labels = c("(A)", "(B)", "(C)", "(D)"),
                          font.label = list(size = 12, face = "plain")
        ),
        bottom = ggpubr::text_grob(
            paste(c("Black lines: Observed value", 
                    "Blue lines: Smoothed valued (LOESS)",
                    "February 27th, 2020 - October 28th, 2020"),
                  collapse = "\n"),
            hjust = 1, x = 1, size = 10
            )
    )

## Exporing plots (Figure 1)
ggsave("output/main/figures/figure1_temporal_admissions_deaths.pdf", 
       plot_comb_volume
       , width = 10, height = 7, dpi = 800)




# Finished
    
