#################################################################
##                      Author: Qier Meng                      ##
#################################################################

#################################################################
##          PA as a biomarker for future MS diagnosis          ##
#################################################################

## Disease history is defined as
## having been diagnosed with it on or before device wearing date

library(data.table) ## fread
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(table1)

## load demographic dataset
df_demog <- fread("/datapath/demogdata.tab", header = TRUE, sep = "\t")
df_demog <- df_demog %>%
  mutate(eid = f.eid) %>%
  dplyr::select(-f.eid) 

## load physical activity data
df_accel <- readRDS("/datapath/padata.rds.rds")
df_new_limvpa <- readRDS("/datapath/newthreshold.rds")
df_accel <- df_accel %>%
  select(eid, TAC_mean, TLAC_mean, Sed_Mins_mean, ASTP_mean,
         SATP_mean, relA_mean) %>%
  mutate(eid = as.numeric(eid)) %>%
  inner_join(df_new_limvpa, by = "eid")

## diagnostic date data
df_date <- fread("/datapath/diagnosticdata.tab", header = TRUE, sep = "\t")
df_date <- df_date %>%
  mutate(eid = f.eid) %>%
  dplyr::select(-f.eid)

## join demographic and pa data; 90010 is start wear time
df_demog_pa <- inner_join(df_demog, df_accel, by = "eid") %>%
  mutate(f.90010.0.0 = as.Date(f.90010.0.0)) ## 93370 3267

# join df_demog_pa and ms diagnosis date and major diseases diagnosis dates
df_demog_pa_date <- inner_join(df_demog_pa,
                               select(df_date, eid, f.131042.0.0,
                                      f.130706.0.0, 
                                      f.130708.0.0, 
                                      f.130710.0.0, 
                                      f.130712.0.0, 
                                      f.130714.0.0, 
                                      f.131180.0.0, 
                                      f.131366.0.0, 
                                      f.131368.0.0, 
                                      f.131362.0.0, 
                                      f.131360.0.0, 
                                      f.131056.0.0, 
                                      f.131306.0.0), 
                               by = "eid") ## 93363 3280 7 subjects dropped

## add age at wear (use years to calculate)
df_demog_pa_date <- df_demog_pa_date %>%
  mutate(mob = make_date(year = f.34.0.0, month = f.52.0.0),
         age_atwear = year(f.90010.0.0) - year(mob)) ## 93363 3282

## df for survival analysis
df_pa_survival <- df_demog_pa_date %>%
  mutate(date.data = as.IDate("2021-04-01")) %>%
  mutate(status = ifelse(is.na(f.131042.0.0), 0, 1)) %>%
  mutate(date.outcome = if_else(is.na(f.131042.0.0), 
                                if_else(is.na(f.40000.0.0), date.data, f.40000.0.0),
                                f.131042.0.0)) %>%
  mutate(time.diff = difftime(date.outcome, f.90010.0.0)) %>%
  filter(time.diff >= 0) ## 92994  3286

## add disease history variables
df_pa_survival <- df_pa_survival %>%
  mutate(diabetes_1 = if_else(is.na(f.130706.0.0), 0,
                              if_else(f.130706.0.0 > f.90010.0.0, 0, 1)),
         diabetes_2 = if_else(is.na(f.130708.0.0), 0,
                              if_else(f.130708.0.0 > f.90010.0.0, 0, 1)),
         diabetes_3 = if_else(is.na(f.130710.0.0), 0,
                              if_else(f.130710.0.0 > f.90010.0.0, 0, 1)),
         diabetes_4 = if_else(is.na(f.130712.0.0), 0,
                              if_else(f.130712.0.0 > f.90010.0.0, 0, 1)),
         diabetes_5 = if_else(is.na(f.130714.0.0), 0,
                              if_else(f.130714.0.0 > f.90010.0.0, 0, 1)),
         stroke_1 = if_else(is.na(f.131180.0.0), 0,
                            if_else(f.131180.0.0 > f.90010.0.0, 0, 1)),
         stroke_2 = if_else(is.na(f.131366.0.0), 0,
                            if_else(f.131366.0.0 > f.90010.0.0, 0, 1)),
         stroke_3 = if_else(is.na(f.131368.0.0), 0,
                            if_else(f.131368.0.0 > f.90010.0.0, 0, 1)),
         stroke_4 = if_else(is.na(f.131362.0.0), 0,
                            if_else(f.131362.0.0 > f.90010.0.0, 0, 1)),
         stroke_5 = if_else(is.na(f.131360.0.0), 0,
                            if_else(f.131360.0.0 > f.90010.0.0, 0, 1)),
         stroke_6 = if_else(is.na(f.131056.0.0), 0,
                            if_else(f.131056.0.0 > f.90010.0.0, 0, 1))) %>%
  # clean up variable names and clean up disease variables
  mutate(sex = if_else(f.31.0.0 == 1, "Male", "Female"),
         bmi = f.21001.0.0,
         white = as.factor(case_when(
           f.21000.0.0 %in% c(1, 1001, 1002, 1003) ~ "1",
           f.21000.0.0 %in% c(-1, -3, NA) ~ as.character(f.21000.0.0),
           TRUE ~ "0")),
         smoking = as.factor(f.20116.0.0),
         drinking = as.factor(f.1558.0.0),
         ill_disability = as.factor(f.2188.0.0),
         diabetes = as.factor(if_else(
           if_any(contains("diabetes"), ~ . == 1), 1, 0)),
         stroke = as.factor(if_else(
           if_any(contains("stroke"), ~ . == 1), 1, 0)),
         chd = as.factor(
           if_else(is.na(f.131306.0.0), 0, 
                   if_else(f.131306.0.0 > f.90010.0.0, 0, 1))),
         cancer = as.factor(
           if_else(is.na(f.40005.0.0), 0, 
                   if_else(f.40005.0.0 > f.90010.0.0, 0, 1))),
         townsend = f.189.0.0
  ) ## 92994  3308

#################################################################
##                  Cleaned data for analyses                  ##
#################################################################

df_survival <- df_pa_survival %>%
  mutate(white_binary = as.factor(case_when(
    f.21000.0.0 %in% c(1, 1001, 1002, 1003) ~ "1",
    f.21000.0.0 %in% c(-1, -3, NA) ~ NA_character_,
    TRUE ~ "0")),
    smoke_binary = as.factor(case_when(
      f.20116.0.0 %in% c(1, 2) ~ "1",
      f.20116.0.0 == 0 ~ "0",
      TRUE ~ NA_character_)),
    drinking_binary = as.factor(case_when(
      f.1558.0.0 %in% c(1, 2, 3) ~ "1",
      f.1558.0.0 %in% c(4, 5, 6) ~ "0",
      TRUE ~ NA_character_))) %>%
  drop_na(age_atwear, sex, bmi, white_binary, 
          smoke_binary, drinking_binary,
          diabetes, stroke, chd, cancer, townsend) ## 92169  3311

# calculate person-year
surv_time <- df_survival$time.diff %>% as.numeric()
sum(surv_time/365.25)
# incidence rate
47/sum(surv_time/365.25)

#################################################################
##                           Table 1                           ##
#################################################################

## before NA removal columns
label(df_pa_survival$age_atwear) <- "Age at accelerometer wearing"
label(df_pa_survival$sex) <- "Sex"
label(df_pa_survival$bmi) <- "BMI"
label(df_pa_survival$white) <- "Race"
label(df_pa_survival$smoking) <- "Smoking status"
label(df_pa_survival$drinking) <- "Frequency of drinking alcohol"
label(df_pa_survival$ill_disability) <- "Long-standing illness, disability or infirmity"
label(df_pa_survival$diabetes) <- "Diabetes"
label(df_pa_survival$stroke) <- "Stroke"
label(df_pa_survival$chd) <- "Coronary heart disease"
label(df_pa_survival$cancer) <- "Cancer"
label(df_pa_survival$townsend) <- "Townsend deprivation index at recruitment"

table1(~ age_atwear + sex + bmi + 
         white + smoking + drinking + ill_disability + 
         diabetes + stroke + chd + cancer + townsend | status, 
       data = df_pa_survival)

## after NA removal columns
label(df_survival$age_atwear) <- "Age at accelerometer wearing"
label(df_survival$sex) <- "Sex"
label(df_survival$bmi) <- "BMI"
label(df_survival$white) <- "Race"
label(df_survival$smoking) <- "Smoking status"
label(df_survival$drinking) <- "Frequency of drinking alcohol"
label(df_survival$ill_disability) <- "Long-standing illness, disability or infirmity"
label(df_survival$diabetes) <- "Diabetes"
label(df_survival$stroke) <- "Stroke"
label(df_survival$chd) <- "Coronary heart disease"
label(df_survival$cancer) <- "Cancer"
label(df_survival$townsend) <- "Townsend deprivation index at recruitment"

table1(~ age_atwear + sex + bmi + 
         white + smoking + drinking + ill_disability + 
         diabetes + stroke + chd + cancer + townsend | status, 
       data = df_survival)

#################################################################
##                          Flowchart                          ##
#################################################################

# events
df_survival %>% filter(status == 1) %>% dim() ## 47 3296
# censored
df_survival %>% filter(status == 0) %>% dim() ## 92122  3296
# deaths
df_survival %>% filter(!is.na(f.40000.0.0) & is.na(f.131042.0.0)) %>% dim() ## 927 3296
# no diagnosis until April 1, 2021
df_survival %>% filter(is.na(f.40000.0.0) & is.na(f.131042.0.0)) %>% dim() ## 91195  3296

#################################################################
##                     Kaplan-Meier curves                     ##
#################################################################

# sex
km_sex <- survfit(Surv(time.diff, status) ~ sex, 
                  data = df_survival)
km_sex_plot <- ggsurvplot(km_sex, data = df_survival, ylim = c(0.999, 1),
                          legend.labs=c("Female", "Male"),
                          title = "Sex",
                          font.main = 14,
                          font.x = 10,
                          font.y = 10,
                          font.tickslab = 8)
km_sex_ggplot <- km_sex_plot$plot +
  xlab("Time (Days)") + 
  theme(legend.position = c(0.1, 0.1), 
        legend.justification = c(0.1, 0.1))

# race
km_race <- survfit(Surv(time.diff, status) ~ white_binary, 
                   data = df_survival)
km_race_plot <- ggsurvplot(km_race, data = df_survival, ylim = c(0.999, 1),
                           legend.labs=c("Non-white", "White"),
                           title = "Race",
                           font.main = 14,
                           font.x = 10,
                           font.y = 10,
                           font.tickslab = 8) 
km_race_ggplot <- km_race_plot$plot +
  xlab("Time (Days)") + 
  theme(legend.position = c(0.1, 0.1), 
        legend.justification = c(0.1, 0.1))


# smoking
km_smoke <- survfit(Surv(time.diff, status) ~ smoke_binary, 
                    data = df_survival)
km_smoke_plot <- ggsurvplot(km_smoke, data = df_survival, ylim = c(0.999, 1),
                            legend.labs=c("Not previous or current", "Previous or current"),
                            title = "Smoking Status",
                            font.main = 14,
                            font.x = 10,
                            font.y = 10,
                            font.tickslab = 8)
km_smoke_ggplot <- km_smoke_plot$plot +
  xlab("Time (Days)") + 
  theme(legend.position = c(0.1, 0.1), 
        legend.justification = c(0.1, 0.1))

# drinking
km_drink <- survfit(Surv(time.diff, status) ~ drinking_binary, 
                    data = df_survival)
km_drink_plot <- ggsurvplot(km_drink, data = df_survival, ylim = c(0.999, 1),
                            legend.labs=c("Less than once a week", "At least once a week"),
                            title = "Drinking Frequency",
                            font.main = 14,
                            font.x = 10,
                            font.y = 10,
                            font.tickslab = 8)
km_drink_ggplot <- km_drink_plot$plot +
  xlab("Time (Days)") + 
  theme(legend.position = c(0.1, 0.1), 
        legend.justification = c(0.1, 0.1))

# demographics
ggarrange(km_sex_ggplot, km_race_ggplot, km_smoke_ggplot, km_drink_ggplot,
          labels = "AUTO",
          ncol=2, nrow=2)

# diabetes
km_diabetes <- survfit(Surv(time.diff, status) ~ diabetes, 
                       data = df_survival)
km_diabetes_plot <- ggsurvplot(km_diabetes, data = df_survival, ylim = c(0.9975, 1),
                               legend.labs=c("No", "Yes"),
                               title = "Diabetes",
                               font.main = 14,
                               font.x = 10,
                               font.y = 10,
                               font.tickslab = 8)
km_diabetes_ggplot <- km_diabetes_plot$plot +
  xlab("Time (Days)") + 
  theme(legend.position = c(0.1, 0.1), 
        legend.justification = c(0.1, 0.1))

# stroke (different y scale)
km_stroke <- survfit(Surv(time.diff, status) ~ stroke, 
                     data = df_survival)
km_stroke_plot <- ggsurvplot(km_stroke, data = df_survival, ylim = c(0.9975, 1),
                             legend.labs=c("No", "Yes"),
                             title = "Stroke",
                             font.main = 14,
                             font.x = 10,
                             font.y = 10,
                             font.tickslab = 8)
km_stroke_ggplot <- km_stroke_plot$plot +
  xlab("Time (Days)") + 
  theme(legend.position = c(0.1, 0.1), 
        legend.justification = c(0.1, 0.1))

# chd
km_chd <- survfit(Surv(time.diff, status) ~ chd, 
                  data = df_survival)
km_chd_plot <- ggsurvplot(km_chd, data = df_survival, ylim = c(0.9975, 1),
                          legend.labs=c("No", "Yes"),
                          title = "Coronary Heart Disease",
                          font.main = 14,
                          font.x = 10,
                          font.y = 10,
                          font.tickslab = 8)
km_chd_ggplot <- km_chd_plot$plot +
  xlab("Time (Days)") + 
  theme(legend.position = c(0.1, 0.1), 
        legend.justification = c(0.1, 0.1))

# cancer
km_cancer <- survfit(Surv(time.diff, status) ~ cancer, 
                     data = df_survival)
km_cancer_plot <- ggsurvplot(km_cancer, data = df_survival, ylim = c(0.9975, 1),
                             legend.labs=c("No", "Yes"),
                             title = "Cancer",
                             font.main = 14,
                             font.x = 10,
                             font.y = 10,
                             font.tickslab = 8)
km_cancer_ggplot <- km_cancer_plot$plot +
  xlab("Time (Days)") + 
  theme(legend.position = c(0.1, 0.1), 
        legend.justification = c(0.1, 0.1))

# disease history
ggarrange(km_diabetes_ggplot, km_stroke_ggplot, km_chd_ggplot, km_cancer_ggplot,
          labels = "AUTO",
          ncol=2, nrow=2)


#################################################################
##                   Single predictor models                   ##
#################################################################

single_vars <- c("age_atwear", "sex", "bmi", "white_binary", "smoke_binary", 
                 "drinking_binary", "diabetes", "stroke", "chd",
                 "cancer", "townsend", 
                 "TAC_mean", "TLAC_mean", "Sed_Mins_mean",
                 "LIPA_Mins_mean", "MVPA_Mins_mean", "ASTP_mean",
                 "SATP_mean", "relA_mean")
single_res <- NULL

for (i in single_vars) {
  formula <- as.formula(paste0("Surv(time.diff, status) ~ ", i))
  coxreg <- coxph(formula, data = df_survival)
  single_res <- bind_rows(single_res,
                          bind_cols(summary(coxreg)$conf.int,
                                    t(as.matrix(summary(coxreg)$concordance))
                          ) %>%
                            mutate(p_value = summary(coxreg)$coef[1, 5],
                                   variable = rownames(summary(coxreg)$conf.int))
  )
}

single_res_clean <- single_res %>% 
  select(variable, `exp(coef)`, `lower .95`, `upper .95`,
         p_value, C, `se(C)`) %>%
  arrange(desc(C))

single_res_clean


#################################################################
##                  Grouped forward selection                  ##
#################################################################

# characteristic predictors
char_vars <- c("age_atwear", "sex", "bmi", "white_binary", "smoke_binary", 
               "drinking_binary", "diabetes", "stroke", "chd",
               "cancer", "townsend")

# PA predictors
pa_vars <- c("TAC_mean", "TLAC_mean", "Sed_Mins_mean",
             "LIPA_Mins_mean", "MVPA_Mins_mean", "ASTP_mean",
             "SATP_mean", "relA_mean")

# forward selection function
forward_selection <- function(variable_list, stoprule, include_vars = NULL) {
  # initialization
  selected_vars <- NULL
  c_track <- NULL
  delta_c_track <- NULL
  pre_c <- if(length(include_vars) == 0) {0} else{
    base_model <- coxph(as.formula(paste0("Surv(time.diff, status) ~ ", 
                                          paste0(include_vars, collapse = " + "))), 
                        data = df_survival)
    summary(base_model)$concordance[1]
  }
  
  repeat {
    vars <- setdiff(variable_list, selected_vars)
    
    if (length(vars) == 0) { break }
    
    # track the models
    track <- NULL
    for (i in vars) {
      formula <- as.formula(paste0("Surv(time.diff, status) ~ ", 
                                   paste0(c(include_vars, selected_vars, i), collapse = " + ")))
      coxreg <- coxph(formula, data = df_survival)
      track <- bind_rows(track, data.frame(variable = i, 
                                           c = summary(coxreg)$concordance[1]))
    }
    
    delta_c <- max(track$c) - pre_c
    
    if (delta_c < stoprule) { break }
    
    selected_vars <- c(selected_vars, track[which.max(track$c), 1])
    c_track <- c(c_track, max(track$c))
    delta_c_track <- c(delta_c_track, delta_c)
    
    pre_c <- max(track$c)
  }
  
  model <- coxph(as.formula(paste0("Surv(time.diff, status) ~ ", 
                                   paste0(c(include_vars, selected_vars), collapse = " + "))), 
                 data = df_survival)
  df_selection <- data.frame("Selected variable" = selected_vars, 
                             "Cumulative C" = c_track,
                             "Delta C" = delta_c_track)
  df_model <- bind_cols(data.frame(variable = rownames(summary(model)$conf.int),
                                   p_value = summary(model)$coef[, 5]),
                        summary(model)$conf.int[, c(1, 3, 4)])
  rownames(df_model) <- NULL
  return(list("Selection Process" = df_selection, 
              "Model Summary" = df_model))
}

# characteristic predictors forward selection
fs_char <- forward_selection(char_vars, 0.01)

char_selected <- c("age_atwear", "stroke", "townsend")

# PA predictors forward selection
fs_pa <- forward_selection(pa_vars, 0.01, char_selected)

# selected model
final_vars <- c(fs_char[[1]] %>% pull(Selected.variable),
                fs_pa[[1]] %>% pull(Selected.variable))

selected_model <- coxph(as.formula(paste0("Surv(time.diff, status) ~ ", 
                                       paste0(final_vars, collapse = " + "))), 
                     data = df_survival)
summary(selected_model)

# final three-predictor model
final_model <- coxph(Surv(time.diff, status) ~ age_atwear + stroke + relA_mean, 
                     data = df_survival)
summary(final_model)

# time gap between initial assessment and device wearing
timegap <- difftime(df_survival$f.90010.0.0, df_survival$f.53.0.0) %>% as.numeric()
summary(timegap)
sd(timegap)

# time to event summary stats
df_survival %>% 
  select(time.diff, status) %>%
  mutate(time.diff.year = as.numeric(time.diff)/365.25) %>%
  group_by(status) %>%
  summarise(mean_year = mean(time.diff.year), sd_year = sd(time.diff.year),
            n = n())

