#################################################################
##                      Author: Qier Meng                      ##
#################################################################

##################################################################
##                 PA among current MS patients                 ##
##################################################################

library(data.table) ## fread
library(tidyverse)
library(lubridate)
library(table1)
library(ggpubr)

## load demographic dataset
df_demog <- fread("/datapath/demogdata.tab", header = TRUE, sep = "\t")
df_demog <- df_demog %>%
  mutate(eid = f.eid) %>%
  dplyr::select(-f.eid)

## load physical activity data
df_accel <- readRDS("/datapath/padata.rds")
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

## drop participants with missing age, sex, or BMI
df_demog_pa_date <- df_demog_pa_date %>%
  drop_na(age_atwear, f.31.0.0, f.21001.0.0) ## 93160  3282

## add disease history variables (ever diagnosed - on or before April 1, 2021)
df_demog_pa_date <- df_demog_pa_date %>%
  mutate(diabetes_1 = if_else(is.na(f.130706.0.0), 0, 1),
         diabetes_2 = if_else(is.na(f.130708.0.0), 0, 1),
         diabetes_3 = if_else(is.na(f.130710.0.0), 0, 1),
         diabetes_4 = if_else(is.na(f.130712.0.0), 0, 1),
         diabetes_5 = if_else(is.na(f.130714.0.0), 0, 1),
         stroke_1 = if_else(is.na(f.131180.0.0), 0, 1),
         stroke_2 = if_else(is.na(f.131366.0.0), 0, 1),
         stroke_3 = if_else(is.na(f.131368.0.0), 0, 1),
         stroke_4 = if_else(is.na(f.131362.0.0), 0, 1),
         stroke_5 = if_else(is.na(f.131360.0.0), 0, 1),
         stroke_6 = if_else(is.na(f.131056.0.0), 0, 1)) %>%
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
         chd = as.factor(if_else(is.na(f.131306.0.0), 0, 1)),
         cancer = as.factor(
           if_else(is.na(f.40005.0.0), 0, 1)),
         townsend = f.189.0.0
  ) ## 93160  3304

## drop participants who had ever (on or before April 1, 2021) been diagnosed with 
## diabetes, stroke, or coronary heart disease

df_casecontrol <- df_demog_pa_date %>% 
  filter(if_all(c(diabetes, stroke, chd), ~ . == 0)) ## 79145  3304

## drop participants diagnosed with MS after accelerometer wearing
df_casecontrol <- df_casecontrol %>%
  mutate(time.diff = difftime(f.131042.0.0, f.90010.0.0)) %>%
  filter(is.na(time.diff) | time.diff <= 0) ## 79108  3305

##################################################################
##                           Matching                           ##
##################################################################

df_ms_cases <- df_casecontrol %>%
  filter(time.diff <= 0) ## 316 3305

df_healthy_pre <- df_casecontrol %>%
  filter(is.na(f.131042.0.0)) ## 78792  3305

# match case control
df_controls <- NULL
df_nomatch <- NULL
n_match <- 30

df_healthy <- df_healthy_pre %>%
  mutate(selected = 0) ## 78792  3306

set.seed(123)
for (i in 1:nrow(df_ms_cases)) {
  temp <- NULL
  temp <- df_healthy %>%
    filter(selected == 0) %>%
    # sex
    filter(f.31.0.0 == df_ms_cases[i,]$f.31.0.0) %>%
    # age at wear
    filter(age_atwear >= df_ms_cases[i,]$age_atwear - 3 &
             age_atwear <= df_ms_cases[i,]$age_atwear + 3) %>%
    # bmi
    filter(f.21001.0.0 >= df_ms_cases[i,]$f.21001.0.0 - 3 &
             f.21001.0.0 <= df_ms_cases[i,]$f.21001.0.0 + 3)
  if (nrow(temp) >= n_match) {
    sampled_controls <- sample_n(temp, n_match)
    df_controls <- bind_rows(df_controls, sampled_controls)
    df_healthy <- df_healthy %>%
      mutate(selected = if_else(eid %in% sampled_controls$eid, 1, selected))
  } else {
    df_nomatch <- bind_rows(df_nomatch, df_ms_cases[i,])
  }
}

# add group variable
df_ms_cases <- df_ms_cases %>%
  mutate(group = "MS") ## 316 3306

df_controls <- df_controls %>%
  mutate(group = "Controls") ## 9480 3307

#################################################################
##                           Table 1                           ##
#################################################################

df_table1 <- bind_rows(df_ms_cases, df_controls)

label(df_table1$age_atwear) <- "Age at accelerometer wearing"
label(df_table1$sex) <- "Sex"
label(df_table1$bmi) <- "BMI"
label(df_table1$white) <- "Race (white)"
label(df_table1$smoking) <- "Smoking status"
label(df_table1$drinking) <- "Frequency of drinking alcohol"
lable(df_table1$cancer) <- "Cancer"
label(df_table1$ill_disability) <- "Long-standing illness, disability or infirmity"
label(df_table1$townsend) <- "Townsend deprivation index at recruitment"

table1(~ age_atwear + sex + bmi + 
         white + smoking + drinking +
         ill_disability + cancer + townsend | group, data = df_table1)

#################################################################
##                           Results                           ##
#################################################################

casecontrol_result <- NULL
acc_vars <- c("TAC_mean", "TLAC_mean", "Sed_Mins_mean",
              "LIPA_Mins_mean", "MVPA_Mins_mean", "ASTP_mean",
              "SATP_mean", "relA_mean")

for (i in acc_vars) {
  test_res <- t.test(pull(df_controls, i),
                     pull(df_ms_cases, i), 
                     alternative = "two.sided")
  casecontrol_result <- casecontrol_result %>%
    bind_rows(data.frame(Variable = i,
                         Control = paste0(round(test_res$estimate[1],2), " (",
                                          round(sd(pull(df_controls, i)),2), ")"),
                         MS = paste0(round(test_res$estimate[2],2), " (",
                                     round(sd(pull(df_ms_cases, i)),2), ")"),
                         `p-value` = if_else(test_res$p.value < 0.001, 
                                             "<0.001", 
                                             as.character(round(test_res$p.value, 3)))
    )
    )
}

## race (%non-white)
race_test <- prop.test(x = c(sum(!(df_controls$f.21000.0.0 %in% c(1, 1001, 1002, 1003, -1, -3, NA))), 
                             sum(!(df_ms_cases$f.21000.0.0 %in% c(1, 1001, 1002, 1003, -1, -3, NA)))), 
                       n = c(nrow(df_controls) - sum(df_controls$f.21000.0.0 %in% c(-1, -3, NA)), 
                             nrow(df_ms_cases) - sum(df_controls$f.21000.0.0 %in% c(-1, -3, NA))), 
                       alternative = "two.sided",
                       correct = FALSE)
casecontrol_result <- casecontrol_result %>%
  bind_rows(data.frame(Variable = "Race",
                       Control = paste0(sum(!(df_controls$f.21000.0.0 %in% c(1, 1001, 1002, 1003, -1, -3, NA))), " (",
                                        round(race_test$estimate[1]*100,2), "%)"),
                       MS = paste0(sum(!(df_ms_cases$f.21000.0.0 %in% c(1, 1001, 1002, 1003, -1, -3, NA))), " (",
                                   round(race_test$estimate[2]*100,2), "%)"),
                       `p-value` = if_else(race_test$p.value < 0.001, 
                                           "<0.001", 
                                           as.character(round(race_test$p.value, 3))
                       )
  )
  )

## smoking (previous or current))
smoke_test <- prop.test(x = c(sum(df_controls$f.20116.0.0 %in% c(1, 2)), 
                              sum(df_ms_cases$f.20116.0.0 %in% c(1, 2))), 
                        n = c(nrow(df_controls) - sum(df_controls$f.20116.0.0 %in% c(-3, NA)), 
                              nrow(df_ms_cases) - sum(df_ms_cases$f.20116.0.0 %in% c(-3, NA))), 
                        alternative = "two.sided",
                        correct = FALSE)
casecontrol_result <- casecontrol_result %>%
  bind_rows(data.frame(Variable = "Smoking",
                       Control = paste0(sum(df_controls$f.20116.0.0 %in% c(1, 2)), " (",
                                        round(smoke_test$estimate[1]*100,2), "%)"),
                       MS = paste0(sum(df_ms_cases$f.20116.0.0 %in% c(1, 2)), " (",
                                   round(smoke_test$estimate[2]*100,2), "%)"),
                       `p-value` = if_else(smoke_test$p.value < 0.001, 
                                           "<0.001", 
                                           as.character(round(smoke_test$p.value, 3))
                       )
  )
  )

## drinking (% at least once a week)
drink_test <- prop.test(x = c(sum(df_controls$f.1558.0.0 %in% c(1, 2, 3)), 
                              sum(df_ms_cases$f.1558.0.0 %in% c(1, 2, 3))), 
                        n = c(nrow(df_controls) - sum(df_controls$f.1558.0.0 %in% c(-3, NA)), 
                              nrow(df_ms_cases) - sum(df_ms_cases$f.1558.0.0 %in% c(-3, NA))), 
                        alternative = "two.sided",
                        correct = FALSE)
casecontrol_result <- casecontrol_result %>%
  bind_rows(data.frame(Variable = "Drinking",
                       Control = paste0(sum(df_controls$f.1558.0.0 %in% c(1, 2, 3)), " (",
                                        round(drink_test$estimate[1]*100,2), "%)"),
                       MS = paste0(sum(df_ms_cases$f.1558.0.0 %in% c(1, 2, 3)), " (",
                                   round(drink_test$estimate[2]*100,2), "%)"),
                       `p-value` = if_else(drink_test$p.value < 0.001, 
                                           "<0.001", 
                                           as.character(round(drink_test$p.value, 3))
                       )
  )
  )

casecontrol_result

# proportions of no MVPA in each group
# ms
sum(df_ms_cases$MVPA_Mins_mean == 0)/length(df_ms_cases$MVPA_Mins_mean)
sum(df_controls$MVPA_Mins_mean == 0)/length(df_controls$MVPA_Mins_mean)

#################################################################
##                        Density Plots                        ##
#################################################################

## TAC
dens_tac <- ggplot(rbind(select(df_ms_cases, group, TAC_mean),
                         select(df_controls, group, TAC_mean)),
                   aes(x = TAC_mean, color = group, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(x = "TA", y = "Density") +
  theme(text=element_text(size=10), legend.title = element_blank())

## TLAC
dens_tlac <- ggplot(rbind(select(df_ms_cases, group, TLAC_mean),
                          select(df_controls, group, TLAC_mean)),
                    aes(x = TLAC_mean, color = group, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(x = "TLA", y = "Density") +
  theme(text=element_text(size=10), legend.title = element_blank())

## Sedentary Minutes
dens_sedm<- ggplot(rbind(select(df_ms_cases, group, Sed_Mins_mean),
                         select(df_controls, group, Sed_Mins_mean)),
                   aes(x = Sed_Mins_mean, color = group, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(x = "ST", y = "Density") +
  theme(text=element_text(size=10), legend.title = element_blank())

## Relative Amplitude
dens_ra<- ggplot(rbind(select(df_ms_cases, group, relA_mean),
                       select(df_controls, group, relA_mean)),
                 aes(x = relA_mean, color = group, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(x = "RA", y = "Density") +
  theme(text=element_text(size=10), legend.title = element_blank())

## LIPA
dens_lipa<- ggplot(rbind(select(df_ms_cases, group, LIPA_Mins_mean),
                         select(df_controls, group, LIPA_Mins_mean)),
                   aes(x = LIPA_Mins_mean, color = group, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(x = "LIPA", y = "Density") +
  theme(text=element_text(size=10), legend.title = element_blank())

## MVPA
dens_mvpa<- ggplot(rbind(select(df_ms_cases, group, MVPA_Mins_mean),
                         select(df_controls, group, MVPA_Mins_mean)),
                   aes(x = MVPA_Mins_mean, color = group, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(x = "MVPA", y = "Density") +
  theme(text=element_text(size=10), legend.title = element_blank())

## ASTP
dens_astp<- ggplot(rbind(select(df_ms_cases, group, ASTP_mean),
                         select(df_controls, group, ASTP_mean)),
                   aes(x = ASTP_mean, color = group, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(x = "ASTP", y = "Density") +
  theme(text=element_text(size=10), legend.title = element_blank())

## SATP
dens_satp<- ggplot(rbind(select(df_ms_cases, group, SATP_mean),
                         select(df_controls, group, SATP_mean)),
                   aes(x = SATP_mean, color = group, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(x = "SATP", y = "Density") +
  theme(text=element_text(size=10), legend.title = element_blank())

ggarrange(dens_tac, dens_tlac, dens_sedm, dens_ra,
          dens_lipa, dens_mvpa, dens_astp, dens_satp,
          labels = "AUTO",
          ncol=2, nrow=4, 
          common.legend = TRUE, legend="bottom")

## disease duration stats for current MS patients
duration_years <- as.numeric(abs(df_ms_cases$time.diff))/365.25

mean(duration_years)
sd(duration_years)
min(duration_years)
max(duration_years)


