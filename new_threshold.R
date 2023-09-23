#################################################################
##                       New MVPA Cutoff                       ##
#################################################################

library(tidyverse)
library(lubridate)
library(refund)

# minute level accelerometry data, rows correspond to a subject-day
accel_mat <- read_rds("/datapath/acceldata.rds")

## perform fpca on the log transformed acceleration data
## we will use this for two purposes:
##   1) to impute the (relatively small) number of "missing" data
##   2) calculate surrogate measures based on PCA 
X        <- as.matrix(accel_mat[,paste0("MIN",1:1440)])
fit_fpca <- fpca.face(log(1+X), knots=50)

## impute missing acceleration data among good days using fpca predicted values
## truncate below by 0 since ENMO is bounded below by 0 by definition
inx_na    <- is.na(X)
X[inx_na] <- pmax(0,exp(fit_fpca$Yhat[inx_na]) + 1)

## df for new daily LIPA and MVPA
df_new_cutoff_daily <- accel_mat %>%
  select(eid, date)

df_new_cutoff_daily$MVPA_Mins <- rowSums(X >= 193)
df_new_cutoff_daily$LIPA_Mins <- rowSums(X < 193 & X >= 30)

## df for new avg LIPA and MVPA
df_new_cutoff <- df_new_cutoff_daily %>%
  group_by(eid) %>%
  summarise(LIPA_Mins_mean = mean(LIPA_Mins),
            MVPA_Mins_mean = mean(MVPA_Mins)) %>%
  mutate(eid = as.numeric(eid))

# saveRDS(df_new_cutoff, file = "/datapath/newthreshold.rds")

