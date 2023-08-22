#THERE IS NO USE OF COMORBIDITIES TO DEVELOP THE MODEL SO TO LESSEN THE COMPUTATIONAL LOAD, COMORBIDITIES WERE NOT CONSIDERED

#LOADING THE REQUIRED LIBRARIES
.libPaths(c("/mnt/bmh01-rds/mrc-multi-outcome/R\ packages/"))
library(dplyr)
library(mstate)

#READING THE DATA
df <- readRDS("multimorbidity_data.rds")

#FILTERING THE DATA
mmd <- select(filter(df, CVD_hist == "0", Diab_t2_hist == "0", CKD_hist == "0"), 
              c("person_id", "care_site_id", "study_dtindex", "study_dtindex_r", "Age", 
                "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD", 
                "Death_t", "Death_c", "Death_NelsonAalen_link", "CKD_hist", 
                "CKD_hist_t", "CKD_ev_c", "CKD_ev_t", "Diab_t2_hist", "Diab_t2_hist_t", 
                "Diab_t2_ev_c", "Diab_t2_ev_t", "dtcens_combdeath_r", "CVD_ev_t", 
                "CVD_ev_c", "CVD_hist", "dtcens_var", "o2", "c2", "o3", "c3", "o4", 
                "c4", "o5", "c5", "o6", "c6", "o7", "c7", "o8", "c8", "o9", "c9"))

#CREATING TRAIN AND TEST DATASETS
set.seed(1)  #to make reproducible
sample <- sample(c(TRUE, FALSE), nrow(mmd), replace=TRUE, prob=c(0.7,0.3))
train <- mmd[sample, ]
test <- mmd[!sample, ]

#TRANSITION MATRIX
tmat <- transMat(x = list(c(2,3,4,9), c(5,6,9), c(5,7,9), c(6,7,9), c(8,9), c(8,9), c(8,9), c(9), c()), 
                 names = c("Healthy", "CVD", "T2D", "CKD", "CVD+T2D", "CVD+CKD", "CKD+T2D", "CVD+T2D+CKD", "Death"))
tmat
mstrain <- msprep(data = train, trans = tmat, time = c(NA, "o2", "o3", "o4", "o5", "o6", "o7", "o8", "o9"), 
                  status = c(NA, "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9"), 
                  keep = c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD"))
events(mstrain)

#CREATION OF TRANSITION-SPECIFIC COVARIATES
covs <- c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD")
mstrain <- expand.covs(mstrain, covs, longnames = FALSE)
mstrain[, c("Tstart", "Tstop", "time")] <- mstrain[, c("Tstart", "Tstop", "time")]/365.25 #to convert days into years

#MODEL BUILDUP
cfull <- coxph(Surv(Tstart, Tstop, status) ~ Age.1 + Age.2 + Age.3 + Age.4 + Age.5 + Age.6 + Age.7 + Age.8 +
                 Age.9 + Age.10 + Age.11 + Age.12 + Age.13 + Age.14 + Age.15 + Age.16 + Age.17 + Age.18 + Age.19 + 
                 Age.20 + gender.1  + gender.2 + gender.3 + gender.4 + gender.5 + gender.6 + gender.7 + gender.8 + 
                 gender.9 + gender.10 + gender.11 + gender.12 + gender.13 + gender.14 + gender.15 + gender.16 + 
                 gender.17 + gender.18 + gender.19 + gender.20 + BMI.1 + BMI.2 + BMI.3 + BMI.4 + BMI.5 + BMI.6 +
                 BMI.7 + BMI.8 + BMI.9 + BMI.10 + BMI.11 + BMI.12 + BMI.13 + BMI.14 + BMI.15 + BMI.16 + BMI.17 + 
                 BMI.18 + BMI.19 + BMI.20 + Cholhdl_ratio.1 + Cholhdl_ratio.2 + Cholhdl_ratio.3 + Cholhdl_ratio.4 +
                 Cholhdl_ratio.5 + Cholhdl_ratio.6 + Cholhdl_ratio.7 + Cholhdl_ratio.8 + Cholhdl_ratio.9 + 
                 Cholhdl_ratio.10 + Cholhdl_ratio.11 + Cholhdl_ratio.12 + Cholhdl_ratio.13 + Cholhdl_ratio.14 + 
                 Cholhdl_ratio.15 + Cholhdl_ratio.16 + Cholhdl_ratio.17 + Cholhdl_ratio.18 + Cholhdl_ratio.19 + 
                 Cholhdl_ratio.20 + Ethnicity61.1 + Ethnicity61.2 + Ethnicity61.3 + Ethnicity61.4 + Ethnicity61.5 +
                 Ethnicity61.6 +Ethnicity61.7 + Ethnicity61.8 + Ethnicity61.9 +Ethnicity61.10 + Ethnicity61.11 +
                 Ethnicity61.12 + Ethnicity61.13 + Ethnicity61.14 + Ethnicity61.15 + Ethnicity61.16 + Ethnicity61.17 +
                 Ethnicity61.18 + Ethnicity61.19 + Ethnicity61.20 + Ethnicity62.1 + Ethnicity62.2 + Ethnicity62.3 +
                 Ethnicity62.4 + Ethnicity62.5 + Ethnicity62.6 + Ethnicity62.7 + Ethnicity62.8 + Ethnicity62.9 +
                 Ethnicity62.10 + Ethnicity62.11 + Ethnicity62.12 + Ethnicity62.13 + Ethnicity62.14 + Ethnicity62.15 +
                 Ethnicity62.16 + Ethnicity62.17 + Ethnicity62.18 + Ethnicity62.19 + Ethnicity62.20 + Ethnicity63.1 +
                 Ethnicity63.2 + Ethnicity63.3 + Ethnicity63.4 + Ethnicity63.5 + Ethnicity63.6 + Ethnicity63.7 +
                 Ethnicity63.8 + Ethnicity63.9 + Ethnicity63.10 + Ethnicity63.11 + Ethnicity63.12 + Ethnicity63.13 +
                 Ethnicity63.14 + Ethnicity63.15 + Ethnicity63.16 + Ethnicity63.17 + Ethnicity63.18 + Ethnicity63.19 +
                 Ethnicity63.20 + Ethnicity64.1 + Ethnicity64.2 + Ethnicity64.3 + Ethnicity64.4 + Ethnicity64.5 +
                 Ethnicity64.6 + Ethnicity64.7 + Ethnicity64.8 + Ethnicity64.9 + Ethnicity64.10 + Ethnicity64.11 +
                 Ethnicity64.12 + Ethnicity64.13 + Ethnicity64.14 + Ethnicity64.15 + Ethnicity64.16 + Ethnicity64.17 +
                 Ethnicity64.18 + Ethnicity64.19 + Ethnicity64.20 + SBP.1 + SBP.2 + SBP.3 + SBP.4 + SBP.5 + SBP.6 +
                 SBP.7 + SBP.8 + SBP.9 + SBP.10 + SBP.11 + SBP.12 + SBP.13 + SBP.14 + SBP.15 + SBP.16 + SBP.17 +
                 SBP.18 + SBP.19 + SBP.20 + Smoking1.1 + Smoking1.2 + Smoking1.3 + Smoking1.4 + Smoking1.5 + 
                 Smoking1.6 + Smoking1.7 + Smoking1.8 + Smoking1.9 + Smoking1.10 + Smoking1.11 + Smoking1.12 + 
                 Smoking1.13 + Smoking1.14 + Smoking1.15 + Smoking1.16 + Smoking1.17 + Smoking1.18 + Smoking1.19 +
                 Smoking1.20 + Smoking2.1 + Smoking2.2 + Smoking2.3 + Smoking2.4 + Smoking2.5 + Smoking2.6 + 
                 Smoking2.7 + Smoking2.8 + Smoking2.9 + Smoking2.10 + Smoking2.11 + Smoking2.12 + Smoking2.13 + 
                 Smoking2.14 + Smoking2.15 + Smoking2.16 + Smoking2.17 + Smoking2.18 + Smoking2.19 + Smoking2.20 + 
                 IMD1.1 + IMD1.2 + IMD1.3 + IMD1.4 + IMD1.5 + IMD1.6 + IMD1.7 + IMD1.8 + IMD1.9 + IMD1.10 + IMD1.11 +
                 IMD1.12 + IMD1.13 + IMD1.14 + IMD1.15 + IMD1.16 + IMD1.17 + IMD1.18 + IMD1.19 + IMD1.20 + IMD2.1 +
                 IMD2.2 + IMD2.3 + IMD2.4 + IMD2.5 + IMD2.6 + IMD2.7 + IMD2.8 + IMD2.9 + IMD2.10 + IMD2.11 + IMD2.12 +
                 IMD2.13 + IMD2.14 + IMD2.15 + IMD2.16 + IMD2.17 + IMD2.18 + IMD2.19 + IMD2.20 + IMD3.1 + IMD3.2 +
                 IMD3.3 + IMD3.4 + IMD3.5 + IMD3.6 + IMD3.7 + IMD3.8 + IMD3.9 + IMD3.10 + IMD3.11 + IMD3.12 + IMD3.13 +
                 IMD3.14 + IMD3.15 + IMD3.16 + IMD3.17 + IMD3.18 + IMD3.19 + IMD3.20 + IMD4.1 + IMD4.2 + IMD4.3 + 
                 IMD4.4 + IMD4.5 + IMD4.6 + IMD4.7 + IMD4.8 + IMD4.9 + IMD4.10 + IMD4.11 + IMD4.12 + IMD4.13 + IMD4.14 +
                 IMD4.15 + IMD4.16 + IMD4.17 + IMD4.18 + IMD4.19 + IMD4.20 + strata(trans), data = mstrain, method = "breslow")
summary(cfull)
#AN ATTEMPT WAS MADE TO BUILD A MODEL WITH ALL THE DUMMY VARIABLES BUT DUE TO COMPUTATIONAL RESTRAINT IT WAS NOT FEASIBLE

#REDUCED RANK MODEL
rr1 <- redrank(Surv(Tstart, Tstop, status) ~ Age + gender + BMI + Cholhdl_ratio + Ethnicity6 + SBP + Smoking + IMD, 
               data = mstrain, R = 1, print.level = 0)
rr1$Alpha
rr1$Gamma
rr1$Beta
