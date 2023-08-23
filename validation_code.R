#THERE WAS NO USE OF COMORBIDITIES TO VALIDATE THE MODEL SO TO LESSEN THE COMPUTATIONAL LOAD, COMORBIDITIES WERE NOT CONSIDERED

#THE REQUIRED LIBRARIES
.libPaths(c("/mnt/bmh01-rds/mrc-multi-outcome/R\ packages/"))
library(dplyr)
library(mstate)
library(ggplot2)

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

#TRAIN AND TEST DATASETS
set.seed(1)  #to make reproducible
sample <- sample(c(TRUE, FALSE), nrow(mmd), replace = TRUE, prob = c(0.7, 0.3))
train <- mmd[sample, ]
test <- mmd[!sample, ]

#TRAIN DATASET TRANSITION PROBABILITIES
tmat <- transMat(x = list(c(2,3,4,9), c(5,6,9), c(5,7,9), c(6,7,9), c(8,9), c(8,9), c(8,9), c(9), c()), 
                 names = c("Healthy", "CVD", "T2D", "CKD", "CVD+T2D", "CVD+CKD", "CKD+T2D", "CVD+T2D+CKD", "Death"))
mstrain <- msprep(data = train, trans = tmat, time = c(NA, "o2", "o3", "o4", "o5", "o6", "o7", "o8", "o9"), 
                  status = c(NA, "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9"), 
                  keep = c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD"))
covs <- c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD")
mstrain <- expand.covs(mstrain, covs, longnames = FALSE)
mstrain[, c("Tstart", "Tstop", "time")] <- mstrain[, c("Tstart", "Tstop", "time")]/365.25 #to convert days into years
c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = mstrain, method = "breslow")
msf0 <- msfit(object = c0, vartype = "aalen", trans = tmat)
pt0 <- probtrans(msf0, predt = 0, method = "aalen")
ord <- c(2, 3, 4, 5, 6, 7, 8, 9, 1)
png(file = "Stacked transition probabilities.png", width = 1080, height = 1080, type = "cairo")
plot(pt0, ord = ord, xlab = "Duration", las = 1)
dev.off()

#MODEL FITTING INTO TEST DATASET 
msf0v <- msfit(object = c0, newdata = test, vartype = "aalen", trans = tmat)

#STACKED TRANSITIONS PROBABILITIES AFTER 5 YEARS (1826 DAYS)
pt0v5 <- probtrans(msf0v, predt = 1826/365.25, method = "aalen")
png(file = "Stacked transitions probabilities after 5 years.png", width = 1080, height = 1080, type = "cairo")
plot(pt0v5, ord = ord, xlab = "Duration", las = 1)
dev.off()

#STACKED TRANSITIONS PROBABILITIES AFTER 10 YEARS (3652 DAYS)
pt0v10 <- probtrans(msf0v, predt = 3652/365.25, method = "aalen")
png(file = "Stacked transitions probabilities after 10 years.png", width = 1080, height = 1080, type = "cairo")
plot(pt0v10, ord = ord, xlab = "Duration", las = 1)
dev.off()
