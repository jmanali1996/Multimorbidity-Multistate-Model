#THERE IS NO USE OF COMORBIDITIES FOR VISUALIZATIONS SO TO LESSEN THE COMPUTATIONAL LOAD, COMORBIDITIES WERE NOT CONSIDERED

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

#STACKED TRANSITIONS PROBABILITIES
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
png(file = "Stacked transitions probabilities.png", width = 1080, height = 1080, type = "cairo")
plot(pt0, ord = ord, xlab = "Duration", las = 1)
dev.off()

#ESTIMATED CUMULATIVE HAZARDS FROM HEALTHY STAGE
tmat_H <- transMat(x = list(c(2,3,4,9), c(), c(), c(), c(), c(), c(), c(), c()), 
                   names = c("Healthy", "CVD", "T2D", "CKD", "CVD+T2D", "CVD+CKD", "CKD+T2D", "CVD+T2D+CKD", "Death"))
tmat_H
mstrain_H <- msprep(data = train, trans = tmat_H, time = c(NA, "o2", "o3", "o4", "o5", "o6", "o7", "o8", "o9"), 
                    status = c(NA, "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9"), 
                    keep = c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD"))
mstrain_H <- expand.covs(mstrain_H, covs, longnames = FALSE)
mstrain_H[, c("Tstart", "Tstop", "time")] <- mstrain_H[, c("Tstart", "Tstop", "time")]/365.25 #to convert days into years
c0_H <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = mstrain_H, method = "breslow")
msf0_H <- msfit(object = c0_H, vartype = "aalen", trans = tmat_H)
png(file = "Estimated cumulative hazard from healthy state.png", width = 1080, height = 1080, type = "cairo")
plot(msf0_H, las = 1, xlab = "Duration", main = "Estimated cumulative hazard of healthy state", use.ggplot = T) 
dev.off()

#ESTIMATED CUMULATIVE HAZARDS FROM CVD STAGE
tmat_CVD <- transMat(x = list(c(4,5,8), c(), c(), c(), c(), c(), c(), c()), 
                     names = c("CVD", "T2D", "CKD", "CVD+T2D", "CVD+CKD", "CKD+T2D", "CVD+T2D+CKD", "Death"))
tmat_CVD
mstrain_CVD <- msprep(data = train, trans = tmat_CVD, time = c("o2", "o3", "o4", "o5", "o6", "o7", "o8", "o9"), 
                      status = c("c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9"), 
                      keep = c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD"))
mstrain_CVD <- expand.covs(mstrain_CVD, covs, longnames = FALSE)
mstrain_CVD[, c("Tstart", "Tstop", "time")] <- mstrain_CVD[, c("Tstart", "Tstop", "time")]/365.25 #to convert days into years
c0_CVD <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = mstrain_CVD, method = "breslow")
msf0_CVD <- msfit(object = c0_CVD, vartype = "aalen", trans = tmat_CVD)
png(file = "Estimated cumulative hazard from CVD state.png", width = 1080, height = 1080, type = "cairo")
plot(msf0_CVD, las = 1, xlab = "Duration", main = "Estimated cumulative hazard of CVD state", use.ggplot = T) 
dev.off()

#ESTIMATED CUMULATIVE HAZARDS FROM T2D STAGE
tmat_T2D <- transMat(x = list(c(3,5,7), c(), c(), c(), c(), c(), c()), 
                     names = c("T2D", "CKD", "CVD+T2D", "CVD+CKD", "CKD+T2D", "CVD+T2D+CKD", "Death"))
tmat_T2D
mstrain_T2D <- msprep(data = train, trans = tmat_T2D, time = c("o3", "o4", "o5", "o6",
                                                               "o7", "o8", "o9"), status = c("c3", "c4", "c5", "c6", "c7",
                                                                                             "c8", "c9"), keep = c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD"))
mstrain_T2D <- expand.covs(mstrain_T2D, covs, longnames = FALSE)
mstrain_T2D[, c("Tstart", "Tstop", "time")] <- mstrain_T2D[, c("Tstart", "Tstop", "time")]/365.25 #to convert days into years
c0_T2D <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = mstrain_T2D, method = "breslow")
msf0_T2D <- msfit(object = c0_T2D, vartype = "aalen", trans = tmat_T2D)
png(file = "Estimated cumulative hazard from T2D state.png", width = 1080, height = 1080, type = "cairo")
plot(msf0_T2D, las = 1, xlab = "Duration", main = "Estimated cumulative hazard of T2D state", use.ggplot = T) 
dev.off()

#ESTIMATED CUMULATIVE HAZARDS FROM CKD STAGE
tmat_CKD <- transMat(x = list(c(3,4,6), c(), c(), c(), c(), c()), 
                     names = c("CKD", "CVD+T2D", "CVD+CKD", "CKD+T2D", "CVD+T2D+CKD", "Death"))
tmat_CKD
mstrain_CKD <- msprep(data = train, trans = tmat_CKD, time = c("o4", "o5", "o6", "o7", "o8", "o9"), 
                      status = c("c4", "c5", "c6", "c7", "c8", "c9"), 
                      keep = c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD"))
mstrain_CKD <- expand.covs(mstrain_CKD, covs, longnames = FALSE)
mstrain_CKD[, c("Tstart", "Tstop", "time")] <- mstrain_CKD[, c("Tstart", "Tstop", "time")]/365.25 #to convert days into years
c0_CKD <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = mstrain_CKD, method = "breslow")
msf0_CKD <- msfit(object = c0_CKD, vartype = "aalen", trans = tmat_CKD)
png(file = "Estimated cumulative hazard from CKD state.png", width = 1080, height = 1080, type = "cairo")
plot(msf0_CKD, las = 1, xlab = "Duration", main = "Estimated cumulative hazard of CKD state", use.ggplot = T) 
dev.off()

#ESTIMATED CUMULATIVE HAZARDS FROM CVD+T2D STAGE
tmat_CD <- transMat(x = list(c(4,5), c(), c(), c(), c()), 
                    names = c("CVD+T2D", "CVD+CKD", "CKD+T2D", "CVD+T2D+CKD", "Death"))
tmat_CD
mstrain_CD <- msprep(data = train, trans = tmat_CD, time = c("o5", "o6", "o7", "o8", "o9"), 
                     status = c("c5", "c6", "c7", "c8", "c9"), 
                     keep = c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD"))
mstrain_CD <- expand.covs(mstrain_CD, covs, longnames = FALSE)
mstrain_CD[, c("Tstart", "Tstop", "time")] <- mstrain_CD[, c("Tstart", "Tstop", "time")]/365.25 #to convert days into years
c0_CD <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = mstrain_CD, method = "breslow")
msf0_CD <- msfit(object = c0_CD, vartype = "aalen", trans = tmat_CD)
png(file = "Estimated cumulative hazard from CVD+T2D state.png", width = 1080, height = 1080, type = "cairo")
plot(msf0_CD, las = 1, xlab = "Duration", main = "Estimated cumulative hazard of CVD+T2D state", use.ggplot = T) 
dev.off()

#ESTIMATED CUMULATIVE HAZARDS FROM CVD+CKD STAGE
tmat_CK <- transMat(x = list(c(3,4), c(), c(), c()), 
                    names = c("CVD+CKD", "CKD+T2D", "CVD+T2D+CKD", "Death"))
tmat_CK
mstrain_CK <- msprep(data = train, trans = tmat_CK, time = c("o6", "o7", "o8", "o9"), 
                     status = c("c6", "c7", "c8", "c9"), 
                     keep = c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD"))
mstrain_CK <- expand.covs(mstrain_CK, covs, longnames = FALSE)
mstrain_CK[, c("Tstart", "Tstop", "time")] <- mstrain_CK[, c("Tstart", "Tstop", "time")]/365.25 #to convert days into years
c0_CK <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = mstrain_CK, method = "breslow")
msf0_CK <- msfit(object = c0_CK, vartype = "aalen", trans = tmat_CK)
png(file = "Estimated cumulative hazard from CVD+CKD state.png", width = 1080, height = 1080, type = "cairo")
plot(msf0_CK, las = 1, xlab = "Duration", main = "Estimated cumulative hazard of CVD+CKD state", use.ggplot = T) 
dev.off()

#ESTIMATED CUMULATIVE HAZARDS FROM CKD+T2D STAGE
tmat_KD <- transMat(x = list(c(2,3), c(), c()), 
                    names = c("CKD+T2D", "CVD+T2D+CKD", "Death"))
tmat_KD
mstrain_KD <- msprep(data = train, trans = tmat_KD, time = c("o7", "o8", "o9"), 
                     status = c("c7", "c8", "c9"),
                     keep = c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD"))
mstrain_KD <- expand.covs(mstrain_KD, covs, longnames = FALSE)
mstrain_KD[, c("Tstart", "Tstop", "time")] <- mstrain_KD[, c("Tstart", "Tstop", "time")]/365.25 #to convert days into years
c0_KD <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = mstrain_KD, method = "breslow")
msf0_KD <- msfit(object = c0_KD, vartype = "aalen", trans = tmat_KD)
png(file = "Estimated cumulative hazard from CKD+T2D state.png", width = 1080, height = 1080, type = "cairo")
plot(msf0_KD, las = 1, xlab = "Duration", main = "Estimated cumulative hazard of CKD+T2D state", use.ggplot = T) 
dev.off()

#ESTIMATED CUMULATIVE HAZARDS FROM CVD+T2D+CKD STAGE
df_CDK <- subset(train, c8 == "1", select = c("person_id", "care_site_id", "study_dtindex", "study_dtindex_r", "c8", "o9", "c9"))
surv_obj <- survfit(Surv(o9, c9) ~ 1, data = df_CDK)
cum_haz <- -log(surv_obj$surv)
png(file = "Estimated cumulative hazard of CVD+T2D+CKD state.png", width = 1080, height = 1080, type = "cairo")
plot(cum_haz, xlab = "Duration", ylab = "Cumulative hazard", main = "Estimated cumulative hazard of CVD+T2D+CKD state", use.ggplot = T) 
dev.off()