#THIS WAS THE PRIMARY APPROACH TO VALIDATION BUT TO COMPUTATIONAL RESTRAINT IT WAS NOT POSSIBLE
#THERE IS NO USE OF COMORBIDITIES TO VALIDATE THE MODEL SO TO LESSEN THE COMPUTATIONAL LOAD, COMORBIDITIES WERE NOT CONSIDERED

#THE REQUIRED LIBRARIES
.libPaths(c("/mnt/bmh01-rds/mrc-multi-outcome/R\ packages/"))
library(dplyr)
library(mstate)
library(ggplot2)
library(calibmsm)

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

#TRAIN DATASET IN MSDATA FORMAT
tmat <- transMat(x = list(c(2,3,4,9), c(5,6,9), c(5,7,9), c(6,7,9), c(8,9), c(8,9), c(8,9), c(9), c()), 
                 names = c("Healthy", "CVD", "T2D", "CKD", "CVD+T2D", "CVD+CKD", "CKD+T2D", "CVD+T2D+CKD", "Death"))
mstrain <- msprep(data = train, trans = tmat, time = c(NA, "o2", "o3", "o4", "o5", "o6", "o7", "o8", "o9"), 
                  status = c(NA, "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9"), 
                  keep = c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD"))
covs <- c("Age", "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD")
mstrain <- expand.covs(mstrain, covs, longnames = FALSE)
mstrain[, c("Tstart", "Tstop", "time")] <- mstrain[, c("Tstart", "Tstop", "time")]/365.25 #to convert days into years

#CREATION OF STACKED DATASET OF PREDICTED TRANSITIONS PROBABILITIES 
s <- 0 #starting time
t.eval <- 1826
#creating data frame to store predicted risks
test0_temp <- vector("list", 9)
for (j in 1:9){
  test0_temp[[j]] <-  data.frame(matrix(NA, ncol = 20, nrow = nrow(mmd)))
  colnames(test0_temp[[j]]) <- c("id", paste("pstate", 1:9, sep = ""), paste("se", 1:9, sep = ""), "j")
}
test0_temp
#running through each patient id
for (id.iter in 1:nrow(mmd)){
  print(paste("id.iter = ", id.iter, Sys.time()))
  cfull <- coxph(Surv(Tstart, Tstop, status) ~ Age + gender + BMI + Cholhdl_ratio + Ethnicity6 + SBP + Smoking + IMD, 
                 data = subset(mstrain, id != id.iter), method = "breslow")
  pat.loc <- which(mstrain$id == id.iter) #getting location of individual in mstrain
  #creating a miniature dataset
  pat.dat <- mstrain[rep(pat.loc[1], 20), 9:16]
  pat.dat$trans <- 1:20
  attr(pat.dat, "trans") <- tmat
  pat.dat <- expand.covs(pat.dat, covs, longnames = FALSE)
  pat.dat$strata <- pat.dat$trans
  msf.pat <- msfit(cfull, pat.dat, trans = tmat) #obtaining cumulative incidence functions for the individual of interest
  pt <- probtrans(msf.pat, predt = s) #generating 5 year transition probabilities at time s
  #function to extract the transition probabilities from state j into each state, after followup time f.time
  extract.tp <- function(tp.object, state, f.time){
    #creating output object
    output.object <- as.numeric(base::subset(tp.object[[state]], time > f.time) |> dplyr::slice(1) |> dplyr::select(-c(time)))
    return(output.object)
  }
  #calculating required transition probabilities and storing in output dataset
  #will be generating risks out of every state j and store in tp.id
  for (j in 1:9){
    test0_temp[[j]][id.iter, ] <- c(id.iter, extract.tp(tp.object = pt, state = j, f.time = t.eval - s), j)
  }
}
test0 <- do.call("rbind", test0_temp) #combining into one dataset
test0 <- test0 |> dplyr::filter(j == 1) #deleting j != 1 at time s = 0 as individuals can't be in that stage 

#CALIBRATION PLOTS WITH EVALUATION PERIOD OF 5 YEARS (1826 DAYS)
eval5 <- 1826
#binary logistic regression models with inverse probability of censoring weights (BLR-IPCW)
dat.calib.blr <- calib_blr(data.mstate = mstrain,
                           data.raw = test0,
                           j = 1,
                           s = 0,
                           t = eval5,
                           tp.pred = test0 %>%
                             dplyr::filter(j == 1) %>%
                             dplyr::select(any_of(paste("pstate", 1:9, sep = ""))),
                           curve.type = "rcs",
                           rcs.nk = 3,
                           w.covs = covs,
                           CI = 95,
                           CI.R.boot = 200)
png(file = "5yrs calib blr.png", width = 1080, height = 1080, type = "cairo")
plot(dat.calib.blr, combine = TRUE, nrow = 3, ncol = 3)
dev.off()
#multinomial logistic regression models with inverse probability of censoring weights (MLR-IPCW)
dat.calib.mlr <- calib_mlr(data.mstate = mstrain,
                           data.raw = test,
                           j = 1,
                           s = 0,
                           t = eval5,
                           tp.pred = test0 %>%
                             dplyr::filter(j == 1) %>%
                             dplyr::select(any_of(paste("pstate", 1:9, sep = ""))),
                           curve.type = "rcs",
                           rcs.nk = 3,
                           w.covs = covs)
png(file = "5yrs calib mlr.png", width = 1080, height = 1080, type = "cairo")
plot(dat.calib.mlr, combine = TRUE, nrow = 3, ncol = 3)
dev.off()
#pseudo-value approach (PV)
dat.calib.pv <- calib_pv(data.mstate = mstrain,
                         data.raw = test,
                         j = 1,
                         s = 0,
                         t = eval5,
                         tp.pred = test0 %>%
                           dplyr::filter(j == 1) %>%
                           dplyr::select(any_of(paste("pstate", 1:9, sep = ""))),
                         curve.type = "rcs",
                         rcs.nk = 3,
                         group.vars = c("Ethnicity6"),
                         n.pctls = 5,
                         CI = 95,
                         CI.type = "parametric")
png(file = "5yrs calib pv.png", width = 1080, height = 1080, type = "cairo")
plot(dat.calib.pv, combine = TRUE, nrow = 3, ncol = 3)
dev.off()

#CALIBRATION PLOTS WITH EVALUATION PERIOD OF 10 YEARS (3652 DAYS)
eval10 <- 3652
#binary logistic regression models with inverse probability of censoring weights (BLR-IPCW)
dat.calib.blr <- calib_blr(data.mstate = mstrain,
                           data.raw = test,
                           j = 1,
                           s = 0,
                           t = eval10,
                           tp.pred = test0 %>%
                             dplyr::filter(j == 1) %>%
                             dplyr::select(any_of(paste("pstate", 1:9, sep = ""))),
                           curve.type = "rcs",
                           rcs.nk = 3,
                           w.covs = covs,
                           CI = 95,
                           CI.R.boot = 200)
png(file = "10yrs calib blr.png", width = 1080, height = 1080, type = "cairo")
plot(dat.calib.blr, combine = TRUE, nrow = 3, ncol = 3)
dev.off()
#multinomial logistic regression models with inverse probability of censoring weights (MLR-IPCW)
dat.calib.mlr <- calib_mlr(data.mstate = mstrain,
                           data.raw = test,
                           j = 1,
                           s = 0,
                           t = eval10,
                           tp.pred = test0 %>%
                             dplyr::filter(j == 1) %>%
                             dplyr::select(any_of(paste("pstate", 1:9, sep = ""))),
                           curve.type = "rcs",
                           rcs.nk = 3,
                           w.covs = covs)
png(file = "10yrs calib mlr.png", width = 1080, height = 1080, type = "cairo")
plot(dat.calib.mlr, combine = TRUE, nrow = 3, ncol = 3)
dev.off()
#pseudo-value approach (PV)
dat.calib.pv <- calib_pv(data.mstate = mstrain,
                         data.raw = test,
                         j = 1,
                         s = 0,
                         t = eval10,
                         tp.pred = test0 %>%
                           dplyr::filter(j == 1) %>%
                           dplyr::select(any_of(paste("pstate", 1:9, sep = ""))),
                         curve.type = "rcs",
                         rcs.nk = 3,
                         group.vars = c("Ethnicity6"),
                         n.pctls = 5,
                         CI = 95,
                         CI.type = "parametric")
png(file = "10yrs calib pv.png", width = 1080, height = 1080, type = "cairo")
plot(dat.calib.pv, combine = TRUE, nrow = 3, ncol = 3)
dev.off()