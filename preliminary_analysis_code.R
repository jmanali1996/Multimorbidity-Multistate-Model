#LOADING THE REQUIRED LIBRARIES
.libPaths(c("/mnt/bmh01-rds/mrc-multi-outcome/R\ packages/"))
library(dplyr)

#READING THE DATA
df <- readRDS("multimorbidity_data.rds")

#FILTERING THE DATA
mmd <- select(filter(df, CVD_hist == "0", Diab_t2_hist == "0", CKD_hist == "0"), 
              c("person_id", "care_site_id", "study_dtindex", "study_dtindex_r", "Age", 
                "gender", "BMI", "Cholhdl_ratio", "Ethnicity6", "SBP", "Smoking", "IMD",
                "Alcohol_misuse", "Eating_disorders", "Asthma", "Anxiety_disorders", 
                "Depression", "Visual_impairment", "Bronchiectasis", "Hepatic_failure", 
                "Viral_hepatitis", "Sinusitis", "COPD", "Dementia", "Diverticular", 
                "Epilepsy", "Hearing_loss", "Hypertension", "IBS", "Intellectual_dis",
                "MS", "Parkinsons", "Perip_vascular", "Psoriasis", "Substance_misuse",
                "RA", "Schizophrenia", "Bipolar", "Thyroid", "Peptic_ulcer", "IBD", 
                "Prostate", "Death_t", "Death_c", "Death_NelsonAalen_link", "CKD_hist", 
                "CKD_hist_t", "CKD_ev_c", "CKD_ev_t", "Diab_t2_hist", "Diab_t2_hist_t", 
                "Diab_t2_ev_c", "Diab_t2_ev_t", "dtcens_combdeath_r", "CVD_ev_t", 
                "CVD_ev_c", "CVD_hist", "dtcens_var", "o2", "c2", "o3", "c3", "o4", 
                "c4", "o5", "c5", "o6", "c6", "o7", "c7", "o8", "c8", "o9", "c9"))

#PRELIMINARY ANALYSIS OF THE COMPLETE DATA
summary(mmd) #to know the data summary at a glance
#to calculate standard deviation
sd(mmd$Age) 
sd(mmd$BMI) 
sd(mmd$Cholhdl_ratio) 
sd(mmd$SBP) 
#to calculate coefficient of variation
sd(mmd$Age)/mean(mmd$Age)*100 
sd(mmd$BMI)/mean(mmd$BMI)*100 
sd(mmd$Cholhdl_ratio)/mean(mmd$Cholhdl_ratio)*100 
sd(mmd$SBP)/mean(mmd$SBP)*100 
#to calculate proportion distribution
round(prop.table(table(mmd$gender))*100, 0) 
round(prop.table(table(mmd$Ethnicity6))*100, 0)
round(prop.table(table(mmd$Smoking))*100, 0)
round(prop.table(table(mmd$IMD))*100, 0)
round(prop.table(table(mmd$Alcohol_misuse))*100, 0)
round(prop.table(table(mmd$Eating_disorders))*100, 0)
round(prop.table(table(mmd$Asthama))*100, 0)
round(prop.table(table(mmd$Anxiety_disorders))*100, 0)
round(prop.table(table(mmd$Depression))*100, 0)
round(prop.table(table(mmd$Visual_impairment))*100, 0)
round(prop.table(table(mmd$Bronchiectasis))*100, 0)
round(prop.table(table(mmd$Hepatic_failure))*100, 0)
round(prop.table(table(mmd$Viral_hepatitis))*100, 0)
round(prop.table(table(mmd$Sinusitis))*100, 0)
round(prop.table(table(mmd$COPD))*100, 0)
round(prop.table(table(mmd$Dementia))*100, 0)
round(prop.table(table(mmd$Diverticular))*100, 0)
round(prop.table(table(mmd$Epilepsy))*100, 0)
round(prop.table(table(mmd$Hearing_loss))*100, 0)
round(prop.table(table(mmd$Hypertension))*100, 0)
round(prop.table(table(mmd$IBS))*100, 0)
round(prop.table(table(mmd$Intellectual_dis))*100, 0)
round(prop.table(table(mmd$MS))*100, 0)
round(prop.table(table(mmd$Parkinsons))*100, 0)
round(prop.table(table(mmd$Perip_vascular))*100, 0)
round(prop.table(table(mmd$Psoriasis))*100, 0)
round(prop.table(table(mmd$Substance_misuse))*100, 0)
round(prop.table(table(mmd$RA))*100, 0)
round(prop.table(table(mmd$Schizophrenia))*100, 0)
round(prop.table(table(mmd$Bipolar))*100, 0)
round(prop.table(table(mmd$Thyroid))*100, 0)
round(prop.table(table(mmd$Peptic_ulcer))*100, 0)
round(prop.table(table(mmd$IBD))*100, 0)
round(prop.table(table(mmd$Prostate))*100, 0)

#CREATING TRAIN AND TEST DATASETS
set.seed(1)  #to make reproducible
sample <- sample(c(TRUE, FALSE), nrow(mmd), replace=TRUE, prob=c(0.7,0.3))
train <- mmd[sample, ]
test <- mmd[!sample, ]

#PRELIMINARY ANALYSIS OF THE TRAINING DATA
summary(train) #to know the data summary at a glance
#to calculate standard deviation
sd(train$Age) 
sd(train$BMI) 
sd(train$Cholhdl_ratio) 
sd(train$SBP) 
#to calculate coefficient of variation
sd(train$Age)/mean(train$Age)*100 
sd(train$BMI)/mean(train$BMI)*100 
sd(train$Cholhdl_ratio)/mean(train$Cholhdl_ratio)*100 
sd(train$SBP)/mean(train$SBP)*100 
#to calculate proportion distribution
round(prop.table(table(train$gender))*100, 0) 
round(prop.table(table(train$Ethnicity6))*100, 0)
round(prop.table(table(train$Smoking))*100, 0)
round(prop.table(table(train$IMD))*100, 0)
round(prop.table(table(train$Alcohol_misuse))*100, 0)
round(prop.table(table(train$Eating_disorders))*100, 0)
round(prop.table(table(train$Asthama))*100, 0)
round(prop.table(table(train$Anxiety_disorders))*100, 0)
round(prop.table(table(train$Depression))*100, 0)
round(prop.table(table(train$Visual_impairment))*100, 0)
round(prop.table(table(train$Bronchiectasis))*100, 0)
round(prop.table(table(train$Hepatic_failure))*100, 0)
round(prop.table(table(train$Viral_hepatitis))*100, 0)
round(prop.table(table(train$Sinusitis))*100, 0)
round(prop.table(table(train$COPD))*100, 0)
round(prop.table(table(train$Dementia))*100, 0)
round(prop.table(table(train$Diverticular))*100, 0)
round(prop.table(table(train$Epilepsy))*100, 0)
round(prop.table(table(train$Hearing_loss))*100, 0)
round(prop.table(table(train$Hypertension))*100, 0)
round(prop.table(table(train$IBS))*100, 0)
round(prop.table(table(train$Intellectual_dis))*100, 0)
round(prop.table(table(train$MS))*100, 0)
round(prop.table(table(train$Parkinsons))*100, 0)
round(prop.table(table(train$Perip_vascular))*100, 0)
round(prop.table(table(train$Psoriasis))*100, 0)
round(prop.table(table(train$Substance_misuse))*100, 0)
round(prop.table(table(train$RA))*100, 0)
round(prop.table(table(train$Schizophrenia))*100, 0)
round(prop.table(table(train$Bipolar))*100, 0)
round(prop.table(table(train$Thyroid))*100, 0)
round(prop.table(table(train$Peptic_ulcer))*100, 0)
round(prop.table(table(train$IBD))*100, 0)
round(prop.table(table(train$Prostate))*100, 0)

#PRELIMINARY ANALYSIS OF THE TESTING DATA
summary(test) #to know the data summary at a glance
#to calculate standard deviation
sd(test$Age) 
sd(test$BMI) 
sd(test$Cholhdl_ratio) 
sd(test$SBP) 
#to calculate coefficient of variation
sd(test$Age)/mean(test$Age)*100 
sd(test$BMI)/mean(test$BMI)*100 
sd(test$Cholhdl_ratio)/mean(test$Cholhdl_ratio)*100 
sd(test$SBP)/mean(test$SBP)*100 
#to calculate proportion distribution
round(prop.table(table(test$gender))*100, 0) 
round(prop.table(table(test$Ethnicity6))*100, 0)
round(prop.table(table(test$Smoking))*100, 0)
round(prop.table(table(test$IMD))*100, 0)
round(prop.table(table(test$Alcohol_misuse))*100, 0)
round(prop.table(table(test$Eating_disorders))*100, 0)
round(prop.table(table(test$Asthama))*100, 0)
round(prop.table(table(test$Anxiety_disorders))*100, 0)
round(prop.table(table(test$Depression))*100, 0)
round(prop.table(table(test$Visual_impairment))*100, 0)
round(prop.table(table(test$Bronchiectasis))*100, 0)
round(prop.table(table(test$Hepatic_failure))*100, 0)
round(prop.table(table(test$Viral_hepatitis))*100, 0)
round(prop.table(table(test$Sinusitis))*100, 0)
round(prop.table(table(test$COPD))*100, 0)
round(prop.table(table(test$Dementia))*100, 0)
round(prop.table(table(test$Diverticular))*100, 0)
round(prop.table(table(test$Epilepsy))*100, 0)
round(prop.table(table(test$Hearing_loss))*100, 0)
round(prop.table(table(test$Hypertension))*100, 0)
round(prop.table(table(test$IBS))*100, 0)
round(prop.table(table(test$Intellectual_dis))*100, 0)
round(prop.table(table(test$MS))*100, 0)
round(prop.table(table(test$Parkinsons))*100, 0)
round(prop.table(table(test$Perip_vascular))*100, 0)
round(prop.table(table(test$Psoriasis))*100, 0)
round(prop.table(table(test$Substance_misuse))*100, 0)
round(prop.table(table(test$RA))*100, 0)
round(prop.table(table(test$Schizophrenia))*100, 0)
round(prop.table(table(test$Bipolar))*100, 0)
round(prop.table(table(test$Thyroid))*100, 0)
round(prop.table(table(test$Peptic_ulcer))*100, 0)
round(prop.table(table(test$IBD))*100, 0)
round(prop.table(table(test$Prostate))*100, 0)