## Load packages
#install.packages("lavaan", dependencies = TRUE)
library(lavaan)
library(dplyr)
library(tidyverse)
library(rstatix)

## Indirect effect: Product of Coefficients
## Sobel test (Delta Method)

## Load in data
## Load in Chronic Pain (Keira's Phenotype) and Depression (Natasha's Phenotype) data
UKB_CP <- read.csv("/Users/hannahcasey/Desktop/PhD/output/UKB/pain/UKBBChronicPainJohnsonEligible.csv")

## Remame ID column, f_eid -> f_eid
UKB_CP <- UKB_CP %>% 
  rename(
    f.eid = n_eid)

## Identify chronic pain individuals (at least one site of chronic pain)
UKB_CP$ChronicPainStatus <- NA
UKB_CP$ChronicPainStatus[UKB_CP$sites > 0] <- 1
UKB_CP$ChronicPainStatus[UKB_CP$sites < 1] <- 0

## Load in MDD derived from Naatasha's Algorithm
UKB_MDD <- read.csv("/Users/hannahcasey/Desktop/PhD/output/UKB/depression/UKBBNatashaPhenotypeDepressionMania.csv")

## Remame ID column, f_eid -> f_eid
UKB_MDD <- UKB_MDD %>% 
  rename(f.eid = f_eid)

## Define MDD (probable recurrent MDD) and no MDD (No probable MDD and probable single MDD)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Add extra column indicating depression status
UKB_MDD$DepressionStatus <- NA
UKB_MDD$DepressionStatus[UKB_MDD$patient_group == "recurrent depression"] <- "Probable Recurrent MDD"
UKB_MDD$DepressionStatus[UKB_MDD$patient_group == "single depression"] <- "Probable Single MDD"
UKB_MDD$DepressionStatus[UKB_MDD$patient_group %in% c("control", "bipolar depression", "unipolar mania")] <- "No Probable MDD"

UKB_MDD$RecurrentMDDStatus <- "No Probable Recurrent MDD"
UKB_MDD$RecurrentMDDStatus[UKB_MDD$DepressionStatus == "Probable Recurrent MDD"] <-  "Probable Recurrent MDD"
freq_table(UKB_MDD$RecurrentMDDStatus) 

## Sort participants into comorbidity groups
## Merge depression and chronic pain data 
UKB_CP_MDD <- full_join(UKB_MDD, UKB_CP[,c("f.eid", "ChronicPainStatus", "sites")], by = "f.eid")

## Create new column indicating which comorbid group participant belong to
UKB_CP_MDD$ComorbidStatus <- NA
UKB_CP_MDD$ComorbidStatus[UKB_CP_MDD$RecurrentMDDStatus == "No Probable Recurrent MDD"  & UKB_CP_MDD$ChronicPainStatus == 0] <- "No Probable Recurrent MDD + No Chronic Pain"
UKB_CP_MDD$ComorbidStatus[UKB_CP_MDD$RecurrentMDDStatus == "Probable Recurrent MDD" & UKB_CP_MDD$ChronicPainStatus == 0] <- "Probable Recurrent MDD + No Chronic Pain"
UKB_CP_MDD$ComorbidStatus[UKB_CP_MDD$RecurrentMDDStatus == "No Probable Recurrent MDD"  & UKB_CP_MDD$ChronicPainStatus == 1] <- "No Probable Recurrent MDD + Chronic Pain"
UKB_CP_MDD$ComorbidStatus[UKB_CP_MDD$RecurrentMDDStatus == "Probable Recurrent MDD" & UKB_CP_MDD$ChronicPainStatus == 1] <- "Probable Recurrent MDD + Chronic Pain"


## Check frequency of each comorbidity group
freq_table(UKB_CP_MDD$ComorbidStatus)

## Load in datasets containing inflammation associated imaging features
UKB_CRP_PRS_DK <- read.csv("~/Desktop/PhD/projects/UKBCRPImagingPRS/resources/UKB_CRP_PRS_DK.csv")

## Load in NMR metabolomics data
UKB_NMR <- read.table("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-10-nmr-metabolomics-ukb48936/NMRMetabolomics.tsv.gz", header = TRUE)

## Extract ID column and columns pretaining to GlycA measurment
UKB_NMR_glycA <- cbind(UKB_NMR[,1], UKB_NMR[grepl("f.23480.", names(UKB_NMR))])

## Remove redundant dataframes 
rm(UKB_NMR)

## Rename ID column
UKB_NMR_glycA <- UKB_NMR_glycA %>%
  rename(f.eid = "UKB_NMR[, 1]",
         glycA_1 = f.23480.0.0,
         glycA_2 = f.23480.1.0)

UKB_NMR_glycA <- UKB_NMR_glycA[!is.na(UKB_NMR_glycA$glycA_1),]


## Join all data in single dataframe
UKB_glycA_CPMDD_imaging_SEM <- full_join(UKB_CP_MDD, UKB_CRP_PRS_DK, by = "f.eid")
UKB_glycA_CPMDD_imaging_SEM <- full_join(UKB_glycA_CPMDD_imaging_SEM, UKB_NMR_glycA, by = "f.eid")

## Only keep essential columns: FIDs, GlycA (1st instance), comorbidity status, MDD status, CP status (binary and ordinal), consistently associated imaging features, covariates (age, age squared, sex, BMI, assessment centre, ICV)
## Imaging features associated with both GlycA and CP+MDD+ in the same direction:
## Temporal Lobe - decrease
## Middletemporal - decrease (f.26735.2.0, f.26836.2.0)
## Thalamus - decrease (f.25011.2.0 , f.25012.2.0)
## Amygdala - decrease (f.25021.2.0, f.25022.2.0)

UKB_glycA_CPMDD_imaging_SEM <- UKB_glycA_CPMDD_imaging_SEM %>%
  select(f.eid, glycA_1, ComorbidStatus, ChronicPainStatus, sites, DepressionStatus, temporal_lobe, f.26735.2.0, f.26836.2.0, 
         f.25011.2.0, f.25012.2.0, f.25021.2.0, f.25022.2.0, age, age_squared, sex, BMI, assessment_centre_first_imaging, ICV) %>%
  rename(chronic_pain_sites = sites,
         middletemporal_left = f.26735.2.0,
         middletemporal_right = f.26836.2.0,
         thalamus_left = f.25011.2.0,
         thalamus_right = f.25012.2.0,
         amygdala_left = f.25021.2.0,
         amygdala_right = f.25022.2.0)


## Recode catagorical data
UKB_glycA_CPMDD_imaging_SEM <- UKB_glycA_CPMDD_imaging_SEM %>% mutate(sex=recode(sex, "Male" = 0, "Female" = 1),
    DepressionStatus=recode(DepressionStatus, "No Probable MDD" = 0, "Probable Recurrent MDD" = 1, "Probable Single MDD" = 0),
    ComorbidStatus_CPMDD =recode(ComorbidStatus, "Probable Recurrent MDD + Chronic Pain" = 1, "No Probable Recurrent MDD + Chronic Pain" = 0, "Probable Recurrent MDD + No Chronic Pain" = 0, "No Probable Recurrent MDD + No Chronic Pain" = 0),
    ComorbidStatus_CPnoMDD =recode(ComorbidStatus, "Probable Recurrent MDD + Chronic Pain" = 0, "No Probable Recurrent MDD + Chronic Pain" = 1, "Probable Recurrent MDD + No Chronic Pain" = 0, "No Probable Recurrent MDD + No Chronic Pain" = 0),
    ComorbidStatus_noCP_MDD =recode(ComorbidStatus,  "Probable Recurrent MDD + Chronic Pain" = 0, "No Probable Recurrent MDD + Chronic Pain" = 0, "Probable Recurrent MDD + No Chronic Pain" = 1, "No Probable Recurrent MDD + No Chronic Pain" = 0),
    ComorbidStatus_noCP_noMDD =recode(ComorbidStatus,  "Probable Recurrent MDD + Chronic Pain" = 0, "No Probable Recurrent MDD + Chronic Pain" = 0, "Probable Recurrent MDD + No Chronic Pain" = 0, "No Probable Recurrent MDD + No Chronic Pain" = 1))


## Carry out mediation analysis
## Independent variable = GlycA
## Dependent = MDD status
## Mediator = Middletemporal


## Remove outliers from total (hemispheres summed) middletemporal volume
## Combine middletemporal 
UKB_glycA_CPMDD_imaging_SEM$middletemporal_total <- UKB_glycA_CPMDD_imaging_SEM$middletemporal_left + UKB_glycA_CPMDD_imaging_SEM$middletemporal_right

## Plot distribution of total middletemporal volume
hist(UKB_glycA_CPMDD_imaging_SEM$middletemporal_total)

## Remove outliers
outliers <- boxplot(UKB_glycA_CPMDD_imaging_SEM$middletemporal_total, plot=FALSE)$out
UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed <- UKB_glycA_CPMDD_imaging_SEM
UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed <- UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed[-which(UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed$middletemporal_total %in% outliers),]
hist(UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed$middletemporal_total)

## Scale neuroimaging variables
#UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed_temp$glycA_1 <- scale(UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed_temp$glycA_1)
UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed$middletemporal_total <- scale(UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed$middletemporal_total)
UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed$ICV <- scale(UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed$ICV)


## Specify model
specmod <- "
# Path c` (direct effect)
ComorbidStatus_CPMDD ~ c*glycA_1
# Path a 
middletemporal_total ~ a*glycA_1 + age + age_squared + BMI + assessment_centre_first_imaging  + ICV + sex

# Path b
ComorbidStatus_CPMDD ~ b*middletemporal_total + age + age_squared + BMI + assessment_centre_first_imaging + ICV + sex

## Indirect effect (a*b): Sobel Test (Delta Method)
ab := a*b
"

## Fit model
fitmod <- sem(specmod, data = UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed)
## Summarize the results
summary(fitmod, fit.measures = TRUE, rsquare = TRUE)


## Resampling Method: Percentile bootstrapping
set.seed(151097)

fitmod2 <- sem(specmod, data = UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed, se = "bootstrap", bootstrap = 5000)
summary(fitmod2, fit.measures = TRUE, rsquare = TRUE)
parameterEstimates(fitmod2, ci = TRUE, level = 0.95, boot.ci.type = "perc")




## Specify model
specmod <- "
# Path c` (direct effect)
ComorbidStatus_CPnoMDD ~ c*glycA_1
# Path a 
middletemporal_total ~ a*glycA_1 + age + age_squared + BMI + assessment_centre_first_imaging  + ICV + sex

# Path b
ComorbidStatus_CPnoMDD ~ b*middletemporal_total + age + age_squared + BMI + assessment_centre_first_imaging + ICV + sex

## Indirect effect (a*b): Sobel Test (Delta Method)
ab := a*b
"
## Fit model
fitmod <- sem(specmod, data = UKB_glycA_CPMDD_imaging_SEM_middletemporal_total_outliers_removed)
## Summarize the results
summary(fitmod, fit.measures = TRUE, rsquare = TRUE)


## Carry out mediation analysis
## Independent variable = GlycA
## Dependent = MDD status
## Mediator = Middletemporal


## Remove outliers from total (hemispheres summed) middletemporal volume
## Combine middletemporal 
UKB_glycA_CPMDD_imaging_SEM$thalamus_total <- UKB_glycA_CPMDD_imaging_SEM$thalamus_left + UKB_glycA_CPMDD_imaging_SEM$thalamus_right

## Plot distribution of total middletemporal volume
hist(UKB_glycA_CPMDD_imaging_SEM$thalamus_total)

## Remove outliers
outliers <- boxplot(UKB_glycA_CPMDD_imaging_SEM$thalamus_total, plot=FALSE)$out
UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed <- UKB_glycA_CPMDD_imaging_SEM
UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed <- UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed[-which(UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed$thalamus_total %in% outliers),]
hist(UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed$thalamus_total)

## Scale neuroimaging variables
#UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed_temp$glycA_1 <- scale(UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed_temp$glycA_1)
UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed$thalamus_total <- scale(UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed$thalamus_total)
UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed$ICV <- scale(UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed$ICV)


## Specify model
specmod <- "
# Path c` (direct effect)
ComorbidStatus_CPMDD ~ c*glycA_1
# Path a 
thalamus_total ~ a*glycA_1 + age + age_squared + BMI + assessment_centre_first_imaging  + ICV + sex

# Path b
ComorbidStatus_CPMDD ~ b*thalamus_total + age + age_squared + BMI + assessment_centre_first_imaging + ICV + sex

## Indirect effect (a*b): Sobel Test (Delta Method)
ab := a*b
"

## Fit model
fitmod <- sem(specmod, data = UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed)
## Summarize the results
summary(fitmod, fit.measures = TRUE, rsquare = TRUE)

## Resampling Method: Percentile bootstrapping
set.seed(151097)

fitmod2 <- sem(specmod, data = UKB_glycA_CPMDD_imaging_SEM_thalamus_total_outliers_removed, se = "bootstrap", bootstrap = 500)
summary(fitmod2, fit.measures = TRUE, rsquare = TRUE)
parameterEstimates(fitmod2, ci = TRUE, level = 0.95, boot.ci.type = "perc")



## Carry out mediation analysis
## Independent variable = GlycA
## Dependent = MDD status
## Mediator = Middletemporal


## Remove outliers from total (hemispheres summed) middletemporal volume
## Combine middletemporal 
UKB_glycA_CPMDD_imaging_SEM$amygdala_total <- UKB_glycA_CPMDD_imaging_SEM$amygdala_left + UKB_glycA_CPMDD_imaging_SEM$amygdala_right

## Plot distribution of total middletemporal volume
hist(UKB_glycA_CPMDD_imaging_SEM$amygdala_total)

## Remove outliers
outliers <- boxplot(UKB_glycA_CPMDD_imaging_SEM$amygdala_total, plot=FALSE)$out
UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed <- UKB_glycA_CPMDD_imaging_SEM
UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed <- UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed[-which(UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$amygdala_total %in% outliers),]
hist(UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$amygdala_total)

## Scale neuroimaging variables
#UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed_temp$glycA_1 <- scale(UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed_temp$glycA_1)
UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$amygdala_total <- scale(UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$amygdala_total)
UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$ICV <- scale(UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$ICV)


## Specify model
specmod <- "
# Path c` (direct effect)
ComorbidStatus_CPMDD ~ c*glycA_1
# Path a 
amygdala_total ~ a*glycA_1 + age + age_squared + BMI + assessment_centre_first_imaging  + ICV

# Path b
ComorbidStatus_CPMDD ~ b*amygdala_total + age + age_squared + BMI + assessment_centre_first_imaging + ICV

## Indirect effect (a*b): Sobel Test (Delta Method)
ab := a*b
"

## Fit model
fitmod <- sem(specmod, data = UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed)
## Summarize the results
summary(fitmod, fit.measures = TRUE, rsquare = TRUE)





## Carry out mediation analysis
## Independent variable = GlycA
## Dependent = MDD status
## Mediator = temporal lobe


## Remove outliers from total (hemispheres summed) middletemporal volume
## Combine middletemporal 
UKB_glycA_CPMDD_imaging_SEM$amygdala_total <- UKB_glycA_CPMDD_imaging_SEM$amygdala_left + UKB_glycA_CPMDD_imaging_SEM$amygdala_right

## Plot distribution of total middletemporal volume
hist(UKB_glycA_CPMDD_imaging_SEM$amygdala_total)

## Remove outliers
outliers <- boxplot(UKB_glycA_CPMDD_imaging_SEM$amygdala_total, plot=FALSE)$out
UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed <- UKB_glycA_CPMDD_imaging_SEM
UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed <- UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed[-which(UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$amygdala_total %in% outliers),]
hist(UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$amygdala_total)

## Scale neuroimaging variables
#UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed_temp$glycA_1 <- scale(UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed_temp$glycA_1)
UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$amygdala_total <- scale(UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$amygdala_total)
UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$ICV <- scale(UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed$ICV)


## Specify model
specmod <- "
# Path c` (direct effect)
ComorbidStatus_CPMDD ~ c*glycA_1
# Path a 
amygdala_total ~ a*glycA_1 + age + age_squared + BMI + assessment_centre_first_imaging  + ICV

# Path b
ComorbidStatus_CPMDD ~ b*amygdala_total + age + age_squared + BMI + assessment_centre_first_imaging + ICV

## Indirect effect (a*b): Sobel Test (Delta Method)
ab := a*b
"

## Fit model
fitmod <- sem(specmod, data = UKB_glycA_CPMDD_imaging_SEM_amygdala_total_outliers_removed)
## Summarize the results
summary(fitmod, fit.measures = TRUE, rsquare = TRUE)

