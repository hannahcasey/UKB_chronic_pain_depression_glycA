## Script Overview
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`



## Load Packages
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
library(foreign)
#install.packages("ggplot2")
library(ggplot2)
library(dplyr)
#library(mosaic)
#install.packages("hrbrthemes")
library(hrbrthemes)
#library(jtools)
#install.packages("rstatix")
library(rstatix)
#install.packages("emmeans")
library(emmeans)
#install.packages("RNOmni")
library(RNOmni)
#install.packages("effects")
#install.packages("hrbrthemes")
library(hrbrthemes)
library(car) ## levenes test
library(tidyverse)

## Load and tidy data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
## Load in NMR metabolomics data
UKB_NMR <- read.table("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-10-nmr-metabolomics-ukb48936/NMRMetabolomics.tsv.gz", header = TRUE)
## f.30710 = C-reactive protein
## Extract ID column and columns pretaining to CRP measurment
UKB_NMR_glycA <- cbind(UKB_NMR[,1], UKB_NMR[grepl("f.23480.", names(UKB_NMR))])

## Remove redundant dataframes 
rm(UKB_NMR)
## Rename ID column
UKB_NMR_glycA <- UKB_NMR_glycA %>%
  rename(f.eid = "UKB_NMR[, 1]",
         glycA_1 = f.23480.0.0,
         glycA_2 = f.23480.1.0)


## Load in Natasha's Phenotype of MDD
UKB_Natasha_MDD_phenotype <- read.csv("~/Desktop/PhD/output/UKB/depression/UKBBNatashaPhenotypeDepressionMania.csv", header = T)
## Rename ID column\
UKB_Natasha_MDD_phenotype <- UKB_Natasha_MDD_phenotype %>%
  rename(f.eid = f_eid)

## Add extra column indicating depression status
UKB_Natasha_MDD_phenotype$DepressionStatus <- NA
UKB_Natasha_MDD_phenotype$DepressionStatus[UKB_Natasha_MDD_phenotype$patient_group == "recurrent depression"] <- "Probable Recurrent MDD"
UKB_Natasha_MDD_phenotype$DepressionStatus[UKB_Natasha_MDD_phenotype$patient_group == "single depression"] <- "Probable Single MDD"
UKB_Natasha_MDD_phenotype$DepressionStatus[UKB_Natasha_MDD_phenotype$patient_group %in% c("control", "bipolar depression", "unipolar mania")] <- "No Probable MDD"



## Define MDD (probable recurrent MDD) and no MDD (No probable MDD and probable single MDD)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
UKB_Natasha_MDD_phenotype$RecurrentMDDStatus <- "No Probable Recurrent MDD"
UKB_Natasha_MDD_phenotype$RecurrentMDDStatus[UKB_Natasha_MDD_phenotype$DepressionStatus == "Probable Recurrent MDD"] <-  "Probable Recurrent MDD"
freq_table(UKB_Natasha_MDD_phenotype$RecurrentMDDStatus)

## Load in chronic pain data (Keira Johnson Phenotype)
UKB_chronic_pain_Johnson_eligible <- read.csv("~/Desktop/PhD/output/UKB/pain/UKBBChronicPainJohnsonEligible.csv")

## Rename ID column\
UKB_chronic_pain_Johnson_eligible <- UKB_chronic_pain_Johnson_eligible %>%
  rename(f.eid = n_eid)

## Identify chronic pain individuals (at least one site of chronic pain)
UKB_chronic_pain_Johnson_eligible$ChronicPainStatus <- NA
UKB_chronic_pain_Johnson_eligible$ChronicPainStatus[UKB_chronic_pain_Johnson_eligible$sites > 0] <- 1
UKB_chronic_pain_Johnson_eligible$ChronicPainStatus[UKB_chronic_pain_Johnson_eligible$sites < 1] <- 0

## Remove chronic pain status of those with widespread chronic pain
UKB_chronic_pain_Johnson_eligible$ChronicPainStatus[UKB_chronic_pain_Johnson_eligible$sites == 8] <- NA


##Quantify chronic pain sites
UKB_chronic_pain_Johnson_eligible$ChronicPainGroup<- NA
UKB_chronic_pain_Johnson_eligible$ChronicPainGroup[UKB_chronic_pain_Johnson_eligible$sites == 0] <- "No Sites"
UKB_chronic_pain_Johnson_eligible$ChronicPainGroup[UKB_chronic_pain_Johnson_eligible$sites > 0 & UKB_chronic_pain_Johnson_eligible$sites < 3] <- "1-2 Sites"
UKB_chronic_pain_Johnson_eligible$ChronicPainGroup[UKB_chronic_pain_Johnson_eligible$sites > 2 & UKB_chronic_pain_Johnson_eligible$sites < 5] <- "3-4 Sites"
UKB_chronic_pain_Johnson_eligible$ChronicPainGroup[UKB_chronic_pain_Johnson_eligible$sites > 4 & UKB_chronic_pain_Johnson_eligible$sites < 8] <- "5-7 Sites"
UKB_chronic_pain_Johnson_eligible$ChronicPainGroup[UKB_chronic_pain_Johnson_eligible$sites == 8] <- "Widespread"


## load in CRP covariate file
UKB_covariates <- read.table("~/Desktop/PhD/projects/UKBChronicPainDepressionCRP/resources/UKB_CRP_covariates_CRP.txt", header = T)

UKB_covariates <- UKB_covariates %>%
  rename(f.eid = n_eid)

## convert covariates to factors form statistical analysis

covars <- names(UKB_covariates)[-1]
UKB_covariates[covars] <- lapply(UKB_covariates[covars], factor)  ## as.factor() could also be used

## Combine data into one dataset
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
UKB_CP_MDD_glycA <- left_join(UKB_NMR_glycA, UKB_Natasha_MDD_phenotype, by = "f.eid")
UKB_CP_MDD_glycA <- left_join(UKB_CP_MDD_glycA, UKB_chronic_pain_Johnson_eligible, by = "f.eid")
UKB_CP_MDD_glycA <- left_join(UKB_CP_MDD_glycA, UKB_covariates, by = "f.eid")

## Tidy chronic pain depression comorbidity data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create new column indicating which comorbid group participant belong to
UKB_CP_MDD_glycA$ComorbidStatus <- NA
UKB_CP_MDD_glycA$ComorbidStatus[UKB_CP_MDD_glycA$RecurrentMDDStatus == "No Probable Recurrent MDD"  & UKB_CP_MDD_glycA$ChronicPainStatus == 0] <- "No Probable Recurrent MDD + No Chronic Pain"
UKB_CP_MDD_glycA$ComorbidStatus[UKB_CP_MDD_glycA$RecurrentMDDStatus == "Probable Recurrent MDD" & UKB_CP_MDD_glycA$ChronicPainStatus == 0] <- "Probable Recurrent MDD + No Chronic Pain"
UKB_CP_MDD_glycA$ComorbidStatus[UKB_CP_MDD_glycA$RecurrentMDDStatus == "No Probable Recurrent MDD"  & UKB_CP_MDD_glycA$ChronicPainStatus == 1] <- "No Probable Recurrent MDD + Chronic Pain"
UKB_CP_MDD_glycA$ComorbidStatus[UKB_CP_MDD_glycA$RecurrentMDDStatus == "Probable Recurrent MDD" & UKB_CP_MDD_glycA$ChronicPainStatus == 1] <- "Probable Recurrent MDD + Chronic Pain"

## Check that all participants have been assigned to a comorbidity group
freq_table(UKB_CP_MDD_glycA$ComorbidStatus, na.rm = TRUE)


## Plot GlycA in chronic pain groups/status and MDD status
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Sort chronic pain groups prior to plotting:
UKB_CP_MDD_glycA$ChronicPainGroup <- factor(UKB_CP_MDD_glycA$ChronicPainGroup  , levels=c("No Sites" , "1-2 Sites", "3-4 Sites", "5-7 Sites", "Widespread"))


## Create temp df with no NAs in chronic pain status and glycA columns
UKB_CP_MDD_glycA_temp <- UKB_CP_MDD_glycA[!is.na(UKB_CP_MDD_glycA$ChronicPainStatus) & !is.na(UKB_CP_MDD_glycA$glycA_1), ]

## Plot glycA levels in chronic pain status
ggplot(UKB_CP_MDD_glycA_temp, aes(x=as.factor(ChronicPainStatus), y= glycA_1, group = ChronicPainStatus)) + 
  geom_boxplot(fill = c("#ffffbf", "#d7191c")) +
  theme_classic(base_size = 20) +
  xlab("Chronic Pain Group") + ylab("GlycA (mmol/l)") +
  scale_x_discrete(labels=c("No Chronic Pain","Chronic Pain")) +
  ggtitle("GlycA in Chronic Pain")

## Create temp df with no NAs in chronic pain status and glycA columns
UKB_CP_MDD_glycA_temp <- UKB_CP_MDD_glycA[!is.na(UKB_CP_MDD_glycA$ChronicPainGroup) & !is.na(UKB_CP_MDD_glycA$glycA_1), ]

## Plot glycA levels in chronic pain site groups
ggplot(UKB_CP_MDD_glycA_temp, aes(x=ChronicPainGroup, y= glycA_1, group = ChronicPainGroup, fill = ChronicPainGroup)) + 
  geom_boxplot(fill = c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000")) +
  theme_classic(base_size = 20) +
  labs(title="GlycA in Chronic Pain Groups",
       x ="Chronic Pain Group", y = "GlycA (mmol/l)")


## Plot GlycA in recurrent MDD status groups
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create temp df with no NAs in chronic pain status and glycA columns
UKB_CP_MDD_glycA_temp <- UKB_CP_MDD_glycA[!is.na(UKB_CP_MDD_glycA$RecurrentMDDStatus) & !is.na(UKB_CP_MDD_glycA$glycA_1), ]


ggplot(UKB_CP_MDD_glycA_temp, aes(x=as.factor(RecurrentMDDStatus), y= glycA_1, group = RecurrentMDDStatus)) + 
  geom_boxplot(fill = c("#ffffbf", "#d7191c")) +
  theme_classic(base_size = 20)  +
  labs(title = "GlycA Levels in Depression", x = "Probable Recurrent MDD Status", y = "GlycA (mmol/l)")

## `Plot GlycA in chronic pain depression comorbidity
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
UKB_CP_MDD_glycA$ComorbidStatus <- factor(UKB_CP_MDD_glycA$ComorbidStatus  , levels=c("No Probable Recurrent MDD + No Chronic Pain",  "Probable Recurrent MDD + No Chronic Pain", "No Probable Recurrent MDD + Chronic Pain",
                                                                                      "Probable Recurrent MDD + Chronic Pain"))
## Create temp df with no NAs in chronic pain status and glycA columns
UKB_CP_MDD_glycA_temp <- UKB_CP_MDD_glycA[!is.na(UKB_CP_MDD_glycA$ComorbidStatus) & !is.na(UKB_CP_MDD_glycA$glycA_1), ]

ggplot(data=UKB_CP_MDD_glycA_temp) +
  (aes(x=as.factor(ComorbidStatus), y=glycA_1)) +
  geom_boxplot(fill = c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33")) +
  theme_minimal(base_size = 15) +
  labs(title = "GlycA levels in Comorbid Chronic Pain and Depression", x = "Comorbidity Group", y = "GlycA (mmol/l)") +
  scale_x_discrete(labels=c("No Probable Recurrent MDD + No Chronic Pain" = "CP-MDD-",
                            "Probable Recurrent MDD + No Chronic Pain" = "CP-MDD+",
                            "No Probable Recurrent MDD + Chronic Pain" = "CP+MDD-",
                            "Probable Recurrent MDD + Chronic Pain" = "CP+MDD+"))

## Statistical analysis: (x ~ covariate + group) - CRP ~ Chronic Pain Status 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create temp df with no NAs in chronic pain status and glycA columns
UKB_CP_MDD_glycA_temp <- UKB_CP_MDD_glycA[!is.na(UKB_CP_MDD_glycA$ChronicPainStatus) & !is.na(UKB_CP_MDD_glycA$glycA_1), ]

## Run linear regression
UKB_CP_MDD_glycA_temp$ChronicPainStatus <- as.factor(UKB_CP_MDD_glycA_temp$ChronicPainStatus)

## run first model: adjusted for BMI, sex, age
model_1 <- lm(glycA_1 ~ ChronicPainStatus + BMI_cat + Sex + Age, data = UKB_CP_MDD_glycA_temp)
summary(model_1)

## run second model: adjusted for all covariates
model_2 <- lm(glycA_1 ~ ChronicPainStatus + BMI_cat + Sex + Alcohol + Smoking + Age + Deprivation_rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, data = UKB_CP_MDD_glycA_temp)
summary(model_2)

layout(matrix(c(1,2,3,4),2,2))
plot(model_1)
plot(model_2)

## Statistical analysis: (x ~ covariate + group) - GlycA ~ Chronic Pain Group
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create temp df with no NAs in chronic pain status and glycA columns
UKB_CP_MDD_glycA_temp <- UKB_CP_MDD_glycA[!is.na(UKB_CP_MDD_glycA$ChronicPainGroup) & !is.na(UKB_CP_MDD_glycA$glycA_1), ]

## run first model: adjusted for BMI, sex, age
anova_cp_staus_glycA_1 <- aov(glycA_1 ~ ChronicPainGroup + BMI_cat + Sex + Age, data = UKB_CP_MDD_glycA_temp)
summary(anova_cp_staus_glycA_1)

## Post-hoc analysis
TukeyHSD(anova_cp_staus_glycA_1, "ChronicPainGroup", conf.level = 0.95)

## run second model: adjusted for all covariates
anova_cp_staus_glycA_2 <- aov(glycA_1 ~ ChronicPainGroup + BMI_cat + Sex + Alcohol + Smoking + Age + Deprivation_rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, data = UKB_CP_MDD_glycA_temp)
summary(anova_cp_staus_glycA_2)

## Post-hoc analysis
TukeyHSD(anova_cp_staus_glycA_2, "ChronicPainGroup", conf.level = 0.95)


## Statistical analysis: (x ~ covariate + group) - CRP ~ Depression Status 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create temp df with no NAs in chronic pain status and glycA columns
UKB_CP_MDD_glycA_temp <- UKB_CP_MDD_glycA[!is.na(UKB_CP_MDD_glycA$RecurrentMDDStatus) & !is.na(UKB_CP_MDD_glycA$glycA_1), ]

## Run linear regression
UKB_CP_MDD_glycA_temp$RecurrentMDDStatus <- as.factor(UKB_CP_MDD_glycA_temp$RecurrentMDDStatus)

## run first model: adjusted for BMI, sex, age
model_1 <- lm(glycA_1 ~ RecurrentMDDStatus + BMI_cat + Sex + Age, data = UKB_CP_MDD_glycA_temp)
summary(model_1)

## run second model: adjusted for all covariates
model_2 <- lm(glycA_1 ~ RecurrentMDDStatus + BMI_cat + Sex + Alcohol + Smoking + Age + Deprivation_rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, data = UKB_CP_MDD_glycA_temp)
summary(model_2)

layout(matrix(c(1,2,3,4),2,2))
plot(model_1)
plot(model_2)


## Statistical analysis: (x ~ covariate + group) - GlycA ~ Comorbidity Group
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create temp df with no NAs in comorbidity group and glycA columns
UKB_CP_MDD_glycA_temp <- UKB_CP_MDD_glycA[!is.na(UKB_CP_MDD_glycA$ComorbidStatus) & !is.na(UKB_CP_MDD_glycA$glycA_1), ]

## run first model: adjusted for BMI, sex, age
anova_cp_staus_glycA_1 <- aov(glycA_1 ~ ComorbidStatus + BMI_cat + Sex + Age, data = UKB_CP_MDD_glycA_temp)
summary(anova_cp_staus_glycA_1)

## Post-hoc analysis
TukeyHSD(anova_cp_staus_glycA_1, "ComorbidStatus", conf.level = 0.95)

## run second model: adjusted for all covariates
anova_cp_staus_glycA_2 <- aov(glycA_1 ~ ComorbidStatus + BMI_cat + Sex + Alcohol + Smoking + Age + Deprivation_rank + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, data = UKB_CP_MDD_glycA_temp)
summary(anova_cp_staus_glycA_2)

## Post-hoc analysis
TukeyHSD(anova_cp_staus_glycA_2, "ComorbidStatus", conf.level = 0.95)