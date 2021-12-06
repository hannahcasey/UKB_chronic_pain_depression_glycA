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



## Load in physical measures data
UKBPhysicalMeasures <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2018-10-phenotypes-ukb24262/PhysicalMeasures.rds")
## Extract ID column and columns pretaining to BMI measurment, waiste measurment and hip measurment
UKBPhysicalMeasuresBMIWH <- cbind(UKBPhysicalMeasures[,1], UKBPhysicalMeasures[grepl("f.21001.", names(UKBPhysicalMeasures))],
                                  UKBPhysicalMeasures[grepl("f.48.", names(UKBPhysicalMeasures))],
                                  UKBPhysicalMeasures[grepl("f.49.", names(UKBPhysicalMeasures))])
## Remove redundant dataframes 
rm(UKBPhysicalMeasures)
## Rename ID column
names(UKBPhysicalMeasuresBMIWH)[names(UKBPhysicalMeasuresBMIWH) == "UKBPhysicalMeasures[, 1]"] <- "n_eid"

## Load in Natasha's Phenotype of MDD
UKBNatashaMDDPhenotype <- read.csv("~/Desktop/PhD/output/UKB/depression/UKBBNatashaPhenotypeDepressionMania.csv", header = T)
## Rename ID column
colnames(UKBNatashaMDDPhenotype)[which(names(UKBNatashaMDDPhenotype) == "f_eid")] <- "n_eid"

## Add extra column indicating depression status
UKBNatashaMDDPhenotype$DepressionStatus <- NA
UKBNatashaMDDPhenotype$DepressionStatus[UKBNatashaMDDPhenotype$patient_group == "recurrent depression"] <- "Probable Recurrent MDD"
UKBNatashaMDDPhenotype$DepressionStatus[UKBNatashaMDDPhenotype$patient_group == "single depression"] <- "Probable Single MDD"
UKBNatashaMDDPhenotype$DepressionStatus[UKBNatashaMDDPhenotype$patient_group %in% c("control", "bipolar depression", "unipolar mania")] <- "No Probable MDD"

## Load in chronic pain data (Keira Johnson Phenotype)
UKBChronicPainJohnsonEligible <- read.csv("~/Desktop/PhD/output/UKB/pain/UKBBChronicPainJohnsonEligible.csv")

## Identify chronic pain individuals (at least one site of chronic pain)
UKBChronicPainJohnsonEligible$ChronicPainStatus <- NA
UKBChronicPainJohnsonEligible$ChronicPainStatus[UKBChronicPainJohnsonEligible$sites > 0] <- 1
UKBChronicPainJohnsonEligible$ChronicPainStatus[UKBChronicPainJohnsonEligible$sites < 1] <- 0

## Remove chronic pain status of those with widespread chronic pain
UKBChronicPainJohnsonEligible$ChronicPainStatus[UKBChronicPainJohnsonEligible$sites == 8] <- NA


##Quantify chronic pain sites
UKBChronicPainJohnsonEligible$ChronicPainGroup<- NA
UKBChronicPainJohnsonEligible$ChronicPainGroup[UKBChronicPainJohnsonEligible$sites == 0] <- "No Sites"
UKBChronicPainJohnsonEligible$ChronicPainGroup[UKBChronicPainJohnsonEligible$sites > 0 & UKBChronicPainJohnsonEligible$sites < 3] <- "1-2 Sites"
UKBChronicPainJohnsonEligible$ChronicPainGroup[UKBChronicPainJohnsonEligible$sites > 2 & UKBChronicPainJohnsonEligible$sites < 5] <- "3-4 Sites"
UKBChronicPainJohnsonEligible$ChronicPainGroup[UKBChronicPainJohnsonEligible$sites > 4 & UKBChronicPainJohnsonEligible$sites < 8] <- "5-7 Sites"
UKBChronicPainJohnsonEligible$ChronicPainGroup[UKBChronicPainJohnsonEligible$sites == 8] <- "Widespread"

## Load in gender data
UKBBaselineStatisics <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2020-09-phenotypes-ukb43743/BaselineCharacteristics.rds")

## f.31 = Sex
## Extract ID column and columns pretaining to CRP measurment
UKBBaselineStatisicsSex <- cbind(UKBBaselineStatisics[,1], UKBBaselineStatisics[grepl("f.31.", names(UKBBaselineStatisics))])

## Rename ID  ans Sex column
names(UKBBaselineStatisicsSex)[names(UKBBaselineStatisicsSex) == "UKBBaselineStatisics[, 1]"] <- "n_eid"
names(UKBBaselineStatisicsSex)[names(UKBBaselineStatisicsSex) == "f.31.0.0"] <- "Sex"

## Remove redundant dataframes 
rm(UKBBaselineStatisics)




