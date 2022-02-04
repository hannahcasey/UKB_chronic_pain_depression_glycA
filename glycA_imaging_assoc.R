## Script overview:
## General coritcal volumes
## Regional cortical volumes
## Individual FA tracts
## Individual MD tracts
## General FA measures
## General MD Measures
## REgional subcortical volumes


## load libraries
library(dplyr)
library(ggplot2)
library(patchwork)
library(nlme)
library(stringr)


## Save Shen's Function
run_model <- function(ls.mod,mod.dat_short,mod.dat_long){
  # define vars
  mod.dep = as.character(ls.mod[1])
  mod.factor = as.character(ls.mod[2])
  mod.covs = as.character(ls.mod[3])
  mod.type = as.character(ls.mod[4])
  
  # run model
  if (mod.type=='lme'){
    # model.expression
    fh_r=1:(nrow(mod.dat_long)/2)
    sh_r=(nrow(mod.dat_long)/2+1):nrow(mod.dat_long)
    
    mod.dat_long[fh_r,mod.dep]=scale(mod.dat_long[fh_r,mod.dep])
    mod.dat_long[sh_r,mod.dep]=scale(mod.dat_long[sh_r,mod.dep])
    if(is.numeric(mod.dat_long[,mod.factor])){
      mod.dat_long[fh_r,mod.factor]=scale(mod.dat_long[fh_r,mod.factor])
      mod.dat_long[sh_r,mod.factor]=scale(mod.dat_long[sh_r,mod.factor])
    }
    
    mod=paste0(mod.dep,'~',mod.covs,'+',mod.factor)
    
    fit=lme(as.formula(as.character(mod)),data=mod.dat_long,
            na.action=na.exclude,random=~1|f.eid,control=lmeControl(opt = "optim"))
    tmp.ci = intervals(fit,which='fixed')$fixed %>% as.data.frame %>% 
      select(Lower_95CI=lower,Upper_95CI=upper) %>% 
      tail(1)
    tmp.res = summary(fit)$tTable %>% 
      as.data.frame %>% 
      select(beta=Value,std=Std.Error,t.value=`t-value`,p.value=`p-value`) %>% 
      tail(1) %>% 
      mutate(mod_name = paste0(mod.dep,'~',mod.factor),tmp.ci) %>% 
      select(mod_name, everything())
    
    
  }else{
    dep.dat=mod.dat_short[,mod.dep]            
    if (length(table(dep.dat))==2){
      mod=paste0(mod.dep,'~',mod.covs,'+scale(',mod.factor,')')
      fit=glm(as.formula(as.character(mod)),data=mod.dat_short,na.action=na.exclude,family = 'binomial')
    }else{
      mod=paste0('scale(',mod.dep,')~',mod.covs,'+scale(',mod.factor,')')
      fit=glm(as.formula(as.character(mod)),data=mod.dat_short,na.action=na.exclude)
    }            
    tmp.ci = confint(fit) %>% as.data.frame %>% 
      select(Lower_95CI=`2.5 %`,Upper_95CI=`97.5 %`) %>% 
      tail(1)
    tmp.res = summary(fit)$coefficients %>% 
      as.data.frame %>% 
      select(beta=Estimate,std=`Std. Error`,t.value=`t value`,p.value=`Pr(>|t|)`) %>% 
      tail(1)%>% 
      mutate(mod_name = paste0(mod.dep,'~',mod.factor),tmp.ci) %>% 
      select(mod_name, everything())
  }
  
  return(tmp.res)
}
## Load and tidy data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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



## Load in phenotypic serum CRP measures
## Load in assay data
UKB_assay_results <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2020-09-phenotypes-ukb43743/AssayResults.rds")
## f.30710 = C-reactive protein
## Extract ID column and columns pretaining to CRP measurment
UKB_assay_results_filtered <- cbind(UKB_assay_results[,1], UKB_assay_results[grepl("f.30710.", names(UKB_assay_results))])

## rename FID and CRP columns
UKB_assay_results_filtered <- UKB_assay_results_filtered %>%
  rename(f.eid = "UKB_assay_results[, 1]",
         CRP = f.30710.0.0)


## Load in covariate information:
## f.21022 = Age at recruitment -
## f.31 = Sex -
## f.21001 = BMI -
## f.54 = assessment centre -
## Volume of EstimatedTotalIntraCranial (whole brain)

## Load in baseline characteristics - sex
UKB_baseline_characteristics <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2020-09-phenotypes-ukb43743/BaselineCharacteristics.rds")

## Extract sex status
UKB_baseline_characteristics_filtered <- cbind(UKB_baseline_characteristics[,1],
                                               UKB_baseline_characteristics[grepl("f.31.0.0", names(UKB_baseline_characteristics))],
                                               UKB_baseline_characteristics[grepl("f.21022.0.0", names(UKB_baseline_characteristics))])

## Rename ID and sex column in new filtered dataframe
UKB_baseline_characteristics_filtered <- UKB_baseline_characteristics_filtered %>%
  rename(f.eid = "UKB_baseline_characteristics[, 1]",
         sex = f.31.0.0,
         age = f.21022.0.0)

## Create new variable, age^2
UKB_baseline_characteristics_filtered$age_squared <- (UKB_baseline_characteristics_filtered$age^2)

## Load in physical measures data - BMI
UKB_physical_measures <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-04-phenotypes-ukb44797/PhysicalMeasures.rds")

## Extract BMI data
UKB_physical_measures_filtered <- cbind(UKB_physical_measures[,1], UKB_physical_measures[grepl("f.21001.", names(UKB_physical_measures))])

## Rename ID and BMI column in new filtered dataframe
UKB_physical_measures_filtered <- UKB_physical_measures_filtered %>%
  rename(f.eid = "UKB_physical_measures[, 1]",
         BMI = f.21001.0.0)

## Load in recruitment data - assessment centre
UKB_recruitment <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-04-phenotypes-ukb44797/Recruitment.rds")
UKB_recruitment_filtered <- cbind(UKB_recruitment[,1],
                                  UKB_recruitment[grepl("f.21003.", names(UKB_recruitment))],
                                  UKB_recruitment[grepl("f.54.", names(UKB_recruitment))])

## Rename ID and 2nd instance assessment centre column in new filtered dataframe
UKB_recruitment_filtered <- UKB_recruitment_filtered %>%
  rename(f.eid = "UKB_recruitment[, 1]",
         assessment_centre_first_imaging = f.54.2.0)



## Load in structural imaging data
UKB_structural_MRI <- readRDS("~/Desktop/PhD/projects/UKBCRPImagingPRS/resources/Imaging.rds")

## Extract total ICV data
UKB_structural_MRI_filtered <- cbind(UKB_structural_MRI[,1],
                                     UKB_structural_MRI[grepl("f.25008.", names(UKB_structural_MRI))],
                                     UKB_structural_MRI[grepl("f.25006.", names(UKB_structural_MRI))],
                                     UKB_structural_MRI[grepl("f.25004.", names(UKB_structural_MRI))])

## Add white matter, grey matter and cerebrospinal fluid volume to calcuate total ICV
UKB_structural_MRI_filtered$ICV <- UKB_structural_MRI_filtered$f.25008.2.0 + UKB_structural_MRI_filtered$f.25006.2.0 + UKB_structural_MRI_filtered$f.25004.2.0


## Rename ID in new filtered dataframe
UKB_structural_MRI_filtered <- UKB_structural_MRI_filtered %>%
  rename(f.eid = "UKB_structural_MRI[, 1]")

## Extract Freesurfer DKT
rx_DK <- number_range(26721, 26922)
rx_FA <- number_range(25488, 25514)
rx_MA <- number_range(25515, 25541)
rx_FA_dMRI <- number_range(25056, 25103)
rx_MD_dMRI <- number_range(25104, 25151)
rx_subcort <- number_range(25011, 25024)

UKB_structural_MRI_DK <- cbind(UKB_structural_MRI[,1],
                               UKB_structural_MRI[grepl(rx_DK, names(UKB_structural_MRI))],
                               UKB_structural_MRI[grepl(rx_FA, names(UKB_structural_MRI))],
                               UKB_structural_MRI[grepl(rx_MA, names(UKB_structural_MRI))],
                               UKB_structural_MRI[grepl(rx_FA_dMRI, names(UKB_structural_MRI))],
                               UKB_structural_MRI[grepl(rx_MD_dMRI, names(UKB_structural_MRI))],
                               UKB_structural_MRI[grepl(rx_subcort, names(UKB_structural_MRI))])

## Rename ID column in new filtered dataframe
UKB_structural_MRI_DK <- UKB_structural_MRI_DK %>%
  rename(f.eid = "UKB_structural_MRI[, 1]")

## remove redudnet dataframes
rm(UKB_baseline_characteristics)
rm(UKB_recruitment)
rm(UKB_physical_measures)

## Combine ID, covariates, CRP PRS and DKT structural measures into single dataframe
UKB_imaging_glycA <- inner_join(UKB_NMR_glycA, UKB_baseline_characteristics_filtered, by = "f.eid")
UKB_imaging_glycA <- inner_join(UKB_imaging_glycA,UKB_assay_results_filtered, by = "f.eid")
UKB_imaging_glycA <- inner_join(UKB_imaging_glycA, UKB_recruitment_filtered, by = "f.eid")
UKB_imaging_glycA <- inner_join(UKB_imaging_glycA, UKB_physical_measures_filtered, by = "f.eid")
UKB_imaging_glycA <- inner_join(UKB_imaging_glycA, UKB_structural_MRI_filtered, by = "f.eid")
UKB_imaging_glycA <- inner_join(UKB_imaging_glycA, UKB_structural_MRI_DK, by = "f.eid")


## Load in UKB DKW field ID key file
UKB_imaging_DKW_key <- read.csv("~/Desktop/PhD/projects/UKBCRPImagingPRS/resources/UKB_DKW_field_IDs.csv", header = F)
UKB_imaging_DKW_key <- UKB_imaging_DKW_key %>%
  rename(feild_ID = V1,
         cortical_volume = V2)

## Load in UKB FA field ID key file
UKB_imaging_FA_key <- read.csv("~/Desktop/PhD/projects/UKBCRPImagingPRS/resources/UKB_FA_field_IDs.csv", header = F)
UKB_imaging_FA_key <- UKB_imaging_FA_key %>%
  rename(feild_ID = V1,
         FA_value = V2)


## Load in UKB FA field ID key file
UKB_imaging_MD_key <- read.csv("~/Desktop/PhD/projects/UKBCRPImagingPRS/resources/UKB_MD_field_IDs.csv", header = F)
UKB_imaging_MD_key <- UKB_imaging_MD_key %>%
  rename(feild_ID = V1,
         MD_value = V2)


## Load in UKB DKW field ID key file
UKB_subcortical_key <- read.csv("~/Desktop/PhD/projects/UKBCRPImagingPRS/resources/UKB_subcortical_field_IDs.csv", header = F)
UKB_subcortical_key <- UKB_subcortical_key %>%
  rename(feild_ID = V1,
         subcortical_volume = V2)

## Calculate lobar measures

## Frontal lobe FIDs:
## Superior frontal gyrus = 26815 + 26916
## rostralmiddlefrontal = 26814 + 26915
## Caudalmiddlefrontal = 26791 + 26892
## Pars orbitalis = 26806 + 26907
## pars triangularis = 26807 + 26908
## pars opercularis = 26805 + 26906
## Frontal pole = 26819 + 26920
## Lateral orbitofrontal = 26799 + 26900
## medial orbitofrontal = 26801 + 26902
## Precentral gyrus = 26811 + 26912
## paracentral cortex = 26804 + 26905

## Sum volume of individual structures in frontal lobe
UKB_imaging_glycA$frontal_lobe <- UKB_imaging_glycA$f.26815.2.0 + UKB_imaging_glycA$f.26916.2.0 + UKB_imaging_glycA$f.26814.2.0 +
  UKB_imaging_glycA$f.26915.2.0 + UKB_imaging_glycA$f.26791.2.0 + UKB_imaging_glycA$f.26892.2.0 + UKB_imaging_glycA$f.26806.2.0 +
  UKB_imaging_glycA$f.26907.2.0 + UKB_imaging_glycA$f.26807.2.0 + UKB_imaging_glycA$f.26908.2.0 + UKB_imaging_glycA$f.26805.2.0 + 
  UKB_imaging_glycA$f.26906.2.0 + UKB_imaging_glycA$f.26819.2.0 + UKB_imaging_glycA$f.26920.2.0 + UKB_imaging_glycA$f.26799.2.0 +
  UKB_imaging_glycA$f.26900.2.0 + UKB_imaging_glycA$f.26801.2.0 + UKB_imaging_glycA$f.26902.2.0 + UKB_imaging_glycA$f.26811.2.0 + 
  UKB_imaging_glycA$f.26912.2.0 + UKB_imaging_glycA$f.26804.2.0 + UKB_imaging_glycA$f.26905.2.0

mean(UKB_imaging_glycA$frontal_lobe, na.rm = T)

## Temporal lobe FIDs:
## Insula= 26821 + 26922
## Superior temporal = 26817 + 26918
## transverse temporal = 26820 + 26921
## banks STS = 	26789 + 26890
## Middle temporal gyrus = 26802 + 26903
## Inferior temporal gyrus = 26796 + 26897
## Fusiform = 26794 + 26895
## parahippocampal = 26803 + 26904
## entorhinal = 26793 + 26894

## Sum volume of individual structures in temporal lobe
UKB_imaging_glycA$temporal_lobe <- UKB_imaging_glycA$f.26821.2.0 + UKB_imaging_glycA$f.26922.2.0 + UKB_imaging_glycA$f.26817.2.0 +
  UKB_imaging_glycA$f.26918.2.0 + UKB_imaging_glycA$f.26820.2.0 + UKB_imaging_glycA$f.26921.2.0 + UKB_imaging_glycA$f.26789.2.0 +
  UKB_imaging_glycA$f.26890.2.0 + UKB_imaging_glycA$f.26802.2.0 + UKB_imaging_glycA$f.26903.2.0 + UKB_imaging_glycA$f.26796.2.0 +
  UKB_imaging_glycA$f.26897.2.0 + UKB_imaging_glycA$f.26794.2.0 + UKB_imaging_glycA$f.26895.2.0 + UKB_imaging_glycA$f.26803.2.0 +
  UKB_imaging_glycA$f.26904.2.0 + UKB_imaging_glycA$f.26793.2.0 + UKB_imaging_glycA$f.26894.2.0

mean(UKB_imaging_glycA$temporal_lobe, na.rm = T)

## Parietal lobe:
## Postcentral gyrus = 26809 + 26910
## paracentral cortex = 26804 + 26905
## Superior parietal cortex = 26816 + 26917
## Inferior parietal cortex = 26795 + 26896
## Supramarginal gyrus = 26818 + 26919
## Precuneus = 26812 + 26913

## Sum volume of individual structures in parietal lobe
UKB_imaging_glycA$parietal_lobe <- UKB_imaging_glycA$f.26809.2.0 + UKB_imaging_glycA$f.26910.2.0 + UKB_imaging_glycA$f.26804.2.0 +
  UKB_imaging_glycA$f.26905.2.0 + UKB_imaging_glycA$f.26816.2.0 + UKB_imaging_glycA$f.26917.2.0 + UKB_imaging_glycA$f.26795.2.0 +
  UKB_imaging_glycA$f.26896.2.0 + UKB_imaging_glycA$f.26818.2.0 + UKB_imaging_glycA$f.26919.2.0 + UKB_imaging_glycA$f.26812.2.0 +
  UKB_imaging_glycA$f.26913.2.0

mean(UKB_imaging_glycA$parietal_lobe, na.rm = T)

## Occipital lobe:
## Lateral occipital cortex = 26798, 26899
## Cuneus = 26792, 26893
## Pericalcarine cortex = 26808, 26909
## Lingual gyrus = 26800, 26901

UKB_imaging_glycA$occipital_lobe <- UKB_imaging_glycA$f.26798.2.0 + UKB_imaging_glycA$f.26899.2.0 + UKB_imaging_glycA$f.26792.2.0 +
  UKB_imaging_glycA$f.26893.2.0 + UKB_imaging_glycA$f.26808.2.0 + UKB_imaging_glycA$f.26909.2.0 + UKB_imaging_glycA$f.26800.2.0 +
  UKB_imaging_glycA$f.26901.2.0

mean(UKB_imaging_glycA$occipital_lobe, na.rm = T)

## Cingulate lobe:
## Rostral ACC = 26813, 26914
## Caudal ACC = 26790, 26891
## Posterior cingulate cortex = 26810 + 26911
## Cingulate isthmus = 26797 + 26898

UKB_imaging_glycA$cingulate_lobe <- UKB_imaging_glycA$f.26813.2.0 + UKB_imaging_glycA$f.26914.2.0 + UKB_imaging_glycA$f.26790.2.0 +
  UKB_imaging_glycA$f.26891.2.0 + UKB_imaging_glycA$f.26810.2.0 + UKB_imaging_glycA$f.26911.2.0 + UKB_imaging_glycA$f.26797.2.0 +
  UKB_imaging_glycA$f.26898.2.0

mean(UKB_imaging_glycA$cingulate_lobe, na.rm = T)

## Sum lobar volumes to get global cortical volume
UKB_imaging_glycA$global_cortical_volume <- UKB_imaging_glycA$frontal_lobe + UKB_imaging_glycA$temporal_lobe + UKB_imaging_glycA$parietal_lobe +
  UKB_imaging_glycA$occipital_lobe +UKB_imaging_glycA$cingulate_lobe 

mean(UKB_imaging_glycA$global_cortical_volume, na.rm = T)

## Save dataframe with CRP PRS, serum CRP, covariates and strucutral data
#write.csv(UKB_imaging_glycA, "~/Desktop/PhD/projects/UKBCRPImagingPRS/resources/UKB_imaging_glycA.csv", quote = F)

UKB_imaging_glycA <- read.csv("~/Desktop/PhD/projects/UKBCRPImagingPRS/resources/UKB_imaging_glycA.csv", header = T)

## Standardise all cortical global and lobar volumes
UKB_imaging_glycA$global_cortical_volume <- scale(UKB_imaging_glycA$global_cortical_volume)
UKB_imaging_glycA$frontal_lobe <- scale(UKB_imaging_glycA$frontal_lobe)
UKB_imaging_glycA$temporal_lobe <- scale(UKB_imaging_glycA$temporal_lobe)
UKB_imaging_glycA$parietal_lobe <- scale(UKB_imaging_glycA$parietal_lobe)
UKB_imaging_glycA$occipital_lobe <- scale(UKB_imaging_glycA$occipital_lobe)
UKB_imaging_glycA$cingulate_lobe <- scale(UKB_imaging_glycA$cingulate_lobe)

## Remove outliers - IQR method
outliers <- boxplot(UKB_imaging_glycA$global_cortical_volume, plot=FALSE)$out
UKB_imaging_glycA_global_cortical_volume_outliers_removed <- UKB_imaging_glycA
UKB_imaging_glycA_global_cortical_volume_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA$global_cortical_volume %in% outliers),]

## Statistical analysis
## Create dataframe to store glm output
glm_glycA_cortical <- data.frame(cortical_volume=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                   p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric())


## CRP PRS ~ global corical volume + covariates
glm1 <- glm(global_cortical_volume  ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV,
            data = UKB_imaging_glycA_global_cortical_volume_outliers_removed)

summary(glm1)

glm_glycA_cortical[1,"cortical_volume"] <- "global_cortical_volume"
glm_glycA_cortical[1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_cortical[1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_cortical[1, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
glm_glycA_cortical[1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]


cortical_lobes <- c("frontal_lobe", "temporal_lobe", "parietal_lobe", "occipital_lobe", "cingulate_lobe")

for (i in 1:5){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_imaging_glycA[,cortical_lobes[i]], plot=FALSE)$out
  UKB_imaging_glycA_lobar_outliers_removed <- UKB_imaging_glycA
  UKB_imaging_glycA_lobar_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA[,cortical_lobes[i]] %in% outliers),]
  
  
  mod <- paste0(cortical_lobes[i], "~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age + age_squared + ICV")
  glm1 <- glm(as.formula(mod),
              data = UKB_imaging_glycA_lobar_outliers_removed)
  
  glm_glycA_cortical[i+1,"cortical_volume"] <- cortical_lobes[i]
  glm_glycA_cortical[i+1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
  glm_glycA_cortical[i+1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
  glm_glycA_cortical[i+1, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
  glm_glycA_cortical[i+1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]
  glm_glycA_cortical[i+1, "p.adjust"] <- p.adjust(coef(summary(glm1))["glycA_1",4], method = "fdr", n = 5)
}

## Plot results
# lock in factor level order
glm_glycA_cortical$cortical_volume  = with(glm_glycA_cortical, reorder(cortical_volume, beta))

p1 <- ggplot(data=glm_glycA_cortical[!glm_glycA_cortical$cortical_volume=="global_cortical_volume",], aes(x=cortical_volume, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Mean (95% CI)") +
  theme_bw() + # use a white background
  ylim(-0.3,0.3) +
  labs(title ="Association of GlycA and General Cortical Volumes") + 
  scale_x_discrete(labels=c("frontal_lobe" = "Frontal Lobe", "temporal_lobe" = "Temporal Lobe","parietal_lobe" = "Parietal Lobe", "occipital_lobe" = "Occipital Lobe", "cingulate_lobe" = "Cingulate Lobe"))

p2 <- ggplot(data=glm_glycA_cortical[glm_glycA_cortical$cortical_volume=="global_cortical_volume",], aes(x=cortical_volume, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Mean (95% CI)") +
  theme_bw() +  # use a white background
  ylim(-0.12,0.12) +
  labs(title ="Association of GlycA and Global Cortical Volume") +
  scale_x_discrete(labels=c("global_cortical_volume" = "Global"))


layout <- c(
  area(t = 1, l = 1, b = 3, r = 4),
  area(t = 4, l = 1, b = 4, r = 4))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_general_cortical_structure.jpg")
p1 / p2 + 
  plot_layout(design = layout)
dev.off()

write.csv(glm_glycA_cortical, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_cortical_general.csv",
          quote = FALSE, row.names = FALSE)

## Association analysis of individual cortical volumes

## Run interaction analysis to see if there is an effect modification between glycA and hemisphere on cortical volumes
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## list of column names for each cortical volume
cortical_volume_FIDs <- paste0("f.", c(26789:26821, 26890:26922), ".2.0")

## Create dataframe to store output
interaction_glycA_hemisphere <- data.frame(Volume=character(),  interaction_p=numeric(), interaction_p_adjust=numeric())

for (i in 1:33){
  
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_imaging_glycA[,cortical_volume_FIDs[i]], plot=FALSE)$out
  UKB_imaging_glycA_cortical_regions_outliers_removed <- UKB_imaging_glycA
  UKB_imaging_glycA_cortical_regions_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA[,cortical_volume_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_imaging_glycA[,cortical_volume_FIDs[i +33]], plot=FALSE)$out
  UKB_imaging_glycA_cortical_regions_outliers_removed <- UKB_imaging_glycA_cortical_regions_outliers_removed[-which(UKB_imaging_glycA_cortical_regions_outliers_removed[,cortical_volume_FIDs[i+33]] %in% outliers),]
  
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_imaging_glycA_cortical_regions_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                            cortical_volume_FIDs[i], cortical_volume_FIDs[i +33])]
  
  ## Convert short dataframe to long format
  UKB_glycA_DK_small_long <- reshape(UKB_glycA_DK_small, 
                                       varying = c(cortical_volume_FIDs[i], cortical_volume_FIDs[i + 33]), 
                                       v.names = "Volume",
                                       timevar = "Region", 
                                       times = c(cortical_volume_FIDs[i], cortical_volume_FIDs[i + 33]), 
                                       direction = "long")
  
  ## Add column indicating hemisphere
  UKB_glycA_DK_small_long$Hemisphere <- NA
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_volume_FIDs[i]] <- "left"
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_volume_FIDs[i + 33]] <- "right"
  
  ## Carry out interaction analysis
  lm1 <- lm(Volume ~ glycA_1 * Hemisphere + sex + age + age + age_squared + BMI + ICV + assessment_centre_first_imaging, data = UKB_glycA_DK_small_long)
  
  interaction_glycA_hemisphere[i,"Volume"] <- (cortical_volume_FIDs)[i]
  interaction_glycA_hemisphere[i,"interaction_p"] <- summary(lm1)[["coefficients"]]["glycA_1:Hemisphereright", "Pr(>|t|)"]
  interaction_glycA_hemisphere[i,"interaction_p_adjust"] <- p.adjust(interaction_glycA_hemisphere[i,"interaction_p"], n = 33, method = "fdr")
}
## No interaction found between glycA_1 and hemisphere on volumes

## Run linear mixed effect model analysis to look at the effect of CRP PRS on repeat measures (hemispheres) of corical volumes
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate dataframe to record P-values, CIs and beta coefficients
lme_cortical_regions_glycA <- data.frame(mod_name=character(), beta=numeric(), std=numeric(), t.value=numeric(), p.value=numeric(),
                                           Lower_95CI=numeric(), Upper_95CI=numeric())
## Generate ls.mod
ls.mod.PRS <- data.frame(Dep ="Volume",
                         Factor =  "glycA_1",
                         Covariates = "sex + age + age_squared + BMI + ICV + assessment_centre_first_imaging + Hemisphere",
                         Model = "lme")

## list of column names for each cortical volume
cortical_volume_FIDs <- paste0("f.", c(26789:26821, 26890:26922), ".2.0")

## Iterate through each cortical volume and run lme 
for (i in 1:33){
  
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_imaging_glycA[,cortical_volume_FIDs[i]], plot=FALSE)$out
  UKB_imaging_glycA_cortical_regions_outliers_removed <- UKB_imaging_glycA
  UKB_imaging_glycA_cortical_regions_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA[,cortical_volume_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_imaging_glycA[,cortical_volume_FIDs[i +33]], plot=FALSE)$out
  UKB_imaging_glycA_cortical_regions_outliers_removed <- UKB_imaging_glycA_cortical_regions_outliers_removed[-which(UKB_imaging_glycA_cortical_regions_outliers_removed[,cortical_volume_FIDs[i+33]] %in% outliers),]
  
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_imaging_glycA_cortical_regions_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                            cortical_volume_FIDs[i], cortical_volume_FIDs[i +33])]
  
  ## Convert short dataframe to long format
  UKB_glycA_DK_small_long <- reshape(UKB_glycA_DK_small, 
                                       varying = c(cortical_volume_FIDs[i], cortical_volume_FIDs[i + 33]), 
                                       v.names = "Volume",
                                       timevar = "Region", 
                                       times = c(cortical_volume_FIDs[i], cortical_volume_FIDs[i + 33]), 
                                       direction = "long")
  
  ## Add column indicating hemisphere
  UKB_glycA_DK_small_long$Hemisphere <- NA
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_volume_FIDs[i]] <- "left"
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_volume_FIDs[i + 33]] <- "right"
  
  ## Generate ls.mode
  
  ## Rum lme and get output (BTEA, SD, T, P, 95% CIs)
  output <- run_model(ls.mod.PRS, UKB_glycA_DK_small, UKB_glycA_DK_small_long)
  ## Add region field ID
  lme_cortical_region_iterate_glycA <- cbind(data_frame(region_field_ID = paste0(cortical_volume_FIDs[i],"_",cortical_volume_FIDs[i + 33])), output)
  ## Adjust P value for multiple comparisons (FDR)
  lme_cortical_region_iterate_glycA$p.value.adjust <- p.adjust(lme_cortical_region_iterate_glycA$p.value, n = 33, method = "fdr")
  ## Add cortical volume name
  cortical_volume <- UKB_imaging_DKW_key$cortical_volume[grepl(substring(lme_cortical_region_iterate_glycA$region_field_ID, 3,7), UKB_imaging_DKW_key$feild_ID)]
  cortical_volume_sub <- sub(("\\(.*"), "", cortical_volume)
  lme_cortical_region_iterate_glycA$cortical_volume <- str_to_title(sub(("Volume of "), "", cortical_volume_sub))
  ## Append to data frame
  lme_cortical_regions_glycA <- rbind(lme_cortical_region_iterate_glycA, lme_cortical_regions_glycA)
}


## Plot results
lme_cortical_regions_glycA$cortical_volume = with(lme_cortical_regions_glycA, reorder(cortical_volume, beta))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_regional_cortical_structures.jpg")

ggplot(data=lme_cortical_regions_glycA, aes(x=cortical_volume, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Label") + ylab("Mean (95% CI)") +
  theme_bw() +
  labs(title = "Association of GlycA and Regional Cortical Volumes")

dev.off()

write.csv(lme_cortical_regions_glycA, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_cortical_regions.csv",
          quote = FALSE, row.names = FALSE)


## Run interaction analysis to see if there is an effect modification between glycAand hemisphere on FA values
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## list of column names for each cortical volume
FA_bilateral_FIDs <- paste0("f.", c(25488:25497, 25500:25503, 25505:25514), ".2.0")
FA_unilateral_FIDs <- paste0("f.", c(25498, 25499, 25504), ".2.0")

## Create dataframe to store output
interaction_glycA_hemisphere_FA <- data.frame(Volume=character(),  interaction_p=numeric(), interaction_p_adjust=numeric())

for (i in seq(1, 23, by=2)){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_imaging_glycA[,FA_bilateral_FIDs[i]], plot=FALSE)$out
  UKB_imaging_glycA_FA_tracts_outliers_removed <- UKB_imaging_glycA
  UKB_imaging_glycA_FA_tracts_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA[,FA_bilateral_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_imaging_glycA[,FA_bilateral_FIDs[i +1]], plot=FALSE)$out
  UKB_imaging_glycA_FA_tracts_outliers_removed <- UKB_imaging_glycA_FA_tracts_outliers_removed[-which(UKB_imaging_glycA_FA_tracts_outliers_removed[,FA_bilateral_FIDs[i+1]] %in% outliers),]
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_imaging_glycA_FA_tracts_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                             FA_bilateral_FIDs[i], FA_bilateral_FIDs[i +1])]
  ## Convert short dataframe to long format
  UKB_glycA_DK_small_long <- reshape(UKB_glycA_DK_small, 
                                     varying = c(FA_bilateral_FIDs[i], FA_bilateral_FIDs[i + 1]), 
                                     v.names = "FA_measure",
                                     timevar = "Region", 
                                     times = c(FA_bilateral_FIDs[i], FA_bilateral_FIDs[i + 1]), 
                                     direction = "long")
  
  ## Add column indicating hemisphere
  UKB_glycA_DK_small_long$Hemisphere <- NA
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == FA_bilateral_FIDs[i]] <- "left"
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == FA_bilateral_FIDs[i + 1]] <- "right"
  
  ## Carry out interaction analysis
  lm1 <- lm(FA_measure ~ glycA_1 * Hemisphere + sex + age + age_squared + BMI + ICV + assessment_centre_first_imaging, data = UKB_glycA_DK_small_long)
  
  interaction_glycA_FA <- data_frame(FA_measure=NA,  interaction_p=NA, interaction_p_adjust=NA)
  interaction_glycA_FA$FA_measure <- FA_bilateral_FIDs[i]
  interaction_glycA_FA$interaction_p <- summary(lm1)[["coefficients"]]["glycA_1:Hemisphereright", "Pr(>|t|)"]
  interaction_glycA_FA$interaction_p_adjust <- p.adjust(interaction_glycA_FA$interaction_p, n = 15, method = "fdr")
  
  interaction_glycA_hemisphere_FA <- rbind(interaction_glycA_hemisphere_FA, interaction_glycA_FA)
}


## Run linear mixed effect model analysis to look at the effect of glycA on repeat measures (hemispheres) of FA values
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate dataframe to record P-values, CIs and beta coefficients
lme_FA_values_glycA <- data.frame(mod_name=character(), beta=numeric(), std=numeric(), t.value=numeric(), p.value=numeric(),
                                  Lower_95CI=numeric(), Upper_95CI=numeric())
## Generate ls.mod
ls.mod.glycA <- data.frame(Dep ="FA",
                           Factor =  "glycA_1",
                           Covariates = "sex + age + age_squared + BMI + ICV + assessment_centre_first_imaging + Hemisphere",
                           Model = "lme")

## Iterate through each bilateral FA tract and run lme 
## Skip first two FA measurements as significant hemisphere interaction was found
for (i in seq(1, 23, by=2)){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_imaging_glycA[,FA_bilateral_FIDs[i]], plot=FALSE)$out
  UKB_imaging_glycA_FA_tracts_outliers_removed <- UKB_imaging_glycA
  UKB_imaging_glycA_FA_tracts_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA[,FA_bilateral_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_imaging_glycA[,FA_bilateral_FIDs[i +1]], plot=FALSE)$out
  UKB_imaging_glycA_FA_tracts_outliers_removed <- UKB_imaging_glycA_FA_tracts_outliers_removed[-which(UKB_imaging_glycA_FA_tracts_outliers_removed[,FA_bilateral_FIDs[i+1]] %in% outliers),]
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_imaging_glycA_FA_tracts_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                             FA_bilateral_FIDs[i + 1], FA_bilateral_FIDs[i])]
  
  ## Convert short dataframe to long format
  UKB_glycA_DK_small_long <- reshape(UKB_glycA_DK_small, 
                                     varying = c(FA_bilateral_FIDs[i], FA_bilateral_FIDs[i + 1]), 
                                     v.names = "FA",
                                     timevar = "Region", 
                                     times = c(FA_bilateral_FIDs[i], FA_bilateral_FIDs[i + 1]), 
                                     direction = "long")
  
  ## Add column indicating hemisphere
  UKB_glycA_DK_small_long$Hemisphere <- NA
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == FA_bilateral_FIDs[i]] <- "left"
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == FA_bilateral_FIDs[i + 1]] <- "right"
  
  ## Generate ls.mode
  
  ## Rum lme and get output (BTEA, SD, T, P, 95% CIs)
  output <- run_model(ls.mod.glycA, UKB_glycA_DK_small, UKB_glycA_DK_small_long)
  ## Add region field ID
  lme_FA_meaures_glycA <- cbind(data_frame(region_field_ID = paste0(FA_bilateral_FIDs[i],"_",FA_bilateral_FIDs[i + 1])), output)
  ## Adjust P value for multiple comparisons (FDR)
  lme_FA_meaures_glycA$p.value.adjust <- p.adjust(lme_FA_meaures_glycA $p.value, n = 15, method = "fdr")
  ## Add cortical volume name
  FA_value <- UKB_imaging_FA_key$FA_value[grepl(substring(lme_FA_meaures_glycA $region_field_ID, 3,7), UKB_imaging_FA_key$feild_ID)]
  FA_value_sub <- sub(("\\(.*"), "", FA_value)
  lme_FA_meaures_glycA$FA_tract <- str_to_title(sub(("Weighted-mean FA in tract "), "", FA_value_sub))
  ## Append to data frame
  lme_FA_values_glycA <- rbind(lme_FA_meaures_glycA, lme_FA_values_glycA)
}


# ## Carry out GLM on each hemispheric measurement of FA in tract acoustic radiation
# for (i in 1:2){
#   
#   ## Scale FA value in unilateral tracts
#   UKB_imaging_glycA[,FA_bilateral_FIDs[i]] <- scale(UKB_imaging_glycA[,FA_bilateral_FIDs[i]])
#   ## Run glm
#   mod <- paste0("glycA_1~ ",FA_bilateral_FIDs[i], "+ sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV")
#   glm1 <- glm(as.formula(mod), data = UKB_imaging_glycA)
#   
#   
#   lme_FA_values_glycA[10 + i,"region_field_ID"] <- FA_bilateral_FIDs[i]
#   lme_FA_values_glycA[10 + i, "beta"] <- summary(glm1)[["coefficients"]][FA_bilateral_FIDs[i], "Estimate"]
#   lme_FA_values_glycA[10 + i, "std"] <- summary(glm1)[["coefficients"]][FA_bilateral_FIDs[i], "Std. Error"]
#   lme_FA_values_glycA[10 + i, "p.value"] <- summary(glm1)[["coefficients"]][FA_bilateral_FIDs[i], "Pr(>|t|)"]
#   lme_FA_values_glycA[10 + i, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)[FA_bilateral_FIDs[i],]
#   lme_FA_values_glycA[10 + i, "p.value.adjust"] <- p.adjust(lme_FA_values_glycA[10 + i, "p.value"], method = "fdr", n = 13)
#   
#   FA_value <- UKB_imaging_FA_key$FA_value[grepl(substring(lme_FA_values_glycA$region_field_ID[10 + i], 3,7), UKB_imaging_FA_key$feild_ID)]
#   lme_FA_values_glycA$FA_tract[10 + i] <- str_to_title(sub(("Weighted-mean FA in tract "), "", FA_value))
#   
#   
# }


## Carry out GLM on unilateral FA tracts

## Create new dataframe to store glm output
lme_FA_values_glycA_unilateral <- data.frame(region_field_ID=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                             Lower_95CI=numeric(), Upper_95CI=numeric()) 

for (i in 1:3){
  
  ## Scale FA value in unilateral tracts
  UKB_imaging_glycA[,FA_unilateral_FIDs[i]] <- scale(UKB_imaging_glycA[,FA_unilateral_FIDs[i]])
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_imaging_glycA[,FA_unilateral_FIDs[i]], plot=FALSE)$out
  UKB_imaging_glycA_FA_unilateral_outliers_removed <- UKB_imaging_glycA
  UKB_imaging_glycA_FA_unilateral_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA[,FA_unilateral_FIDs[i]] %in% outliers),]
  
  ## Run glm
  mod <- paste0("glycA_1~ ",FA_unilateral_FIDs[i], "+ sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV")
  glm1 <- glm(as.formula(mod), data = UKB_imaging_glycA)
  
  
  lme_FA_values_glycA_unilateral[i,"region_field_ID"] <- FA_unilateral_FIDs[i]
  lme_FA_values_glycA_unilateral[i, "beta"] <- summary(glm1)[["coefficients"]][FA_unilateral_FIDs[i], "Estimate"]
  lme_FA_values_glycA_unilateral[i, "std"] <- summary(glm1)[["coefficients"]][FA_unilateral_FIDs[i], "Std. Error"]
  lme_FA_values_glycA_unilateral[i, "p.value"] <- summary(glm1)[["coefficients"]][FA_unilateral_FIDs[i], "Pr(>|t|)"]
  lme_FA_values_glycA_unilateral[i, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)[FA_unilateral_FIDs[i],]
  lme_FA_values_glycA_unilateral[i, "p.value.adjust"] <- p.adjust(lme_FA_values_glycA_unilateral[i, "p.value"], method = "fdr", n = 15)
  
  FA_value <- UKB_imaging_FA_key$FA_value[grepl(substring(lme_FA_values_glycA_unilateral$region_field_ID[i], 3,7), UKB_imaging_FA_key$feild_ID)]
  lme_FA_values_glycA_unilateral$FA_tract[i] <- str_to_title(sub(("Weighted-mean FA in tract "), "", FA_value))
}

## Plot results

## Plot effect of CRP PRS on bilateral FA tracts
lme_FA_values_glycA$FA_tract = with(lme_FA_values_glycA, reorder(FA_tract, beta))

p1 <- ggplot(data=lme_FA_values_glycA, aes(x=FA_tract, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Bilateral FA Tracts") + ylab("Mean (95% CI)") +
  theme_bw() +
  ylim(-0.06, 0.06) +
  labs(title ="Association of GlycA and FA Bilateral Tracts")


write.csv(lme_FA_values_glycA, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_FA_bilateral_tracts.csv",
          quote = FALSE, row.names = FALSE)

## Plot effect of CRP PRS on unilateral FA tracts
lme_FA_values_glycA_unilateral$FA_tract = with(lme_FA_values_glycA_unilateral, reorder(FA_tract, beta))

p2 <- ggplot(data=lme_FA_values_glycA_unilateral, aes(x=FA_tract, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Unilateral FA Tracts") + ylab("Mean (95% CI)") +
  theme_bw() +
  ylim(-5e-03, 5e-03) +
  labs(title ="Association of GlycA and FA Unilateral Tracts")

layout <- c(
  area(t = 1, l = 1, b = 3, r = 4),
  area(t = 4, l = 1, b = 4, r = 4))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_FA_tracts.jpg")

p1 / p2 + 
  plot_layout(design = layout)

dev.off()

write.csv(lme_FA_values_glycA_unilateral, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_FA_unilateral_tracts.csv",
          quote = FALSE, row.names = FALSE)

## Run interaction analysis to see if there is an effect modification between glycA and hemisphere on MD values
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## list of column names for each cortical volume
MD_bilateral_FIDs <- paste0("f.", c(25515:25524, 25527:25530, 25532:25541), ".2.0")
MD_unilateral_FIDs <- paste0("f.", c(25525, 25526, 25531), ".2.0")

## Create dataframe to store output
interaction_glycA_hemisphere_MD <- data.frame(MD_measurement=character(),  interaction_p=numeric(), interaction_p_adjust=numeric())

for (i in seq(1, 23, by=2)){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_imaging_glycA[,MD_bilateral_FIDs[i]], plot=FALSE)$out
  UKB_imaging_glycA_MD_tracts_outliers_removed <- UKB_imaging_glycA
  UKB_imaging_glycA_MD_tracts_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA[,MD_bilateral_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_imaging_glycA[,MD_bilateral_FIDs[i +1]], plot=FALSE)$out
  UKB_imaging_glycA_MD_tracts_outliers_removed <- UKB_imaging_glycA_MD_tracts_outliers_removed[-which(UKB_imaging_glycA_MD_tracts_outliers_removed[,MD_bilateral_FIDs[i+1]] %in% outliers),]
  
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_imaging_glycA_MD_tracts_outliers_removed[,c("f.eid", "glycA_1",  "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                             MD_bilateral_FIDs[i] , MD_bilateral_FIDs[i +1])]
  
  
  ## Convert short dataframe to long format
  UKB_glycA_DK_small_long <- reshape(UKB_glycA_DK_small, 
                                     varying = c(MD_bilateral_FIDs[i], MD_bilateral_FIDs[i + 1]), 
                                     v.names = "MD_measurement",
                                     timevar = "Region", 
                                     times = c(MD_bilateral_FIDs[i], MD_bilateral_FIDs[i + 1]), 
                                     direction = "long")
  
  ## Add column indicating hemisphere
  UKB_glycA_DK_small_long$Hemisphere <- NA
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == MD_bilateral_FIDs[i]] <- "left"
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == MD_bilateral_FIDs[i + 1]] <- "right"
  
  ## Carry out interaction analysis
  lm1 <- lm(MD_measurement ~ glycA_1 * Hemisphere + sex + age + age_squared + BMI + ICV + assessment_centre_first_imaging, data = UKB_glycA_DK_small_long)
  
  interaction_glycA_MD <- data_frame(MD_measurement=NA,  interaction_p=NA, interaction_p_adjust=NA)
  interaction_glycA_MD$MD_measurement <- MD_bilateral_FIDs[i]
  interaction_glycA_MD$interaction_p <- summary(lm1)[["coefficients"]]["glycA_1:Hemisphereright", "Pr(>|t|)"]
  interaction_glycA_MD$interaction_p_adjust <- p.adjust(interaction_glycA_MD$interaction_p, n = 15, method = "fdr")
  
  interaction_glycA_hemisphere_MD <- rbind(interaction_glycA_hemisphere_MD, interaction_glycA_MD)
}
## No interaction found between glycA_1 and hemisphere on MD

## Run linear mixed effect model analysis to look at the effect of glycA on repeat measures (hemispheres) of MD values
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate dataframe to record P-values, CIs and beta coefficients
lme_MD_values_glycA <- data.frame(mod_name=character(), beta=numeric(), std=numeric(), t.value=numeric(), p.value=numeric(),
                                  Lower_95CI=numeric(), Upper_95CI=numeric())
## Generate ls.mod
ls.mod.glycA <- data.frame(Dep ="MD",
                           Factor =  "glycA_1",
                           Covariates = "sex + age + age_squared + BMI + ICV + assessment_centre_first_imaging + Hemisphere",
                           Model = "lme")

## Iterate through each bilateral MD tract and run lme 
for (i in seq(1, 23, by=2)){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_imaging_glycA[,MD_bilateral_FIDs[i]], plot=FALSE)$out
  UKB_imaging_glycA_MD_tracts_outliers_removed <- UKB_imaging_glycA
  UKB_imaging_glycA_MD_tracts_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA[,MD_bilateral_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_imaging_glycA[,MD_bilateral_FIDs[i +1]], plot=FALSE)$out
  UKB_imaging_glycA_MD_tracts_outliers_removed <- UKB_imaging_glycA_MD_tracts_outliers_removed[-which(UKB_imaging_glycA_MD_tracts_outliers_removed[,MD_bilateral_FIDs[i+1]] %in% outliers),]
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_imaging_glycA_MD_tracts_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                             MD_bilateral_FIDs[i + 1], MD_bilateral_FIDs[i])]
  
  ## Convert short dataframe to long format
  UKB_glycA_DK_small_long <- reshape(UKB_glycA_DK_small, 
                                     varying = c(MD_bilateral_FIDs[i], MD_bilateral_FIDs[i + 1]), 
                                     v.names = "MD",
                                     timevar = "Region", 
                                     times = c(MD_bilateral_FIDs[i], MD_bilateral_FIDs[i + 1]), 
                                     direction = "long")
  
  ## Add column indicating hemisphere
  UKB_glycA_DK_small_long$Hemisphere <- NA
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == MD_bilateral_FIDs[i]] <- "left"
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == MD_bilateral_FIDs[i + 1]] <- "right"
  
  ## Generate ls.mode
  
  ## Rum lme and get output (BTEA, SD, T, P, 95% CIs)
  output <- run_model(ls.mod.glycA, UKB_glycA_DK_small, UKB_glycA_DK_small_long)
  ## Add region field ID
  lme_MD_meaures_glycA <- cbind(data_frame(region_field_ID = paste0(MD_bilateral_FIDs[i],"_",MD_bilateral_FIDs[i + 1])), output)
  ## Adjust P value for multiple comparisons (FDR)
  lme_MD_meaures_glycA$p.value.adjust <- p.adjust(lme_MD_meaures_glycA$p.value, n = 12, method = "fdr")
  ## Add cortical volume name
  MD_value <- UKB_imaging_MD_key$MD_value[grepl(substring(lme_MD_meaures_glycA$region_field_ID, 3,7), UKB_imaging_MD_key$feild_ID)]
  MD_value_sub <- sub(("\\(.*"), "", MD_value)
  lme_MD_meaures_glycA$MD_tract <- str_to_title(sub(("Weighted-mean MD in tract "), "", MD_value_sub))
  ## Append to data frame
  lme_MD_values_glycA <- rbind(lme_MD_meaures_glycA, lme_MD_values_glycA)
}


## Carry out GLM on unilateral MD tracts

## Create new dataframe to store glm output
lme_MD_values_glycA_unilateral <- data.frame(region_field_ID=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                             Lower_95CI=numeric(), Upper_95CI=numeric()) 

for (i in 1:3){
  
  ## Scale MD value in unilateral tracts
  UKB_imaging_glycA[,MD_unilateral_FIDs[i]] <- scale(UKB_imaging_glycA[,MD_unilateral_FIDs[i]])
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_imaging_glycA[,MD_unilateral_FIDs[i]], plot=FALSE)$out
  UKB_imaging_glycA_MD_unilateral_outliers_removed <- UKB_imaging_glycA
  UKB_imaging_glycA_MD_unilateral_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA[,MD_unilateral_FIDs[i]] %in% outliers),]
  
  ## Run glm
  mod <- paste0(MD_unilateral_FIDs[i], " ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV")
  glm1 <- glm(as.formula(mod), data = UKB_imaging_glycA_MD_unilateral_outliers_removed)
  
  
  lme_MD_values_glycA_unilateral[i,"region_field_ID"] <- MD_unilateral_FIDs[i]
  lme_MD_values_glycA_unilateral[i, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
  lme_MD_values_glycA_unilateral[i, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
  lme_MD_values_glycA_unilateral[i, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
  lme_MD_values_glycA_unilateral[i, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]
  lme_MD_values_glycA_unilateral[i, "p.value.adjust"] <- p.adjust(lme_MD_values_glycA_unilateral[i, "p.value"], method = "fdr", n = 15)
  
  MD_value <- UKB_imaging_MD_key$MD_value[grepl(substring(lme_MD_values_glycA_unilateral$region_field_ID[i], 3,7), UKB_imaging_MD_key$feild_ID)]
  lme_MD_values_glycA_unilateral$MD_tract[i] <- str_to_title(sub(("Weighted-mean MD in tract "), "", MD_value))
}

## Plot results

## Plot effect of CRP PRS on bilateral FA tracts
lme_MD_values_glycA$MD_tract = with(lme_MD_values_glycA, reorder(MD_tract, beta))

p1 <- ggplot(data=lme_MD_values_glycA, aes(x=MD_tract, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Bilateral MD Tracts") + ylab("Mean (95% CI)") +
  theme_bw() +
  ylim(-0.07, 0.07) +
  labs(title ="Association of GlycA and MD Bilateral Tracts")

write.csv(lme_MD_values_glycA, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_MD_bilateral_tracts.csv",
          quote = FALSE, row.names = FALSE)

## Plot effect of CRP PRS on unilateral FA tracts
lme_MD_values_glycA_unilateral$MD_tract = with(lme_MD_values_glycA_unilateral, reorder(MD_tract, beta))

p2 <- ggplot(data=lme_MD_values_glycA_unilateral, aes(x=MD_tract, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Unilateral MD Tracts") + ylab("Mean (95% CI)") +
  theme_bw() +
  ylim(-0.4, 0.4) +
  labs(title ="Association of GlycA and MD Unilateral Tracts")

layout <- c(
  area(t = 1, l = 1, b = 3, r = 4),
  area(t = 4, l = 1, b = 4, r = 4))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_MD_tracts.jpg")

p1 / p2 + 
  plot_layout(design = layout)

dev.off()

write.csv(lme_MD_values_glycA_unilateral, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_MD_unilateral_tracts.csv",
          quote = FALSE, row.names = FALSE)


## Separate DTI into General DTI Measures:
## MD
## Association/Commissural fibres

# Cingulate gyrus part of cingulum (left) - f.25519
# Cingulate gyrus part of cingulum (right) - f.25520
# Inferior fronto-occipital fasciculus (left) - f.25527
# Inferior fronto-occipital fasciculus (right) - f.25528
# Inferior longitudinal fasciculus (left) - f.25529
# Inferior longitudinal fasciculus (right) - f.25530
# Parahippocampal part of cingulum (left) - f.25521
# Parahippocampal part of cingulum (right) - f.25522
# Superior longitudinal fasciculus (left) - f.25536
# Superior longitudinal fasciculus (right) - f.25537
# Uncinate fasciculus (left) - f.25540
# Uncinate fasciculus (right) - f.25541
# Forceps major - f.25525
# Forceps minor - f.25526


association_fibres_MD <- c("f.25519.2.0", "f.25520.2.0", "f.25527.2.0", "f.25528.2.0", "f.25529.2.0", "f.25530.2.0",
                           "f.25521.2.0", "f.25522.2.0", "f.25536.2.0", "f.25537.2.0", "f.25540.2.0", "f.25541.2.0",
                           "f.25525.2.0", "f.25526.2.0")


## Thalamic radiations

# Anterior thalamic radiation (left)  - f.25517
# Anterior thalamic radiation (right) - f.25518
# Posterior thalamic radiation (left) - f.25534
# Posterior thalamic radiation (right) - f.25535
# Superior thalamic radiation (left) - f.25538
# Superior thalamic radiation (right) - f.25539

thalamic_radiations_MD <- c("f.25517.2.0", "f.25518.2.0", "f.25534.2.0", "f.25535.2.0", "f.25538.2.0", "f.25539.2.0")


## Projection fibres

# Acoustic radiation (left) - f.25516
# Acoustic radiation (right) - f.25517
# Corticospinal tract (left) - f.25523
# Corticospinal tract (right) - f.25524
# Medial lemniscus (left) - f.25532
# Medial lemniscus (right) - f.25533
# Middle cerebellar peduncle - f.25531



projection_fibres_MD <- c("f.25516.2.0", "f.25517.2.0", "f.25523.2.0", "f.25524.2.0", "f.25532.2.0", "f.25533.2.0", "f.25531.2.0")


## Global_MD
global_MD <- c(association_fibres_MD, thalamic_radiations_MD, projection_fibres_MD)


## Create df to store glm output 
glm_glycA_MD_general <- data.frame(general_DTI_MD_measure=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                     p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric())

#### GLOBAL MD 
## Remove participants with NA in any of the MD tracts needed to calculate global MD measures
UKB_imaging_glycA_global_MD <- UKB_imaging_glycA[rowSums(is.na(UKB_imaging_glycA[,global_MD]))!=length(global_MD),]
## Perform PCA
PC_res_MD_global <- prcomp(UKB_imaging_glycA_global_MD[,global_MD],scale = TRUE)
## Add 1st PC to df
UKB_imaging_glycA_global_MD <- cbind(UKB_imaging_glycA_global_MD, (PC_res_MD_global$x[,1]))
## Rename 1st PC column
UKB_imaging_glycA_global_MD_PC <- UKB_imaging_glycA_global_MD %>%
  rename(global_MD= "(PC_res_MD_global$x[, 1])")

## Remove outliers - IQR method
outliers <- boxplot(UKB_imaging_glycA_global_MD_PC$global_MD, plot=FALSE)$out
UKB_imaging_glycA_global_MD_PC_outliers_removed <- UKB_imaging_glycA_global_MD_PC
UKB_imaging_glycA_global_MD_PC_outliers_removed <- UKB_imaging_glycA_global_MD_PC[-which(UKB_imaging_glycA_global_MD_PC$global_MD %in% outliers),]

## Carry out glm
glm1 <- glm(global_MD ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV,
            data = UKB_imaging_glycA_global_MD_PC_outliers_removed)
summary(glm1)

glm_glycA_MD_general[1,"general_DTI_MD_measure"] <- "global_MD"
glm_glycA_MD_general[1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_MD_general[1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_MD_general[1, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
glm_glycA_MD_general[1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]


## Create data frames without NAs present in MD tracts needed to calculate general MD measurents
UKB_glycA_association_fibres_MD <- UKB_imaging_glycA[rowSums(is.na(UKB_imaging_glycA[,association_fibres_MD]))!=length(association_fibres_MD),]
UKB_glycA_thalamic_radiations_MD<- UKB_imaging_glycA[rowSums(is.na(UKB_imaging_glycA[,thalamic_radiations_MD]))!=length(thalamic_radiations_MD),]
UKB_glycA_projection_fibres_MD <- UKB_imaging_glycA[rowSums(is.na(UKB_imaging_glycA[,projection_fibres_MD]))!=length(projection_fibres_MD),]

## Run PCA in each general MD measure
PC_res_MD_association_fibres <- prcomp(UKB_glycA_association_fibres_MD[,association_fibres_MD],scale = TRUE)
PC_res_MD_thalamic_radiations <- prcomp(UKB_glycA_thalamic_radiations_MD[,thalamic_radiations_MD],scale = TRUE)
PC_res_MD_projection_fibres <- prcomp(UKB_glycA_projection_fibres_MD[,projection_fibres_MD],scale = TRUE)

## Add 1st PC to df
UKB_glycA_assoc_MD_PC <- cbind(UKB_glycA_association_fibres_MD, PC_res_MD_association_fibres$x[,1])
UKB_glycA_tr_MD_PC <- cbind(UKB_glycA_thalamic_radiations_MD, PC_res_MD_thalamic_radiations$x[,1])
UKB_glycA_projec_MD_PC <- cbind(UKB_glycA_projection_fibres_MD, PC_res_MD_projection_fibres$x[,1])

## Rename 1st PC column
UKB_glycA_assoc_MD_PC <- UKB_glycA_assoc_MD_PC %>%
  rename(association_fibres_MD= "PC_res_MD_association_fibres$x[, 1]")
UKB_glycA_tr_MD_PC <- UKB_glycA_tr_MD_PC %>%
  rename(thalamic_radiations_MD= "PC_res_MD_thalamic_radiations$x[, 1]")
UKB_glycA_projec_MD_PC <- UKB_glycA_projec_MD_PC %>%
  rename(projection_fibres_MD= "PC_res_MD_projection_fibres$x[, 1]")

## Remove outliers - IQR method
outliers <- boxplot(UKB_glycA_assoc_MD_PC$association_fibres_MD, plot=FALSE)$out
UKB_glycA_assoc_MD_PC_outliers_removed <- UKB_glycA_assoc_MD_PC
UKB_glycA_assoc_MD_PC_outliers_removed <- UKB_glycA_assoc_MD_PC[-which(UKB_glycA_assoc_MD_PC$association_fibres_MD %in% outliers),]

outliers <- boxplot(UKB_glycA_tr_MD_PC$thalamic_radiations_MD, plot=FALSE)$out
UKB_glycA_tr_MD_PC_outliers_removed <- UKB_glycA_tr_MD_PC
UKB_glycA_tr_MD_PC_outliers_removed <- UKB_glycA_tr_MD_PC[-which(UKB_glycA_tr_MD_PC$thalamic_radiations_MD %in% outliers),]

outliers <- boxplot(UKB_glycA_projec_MD_PC$projection_fibres_MD, plot=FALSE)$out
UKB_glycA_projec_MD_PC_outliers_removed <- UKB_glycA_projec_MD_PC
UKB_glycA_projec_MD_PC_outliers_removed <- UKB_glycA_projec_MD_PC[-which(UKB_glycA_projec_MD_PC$projection_fibres_MD %in% outliers),]


## run GLM
glm2 <- glm(association_fibres_MD ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_glycA_assoc_MD_PC_outliers_removed)
summary(glm2)

glm3 <- glm(thalamic_radiations_MD ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_glycA_tr_MD_PC_outliers_removed)
summary(glm3)

glm4 <- glm(projection_fibres_MD ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_glycA_projec_MD_PC_outliers_removed)
summary(glm4)


glm_glycA_MD_general[2,"general_DTI_MD_measure"] <- "association_fibres_MD"
glm_glycA_MD_general[2, "beta"] <- summary(glm2)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_MD_general[2, "std"] <- summary(glm2)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_MD_general[2, "p.value"] <- summary(glm2)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
glm_glycA_MD_general[2, c("Lower_95CI", "Upper_95CI")] <- confint(glm2)["glycA_1",]
glm_glycA_MD_general[2, "p.adjust"] <- p.adjust(coef(summary(glm2))["glycA_1",4], method = "fdr", n = 3)


glm_glycA_MD_general[3,"general_DTI_MD_measure"] <- "thalamic_radiations_MD"
glm_glycA_MD_general[3, "beta"] <- summary(glm3)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_MD_general[3, "std"] <- summary(glm3)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_MD_general[3, "p.value"] <- summary(glm3)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
glm_glycA_MD_general[3, c("Lower_95CI", "Upper_95CI")] <- confint(glm3)["glycA_1",]
glm_glycA_MD_general[3, "p.adjust"] <- p.adjust(coef(summary(glm3))["glycA_1",4], method = "fdr", n = 3)


glm_glycA_MD_general[4,"general_DTI_MD_measure"] <- "projection_fibres_MD"
glm_glycA_MD_general[4, "beta"] <- summary(glm4)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_MD_general[4, "std"] <- summary(glm4)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_MD_general[4, "p.value"] <- summary(glm4)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
glm_glycA_MD_general[4, c("Lower_95CI", "Upper_95CI")] <- confint(glm4)["glycA_1",]
glm_glycA_MD_general[4, "p.adjust"] <- p.adjust(coef(summary(glm4))["glycA_1",4], method = "fdr", n = 3)


## Plot results
# lock in factor level order
glm_glycA_MD_general$general_DTI_MD_measure  = with(glm_glycA_MD_general, reorder(general_DTI_MD_measure, beta))

p1 <- ggplot(data=glm_glycA_MD_general[!glm_glycA_MD_general$general_DTI_MD_measure=="global_MD",], aes(x=general_DTI_MD_measure, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Mean (95% CI)") +
  theme_bw() + # use a white background
  ylim(-0.7,0.7) +
  labs(title ="Association of General MD Measures and GlycA") +
  scale_x_discrete(labels=c("association_fibres_MD" = "Association Fibres", "thalamic_radiations_MD" = "Thalamic Radiations",
                            "projection_fibres_MD" = "Projection Fibres"))

p2 <- ggplot(data=glm_glycA_MD_general[glm_glycA_MD_general$general_DTI_MD_measure=="global_MD",], aes(x=general_DTI_MD_measure, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Mean (95% CI)") +
  theme_bw() +  # use a white background
  ylim(-0.7,0.7) +
  labs(title ="Association of Global MD and GlycA") +
  scale_x_discrete(labels=c("global_MD" = "Global"))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_MD_general.jpg")
p1/p2
dev.off()

write.csv(glm_glycA_MD_general, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_MD_general.csv",
          quote = FALSE, row.names = FALSE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Separate STI into General DTI Measures:
## FA
## Association/Commissural fibres

# Cingulate gyrus part of cingulum (left) - f.25492
# Cingulate gyrus part of cingulum (right) - f.25493
# Inferior fronto-occipital fasciculus (left) - f.25500
# Inferior fronto-occipital fasciculus (right) - f.25501
# Inferior longitudinal fasciculus (left) - f.25502
# Inferior longitudinal fasciculus (right) - f.25503
# Parahippocampal part of cingulum (left) - f.25494
# Parahippocampal part of cingulum (right) - f.25495
# Superior longitudinal fasciculus (left) - f.25509
# Superior longitudinal fasciculus (right) - f.25510
# Uncinate fasciculus (left) - f.25513
# Uncinate fasciculus (right) - f.25514
# Forceps major - f.25498
# Forceps minor - f.25499


association_fibres_FA <- c("f.25492.2.0", "f.25493.2.0", "f.25500.2.0", "f.25501.2.0", "f.25502.2.0", "f.25503.2.0",
                           "f.25494.2.0", "f.25495.2.0", "f.25509.2.0", "f.25510.2.0", "f.25513.2.0", "f.25514.2.0",
                           "f.25498.2.0", "f.25499.2.0")


## Thalamic radiations

# Anterior thalamic radiation (left)  - f.25490
# Anterior thalamic radiation (right) - f.25491
# Posterior thalamic radiation (left) - f.25507
# Posterior thalamic radiation (right) - f.25508
# Superior thalamic radiation (left) - f.25511
# Superior thalamic radiation (right) - f.25512

thalamic_radiations_FA <- c("f.25490.2.0", "f.25491.2.0", "f.25507.2.0", "f.25508.2.0", "f.25511.2.0", "f.25512.2.0")


## Projection fibres

# Acoustic radiation (left) - f.25488
# Acoustic radiation (right) - f.25489
# Corticospinal tract (left) - f.25496
# Corticospinal tract (right) - f.25497
# Medial lemniscus (left) - f.25505
# Medial lemniscus (right) - f.25506
# Middle cerebellar peduncle - f.25504



projection_fibres_FA <- c("f.25488.2.0", "f.25489.2.0", "f.25496.2.0", "f.25497.2.0", "f.25505.2.0", "f.25506.2.0", "f.25504.2.0")


## Global_FA
global_FA <- c(association_fibres_FA, thalamic_radiations_FA, projection_fibres_FA)


## Create df to store glm output 
glm_glycA_FA_general <- data.frame(general_DTI_FA_measure=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                   p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric())

#### GLOBAL FA 
## Remove participants with NA in any of the FA tracts needed to calculate global FA measures
UKB_imaging_glycA_global_FA <- UKB_imaging_glycA[rowSums(is.na(UKB_imaging_glycA[,global_FA]))!=length(global_FA),]
## Perform PCA
PC_res_FA_global <- prcomp(UKB_imaging_glycA_global_FA[,global_FA],scale = TRUE)
## Add 1st PC to df
UKB_imaging_glycA_global_FA <- cbind(UKB_imaging_glycA_global_FA, (PC_res_FA_global$x[,1] * -1))
## Rename 1st PC column
UKB_imaging_glycA_global_FA_PC <- UKB_imaging_glycA_global_FA %>%
  rename(global_FA= "(PC_res_FA_global$x[, 1] * -1)")

## Remove outliers - IQR method
outliers <- boxplot(UKB_imaging_glycA_global_FA_PC$global_FA, plot=FALSE)$out
UKB_imaging_glycA_global_FA_PC_outliers_removed <- UKB_imaging_glycA_global_FA_PC
UKB_imaging_glycA_global_FA_PC_outliers_removed <- UKB_imaging_glycA_global_FA_PC[-which(UKB_imaging_glycA_global_FA_PC$global_FA %in% outliers),]


## Carry out glm
glm1 <- glm(global_FA ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_imaging_glycA_global_FA_PC_outliers_removed)
summary(glm1)

glm_glycA_FA_general[1,"general_DTI_FA_measure"] <- "global_FA"
glm_glycA_FA_general[1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_FA_general[1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_FA_general[1, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
glm_glycA_FA_general[1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]


## Create data frames without NAs present in FA tracts needed to calculate general FA measurents
UKB_glycA_association_fibres_FA <- UKB_imaging_glycA[rowSums(is.na(UKB_imaging_glycA[,association_fibres_FA]))!=length(association_fibres_FA),]
UKB_glycA_thalamic_radiations_FA<- UKB_imaging_glycA[rowSums(is.na(UKB_imaging_glycA[,thalamic_radiations_FA]))!=length(thalamic_radiations_FA),]
UKB_glycA_projection_fibres_FA <- UKB_imaging_glycA[rowSums(is.na(UKB_imaging_glycA[,projection_fibres_FA]))!=length(projection_fibres_FA),]

## Run PCA in each general FA measure
PC_res_FA_association_fibres <- prcomp(UKB_glycA_association_fibres_FA[,association_fibres_FA],scale = TRUE)
PC_res_FA_thalamic_radiations <- prcomp(UKB_glycA_thalamic_radiations_FA[,thalamic_radiations_FA],scale = TRUE)
PC_res_FA_projection_fibres <- prcomp(UKB_glycA_projection_fibres_FA[,projection_fibres_FA],scale = TRUE)

## Add 1st PC to df
UKB_glycA_assoc_FA_PC <- cbind(UKB_glycA_association_fibres_FA, PC_res_FA_association_fibres$x[,1] * - 1)
UKB_glycA_tr_FA_PC <- cbind(UKB_glycA_thalamic_radiations_FA, PC_res_FA_thalamic_radiations$x[,1] * - 1)
UKB_glycA_projec_FA_PC <- cbind(UKB_glycA_projection_fibres_FA, PC_res_FA_projection_fibres$x[,1] * - 1)

## Rename 1st PC column
UKB_glycA_assoc_FA_PC <- UKB_glycA_assoc_FA_PC %>%
  rename(association_fibres_FA= "PC_res_FA_association_fibres$x[, 1] * -1")
UKB_glycA_tr_FA_PC <- UKB_glycA_tr_FA_PC %>%
  rename(thalamic_radiations_FA= "PC_res_FA_thalamic_radiations$x[, 1] * -1")
UKB_glycA_projec_FA_PC <- UKB_glycA_projec_FA_PC %>%
  rename(projection_fibres_FA= "PC_res_FA_projection_fibres$x[, 1] * -1")


## Remove outliers - IQR method
outliers <- boxplot(UKB_glycA_assoc_FA_PC$association_fibres_FA, plot=FALSE)$out
UKB_CRP_PRS_association_fibres_FA_PC_outliers_removed <- UKB_glycA_assoc_FA_PC
UKB_CRP_PRS_association_fibres_FA_PC_outliers_removed <- UKB_glycA_assoc_FA_PC[-which(UKB_glycA_assoc_FA_PC$association_fibres_FA %in% outliers),]

outliers <- boxplot(UKB_glycA_tr_FA_PC$thalamic_radiations_FA, plot=FALSE)$out
UKB_CRP_PRS_thalamic_radiations_FA_PC_outliers_removed <- UKB_glycA_tr_FA_PC
UKB_CRP_PRS_thalamic_radiations_FA_PC_outliers_removed <- UKB_glycA_tr_FA_PC[-which(UKB_glycA_tr_FA_PC$thalamic_radiations_FA %in% outliers),]

outliers <- boxplot(UKB_glycA_projec_FA_PC$projection_fibres_FA, plot=FALSE)$out
UKB_CRP_PRS_projection_fibres_FA_PC_outliers_removed <- UKB_glycA_projec_FA_PC
UKB_CRP_PRS_projection_fibres_FA_PC_outliers_removed <- UKB_glycA_projec_FA_PC[-which(UKB_glycA_projec_FA_PC$projection_fibres_FA %in% outliers),]


## run GLM
glm2 <- glm(association_fibres_FA ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_CRP_PRS_association_fibres_FA_PC_outliers_removed)
summary(glm2)

glm3 <- glm(thalamic_radiations_FA ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_CRP_PRS_thalamic_radiations_FA_PC_outliers_removed)
summary(glm3)

glm4 <- glm(projection_fibres_FA ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_CRP_PRS_projection_fibres_FA_PC_outliers_removed)
summary(glm4)


glm_glycA_FA_general[2,"general_DTI_FA_measure"] <- "association_fibres_FA"
glm_glycA_FA_general[2, "beta"] <- summary(glm2)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_FA_general[2, "std"] <- summary(glm2)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_FA_general[2, "p.value"] <- summary(glm2)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
glm_glycA_FA_general[2, c("Lower_95CI", "Upper_95CI")] <- confint(glm2)["glycA_1",]
glm_glycA_FA_general[2, "p.adjust"] <- p.adjust(coef(summary(glm2))["glycA_1",4], method = "fdr", n = 3)


glm_glycA_FA_general[3,"general_DTI_FA_measure"] <- "thalamic_radiations_FA"
glm_glycA_FA_general[3, "beta"] <- summary(glm3)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_FA_general[3, "std"] <- summary(glm3)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_FA_general[3, "p.value"] <- summary(glm3)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
glm_glycA_FA_general[3, c("Lower_95CI", "Upper_95CI")] <- confint(glm3)["glycA_1",]
glm_glycA_FA_general[3, "p.adjust"] <- p.adjust(coef(summary(glm3))["glycA_1",4], method = "fdr", n = 3)


glm_glycA_FA_general[4,"general_DTI_FA_measure"] <- "projection_fibres_FA"
glm_glycA_FA_general[4, "beta"] <- summary(glm4)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_FA_general[4, "std"] <- summary(glm4)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_FA_general[4, "p.value"] <- summary(glm4)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
glm_glycA_FA_general[4, c("Lower_95CI", "Upper_95CI")] <- confint(glm4)["glycA_1",]
glm_glycA_FA_general[4, "p.adjust"] <- p.adjust(coef(summary(glm4))["glycA_1",4], method = "fdr", n = 3)


## Plot results
# lock in factor level order
glm_glycA_FA_general$general_DTI_FA_measure  = with(glm_glycA_FA_general, reorder(general_DTI_FA_measure, beta))

p1 <- ggplot(data=glm_glycA_FA_general[!glm_glycA_FA_general$general_DTI_FA_measure=="global_FA",], aes(x=general_DTI_FA_measure, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Mean (95% CI)") +
  theme_bw() + # use a white background
  ylim(-0.9,0.9) +
  labs(title ="General FA Measures") +
  labs(title ="Association of General FA Measures and GlycA") +
  scale_x_discrete(labels=c("association_fibres_FA" = "Association Fibres", "thalamic_radiations_FA" = "Thalamic Radiations",
                            "projection_fibres_FA" = "Projection Fibres"))

p2 <- ggplot(data=glm_glycA_FA_general[glm_glycA_FA_general$general_DTI_FA_measure=="global_FA",], aes(x=general_DTI_FA_measure, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Mean (95% CI)") +
  theme_bw() +  # use a white background
  ylim(-1.2,1.2) +
  labs(title ="Global FA") +
  labs(title ="Association of Global FA and GlycA")  +
  scale_x_discrete(labels=c("global_FA" = "Global"))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_FA_general.jpg")
p1/p2
dev.off()

write.csv(glm_glycA_FA_general, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_FA_general.csv",
          quote = FALSE, row.names = FALSE)


## Run interaction analysis to see if there is an effect modification between GlycA and hemisphere on subcortical volumes
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create dataframe to store output
interaction_glycA_hemisphere <- data.frame(Volume=character(),  interaction_p=numeric(), interaction_p_adjust=numeric())

## list of column names for each cortical volume
sub_cortical_volume_FIDs <- paste0("f.", c(25011:25024), ".2.0")

for (i in seq(1,14,2)){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_imaging_glycA[,sub_cortical_volume_FIDs[i]], plot=FALSE)$out
  UKB_imaging_glycA_subcortical_regions_outliers_removed <- UKB_imaging_glycA
  UKB_imaging_glycA_subcortical_regions_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA[,sub_cortical_volume_FIDs[i]] %in% outliers),]
  
  outliers <- boxplot(UKB_imaging_glycA[,sub_cortical_volume_FIDs[i + 1]], plot=FALSE)$out
  UKB_imaging_glycA_subcortical_regions_outliers_removed <- UKB_imaging_glycA_subcortical_regions_outliers_removed[-which(UKB_imaging_glycA_subcortical_regions_outliers_removed[,sub_cortical_volume_FIDs[i + 1]] %in% outliers),]
  
  
  ## Select essential columns
  UKB_imaging_glycA_small <- UKB_imaging_glycA_subcortical_regions_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                                                              sub_cortical_volume_FIDs[i], sub_cortical_volume_FIDs[i + 1])]
  
  ## Convert short dataframe to long format
  UKB_imaging_glycA_small_long <- reshape(UKB_imaging_glycA_small, 
                                       varying = c(sub_cortical_volume_FIDs[i], sub_cortical_volume_FIDs[i + 1]), 
                                       v.names = "Volume",
                                       timevar = "Region", 
                                       times = c(sub_cortical_volume_FIDs[i], sub_cortical_volume_FIDs[i + 1]), 
                                       direction = "long")
  
  ## Add column indicating hemisphere
  UKB_imaging_glycA_small_long$Hemisphere <- NA
  UKB_imaging_glycA_small_long$Hemisphere[UKB_imaging_glycA_small_long$Region == sub_cortical_volume_FIDs[i]] <- "left"
  UKB_imaging_glycA_small_long$Hemisphere[UKB_imaging_glycA_small_long$Region == sub_cortical_volume_FIDs[i + 1]] <- "right"
  
  ## Carry out interaction analysis
  lm1 <- lm(Volume ~ glycA_1 * Hemisphere + sex + age + age_squared + BMI + ICV + assessment_centre_first_imaging, data = UKB_imaging_glycA_small_long)
  
  interaction_glycA_hemisphere[i,"Volume"] <- (sub_cortical_volume_FIDs)[i]
  interaction_glycA_hemisphere[i,"interaction_p"] <- summary(lm1)[["coefficients"]]["glycA_1:Hemisphereright", "Pr(>|t|)"]
  
  interaction_glycA_hemisphere[i,"interaction_p_adjust"] <- p.adjust(interaction_glycA_hemisphere[i,"interaction_p"], n = 7, method = "fdr")
}

interaction_glycA_hemisphere
## No interaction found between glycA and hemisphere on volumes

## Run linear mixed effect model analysis to look at the effect of GlycA on repeat measures (hemispheres) of subcorical volumes
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate dataframe to record P-values, CIs and beta coefficients
lme_subcortical_regions_GlycA<- data.frame(mod_name=character(), beta=numeric(), std=numeric(), t.value=numeric(), p.value=numeric(),
                                              Lower_95CI=numeric(), Upper_95CI=numeric())
## Generate ls.mod
ls.mod.PRS <- data.frame(Dep ="Volume",
                         Factor =  "glycA_1",
                         Covariates = "sex + age  + age_squared + BMI + ICV + assessment_centre_first_imaging + Hemisphere",
                         Model = "lme")


## Iterate through each cortical volume and run lme 
for (i in seq(1,14,2)){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_imaging_glycA[,sub_cortical_volume_FIDs[i]], plot=FALSE)$out
  UKB_imaging_glycA_subcortical_regions_outliers_removed <- UKB_imaging_glycA
  UKB_imaging_glycA_subcortical_regions_outliers_removed <- UKB_imaging_glycA[-which(UKB_imaging_glycA[,sub_cortical_volume_FIDs[i]] %in% outliers),]
  
  outliers <- boxplot(UKB_imaging_glycA[,sub_cortical_volume_FIDs[i + 1]], plot=FALSE)$out
  UKB_imaging_glycA_subcortical_regions_outliers_removed <- UKB_imaging_glycA_subcortical_regions_outliers_removed[-which(UKB_imaging_glycA_subcortical_regions_outliers_removed[,sub_cortical_volume_FIDs[i + 1]] %in% outliers),]
  
  
  ## Select essential columns
  UKB_imaging_glycA_small <- UKB_imaging_glycA_subcortical_regions_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                                                                       sub_cortical_volume_FIDs[i], sub_cortical_volume_FIDs[i + 1])]
  
  ## Convert short dataframe to long format
  UKB_imaging_glycA_small_long <- reshape(UKB_imaging_glycA_small, 
                                          varying = c(sub_cortical_volume_FIDs[i], sub_cortical_volume_FIDs[i + 1]), 
                                          v.names = "Volume",
                                          timevar = "Region", 
                                          times = c(sub_cortical_volume_FIDs[i], sub_cortical_volume_FIDs[i + 1]), 
                                          direction = "long")
  
  ## Add column indicating hemisphere
  UKB_imaging_glycA_small_long$Hemisphere <- NA
  UKB_imaging_glycA_small_long$Hemisphere[UKB_imaging_glycA_small_long$Region == sub_cortical_volume_FIDs[i]] <- "left"
  UKB_imaging_glycA_small_long$Hemisphere[UKB_imaging_glycA_small_long$Region == sub_cortical_volume_FIDs[i + 1]] <- "right"
  
  ## Generate ls.mode
  
  ## Rum lme and get output (BTEA, SD, T, P, 95% CIs)
  try(output <- run_model(ls.mod.PRS, UKB_imaging_glycA_small, UKB_imaging_glycA_small_long))
  ## Add region field ID
  lme_subcortical_regions <- cbind(data_frame(region_field_ID = paste0(sub_cortical_volume_FIDs[i],"_",sub_cortical_volume_FIDs[i + 1])), output)
  ## Adjust P value for multiple comparisons (FDR)
  lme_subcortical_regions$p.value.adjust <- p.adjust(lme_subcortical_regions$p.value, n = 7, method = "fdr")
  ## Add cortical volume name
  subcortical_volume <- UKB_subcortical_key$subcortical_volume[grepl(substring(lme_subcortical_regions$region_field_ID, 3,7), UKB_subcortical_key$feild_ID)]
  subcortical_volume_sub <- sub(("\\(.*"), "", subcortical_volume)
  lme_subcortical_regions$subcortical_volume <- str_to_title(sub(("Volume of "), "", subcortical_volume_sub))
  ## Append to data frame
  lme_subcortical_regions_GlycA <- rbind(lme_subcortical_regions, lme_subcortical_regions_GlycA)
}


## Plot effects of CRP PRS on regional sub cortical volumes
lme_subcortical_regions_GlycA$subcortical_volume = with(lme_subcortical_regions_GlycA, reorder(subcortical_volume, beta))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_subcortical_volume.jpg")

ggplot(data=lme_subcortical_regions_GlycA, aes(x=subcortical_volume, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Bilateral Subcortical Volumes") + ylab("Mean Standardized Effect Size (95% CI)") +
  theme_bw() +
  labs(title ="Association of CRP PRS and Subcortical Volumes")

dev.off()

write.csv(lme_subcortical_regions_GlycA, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_subcortical.csv", quote = F, row.names = F)


