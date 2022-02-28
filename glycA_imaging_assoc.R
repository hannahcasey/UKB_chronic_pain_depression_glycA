## Script overview:
## General coritcal volumes
## Regional cortical volumes
## Individual FA tracts
## Individual MD tracts
## General FA measures
## General MD Measures
## REgional subcortical volumes


## Load libraries
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(nlme)

## Set working directory
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
work_dir <- "~/Desktop/PhD/projects/UKB_inflammation_imaging/"
setwd(work_dir)

## Load in LME function developed by Shen
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("func/lme_Shen.R")


## Load in data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load dataframe containing markers of inflammation, imaging features of interest and covariates
UKB_inflammation_imaging_covariates <- read.csv("~/Desktop/PhD/projects/UKB_inflammation_imaging/resources/UKB_inflammation_imaging_covariates.csv", header = TRUE)

## Load in UKB DKW field ID key file
UKB_cortical_region_key <- read.csv("~/Desktop/PhD/projects/UKB_inflammation_imaging/resources/UKB_cortical_region_field_IDs.csv", header = F)
UKB_cortical_region_key <- UKB_cortical_region_key %>%
  rename(feild_ID = V1,
         cortical_volume = V2)

## Load in UKB FA field ID key file
UKB_FA_key <- read.csv("~/Desktop/PhD/projects/UKB_inflammation_imaging/resources/UKB_FA_field_IDs.csv", header = F)
UKB_FA_key <- UKB_FA_key %>%
  rename(feild_ID = V1,
         FA_value = V2)


## Load in UKB FA field ID key file
UKB_MD_key <- read.csv("~/Desktop/PhD/projects/UKB_inflammation_imaging/resources/UKB_MD_field_IDs.csv", header = F)
UKB_MD_key <- UKB_MD_key %>%
  rename(feild_ID = V1,
         MD_value = V2)


## Load in UKB FIRST subcortical field ID key file
UKB_subcortical_key <- read.csv("~/Desktop/PhD/projects/UKB_inflammation_imaging/resources/UKB_subcortical_field_IDs.csv", header = F)
UKB_subcortical_key <- UKB_subcortical_key %>%
  rename(feild_ID = V1,
         subcortical_volume = V2)

## Statistical analysis: general cortical volume ~ glycA
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Remove outliers - IQR method
outliers <- boxplot(UKB_inflammation_imaging_covariates$global_cortical_volume, plot=FALSE)$out
UKB_inflammation_imaging_covariates_global_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates
UKB_inflammation_imaging_covariates_global_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates$global_cortical_volume %in% outliers),]

## Statistical analysis
## Create dataframe to store glm output
glm_glycA_cortical_volume <- data.frame(cortical_volume=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                   p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric(), n=numeric() )

## CRP PRS ~ global corical volume + covariates
glm1 <- glm(global_cortical_volume  ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV,
            data = UKB_inflammation_imaging_covariates_global_cortical_volume_outliers_removed)

summary(glm1)

glm_glycA_cortical_volume[1,"cortical_volume"] <- "global_cortical_volume"
glm_glycA_cortical_volume[1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_cortical_volume[1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_cortical_volume[1, "p.value"] <- glm_glycA_cortical_volume[1, "p.value"]
glm_glycA_cortical_volume[1, "p.adjust"] <-glm_glycA_cortical_volume[1, "p.value"]
glm_glycA_cortical_volume[1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]
glm_glycA_cortical_volume[1, c("n")] <- nobs(glm1)

cortical_lobe_volumes <- c("frontal_lobe_volume", "temporal_lobe_volume", "parietal_lobe_volume", "occipital_lobe_volume", "cingulate_lobe_volume")

for (i in 1:5){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_lobe_volumes[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_lobar_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_lobar_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,cortical_lobe_volumes[i]] %in% outliers),]
  
  
  mod <- paste0(cortical_lobe_volumes[i], "~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age + age_squared + ICV")
  glm1 <- glm(as.formula(mod),
              data = UKB_inflammation_imaging_covariates_lobar_outliers_removed)
  
  glm_glycA_cortical_volume[i+1,"cortical_volume"] <- cortical_lobe_volumes[i]
  glm_glycA_cortical_volume[i+1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
  glm_glycA_cortical_volume[i+1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
  glm_glycA_cortical_volume[i+1, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
  glm_glycA_cortical_volume[i+1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]
  glm_glycA_cortical_volume[i+1, "p.adjust"] <- p.adjust(coef(summary(glm1))["glycA_1",4], method = "fdr", n = 5)
  glm_glycA_cortical_volume[i+1, c("n")] <- nobs(glm1)
  
  
}

## Plot results
# lock in factor level order
glm_glycA_cortical_volume$cortical_volume  = with(glm_glycA_cortical_volume, reorder(cortical_volume, beta))

p1 <- ggplot(data=glm_glycA_cortical_volume[!glm_glycA_cortical_volume$cortical_volume=="global_cortical_volume",], aes(x=cortical_volume, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Mean (95% CI)") +
  theme_bw() + # use a white background
  ylim(-0.3,0.3) +
  labs(title ="Association of GlycA and General Cortical Volumes") + 
  scale_x_discrete(labels=c("frontal_lobe" = "Frontal Lobe", "temporal_lobe" = "Temporal Lobe","parietal_lobe" = "Parietal Lobe", "occipital_lobe" = "Occipital Lobe", "cingulate_lobe" = "Cingulate Lobe"))

p2 <- ggplot(data=glm_glycA_cortical_volume[glm_glycA_cortical_volume$cortical_volume=="global_cortical_volume",], aes(x=cortical_volume, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
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

write.csv(glm_glycA_cortical_volume, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_cortical_volume_general.csv",
          quote = FALSE, row.names = FALSE)

## Statistical analysis: general cortical area ~ glycA
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Remove outliers - IQR method
outliers <- boxplot(UKB_inflammation_imaging_covariates$global_cortical_area, plot=FALSE)$out
UKB_inflammation_imaging_covariates_global_cortical_area_outliers_removed <- UKB_inflammation_imaging_covariates
UKB_inflammation_imaging_covariates_global_cortical_area_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates$global_cortical_area %in% outliers),]

## Statistical analysis
## Create dataframe to store glm output
glm_glycA_cortical_area <- data.frame(cortical_area=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                      p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric(), n=numeric())


## CRP PRS ~ global corical area + covariates
glm1 <- glm(global_cortical_area  ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV,
            data = UKB_inflammation_imaging_covariates_global_cortical_area_outliers_removed)

summary(glm1)

glm_glycA_cortical_area[1,"cortical_area"] <- "global_cortical_area"
glm_glycA_cortical_area[1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_cortical_area[1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_cortical_area[1, "p.value"] <- glm_glycA_cortical_area[1, "p.value"]
glm_glycA_cortical_area[1, "p.adjust"] <- glm_glycA_cortical_area[1, "p.value"] 
glm_glycA_cortical_area[1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]
glm_glycA_cortical_area[1, c("n")] <- nobs(glm1)


cortical_lobes <- c("frontal_lobe_area", "temporal_lobe_area", "parietal_lobe_area", "occipital_lobe_area", "cingulate_lobe_area")

for (i in 1:5){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_lobes[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_lobar_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_lobar_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,cortical_lobes[i]] %in% outliers),]
  
  
  mod <- paste0(cortical_lobes[i], "~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age + age_squared + ICV")
  glm1 <- glm(as.formula(mod),
              data = UKB_inflammation_imaging_covariates_lobar_outliers_removed)
  
  glm_glycA_cortical_area[i+1,"cortical_area"] <- cortical_lobes[i]
  glm_glycA_cortical_area[i+1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
  glm_glycA_cortical_area[i+1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
  glm_glycA_cortical_area[i+1, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
  glm_glycA_cortical_area[i+1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]
  glm_glycA_cortical_area[i+1, "p.adjust"] <- p.adjust(coef(summary(glm1))["glycA_1",4], method = "fdr", n = 5)
  glm_glycA_cortical_area[i+1, "n"] <-nobs(glm1)
  
}

## Plot results
# lock in factor level order
glm_glycA_cortical_area$cortical_area  = with(glm_glycA_cortical_area, reorder(cortical_area, beta))

p1 <- ggplot(data=glm_glycA_cortical_area[!glm_glycA_cortical_area$cortical_area=="global_cortical_area",], aes(x=cortical_area, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Mean (95% CI)") +
  theme_bw() + # use a white background
  ylim(-0.3,0.3) +
  labs(title ="Association of GlycA and General Cortical areas") + 
  scale_x_discrete(labels=c("frontal_lobe" = "Frontal Lobe", "temporal_lobe" = "Temporal Lobe","parietal_lobe" = "Parietal Lobe", "occipital_lobe" = "Occipital Lobe", "cingulate_lobe" = "Cingulate Lobe"))

p2 <- ggplot(data=glm_glycA_cortical_area[glm_glycA_cortical_area$cortical_area=="global_cortical_area",], aes(x=cortical_area, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Mean (95% CI)") +
  theme_bw() +  # use a white background
  ylim(-0.12,0.12) +
  labs(title ="Association of GlycA and Global Cortical area") +
  scale_x_discrete(labels=c("global_cortical_area" = "Global"))


layout <- c(
  area(t = 1, l = 1, b = 3, r = 4),
  area(t = 4, l = 1, b = 4, r = 4))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_general_cortical_area.jpg")
p1 / p2 + 
  plot_layout(design = layout)
dev.off()

write.csv(glm_glycA_cortical_area, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_cortical_area_general.csv",
          quote = FALSE, row.names = FALSE)

## Statistical analysis: general cortical thickness ~ glycA
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Remove outliers - IQR method
outliers <- boxplot(UKB_inflammation_imaging_covariates$global_cortical_thickness, plot=FALSE)$out
UKB_inflammation_imaging_covariates_global_cortical_thickness_outliers_removed <- UKB_inflammation_imaging_covariates
UKB_inflammation_imaging_covariates_global_cortical_thickness_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates$global_cortical_thickness %in% outliers),]

## Statistical analysis
## Create dataframe to store glm output
glm_glycA_cortical_thickness <- data.frame(cortical_thickness=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                           p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric())


## CRP PRS ~ global corical thickness + covariates
glm1 <- glm(global_cortical_thickness  ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV,
            data = UKB_inflammation_imaging_covariates_global_cortical_thickness_outliers_removed)

summary(glm1)

glm_glycA_cortical_thickness[1,"cortical_thickness"] <- "global_cortical_thickness"
glm_glycA_cortical_thickness[1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_cortical_thickness[1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_cortical_thickness[1, "p.value"] <- glm_glycA_cortical_thickness[1, "p.value"]
glm_glycA_cortical_thickness[1, "p.adjust"] <- glm_glycA_cortical_thickness[1, "p.value"]
glm_glycA_cortical_thickness[1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]
glm_glycA_cortical_thickness[1, c("n")] <- nobs(glm1)


cortical_lobes <- c("frontal_lobe_thickness", "temporal_lobe_thickness", "parietal_lobe_thickness", "occipital_lobe_thickness", "cingulate_lobe_thickness")

for (i in 1:5){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_lobes[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_lobar_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_lobar_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,cortical_lobes[i]] %in% outliers),]
  
  
  mod <- paste0(cortical_lobes[i], "~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age + age_squared + ICV")
  glm1 <- glm(as.formula(mod),
              data = UKB_inflammation_imaging_covariates_lobar_outliers_removed)
  
  glm_glycA_cortical_thickness[i+1,"cortical_thickness"] <- cortical_lobes[i]
  glm_glycA_cortical_thickness[i+1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
  glm_glycA_cortical_thickness[i+1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
  glm_glycA_cortical_thickness[i+1, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
  glm_glycA_cortical_thickness[i+1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]
  glm_glycA_cortical_thickness[i+1, "p.adjust"] <- p.adjust(coef(summary(glm1))["glycA_1",4], method = "fdr", n = 5)
  glm_glycA_cortical_thickness[i+1, "n"] <- nobs(glm1)
  
}

## Plot results
# lock in factor level order
glm_glycA_cortical_thickness$cortical_thickness  = with(glm_glycA_cortical_thickness, reorder(cortical_thickness, beta))

p1 <- ggplot(data=glm_glycA_cortical_thickness[!glm_glycA_cortical_thickness$cortical_thickness=="global_cortical_thickness",], aes(x=cortical_thickness, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Mean (95% CI)") +
  theme_bw() + # use a white background
  ylim(-0.6,0.6) +
  labs(title ="Association of GlycA and General Cortical thicknesss") + 
  scale_x_discrete(labels=c("frontal_lobe" = "Frontal Lobe", "temporal_lobe" = "Temporal Lobe","parietal_lobe" = "Parietal Lobe", "occipital_lobe" = "Occipital Lobe", "cingulate_lobe" = "Cingulate Lobe"))

p2 <- ggplot(data=glm_glycA_cortical_thickness[glm_glycA_cortical_thickness$cortical_thickness=="global_cortical_thickness",], aes(x=cortical_thickness, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Mean (95% CI)") +
  theme_bw() +  # use a white background
  ylim(-0.12,0.12) +
  labs(title ="Association of GlycA and Global Cortical thickness") +
  scale_x_discrete(labels=c("global_cortical_thickness" = "Global"))


layout <- c(
  area(t = 1, l = 1, b = 3, r = 4),
  area(t = 4, l = 1, b = 4, r = 4))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_general_cortical_thickness.jpg")
p1 / p2 + 
  plot_layout(design = layout)
dev.off()

write.csv(glm_glycA_cortical_thickness, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_cortical_thickness_general.csv",
          quote = FALSE, row.names = FALSE)


## Association analysis of individual cortical volumes
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Run interaction analysis to see if there is an effect modification between glycA and hemisphere on cortical volumes
cortical_volume_FIDs <- paste0("f.", c(26789:26821, 26890:26922), ".2.0")

## Create dataframe to store output
interaction_glycA_hemisphere <- data.frame(Volume=character(),  interaction_p=numeric(), interaction_p_adjust=numeric())

for (i in 1:33){
  
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_volume_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,cortical_volume_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_volume_FIDs[i +33]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[-which(UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[,cortical_volume_FIDs[i+33]] %in% outliers),]
  
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
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


## Run linear mixed effect model analysis to look at the effect of GlycA on repeat measures (hemispheres) of cortical volumes
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate dataframe to record P-values, CIs and beta coefficients

lme_cortical_volume_glycA <- data.frame(mod_name=character(), beta=numeric(), std=numeric(), t.value=numeric(), p.value=numeric(),
                                         Lower_95CI=numeric(), Upper_95CI=numeric())
## Generate ls.mod
ls.mod.PRS <- data.frame(Dep ="Volume",
                         Factor =  "glycA_1",
                         Covariates = "sex + age + age_squared + BMI + ICV + assessment_centre_first_imaging + Hemisphere",
                         Model = "lme")

## Iterate through each cortical volume and run lme 
for (i in 1:33){
  
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_volume_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,cortical_volume_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_volume_FIDs[i +33]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[-which(UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[,cortical_volume_FIDs[i+33]] %in% outliers),]
  
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
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
  cortical_volume <- UKB_cortical_region_key$cortical_volume[grepl(substring(lme_cortical_region_iterate_glycA$region_field_ID, 3,7), UKB_cortical_region_key$feild_ID)]
  cortical_volume_sub <- sub(("\\(.*"), "", cortical_volume)
  lme_cortical_region_iterate_glycA$cortical_volume <- str_to_title(sub(("Volume of "), "", cortical_volume_sub))
  ## Append to data frame
  lme_cortical_volume_glycA <- rbind(lme_cortical_region_iterate_glycA, lme_cortical_volume_glycA)
}


## Plot results
lme_cortical_volume_glycA$cortical_volume = with(lme_cortical_volume_glycA, reorder(cortical_volume, beta))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_regional_cortical_structures.jpg")

ggplot(data=lme_cortical_volume_glycA, aes(x=cortical_volume, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Label") + ylab("Mean (95% CI)") +
  theme_bw() +
  labs(title = "Association of GlycA and Regional Cortical Volumes")

dev.off()

write.csv(lme_cortical_volume_glycA, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/lme_cortical_volume.csv",
          quote = FALSE, row.names = FALSE)

## Association analysis of individual cortical areas
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Run interaction analysis to see if there is an effect modification between glycA and hemisphere on cortical areas
cortical_area_FIDs <- paste0("f.", c(26722:26754, 26823:26855), ".2.0")

## Create dataframe to store output
interaction_glycA_hemisphere <- data.frame(area=character(),  interaction_p=numeric(), interaction_p_adjust=numeric())

for (i in 1:33){

  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_area_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,cortical_area_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_area_FIDs[i +33]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[-which(UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[,cortical_area_FIDs[i+33]] %in% outliers),]
  
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                                                               cortical_area_FIDs[i], cortical_area_FIDs[i +33])]
  
  ## Convert short dataframe to long format
  UKB_glycA_DK_small_long <- reshape(UKB_glycA_DK_small, 
                                     varying = c(cortical_area_FIDs[i], cortical_area_FIDs[i + 33]), 
                                     v.names = "area",
                                     timevar = "Region", 
                                     times = c(cortical_area_FIDs[i], cortical_area_FIDs[i + 33]), 
                                     direction = "long")
  
  ## Add column indicating hemisphere
  UKB_glycA_DK_small_long$Hemisphere <- NA
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_area_FIDs[i]] <- "left"
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_area_FIDs[i + 33]] <- "right"
  
  ## Carry out interaction analysis
  lm1 <- lm(area ~ glycA_1 * Hemisphere + sex + age + age + age_squared + BMI + ICV + assessment_centre_first_imaging, data = UKB_glycA_DK_small_long)
  
  interaction_glycA_hemisphere[i,"area"] <- (cortical_area_FIDs)[i]
  interaction_glycA_hemisphere[i,"interaction_p"] <- summary(lm1)[["coefficients"]]["glycA_1:Hemisphereright", "Pr(>|t|)"]
  interaction_glycA_hemisphere[i,"interaction_p_adjust"] <- p.adjust(interaction_glycA_hemisphere[i,"interaction_p"], n = 33, method = "fdr")
}
## No hemisphere interaction found between glycA_1 and cortical areas

## Run lme analysis to look at the effects of CRP PRS on cortical areas
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate dataframe to record P-values, CIs and beta coefficients
lme_cortical_area_glycA <- data.frame(mod_name=character(), beta=numeric(), std=numeric(), t.value=numeric(), p.value=numeric(),
                                      Lower_95CI=numeric(), Upper_95CI=numeric())

## Generate ls.mod
ls.mod.glycA.area <- data.frame(Dep ="area",
                         Factor =  "glycA_1",
                         Covariates = "sex + age + age_squared + BMI + ICV + assessment_centre_first_imaging + Hemisphere",
                         Model = "lme")


## Iterate through each cortical area and run glm 
for (i in 1:33){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_area_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_area_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_cortical_area_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,cortical_area_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_area_FIDs[i +33]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_area_outliers_removed <- UKB_inflammation_imaging_covariates_cortical_area_outliers_removed[-which(UKB_inflammation_imaging_covariates_cortical_area_outliers_removed[,cortical_area_FIDs[i+33]] %in% outliers),]
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_inflammation_imaging_covariates_cortical_area_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                                                               cortical_area_FIDs[i], cortical_area_FIDs[i +33])]
  
  ## Convert short dataframe to long format
  UKB_glycA_DK_small_long <- reshape(UKB_glycA_DK_small, 
                                     varying = c(cortical_area_FIDs[i], cortical_area_FIDs[i + 33]), 
                                     v.names = "area",
                                     timevar = "Region", 
                                     times = c(cortical_area_FIDs[i], cortical_area_FIDs[i + 33]), 
                                     direction = "long")
  
  ## Add column indicating hemisphere
  UKB_glycA_DK_small_long$Hemisphere <- NA
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_area_FIDs[i]] <- "left"
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_area_FIDs[i + 33]] <- "right"
  
  ## Generate ls.mode
  
  ## Rum lme and get output (BTEA, SD, T, P, 95% CIs)
  output <- try(run_model(ls.mod.glycA.area, UKB_glycA_DK_small, UKB_glycA_DK_small_long))
  ## Add region field ID
  lme_cortical_region_iterate_glycA <- cbind(data_frame(region_field_ID = paste0(cortical_area_FIDs[i],"_",cortical_area_FIDs[i + 33])), output)
  ## Adjust P value for multiple comparisons (FDR)
  lme_cortical_region_iterate_glycA$p.value.adjust <- p.adjust(lme_cortical_region_iterate_glycA$p.value, n = 33, method = "fdr")
  ## Add cortical area name
  cortical_area <- UKB_cortical_region_key$cortical_volume[grepl(substring(lme_cortical_region_iterate_glycA$region_field_ID, 3,7), UKB_cortical_region_key$feild_ID)]
  cortical_area_sub <- sub(("\\(.*"), "", cortical_area)
  lme_cortical_region_iterate_glycA$cortical_area <- str_to_title(sub(("Area of "), "", cortical_area_sub))
  ## Append to data frame
  lme_cortical_area_glycA <- rbind(lme_cortical_region_iterate_glycA, lme_cortical_area_glycA)
}

## Plot results
lme_cortical_area_glycA$cortical_area = with(lme_cortical_area_glycA, reorder(cortical_area, beta))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_regional_cortical_structures.jpg")

ggplot(data=lme_cortical_area_glycA, aes(x=cortical_area, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Label") + ylab("Mean (95% CI)") +
  theme_bw() +
  labs(title = "Association of GlycA and Regional Cortical areas")

dev.off()

write.csv(lme_cortical_area_glycA, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/lme_cortical_area_glycA.csv",
          quote = FALSE, row.names = FALSE)



## Association analysis of individual cortical thicknesss
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Run interaction analysis to see if there is an effect modification between glycA and hemisphere on cortical thicknesss
cortical_thickness_FIDs <- paste0("f.", c(26756:26788, 26857:26889), ".2.0")

## Create dataframe to store output
interaction_glycA_hemisphere <- data.frame(thickness=character(),  interaction_p=numeric(), interaction_p_adjust=numeric())

for (i in 1:33){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_thickness_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,cortical_thickness_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_thickness_FIDs[i +33]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[-which(UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[,cortical_thickness_FIDs[i+33]] %in% outliers),]
  
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_inflammation_imaging_covariates_cortical_volume_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                                                              cortical_thickness_FIDs[i], cortical_thickness_FIDs[i +33])]
  
  ## Convert short dataframe to long format
  UKB_glycA_DK_small_long <- reshape(UKB_glycA_DK_small, 
                                     varying = c(cortical_thickness_FIDs[i], cortical_thickness_FIDs[i + 33]), 
                                     v.names = "thickness",
                                     timevar = "Region", 
                                     times = c(cortical_thickness_FIDs[i], cortical_thickness_FIDs[i + 33]), 
                                     direction = "long")
  
  ## Add column indicating hemisphere
  UKB_glycA_DK_small_long$Hemisphere <- NA
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_thickness_FIDs[i]] <- "left"
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_thickness_FIDs[i + 33]] <- "right"
  
  ## Carry out interaction analysis
  lm1 <- lm(thickness ~ glycA_1 * Hemisphere + sex + age + age + age_squared + BMI + ICV + assessment_centre_first_imaging, data = UKB_glycA_DK_small_long)
  
  interaction_glycA_hemisphere[i,"thickness"] <- (cortical_thickness_FIDs)[i]
  interaction_glycA_hemisphere[i,"interaction_p"] <- summary(lm1)[["coefficients"]]["glycA_1:Hemisphereright", "Pr(>|t|)"]
  interaction_glycA_hemisphere[i,"interaction_p_adjust"] <- p.adjust(interaction_glycA_hemisphere[i,"interaction_p"], n = 33, method = "fdr")
}
## No hemisphere interaction found between glycA_1 and thicknesss

## Run lme analysis to look at the effects of CRP PRS on cortical thickness
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate dataframe to record P-values, CIs and beta coefficients
lme_cortical_thickness_glycA <- data.frame(mod_name=character(), beta=numeric(), std=numeric(), t.value=numeric(), p.value=numeric(),
                                      Lower_95CI=numeric(), Upper_95CI=numeric())

## Generate ls.mod
ls.mod.glycA.thickness <- data.frame(Dep ="thickness",
                                Factor =  "glycA_1",
                                Covariates = "sex + age + age_squared + BMI + ICV + assessment_centre_first_imaging + Hemisphere",
                                Model = "lme")


## Iterate through each cortical thickness and run glm 
for (i in 1:33){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_thickness_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_thickness_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_cortical_thickness_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,cortical_thickness_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,cortical_thickness_FIDs[i +33]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_cortical_thickness_outliers_removed <- UKB_inflammation_imaging_covariates_cortical_thickness_outliers_removed[-which(UKB_inflammation_imaging_covariates_cortical_thickness_outliers_removed[,cortical_thickness_FIDs[i+33]] %in% outliers),]
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_inflammation_imaging_covariates_cortical_thickness_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                                                            cortical_thickness_FIDs[i], cortical_thickness_FIDs[i +33])]
  
  ## Convert short dataframe to long format
  UKB_glycA_DK_small_long <- reshape(UKB_glycA_DK_small, 
                                     varying = c(cortical_thickness_FIDs[i], cortical_thickness_FIDs[i + 33]), 
                                     v.names = "thickness",
                                     timevar = "Region", 
                                     times = c(cortical_thickness_FIDs[i], cortical_thickness_FIDs[i + 33]), 
                                     direction = "long")
  
  ## Add column indicating hemisphere
  UKB_glycA_DK_small_long$Hemisphere <- NA
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_thickness_FIDs[i]] <- "left"
  UKB_glycA_DK_small_long$Hemisphere[UKB_glycA_DK_small_long$Region == cortical_thickness_FIDs[i + 33]] <- "right"
  
  ## Generate ls.mode
  
  ## Rum lme and get output (BTEA, SD, T, P, 95% CIs)
  output <- try(run_model(ls.mod.glycA.thickness, UKB_glycA_DK_small, UKB_glycA_DK_small_long))
  ## Add region field ID
  lme_cortical_region_iterate_glycA <- cbind(data_frame(region_field_ID = paste0(cortical_thickness_FIDs[i],"_",cortical_thickness_FIDs[i + 33])), output)
  ## Adjust P value for multiple comparisons (FDR)
  lme_cortical_region_iterate_glycA$p.value.adjust <- p.adjust(lme_cortical_region_iterate_glycA$p.value, n = 33, method = "fdr")
  ## Add cortical thickness name
  cortical_thickness <- UKB_cortical_region_key$cortical_volume[grepl(substring(lme_cortical_region_iterate_glycA$region_field_ID, 3,7), UKB_cortical_region_key$feild_ID)]
  cortical_thickness_sub <- sub(("\\(.*"), "", cortical_thickness)
  lme_cortical_region_iterate_glycA$cortical_thickness <- str_to_title(sub(("Mean thickness of "), "", cortical_thickness_sub))
  ## Append to data frame
  lme_cortical_thickness_glycA <- rbind(lme_cortical_region_iterate_glycA, lme_cortical_thickness_glycA)
}


## Plot results
lme_cortical_thickness_glycA$cortical_thickness = with(lme_cortical_thickness_glycA, reorder(cortical_thickness, beta))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_regional_cortical_structures.jpg")

ggplot(data=lme_cortical_thickness_glycA, aes(x=cortical_thickness, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Label") + ylab("Mean (95% CI)") +
  theme_bw() +
  labs(title = "Association of GlycA and Regional Cortical thicknesss")

dev.off()

write.csv(lme_cortical_thickness_glycA, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/lme_cortical_thickness_glycA.csv",
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
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,FA_bilateral_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,FA_bilateral_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,FA_bilateral_FIDs[i +1]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed <- UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed[-which(UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed[,FA_bilateral_FIDs[i+1]] %in% outliers),]
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
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

## Hemisphere interaction between association of acoustic radiation(25488 and 25489) and GlycA

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
## Skip first two FA measurements as significant hemisphere interaction was found in acoustic radiation(25488 and 25489), i.e. first two field IDs in array
for (i in seq(3, 23, by=2)){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,FA_bilateral_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,FA_bilateral_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,FA_bilateral_FIDs[i +1]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed <- UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed[-which(UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed[,FA_bilateral_FIDs[i+1]] %in% outliers),]
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
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
  FA_value <- UKB_FA_key$FA_value[grepl(substring(lme_FA_meaures_glycA $region_field_ID, 3,7), UKB_FA_key$feild_ID)]
  FA_value_sub <- sub(("\\(.*"), "", FA_value)
  lme_FA_meaures_glycA$FA_tract <- str_to_title(sub(("Weighted-mean FA in tract "), "", FA_value_sub))
  ## Append to data frame
  lme_FA_values_glycA <- rbind(lme_FA_meaures_glycA, lme_FA_values_glycA)
}


# Carry out GLM on each hemispheric measurement of FA in tract acoustic radiation

glm_FA_values_glycA <- data.frame(region_field_ID=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                             Lower_95CI=numeric(), Upper_95CI=numeric())

for (i in 1:2){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,FA_bilateral_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,FA_bilateral_FIDs[i]] %in% outliers),]

  ## Scale FA value in unilateral tracts
  UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed[,FA_bilateral_FIDs[i]] <- scale(UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed[,FA_bilateral_FIDs[i]])
  ## Run glm
  mod <- paste0(FA_bilateral_FIDs[i], "~ glycA_1", "+ sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV")
  glm1 <- glm(as.formula(mod), data = UKB_inflammation_imaging_covariates_FA_tracts_outliers_removed)


  glm_FA_values_glycA[i,"region_field_ID"] <- FA_bilateral_FIDs[i]
  glm_FA_values_glycA[i, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
  glm_FA_values_glycA[i, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
  glm_FA_values_glycA[i, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
  glm_FA_values_glycA[i, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]
  glm_FA_values_glycA[i, "p.value.adjust"] <- p.adjust(glm_FA_values_glycA[i, "p.value"], method = "fdr", n = 16)

  FA_value <- UKB_FA_key$FA_value[grepl(substring(glm_FA_values_glycA$region_field_ID[i], 3,7), UKB_FA_key$feild_ID)]
  glm_FA_values_glycA$FA_tract[i] <- str_to_title(sub(("Weighted-mean FA in tract "), "", FA_value))
}


## Carry out GLM on unilateral FA tracts

## Create new dataframe to store glm output
glm_FA_values_glycA_unilateral <- data.frame(region_field_ID=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                             Lower_95CI=numeric(), Upper_95CI=numeric()) 

for (i in 1:3){
  
  ## Scale FA value in unilateral tracts
  UKB_inflammation_imaging_covariates[,FA_unilateral_FIDs[i]] <- scale(UKB_inflammation_imaging_covariates[,FA_unilateral_FIDs[i]])
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,FA_unilateral_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_FA_unilateral_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_FA_unilateral_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,FA_unilateral_FIDs[i]] %in% outliers),]
  
  ## Run glm
  mod <- paste0("glycA_1~ ",FA_unilateral_FIDs[i], "+ sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV")
  glm1 <- glm(as.formula(mod), data = UKB_inflammation_imaging_covariates)
  
  
  glm_FA_values_glycA_unilateral[i,"region_field_ID"] <- FA_unilateral_FIDs[i]
  glm_FA_values_glycA_unilateral[i, "beta"] <- summary(glm1)[["coefficients"]][FA_unilateral_FIDs[i], "Estimate"]
  glm_FA_values_glycA_unilateral[i, "std"] <- summary(glm1)[["coefficients"]][FA_unilateral_FIDs[i], "Std. Error"]
  glm_FA_values_glycA_unilateral[i, "p.value"] <- summary(glm1)[["coefficients"]][FA_unilateral_FIDs[i], "Pr(>|t|)"]
  glm_FA_values_glycA_unilateral[i, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)[FA_unilateral_FIDs[i],]
  glm_FA_values_glycA_unilateral[i, "p.value.adjust"] <- p.adjust(glm_FA_values_glycA_unilateral[i, "p.value"], method = "fdr", n = 16)
  
  FA_value <- UKB_FA_key$FA_value[grepl(substring(glm_FA_values_glycA_unilateral$region_field_ID[i], 3,7), UKB_FA_key$feild_ID)]
  glm_FA_values_glycA_unilateral$FA_tract[i] <- str_to_title(sub(("Weighted-mean FA in tract "), "", FA_value))
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

## Plot effect of CRP PRS on bilateral FA tracts where there is a hemisphere interaction
glm_FA_values_glycA$FA_tract = with(glm_FA_values_glycA, reorder(FA_tract, beta))

p2 <- ggplot(data=glm_FA_values_glycA, aes(x=FA_tract, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Bilateral FA Tracts") + ylab("Mean (95% CI)") +
  theme_bw() +
  ylim(-0.6, 0.6) +
  labs(title ="Association of GlycA and FA Bilateral Tracts (Hemisphere Interaction)")


write.csv(glm_FA_values_glycA, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_FA_bilateral_tracts_hemisphere_interaction.csv",
          quote = FALSE, row.names = FALSE)


## Plot effect of CRP PRS on unilateral FA tracts
glm_FA_values_glycA_unilateral$FA_tract = with(glm_FA_values_glycA_unilateral, reorder(FA_tract, beta))

p3 <- ggplot(data=glm_FA_values_glycA_unilateral, aes(x=FA_tract, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Unilateral FA Tracts") + ylab("Mean (95% CI)") +
  theme_bw() +
  ylim(-5e-03, 5e-03) +
  labs(title ="Association of GlycA and FA Unilateral Tracts")

layout <- c(
  area(t = 1, l = 1, b = 2, r = 4),
  area(t = 3, l = 1, b = 3, r = 4),
  area(t = 4, l = 1, b = 4, r = 4))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_FA_tracts.jpg")

p1 / p2 / p3 + 
  plot_layout(design = layout)

dev.off()

write.csv(glm_FA_values_glycA_unilateral, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_FA_unilateral_tracts.csv",
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
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,MD_bilateral_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,MD_bilateral_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,MD_bilateral_FIDs[i +1]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed <- UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed[-which(UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed[,MD_bilateral_FIDs[i+1]] %in% outliers),]
  
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed[,c("f.eid", "glycA_1",  "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
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
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,MD_bilateral_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,MD_bilateral_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,MD_bilateral_FIDs[i +1]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed <- UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed[-which(UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed[,MD_bilateral_FIDs[i+1]] %in% outliers),]
  
  ## Select essential columns
  UKB_glycA_DK_small <- UKB_inflammation_imaging_covariates_MD_tracts_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
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
  MD_value <- UKB_MD_key$MD_value[grepl(substring(lme_MD_meaures_glycA$region_field_ID, 3,7), UKB_MD_key$feild_ID)]
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
  UKB_inflammation_imaging_covariates[,MD_unilateral_FIDs[i]] <- scale(UKB_inflammation_imaging_covariates[,MD_unilateral_FIDs[i]])
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,MD_unilateral_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_MD_unilateral_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_MD_unilateral_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,MD_unilateral_FIDs[i]] %in% outliers),]
  
  ## Run glm
  mod <- paste0(MD_unilateral_FIDs[i], " ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV")
  glm1 <- glm(as.formula(mod), data = UKB_inflammation_imaging_covariates_MD_unilateral_outliers_removed)
  
  
  lme_MD_values_glycA_unilateral[i,"region_field_ID"] <- MD_unilateral_FIDs[i]
  lme_MD_values_glycA_unilateral[i, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
  lme_MD_values_glycA_unilateral[i, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
  lme_MD_values_glycA_unilateral[i, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
  lme_MD_values_glycA_unilateral[i, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]
  lme_MD_values_glycA_unilateral[i, "p.value.adjust"] <- p.adjust(lme_MD_values_glycA_unilateral[i, "p.value"], method = "fdr", n = 15)
  
  MD_value <- UKB_MD_key$MD_value[grepl(substring(lme_MD_values_glycA_unilateral$region_field_ID[i], 3,7), UKB_MD_key$feild_ID)]
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



## Run glm to look at the effect of glycA on general measures of MD
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Global DM
## Association fibre MD (multiple by -1 as 1st principle component was flipped)
## Thalamic radiation fibre MD
## Projection fibre MD

## Create df to store glm output 
glm_glycA_MD_general <- data.frame(general_DTI_MD_measure=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                     p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric())


## Remove outliers - IQR method
outliers <- boxplot(UKB_inflammation_imaging_covariates$global_MD, plot=FALSE)$out
UKB_inflammation_imaging_covariates_global_MD_PC_outliers_removed <- UKB_inflammation_imaging_covariates
UKB_inflammation_imaging_covariates_global_MD_PC_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates$global_MD %in% outliers),]

outliers <- boxplot(UKB_inflammation_imaging_covariates$association_fibres_MD, plot=FALSE)$out
UKB_glycA_assoc_MD_PC_outliers_removed <- UKB_inflammation_imaging_covariates
UKB_glycA_assoc_MD_PC_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates$association_fibres_MD %in% outliers),]

outliers <- boxplot(UKB_inflammation_imaging_covariates$thalamic_radiations_MD, plot=FALSE)$out
UKB_glycA_tr_MD_PC_outliers_removed <- UKB_inflammation_imaging_covariates
UKB_glycA_tr_MD_PC_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates$thalamic_radiations_MD %in% outliers),]

outliers <- boxplot(UKB_inflammation_imaging_covariates$projection_fibres_MD, plot=FALSE)$out
UKB_glycA_projec_MD_PC_outliers_removed <- UKB_inflammation_imaging_covariates
UKB_glycA_projec_MD_PC_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates$projection_fibres_MD %in% outliers),]


## run GLM
glm1 <- glm(global_MD ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV,
            data = UKB_inflammation_imaging_covariates_global_MD_PC_outliers_removed)
summary(glm1)

glm2 <- glm((association_fibres_MD*-1) ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_glycA_assoc_MD_PC_outliers_removed)
summary(glm2)

glm3 <- glm(thalamic_radiations_MD ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_glycA_tr_MD_PC_outliers_removed)
summary(glm3)

glm4 <- glm(projection_fibres_MD ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_glycA_projec_MD_PC_outliers_removed)
summary(glm4)




glm_glycA_MD_general[1,"general_DTI_MD_measure"] <- "global_MD"
glm_glycA_MD_general[1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_MD_general[1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_MD_general[1, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
glm_glycA_MD_general[1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]

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


## Run glm to look at the effect of glycA on general measures of FA
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Global FA (multiple by -1 as 1st principle component was flipped)
## Association fibre FA (multiple by -1 as 1st principle component was flipped)
## Thalamic radiation fibre FA (multiple by -1 as 1st principle component was flipped)
## Projection fibre FA (multiple by -1 as 1st principle component was flipped)

## Create df to store glm output 
glm_glycA_FA_general <- data.frame(general_DTI_FA_measure=character(), beta=numeric(), std=numeric(), p.value=numeric(),
                                   p.adjust=numeric(), Lower_95CI=numeric(), Upper_95CI=numeric())


## Remove outliers - IQR method
outliers <- boxplot(UKB_inflammation_imaging_covariates$global_FA, plot=FALSE)$out
UKB_inflammation_imaging_covariates_global_FA_PC_outliers_removed <- UKB_inflammation_imaging_covariates
UKB_inflammation_imaging_covariates_global_FA_PC_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates$global_FA %in% outliers),]

outliers <- boxplot(UKB_inflammation_imaging_covariates$association_fibres_FA, plot=FALSE)$out
UKB_glycA_assoc_FA_PC_outliers_removed <- UKB_inflammation_imaging_covariates
UKB_glycA_assoc_FA_PC_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates$association_fibres_FA %in% outliers),]

outliers <- boxplot(UKB_inflammation_imaging_covariates$thalamic_radiations_FA, plot=FALSE)$out
UKB_glycA_tr_FA_PC_outliers_removed <- UKB_inflammation_imaging_covariates
UKB_glycA_tr_FA_PC_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates$thalamic_radiations_FA %in% outliers),]

outliers <- boxplot(UKB_inflammation_imaging_covariates$projection_fibres_FA, plot=FALSE)$out
UKB_glycA_projec_FA_PC_outliers_removed <- UKB_inflammation_imaging_covariates
UKB_glycA_projec_FA_PC_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates$projection_fibres_FA %in% outliers),]


## run GLM
glm1 <- glm((global_FA*-1) ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + age_squared + ICV,
            data = UKB_inflammation_imaging_covariates_global_FA_PC_outliers_removed)
summary(glm1)

glm2 <- glm((association_fibres_FA * -1) ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_glycA_assoc_FA_PC_outliers_removed)
summary(glm2)

glm3 <- glm((thalamic_radiations_FA*-1) ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_glycA_tr_FA_PC_outliers_removed)
summary(glm3)

glm4 <- glm((projection_fibres_FA*-1) ~ glycA_1 + sex + BMI + assessment_centre_first_imaging + age + ICV,
            data = UKB_glycA_projec_FA_PC_outliers_removed)
summary(glm4)




glm_glycA_FA_general[1,"general_DTI_FA_measure"] <- "global_FA"
glm_glycA_FA_general[1, "beta"] <- summary(glm1)[["coefficients"]]["glycA_1", "Estimate"]
glm_glycA_FA_general[1, "std"] <- summary(glm1)[["coefficients"]]["glycA_1", "Std. Error"]
glm_glycA_FA_general[1, "p.value"] <- summary(glm1)[["coefficients"]]["glycA_1", "Pr(>|t|)"]
glm_glycA_FA_general[1, c("Lower_95CI", "Upper_95CI")] <- confint(glm1)["glycA_1",]

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
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,sub_cortical_volume_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,sub_cortical_volume_FIDs[i]] %in% outliers),]
  
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,sub_cortical_volume_FIDs[i + 1]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed[-which(UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed[,sub_cortical_volume_FIDs[i + 1]] %in% outliers),]
  
  
  ## Select essential columns
  UKB_inflammation_imaging_covariates_small <- UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed[,c("f.eid", "glycA_1", "sex", "age", "age_squared", "BMI", "ICV", "assessment_centre_first_imaging",
                                                                              sub_cortical_volume_FIDs[i], sub_cortical_volume_FIDs[i + 1])]
  
  ## Convert short dataframe to long format
  UKB_inflammation_imaging_covariates_small_long <- reshape(UKB_inflammation_imaging_covariates_small, 
                                       varying = c(sub_cortical_volume_FIDs[i], sub_cortical_volume_FIDs[i + 1]), 
                                       v.names = "Volume",
                                       timevar = "Region", 
                                       times = c(sub_cortical_volume_FIDs[i], sub_cortical_volume_FIDs[i + 1]), 
                                       direction = "long")
  
  ## Add column indicating hemisphere
  UKB_inflammation_imaging_covariates_small_long$Hemisphere <- NA
  UKB_inflammation_imaging_covariates_small_long$Hemisphere[UKB_inflammation_imaging_covariates_small_long$Region == sub_cortical_volume_FIDs[i]] <- "left"
  UKB_inflammation_imaging_covariates_small_long$Hemisphere[UKB_inflammation_imaging_covariates_small_long$Region == sub_cortical_volume_FIDs[i + 1]] <- "right"
  
  ## Carry out interaction analysis
  lm1 <- lm(Volume ~ glycA_1 * Hemisphere + sex + age + age_squared + BMI + ICV + assessment_centre_first_imaging, data = UKB_inflammation_imaging_covariates_small_long)
  
  interaction_glycA_hemisphere[i,"Volume"] <- (sub_cortical_volume_FIDs)[i]
  interaction_glycA_hemisphere[i,"interaction_p"] <- summary(lm1)[["coefficients"]]["glycA_1:Hemisphereright", "Pr(>|t|)"]
  
  interaction_glycA_hemisphere[i,"interaction_p_adjust"] <- p.adjust(interaction_glycA_hemisphere[i,"interaction_p"], n = 7, method = "fdr")
}

interaction_glycA_hemisphere
## No interaction found between glycA and hemisphere on volumes

## Run generalized linear model analysis to look at the effect of GlycA on subcorical volumes
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate dataframe to record P-values, CIs and beta coefficients
glm_subcortical_volume_GlycA<- data.frame(mod_name=character(), beta=numeric(), std=numeric(), t.value=numeric(), p.value=numeric(),
                                              Lower_95CI=numeric(), Upper_95CI=numeric())

## Iterate through each cortical volume and run glm 
for (i in seq(1,14,2)){
  
  ## Remove outliers - IQR method
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,sub_cortical_volume_FIDs[i]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates
  UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates[-which(UKB_inflammation_imaging_covariates[,sub_cortical_volume_FIDs[i]] %in% outliers),]
  outliers <- boxplot(UKB_inflammation_imaging_covariates[,sub_cortical_volume_FIDs[i + 1]], plot=FALSE)$out
  UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed <- UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed[-which(UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed[,sub_cortical_volume_FIDs[i + 1]] %in% outliers),]
  
  ## Sum hemispheres
  UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed$total_volume <- rowSums(UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed[,c(sub_cortical_volume_FIDs[i],sub_cortical_volume_FIDs[i + 1])], na.rm = FALSE)
  
  ## Scale summed volume
  UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed$total_volume <- scale(UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed$total_volume)
  
  ## Rum glm and get output (BTEA, SD, T, P, 95% CIs)
  fit <- glm(total_volume ~ glycA_1 + sex + age + age_squared + BMI + ICV + assessment_centre_first_imaging, data = UKB_inflammation_imaging_covariates_subcortical_volume_outliers_removed)
  
  tmp.ci = confint(fit) %>% as.data.frame %>% 
    select(Lower_95CI=`2.5 %`,Upper_95CI=`97.5 %`) %>% 
    head(2) %>%
    tail(1)
  
  tmp.res = summary(fit)$coefficients %>% 
    as.data.frame %>% 
    select(beta=Estimate,std=`Std. Error`,t.value=`t value`,p.value=`Pr(>|t|)`) %>% 
    head(2) %>%
    tail(1) %>%
    mutate(mod_name = paste0("Cortical Volume",'~',"GlycA"),tmp.ci) %>% 
    select(mod_name, everything())
  
  ## Add region field ID
  glm_subcortical_volume_GlycA_temp <- cbind(data_frame(region_field_ID = paste0(sub_cortical_volume_FIDs[i],"_",sub_cortical_volume_FIDs[i + 1])), tmp.res)
  ## Adjust P value for multiple comparisons (FDR)
  glm_subcortical_volume_GlycA_temp$p.value.adjust <- p.adjust(glm_subcortical_volume_GlycA_temp$p.value, n = 7, method = "fdr")
  ## Add cortical volume name
  subcortical_volume <- UKB_subcortical_key$subcortical_volume[grepl(substring(glm_subcortical_volume_GlycA_temp$region_field_ID, 3,7), UKB_subcortical_key$feild_ID)]
  subcortical_volume_sub <- sub(("\\(.*"), "", subcortical_volume)
  glm_subcortical_volume_GlycA_temp$subcortical_volume <- str_to_title(sub(("Volume of "), "", subcortical_volume_sub))
  ## Append to data frame
  glm_subcortical_volume_GlycA <- rbind(glm_subcortical_volume_GlycA_temp, glm_subcortical_volume_GlycA)
}


## Plot effects of CRP PRS on regional sub cortical volumes
glm_subcortical_volume_GlycA$subcortical_volume = with(glm_subcortical_volume_GlycA, reorder(subcortical_volume, beta))

jpeg("~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/association_glycA_subcortical_volume.jpg")

ggplot(data=glm_subcortical_volume_GlycA, aes(x=subcortical_volume, y=beta, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Bilateral Subcortical Volumes") + ylab("Mean Standardized Effect Size (95% CI)") +
  theme_bw() +
  labs(title ="Association of CRP PRS and Subcortical Volumes")

dev.off()

write.csv(glm_subcortical_volume_GlycA, "~/Desktop/PhD/projects/UKBChronicPainDepressionGlycA/output/imaging_association_results/glycA_subcortical.csv", quote = F, row.names = F)


