#Load appropriate libraries
library(readxl)
library(openxlsx)
library(survival)
library(survminer)
library(tidyverse)

#Set working directory 
setwd("/Users/pettya/Desktop/R Analysis")


#Kaplan-Meier analysis on RT, No RT and TP53 WT, VAF Low, VAF High (no multiples) (Figures 1 C,D)
##Read in excel file with all TCGA patient data
endometrial_all <- read.xlsx("Final TCGA with p53 and RT annotation - R.xlsx", 1)

##Create new categorical "Progressed" and "VAF" columns
endometrial_all$Progressed <- ifelse(endometrial_all$Progression.Free.Status == "0:CENSORED", 0, 1)
endometrial_all$VAF[endometrial_all$Allele.Freq..T.<= 0.50] <- "low"
endometrial_all$VAF[endometrial_all$Allele.Freq..T.> 0.50] <- "high"
endometrial_all$VAF[endometrial_all$P53.Mut== 0] <- "WT"

##Create two separate tables for No RT and RT patients
endometrial_all_RT <- endometrial_all[which(endometrial_all$Radiation.Therapy == "Yes"),]
endometrial_all_NoRT <- endometrial_all[which(endometrial_all$Radiation.Therapy == "No"),]

##Order tables by Progression-Free Survival time
endometrial_all_RT <- endometrial_all_RT[order(endometrial_all_RT$Progress.Free.Survival..Months.),]
endometrial_all_NoRT <- endometrial_all_NoRT[order(endometrial_all_NoRT$Progress.Free.Survival..Months.),]

##Convert VAF column to ordered factor
endometrial_all_RT$VAF <- factor(endometrial_all_RT$VAF, levels= c("WT","low","high"))
endometrial_all_NoRT$VAF <- factor(endometrial_all_NoRT$VAF, levels= c("WT","low","high"))

##Create survival objects, fit curves based on VAF, and generate Kaplan-Meier plots for both No RT and RT data
endo_surv_RT <- Surv(time = endometrial_all_RT$Progress.Free.Survival..Months., event = endometrial_all_RT$Progressed)
endo_fit_RT <- surv_fit(endo_surv_RT ~ VAF, data = endometrial_all_RT)
KMRT <- ggsurvplot(endo_fit_RT, data = endometrial_all_RT, palette = alpha(c("royalblue","springgreen4","firebrick2"), 0.6), risk.table = TRUE, 
                   risk.table.height = 0.30, xlab = "Time in Months", xlim = c(0,150), legend.labs = c("WT","VAF Low","VAF High"), 
                   title = "Progression-Free Survival (RT)")
endo_surv_NoRT <- Surv(time = endometrial_all_NoRT$Progress.Free.Survival..Months., event = endometrial_all_NoRT$Progressed)
endo_fit_NoRT <- surv_fit(endo_surv_NoRT ~ VAF, data = endometrial_all_NoRT)
KMNoRT <- ggsurvplot(endo_fit_NoRT, data = endometrial_all_NoRT, palette = alpha(c("royalblue","springgreen4","firebrick2"), 0.6), risk.table = TRUE, 
                     risk.table.height = 0.30, xlab = "Time in Months", xlim = c(0,150), break.x.by = 50, legend.labs = c("WT","VAF Low","VAF High"), 
                     title = "Progression-Free Survival (No RT)")

##Print graphs as jpeg files
jpeg("KM, RT, VAF, no chemo or IV.jpeg", width= 1000, height= 720, units= "px", res= 150)
KMRT 
dev.off()

jpeg("KM, No RT, VAF, no chemo or IV.jpeg", width= 1000, height= 720, units= "px", res= 150)
KMNoRT 
dev.off()



#Cox Hazards Analysis on RT, No RT and TP53 WT, VAF Low, VAF High (no multiples) (Supp Figure 1)
##If using VAF Low as the reference, use below script for reordering VAF
endometrial_all_RT$VAF <- factor(endometrial_all_RT$VAF, levels= c("low","WT","high"))
endometrial_all_NoRT$VAF <- factor(endometrial_all_NoRT$VAF, levels= c("low","WT","high"))

##Generate Cox analysis objects for No RT and RT based on VAF using previously created survival objects
endo_cox_RT <- coxph(endo_surv_RT ~ VAF, data = endometrial_all_RT)
endo_cox_NoRT <- coxph(endo_surv_NoRT ~ VAF, data = endometrial_all_NoRT)

##Generate forest plots and print plots as jpeg files
jpeg("Cox, RT, VAF, ref WT, no chemo or IV.jpeg", width= 720, height= 1000, units= "px", res= 150)
CoxRT <- ggforest(endo_cox_RT, data = endometrial_all_RT, main = "Progression-Free Survival (RT)")
CoxRT
dev.off()

jpeg("Cox, No RT, VAF, ref WT, no chemo or IV.jpeg", width= 720, height= 1000, units= "px", res= 150)
CoxNoRT <- ggforest(endo_cox_NoRT, data = endometrial_all_NoRT, main = "Progression-Free Survival (No RT)")
CoxNoRT
dev.off()



#Kaplan-Meier Analysis on RT, No RT and TP53 WT, Mutant (no multiples) (Figure 1B)
##Convert P53.Mut column to ordered factor
endometrial_all_NoRT$P53.Mut <- factor(endometrial_all_NoRT$P53.Mut, levels = c(0,1), labels= c("p53 WT","p53 Mutant"))
endometrial_all_RT$P53.Mut <- factor(endometrial_all_RT$P53.Mut, levels = c(0,1), labels= c("p53 WT","p53 Mutant"))

##Fit survival curves based on P53.Mut status using previously created survival objects and generate Kaplan-Meier plots for both No RT and RT 
endo_fit_RT_Mut <- surv_fit(endo_surv_RT ~ P53.Mut, data = endometrial_all_RT)
KMRT_Mut <- ggsurvplot(endo_fit_RT_Mut, data = endometrial_all_RT, palette = alpha(c("royalblue", "firebrick2"), 0.6), risk.table = TRUE, 
           risk.table.height = 0.30, xlab = "Time in Months", xlim = c(0,150), legend.labs = c("p53 WT","p53 Mutant"), 
           title = "Progression-Free Survival (RT)")

endo_fit_NoRT_Mut <- surv_fit(endo_surv_NoRT ~ P53.Mut, data = endometrial_all_NoRT)
KMNoRT_Mut <- ggsurvplot(endo_fit_NoRT_Mut, data = endometrial_all_NoRT, palette = alpha(c("royalblue", "firebrick2"), 0.6), risk.table = TRUE, 
           risk.table.height = 0.30, xlab = "Time in Months", xlim = c(0,150), break.x.by = 50, legend.labs = c("p53 WT","p53 Mutant"), 
           title = "Progression-Free Survival (No RT)")

##Print graphs as jpeg files
jpeg("KM, RT, p53 status, no chemo or IV.jpeg", width= 1000, height= 720, units= "px", res= 150)
KMRT_Mut 
dev.off()

jpeg("KM, No RT, p53 status, no chemo or IV.jpeg", width= 1000, height= 720, units= "px", res= 150)
KMNoRT_Mut 
dev.off()



#Cox Hazards Analysis on RT, No RT and TP53 WT, Mutant (no multiples) (Figure 1 stats)
##Generate Cox analysis objects for No RT and RT based on P53.Mut status using previously created survival objects
endo_cox_RT_Mut <- coxph(endo_surv_RT ~ P53.Mut, data = endometrial_all_RT)

endo_cox_NoRT_Mut <- coxph(endo_surv_NoRT ~ P53.Mut, data = endometrial_all_NoRT)

##Generate forest plots and print plots as jpeg files
jpeg("Cox, RT, p53 status, no chemo or IV.jpeg", width= 720, height= 1000, units= "px", res= 150)
CoxRT_Mut <- ggforest(endo_cox_RT_Mut, data = endometrial_all_RT, main = "Progression-Free Survival (RT)")
CoxRT_Mut
dev.off()

jpeg("Cox, No RT, p53 status, no chemo or IV.jpeg", width= 720, height= 1000, units= "px", res= 150)
CoxNoRT_Mut <- ggforest(endo_cox_NoRT_Mut, data = endometrial_all_NoRT, main = "Progression-Free Survival (No RT)")
CoxNoRT_Mut
dev.off()





