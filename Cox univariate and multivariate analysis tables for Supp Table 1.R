#Load required packages
library(readxl)
library(survival)
library(survminer)
library(gtsummary)
library(tidyverse)
library(webshot)

#Set working directory
setwd("/Users/pettya/Desktop/R Analysis")


#Summary tables for univariate and multivariate Cox hazards analysis of Progression-Free Survival (Supp Figure 1B,C)
##Read in excel file with TCGA patient data table and create new Progression column
endometrial <- read_xlsx("Final TCGA with p53 and RT annotation - R.xlsx", 1)
endometrial$Progression.Free.Status <- ifelse(endometrial$Progression.Free.Status == "0:CENSORED", 0, 1)
names(endometrial)[21] <- "Progression"

##Create new categorical p53 column with VAF values based on Allelic Frequency data
endometrial$p53[endometrial$Allele.Freq..T.<= 0.50] <- "VAF Low"
endometrial$p53[endometrial$Allele.Freq..T.> 0.50] <- "VAF High"
endometrial$p53[endometrial$`P53 Mut`== 0] <- "WT"

##Redo Diagnosis Age column to contain 3 age range categories and convert to ordered factor
endometrial$Diagnosis.Age[endometrial$Diagnosis.Age < 50] <- "< 50"
endometrial$Diagnosis.Age[endometrial$Diagnosis.Age >= 50 & endometrial$Diagnosis.Age <= 70] <- "50-70"
endometrial$Diagnosis.Age[endometrial$Diagnosis.Age > 70] <- "> 70"
names(endometrial)[6] <- "Age"
endometrial$Age <- factor(endometrial$Age, levels= c("< 50","50-70","> 70"))

##Create new sCNA column with categorical values based on Fraction Genome Altered
endometrial$sCNA <- ifelse(endometrial$Fraction.Genome.Altered <= 0.20, "Low", "High")

##Change Grade column values
names(endometrial)[15] <- "Grade"
endometrial$Grade[endometrial$Grade == "G1"] <- "1"
endometrial$Grade[endometrial$Grade == "G2"] <- "2"
endometrial$Grade[endometrial$Grade == "G3"  | endometrial$Grade == "High Grade"] <- "3"

##Convert p53 and sCNA columns to ordered factors
endometrial$p53 <- factor(endometrial$p53, levels= c("WT","VAF Low","VAF High"))
endometrial$sCNA <- factor(endometrial$sCNA, levels= c("Low","High"))

##Order table rows by Progression-Free Survival values
endometrial <- endometrial[order(endometrial$Progress.Free.Survival..Months.),]

##Change Radiation and Stage column names
names(endometrial)[22] <- "Radiotherapy"
names(endometrial)[14] <- "Stage"

##Use gtsummary function tbl_uvregression to generate summary table of univariate Cox hazards analysis results (Supp Fugure 1B)
model <- endometrial %>% 
  select(sCNA, Age, p53, Radiotherapy, Stage, Grade, Progress.Free.Survival..Months., Progression) %>% 
  tbl_uvregression(method = coxph, y= Surv(time = Progress.Free.Survival..Months., event = Progression), exponentiate = TRUE, hide_n = TRUE) %>%
  bold_labels() %>% 
  italicize_levels() %>% 
  modify_spanning_header(everything() ~ "Progression-Free Survival")

##Save table as png file
gt::gtsave(as_gt(model), file = "Cox table vs all.png")


##Use gtsummary function tbl_regression to generate summary table of multivariate Cox hazards analysis results (Supp Fugure 1C)
model_mult <- coxph(Surv(time = endometrial$Progress.Free.Survival..Months., event = endometrial$Progression) ~ sCNA + Stage + Grade + p53, data = endometrial) %>% 
  tbl_regression(exponentiate = TRUE) %>%
  bold_labels() %>% 
  italicize_levels() %>% 
  modify_spanning_header(everything() ~ "Progression-Free Survival")

##Save table as png file
gt::gtsave(as_gt(model_mult), file = "Cox table multivariable.png")


