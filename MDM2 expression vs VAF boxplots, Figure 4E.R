#Load appropriate libraries
library(readxl)
library(openxlsx)
library(RColorBrewer)
library(tidyverse)
library(reshape2)

#Set working directory 
setwd("/Users/pettya/Desktop/R Analysis")

#Create new TCGA patient table with MDM2 expression data
##Read in excel file with all TCGA patient data
endometrial <- read_xlsx("Final TCGA with p53 and RT annotation - R.xlsx", 1)

##Read in text file with TCGA patient MDM2 expression data and only keep patient ids and MDM2 expression
MDM2 <- read.delim("MDM2__mRNA_expression_z-scores_relative_to_normal_samples_(log_RNA_Seq_V2_RSEM).txt", header = TRUE, sep = "\t")
MDM2 <- MDM2[,c(2,4)]

##Merge full TCGA table with MDM2 expression table
endometrial_MDM2 <- merge(endometrial, MDM2, by= "Patient.ID")
names(endometrial_MDM2)[52] <- "MDM2_Expression"

##Create categorical VAF column from Allelic Frequency data and convert to ordered factor
endometrial_MDM2$VAF[endometrial_MDM2$Allele.Freq..T. <= 0.5] <- "Low"
endometrial_MDM2$VAF[endometrial_MDM2$Allele.Freq..T. > 0.5] <- "High"
endometrial_MDM2$VAF[endometrial_MDM2$`P53 Mut` == 0] <- "WT"
endometrial_MDM2$VAF <- factor(endometrial_MDM2$VAF, levels = c("WT", "Low", "High"))



#One-way ANOVA and Tukey's Multiple Comparisons stats on MDM2 expression versus VAF
MDM2_VAF <- aov(MDM2_Expression ~ VAF, endometrial_MDM2)
MDM2_MC <- TukeyHSD(MDM2_VAF)
summary(MDM2_VAF)
MDM2_MC



#Boxplot of MDM2 expression versus TP53 VAF (Figure 4E)
##Create colored boxplot with black dots
qb <- ggplot(endometrial_MDM2,  aes(x= VAF, y= MDM2_Expression)) + theme_classic() +
  labs(x= NULL, y= "MDM2 Expression") + theme(axis.title.y = element_text(size= 18, face= "bold")) + 
  theme(axis.text.x = element_text(size = 18, face= "bold"), axis.text.y = element_text(size = 16, face= "bold")) +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", axis.line = element_line(size= 1))
qb <- qb + geom_boxplot(aes(fill= VAF), width= 0.6) + 
      geom_jitter(width= 0.05) + scale_fill_manual(values= alpha(c("royalblue", "springgreen4", "firebrick2"), 0.6))

##Add significance stars to bars
qb <- qb + annotate("text", x= 2:3, y= 4.5, label= "***", size= 8, colour= "brown")  

##Print graph as jpeg file
jpeg("VAF vs MDM2 expression boxplot, box color.jpeg", width= 720, height= 720, units= "px", res= 150)
print(qb)

dev.off()

##Create gray boxplot with colored dots
qb <- ggplot(endometrial_MDM2,  aes(x= VAF, y= MDM2_Expression)) + theme_classic() +
  labs(x= NULL, y= "MDM2 Expression") + theme(axis.title.y = element_text(size= 18, face= "bold")) + 
  theme(axis.text.x = element_text(size = 18, face= "bold"), axis.text.y = element_text(size = 16, face= "bold")) +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", axis.line = element_line(size= 1))
qb <- qb + geom_boxplot(fill= "grey", width= 0.6, outlier.shape = NA) + 
  geom_jitter(aes(color= VAF), width= 0.05) + scale_color_manual(values= alpha(c("royalblue", "springgreen4", "firebrick2"), 0.6))
qb <- qb + annotate("text", x= 2:3, y= 4.5, label= "***", size= 8, colour= "brown")  

##Print graph as jpeg file
jpeg("VAF vs MDM2 expression boxplot, point color, 6-23-22.jpeg", width= 720, height= 720, units= "px", res= 150)
print(qb)

dev.off()





