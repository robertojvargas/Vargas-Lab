##Load required libraries
library(readxl)
library(openxlsx)
library(RColorBrewer)
library(tidyverse)

##Set the working directory
setwd("/Users/pettya/Desktop/R Analysis")

##Read in filtered (no Stage IV or chemo) and compiled endometrial TP53 TCGA dataset (No RT and RT)
endometrial <- read_xlsx("Final TCGA p53 patients, no chemo or Stage IV, 2-22-21.xlsx", 1)

##Change RT variable and convert to factor
endometrial$Radiation.Therapy <- ifelse(endometrial$Radiation.Therapy == "Yes", "RT", "No RT")
names(endometrial)[23] <- "RT"
endometrial$RT <- factor(endometrial$RT)

##Remove TP53 WT patients and convert Allelic Frequency to numeric
endometrial <- endometrial[which(endometrial$Allele.Freq..T. != "WT"),]
endometrial$Allele.Freq..T. <- as.numeric(endometrial$Allele.Freq..T.)

##Violin plot of VAF vs RT
pVAF <- ggplot(endometrial, aes(x= RT, y= Allele.Freq..T.))
pVAF <- pVAF + geom_violin(aes(fill= RT), draw_quantiles = c(0.5), scale= "count", adjust= 0.3) + theme_classic() + labs(title= NULL, x= NULL, y= "VAF")
pVAF <- pVAF + geom_jitter(width= 0.05)  
pVAF <- pVAF + theme(panel.grid.minor.y = element_blank()) +
  theme(legend.position = "none") + scale_fill_manual(values= alpha(c("royalblue", "firebrick2"), 0.6)) +
  theme(axis.title.y = element_text(size= 18, face= "bold"), axis.line= element_line(size= 1)) + 
  theme(axis.text.x = element_text(size = 18, face= "bold"), axis.text.y = element_text(size = 16, face= "bold")) 

##Save violin plot to jpeg file
jpeg("VAF vs RT violin, TCGA, white background, 6-21-22.jpeg", width= 720, height= 720, units= "px", res= 150)
print(pVAF)
dev.off()





##Boxplot of VAF vs RT
jpeg("VAF vs RT, TCGA, gray, white background, 5-26-22.jpeg", width= 720, height= 720, units= "px", res= 150)
pVAF <- ggplot(endometrial, aes(x= RT, y= Allele.Freq..T.))
pVAF <- pVAF + geom_boxplot(fill= "gray") + theme_classic() + labs(title= NULL, x= NULL, y= "VAF")
pVAF <- pVAF + theme(panel.grid.minor.y = element_blank()) +
        theme(legend.position = "none")
pVAF <- pVAF + geom_jitter()
pVAF

dev.off()


##Boxplot of sCNA vs RT
pCN <- qplot(RT, Fraction.Genome.Altered, data= endometrial, geom= c("boxplot", "jitter"), 
              fill= RT, xlab = "RT", ylab = "sCNA") + scale_fill_manual(values= c("Blue", "Red"))
pCN <- pCN + theme_classic()
pCN <- pCN + theme(axis.text.x = element_blank()) + 
  theme(panel.grid.minor.y = element_blank()) + theme(legend.title = element_blank()) +
  theme(axis.title.y = element_text(size= 12)) + theme(axis.title.x = element_text(size= 11)) + 
  theme(legend.text = element_text(size= 12))
pCN <- pCN + scale_y_continuous(limits = c(0, 1.00))
  


  



