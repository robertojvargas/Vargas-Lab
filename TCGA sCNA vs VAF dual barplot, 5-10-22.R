##Load required packages
library(readxl)
library(openxlsx)
library(RColorBrewer)
library(tidyverse)
library(grid)

##Set the working directory
setwd("/Users/pettya/Desktop/R Analysis")

##Read in filtered (no Stage IV or chemo) and compiled endometrial TP53 TCGA dataset (No RT and RT)
endometrial <- read_xlsx("Final TCGA p53 patients, no chemo or Stage IV, 2-22-21.xlsx", 1)

##Order by increasing Fraction Genome Altered (sCNA Fraction)
endometrial_order <- endometrial[order(endometrial$Fraction.Genome.Altered),]
##Convert Allelic Frequency to numeric, create new variable VAF, and convert VAF to factor
endometrial_order$Allele.Freq..T.[endometrial_order$Allele.Freq..T. == "WT"] <- 0
endometrial_order$Allele.Freq..T. <- as.numeric(endometrial_order$Allele.Freq..T.)
endometrial_order <- within(endometrial_order, {
                            VAF <- NA
                            VAF[Allele.Freq..T. <= 0.50 ] <- "Low"
                            VAF[Allele.Freq..T. > 0.50] <- "High"
})
endometrial_order$VAF <- factor(endometrial_order$VAF, levels= c("Low", "High"))

##Create new variabe x with patients numbered 1-351
x <- 1:nrow(endometrial_order)

##Create sCNA Fraction bar plot of patients
J <- ggplot(endometrial_order, aes(x= x, y= Fraction.Genome.Altered)) + theme_classic() +
     labs(title= "TCGA Distribution of Fraction sCNA and TP53 Variant Allelic Frequency", x= NULL, y= "sCNA Fraction") +
     theme(axis.ticks.x= element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size= 16, face= "bold"), axis.title.y = element_text(size= 18, face= "bold")) +
     theme(axis.line = element_line(size= 1))
J <- J + geom_bar(stat = "identity", width= 0.4) + 
        scale_y_continuous(breaks= c(0.2, 1)) + 
        theme(panel.grid.major.y = element_line(size= 1, color = alpha("royalblue", 0.6), linetype = "dashed"), plot.margin = unit(c(0,0,0,0), "mm"))

##Create inverted VAF barplot of patients with same order
J2 <- ggplot(endometrial_order, aes(x= x, y= Allele.Freq..T.)) + theme_classic() +
  labs(title= NULL, x= NULL, y= "TP53 VAF") +
  theme(axis.ticks.x= element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank(), axis.text.y = element_text(size= 16, face= "bold"), axis.title.y = element_text(size= 18, face= "bold")) +
  theme(axis.line.y= element_line(size= 1))
J2 <- J2 + geom_bar(aes(fill= VAF), stat = "identity", width= 0.4) + 
  scale_y_continuous(breaks= c(0.5), trans = "reverse") + 
  theme(panel.grid.major.y = element_line(size= 1, color = alpha("royalblue", 0.6), linetype = "dashed"), plot.margin = unit(c(0,0,0,0), "mm")) + 
  theme(legend.position = c(0.3, 0.1), legend.direction = "horizontal", legend.text = element_text(size= 14), legend.title = element_text(size= 16)) + 
  scale_fill_manual(values= alpha(c("springgreen4", "firebrick2"), 0.6))

##Combine graphs
grid.newpage()
print(J2, vp = viewport(x = 0.5, y = 0.37, width = 0.8, height = 0.4))
print(J, vp = viewport(x = 0.5, y = 0.75, width = 0.8, height = 0.4))

##Combine graphs and print to jpeg file
jpeg("Patient sCNA vs TP53 VAF, 6-24-22.jpeg", height= 500, width= 1300, units= "px", res= 120)
grid.newpage()
print(J2, vp = viewport(x = 0.45, y = 0.37, width = 0.8, height = 0.4))
print(J, vp = viewport(x = 0.45, y = 0.75, width = 0.8, height = 0.4))

dev.off()




