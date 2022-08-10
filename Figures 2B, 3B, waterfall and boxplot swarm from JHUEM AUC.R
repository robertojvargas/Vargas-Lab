#Load packages
library(readxl)
library(openxlsx)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(ggbeeswarm)

#Set the working directory
setwd("/Users/pettya/Desktop/R Analysis")



#Waterfall plot of Log2 AUC/mCherry for JHUEM1 and JHUEM2 (Figure 3B)
##Read in excel file with JHUEM1 and JHUEM2 variant Log2 AUC/mCherry data in long format
JHLog <- read_xlsx("JHUEM1, 2 NTC + variants, Log2 AUC data, 11-30-21.xlsx", 1)

##Pivot data table and aggregrate with mean Log2 AUC/mCherry of each variant for each cell line
JHLog <- melt(JHLog, id= c('Cell', 'Variant'))
JHLog <- dcast(JHLog, Cell ~ Variant, mean)

##Pivot data again to long form  
JHLog <- melt(JHLog, id= c("Cell"))
names(JHLog)[2] <- "Variant"

##Order rows by decreasing mean Log2 AUC/mCherry value and number rows from 1-24
JHLogorder <- JHLog[order(JHLog$value, decreasing = TRUE),]
x <- 1:nrow(JHLogorder)

##Convert Variant column to a factor and set order of variants
JHLogorder$Variant <- factor(JHLogorder$Variant, 
                             levels= c("R248Q", "R248W", "R273C", "R273H", "Y220C", "WT"))

##Generate waterfall plot with x axis set as row number
J <- ggplot(JHLogorder, aes(x= x, y= value)) + theme_classic() +
    labs(title= NULL, 
          x= NULL, y= expression(paste("Log2 ", Delta, "AUC"))) +
    theme(axis.ticks.x= element_blank(), axis.text.x = element_blank()) +
    theme(axis.text.y = element_text(size= 18, face= "bold"), axis.title.y = element_text(size= 22, face= "bold")) +
    theme(legend.text = element_text(size= 16), legend.title = element_text(size= 18), plot.title = element_text(size= 18), axis.line = element_line(size= 1))

J <- J + geom_bar(aes(fill= Variant, color= Variant), stat = "identity", width= 0.7, 
    position= position_dodge(width = 0.4)) + scale_fill_brewer(palette = "Accent") +
    scale_color_brewer(palette = "Accent")

##Save waterfall plot as jpeg file
jpeg("JHUEM Waterfall plot, Accent, no annotation.jpeg", width= 1100, height= 1300, units= "px", res= 150)
print(J)

dev.off()
   


#Boxplot and swarmplot for JHUEM1, JHUEM2 AUC data (Figure 2B)
##Read in excel file with JHUEM1 and JHUEM2 variant AUC data in long format
JHLog <- read_xlsx("JHUEM1, 2 NTC + variants, Log2 AUC data, 11-30-21.xlsx", 3)

##Delete Cell Line column and rename new first column
JHLog <- JHLog[,-1]
names(JHLog)[1] <- "KO"

##Pivot table to long form and convert KO column to factor
JHLog <- melt(JHLog, id= c("KO"))
JHLog$KO <- factor(JHLog$KO, levels= c("KO 6.1", "KO 5.1", "NTC"))

##Generate boxplot and swarmplot of AUC data
JSw <- ggplot(JHLog, aes(x= value, y= KO)) + 
      theme_classic() + scale_x_continuous(limits= c(0, 6), breaks= seq(0, 6, 1)) +
      labs(title= NULL, x= "Mean AUC", y= NULL) +
      theme(axis.text.x= element_text(face= "bold", size= 16), axis.title.x= element_text(size= 18, face= "bold"), axis.text.y = element_text(face = "bold", size= 18)) +
      theme(legend.position= "none", axis.line = element_line(size= 1))
JSw <- JSw + geom_boxplot(fill= "grey", width= 0.3)
JSw <- JSw + geom_beeswarm(aes(color= KO), priority = "density", cex= 3) + scale_color_manual(values = alpha(c("firebrick2", "springgreen4", "royalblue"), 0.6))

##Annotate with significance stars from Kruskal-Wallis test and Dunn's Multiple Comparisons from Prism
JSw <- JSw + annotate("text", x= 6, y= 2, label= "***", size= 8, colour= "brown") + 
  theme(panel.grid.minor.y = element_blank())
JSw <- JSw + annotate("text", x= 6, y= 1, label= "**", size= 8, colour= "brown") + 
  theme(panel.grid.minor.y = element_blank())

##Save plot as jpeg file
jpeg("JHUEM AUC boxplot swarm.jpeg", width= 1000, height= 720, units= "px", res= 150)
print(JSw)

dev.off()








