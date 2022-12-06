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

##Group cell lines by variant and calculate Mean, SD, and SEM
Grouped <- JHLog %>%
  group_by(variable) %>%
  summarise(Mean= mean(value), SD= sd(value))

Grouped <- Grouped %>%
  mutate(SEM= SD/(sqrt(4)))


##Order rows by decreasing mean Log2 AUC/mCherry value and number rows from 1-6
JHLogorder <- Grouped[order(Grouped$Mean, decreasing = TRUE),]
x <- 1:nrow(JHLogorder)
names(JHLogorder)[1] <- "Variant"

##Convert Variant column to a factor and set order of variants
JHLogorder$Variant <- factor(JHLogorder$Variant, 
                             levels= c("R248Q", "R248W", "R273C", "Y220C", "R273H", "WT"))

##Generate waterfall plot with x axis set as row number
J <- ggplot(JHLogorder, aes(x= x, y= Mean)) + theme_classic() +
  labs(title= NULL, 
       x= "TP53 Allele Co-Expressed with Wild Type", y= expression(paste("Change in Radiation Response (Mean ", Delta, "AUC)"))) +
  theme(axis.ticks.x= element_blank(), axis.text.x = element_blank(), axis.title.x = element_text(size= 22, face= "bold")) +
  theme(axis.text.y = element_text(size= 18, face= "bold"), axis.title.y = element_text(size= 24, face= "bold")) +
  theme(legend.text = element_text(size= 16), legend.title = element_text(size= 18), plot.title = element_text(size= 18), axis.line = element_line(size= 1))

J <- J + geom_bar(aes(fill= Variant, color= Variant), stat = "identity", width= 0.7, 
                  position= position_dodge(width = 0.4)) + scale_fill_brewer(palette = "RdYlBu") +
  scale_color_brewer(palette = "RdYlBu")

J <- J + geom_errorbar(aes(ymin= Mean - SEM, ymax= Mean + SEM), width= 0.3)

##Save waterfall plot as jpeg file
jpeg("JHUEM Mean Waterfall plot, RdYlBu, no annotation.jpeg", width= 1100, height= 1300, units= "px", res= 150)
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
JHLog$KO <- factor(JHLog$KO, levels= c("NTC", "KO 5.1", "KO 6.1"))

JSw <- ggplot(JHLog, aes(x= KO, y= value)) + 
  theme_classic() + scale_y_continuous(limits= c(0, 5), breaks= seq(0, 5, 1)) +
  scale_x_discrete(expand= c(0.5,0.5)) + labs(title= NULL, x= NULL, y= "Radiation Response (Mean AUC)") +
  theme(axis.text.y= element_text(face= "bold", size= 16), axis.title.y= element_text(size= 18, face= "bold"), axis.text.x = element_text(face = "bold", size= 18)) +
  theme(legend.position= "none", axis.line = element_line(size= 1))
JSw <- JSw + geom_boxplot(fill= "grey", width= 0.4)
JSw <- JSw + geom_beeswarm(aes(color= KO), priority = "density", cex= 3) + scale_color_manual(values = alpha(c("royalblue", "springgreen4", "firebrick2"), 0.6))

##Annotate with significance stars from Kruskal-Wallis test and Dunn's Multiple Comparisons from Prism
JSw <- JSw + annotate("text", x= 2, y= 6, label= "***", size= 8, colour= "brown") + 
  theme(panel.grid.minor.y = element_blank())
JSw <- JSw + annotate("text", x= 3, y= 6, label= "**", size= 8, colour= "brown") + 
  theme(panel.grid.minor.y = element_blank())

##Save plot as jpeg file
jpeg("JHUEM AUC boxplot swarm.jpeg", width= 1000, height= 720, units= "px", res= 150)
print(JSw)

dev.off()




















