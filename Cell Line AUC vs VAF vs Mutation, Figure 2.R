#Load required packages
library(readxl)
library(openxlsx)
library(RColorBrewer)
library(tidyverse)
library(grid)

##Set the working directory
setwd("/Users/pettya/Desktop/R Analysis")

endometrial <- read_xlsx("FSCNA versus AUC Cells -- Annotated for AF.xlsx", 1)
endometrial <- endometrial[, c(2, 8, 14:16)]
endometrial_order <- endometrial[order(endometrial$AUC),]
names(endometrial_order)[4] <- "Mutation_Type"
names(endometrial_order)[5] <- "Protein_Change"
names(endometrial_order)[1] <- "Cell_Line"
endometrial_order$Cell_Line[endometrial_order$Cell_Line == "IshikawaHeraklio02ER"] <- "Ishikawa"
endometrial_order$Mutation_Type[endometrial_order$Mutation_Type == "X" | endometrial_order$Mutation_Type == "x"] <- "WT"
endometrial_order$Mutation_Type[endometrial_order$Mutation_Type == "Missense/Nonsense" |
                                endometrial_order$Mutation_Type == "Missense/Fsdel" | endometrial_order$Mutation_Type == "Mis/Mis" | 
                                endometrial_order$Mutation_Type == "Misense"] <- "Missense"
endometrial_order$Mutation_Type[endometrial_order$Mutation_Type == "Non/Miss" | 
                                endometrial_order$Mutation_Type == "Nonsense" | 
                                endometrial_order$Mutation_Type == "Splice Site" | 
                                endometrial_order$Mutation_Type == "Deletion" ] <- "Other"


endometrial_order$AF[endometrial_order$AF == "NA"] <- NA
endometrial_order <- drop_na(endometrial_order, AF)
endometrial_order$AF[endometrial_order$AF == "0.53/0.47"] <- "0.53"
endometrial_order$AF[endometrial_order$AF == "0.43/0.88"] <- "0.88"                     
endometrial_order$AF[endometrial_order$AF == "0.5/0.5"] <- "0.5"
endometrial_order$AF[endometrial_order$AF == "0.53/0.5"] <- "0.53"                     
endometrial_order$AF[endometrial_order$AF == "0.61/0.34"] <- "0.61"                     
endometrial_order$AF <- as.numeric(endometrial_order$AF)                    
endometrial_order <- within(endometrial_order, {
                            VAF <- NA
                            VAF[AF == 0.00] <- "WT"
                            VAF[AF <= 0.50 & AF != 0.00] <- "Low"
                            VAF[AF > 0.50] <- "High"
})                     
endometrial_order$VAF <- factor(endometrial_order$VAF, levels= c("WT", "Low", "High"))
endometrial_order$Mutation_Type <- factor(endometrial_order$Mutation_Type, levels= c("WT", "Other", "Missense"))
endometrial_order$AUC <- scale(endometrial_order$AUC)  
endometrial_order[,7] <- "AUC"
endometrial_order[,8] <- "VAF"
endometrial_order[,9] <- "Mutation"
minAUC <- min(endometrial_order$AUC)
maxAUC <- max(endometrial_order$AUC)
##AUC heatmap with Red and Blue
AUC <- ggplot(endometrial_order, aes(x= reorder(Cell_Line, -AUC), y= ...7)) +
  geom_tile(aes(fill= AUC)) + scale_fill_gradientn(limits= c(minAUC, maxAUC), colours= c("blue", "red")) +
  theme_classic() + theme(axis.title.y = element_blank(), axis.title.x= element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) + theme(legend.title = element_blank(), legend.direction = "horizontal") + theme(plot.margin = unit(c(0,0,0,2), "mm"))

##AUC heatmap with colorblind friendly colors
AUC <- ggplot(endometrial_order, aes(x= reorder(Cell_Line, -AUC), y= ...7)) +
  geom_tile(aes(fill= AUC)) + scale_fill_gradientn(limits= c(minAUC, maxAUC), colours= c("#2c7bb6", "#d7191c")) +
  theme_classic() + theme(axis.title.y = element_blank(), axis.title.x= element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size= 16, face= "bold")) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) + theme(legend.title = element_blank(), legend.direction = "horizontal", legend.text = element_text(size= 16)) + theme(plot.margin = unit(c(0,0,0,2), "mm"))

##VAF heatmap with Blue, Green, Red colors  
VAF <- ggplot(endometrial_order, aes(x= reorder(Cell_Line, -AUC), y= ...8)) + 
      geom_tile(aes(fill= VAF)) + scale_fill_manual(values= c("blue", "green", "red")) + theme_classic() + theme(axis.title.y = element_blank(), axis.title.x= element_blank(), axis.ticks.y = element_blank()) +
      theme(axis.text.x= element_blank(), axis.ticks.x= element_blank()) + theme(axis.line.x = element_blank(), axis.line.y = element_blank()) + theme(legend.title = element_blank(), legend.direction = "horizontal") + theme(plot.margin = unit(c(0,0,0,2), "mm"))

##VAF heatmap with colorblind friendly colors
VAF <- ggplot(endometrial_order, aes(x= reorder(Cell_Line, -AUC), y= ...8)) + 
  geom_tile(aes(fill= VAF)) + scale_fill_manual(values= c("#2c7bb6", "#fdae61", "#d7191c")) + theme_classic() + theme(axis.title.y = element_blank(), axis.title.x= element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size= 16, face= "bold")) +
  theme(axis.text.x= element_blank(), axis.ticks.x= element_blank()) + theme(axis.line.x = element_blank(), axis.line.y = element_blank()) + theme(legend.title = element_blank(), legend.direction = "horizontal", legend.text = element_text(size= 16)) + theme(plot.margin = unit(c(0,0,0,2), "mm"))

VAF2 <- ggplot(endometrial_order, aes(x= reorder(Cell_Line, -AUC), y= ...8)) + 
  geom_tile(aes(fill= VAF)) + scale_fill_manual(values= c("blue", "green", "red")) + theme_classic() + theme(axis.title.y = element_blank(), axis.text.y= element_blank(), axis.title.x= element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.text.x = element_text(angle= 90, size= 12, face= "bold", hjust= 0.5, vjust= 0.5)) + theme(axis.line.x = element_blank(), axis.line.y = element_blank()) + theme(legend.position = "none") + theme(plot.margin = unit(c(0,0,0,2), "mm"))

##Mutation heatmap with Blue, Red, Orange colors 
Mutation <- ggplot(endometrial_order, aes(x= reorder(Cell_Line, -AUC), y= ...9)) + 
  geom_tile(aes(fill= Mutation_Type)) + scale_fill_manual(values= c("blue", "red", "orange")) + theme_classic() + theme(axis.title.y = element_blank(), axis.title.x= element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) + theme(legend.title = element_blank(), legend.direction = "horizontal", legend.text = element_text(size= 8)) + theme(plot.margin = unit(c(0,0,0,2), "mm"))

##Mutation heatmap with colorblind friendly colors
Mutation <- ggplot(endometrial_order, aes(x= reorder(Cell_Line, -AUC), y= ...9)) + 
  geom_tile(aes(fill= Mutation_Type)) + scale_fill_manual(values= c("#2c7bb6", "#fdae61", "#d7191c")) + theme_classic() + theme(axis.title.y = element_blank(), axis.title.x= element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank()) + theme(axis.text.y= element_text(size= 16, face= "bold"), legend.title = element_blank(), legend.direction = "horizontal", legend.text = element_text(size= 16)) + theme(plot.margin = unit(c(0,0,0,2), "mm"))

grid.newpage()
print(VAF, vp = viewport(x = 0.522, y = 0.77, width = 0.99, height = 0.05))
print(VAF2, vp = viewport(x = 0.363, y = 0.68, width = 0.54, height = 0.05))
print(AUC, vp = viewport(x = 0.44, y = 0.89, width = 0.836, height = 0.05))
print(Mutation, vp = viewport(x = 0.494, y = 0.83, width = 0.946, height = 0.05))


jpeg("Cell line heatmap, AUC vs Mutation vs VAF, 6-29-22.jpeg", width= 1200, units= "px", res= 120)
grid.newpage()
print(VAF, vp = viewport(x = 0.4875, y = 0.63, width = 0.836, height = 0.12))
print(VAF2, vp = viewport(x = 0.376, y = 0.46, width = 0.52, height = 0.12))
print(AUC, vp = viewport(x = 0.43, y = 0.89, width = 0.73, height = 0.12))
print(Mutation, vp = viewport(x = 0.495, y = 0.76, width = 0.947, height = 0.12))

dev.off()

