##Load necessary libraries
library("reshape2")
library("openxlsx")
library("readxl")
library("tidyverse")
library("BIGL")
##Make sure files are read from R Analysis folder
setwd("/Users/pettya/Desktop/R Analysis")


#Analyzing AMG-232 + radiation combination experiments

##Data prep for synergy analyses
###Make list of sheetnames from Excel file and loop through to create single dataframe for each sheet 
sheetnumber <- seq(1:9)
sheet_names <- c()
for (i in 1:length(sheetnumber)){
  sheet_names <- append(sheet_names, (paste0("Sheet",i)))
}
AUC <- list()

###Read in AMG-232 + Radiation data file
for (i in 1:length(sheet_names)){
  AUC[[i]] <-read_xlsx("JHUEM2, Hec108, Hec1B -AMG232 + Rad, AP46, 5-7-21.xlsx", sheet_names[i])
}
###Extract list of radiation doses and reverse order for use later
Rad <- as.character(rev(unique(AUC[[1]]$Radiation)))

###Loop through each dataframe, take the mean of 6 replicates, and restructure, making cell names into rownames 
for (i in 1:length(AUC)){
  AUC[[i]] <- melt(AUC[[i]], id = c("Dose", "Radiation"))
  AUC[[i]] <- dcast(AUC[[i]], Dose ~ Radiation, mean)
  row.names(AUC[[i]]) <- AUC[[i]]$Dose
  AUC[[i]] <- AUC[[i]][,-1]
}

###Loop through data tables and normalize each response value to Untreated sample
for (i in sheetnumber){
  AUC[[i]] <- AUC[[i]]/AUC[[i]]["0", "0"]
}

###For JHUEM2, Hec108 and Radiation doses of 2, 4, 8 Gy, divide fraction surviving by 2 because plated at twice the density
for (i in Rad){
  if (i == "2" | i == "4" | i == "8"){
    AUC[[1]][i] <- (AUC[[1]][i])/2
    AUC[[2]][i] <- (AUC[[2]][i])/2
    AUC[[3]][i] <- (AUC[[3]][i])/2
    AUC[[4]][i] <- (AUC[[4]][i])/2
    AUC[[5]][i] <- (AUC[[5]][i])/2
    AUC[[6]][i] <- (AUC[[6]][i])/2
  }
  else {
    AUC[[1]][i] <- AUC[[1]][i]
    AUC[[2]][i] <- AUC[[2]][i]
    AUC[[3]][i] <- AUC[[3]][i]
    AUC[[4]][i] <- AUC[[4]][i]
    AUC[[5]][i] <- AUC[[5]][i]
    AUC[[6]][i] <- AUC[[6]][i]
  }
}

###Transpose each table and extract to separate dataframe
AUC_1 <- as.data.frame(t(AUC[[1]]))
AUC_2 <- as.data.frame(t(AUC[[2]]))
AUC_3 <- as.data.frame(t(AUC[[3]]))
AUC_4 <- as.data.frame(t(AUC[[4]]))
AUC_5 <- as.data.frame(t(AUC[[5]]))
AUC_6 <- as.data.frame(t(AUC[[6]]))
AUC_7 <- as.data.frame(t(AUC[[7]]))
AUC_8 <- as.data.frame(t(AUC[[8]]))
AUC_9 <- as.data.frame(t(AUC[[9]]))

###Create new Radiation column for Excel labels
AUC_1$Radiation <- rownames(AUC_1)
AUC_2$Radiation <- rownames(AUC_2)
AUC_3$Radiation <- rownames(AUC_3)
AUC_4$Radiation <- rownames(AUC_4)
AUC_5$Radiation <- rownames(AUC_5)
AUC_6$Radiation <- rownames(AUC_6)
AUC_7$Radiation <- rownames(AUC_7)
AUC_8$Radiation <- rownames(AUC_8)
AUC_9$Radiation <- rownames(AUC_9)

###Pivot tables to long form and combine replicates for cell lines
AUC_1 <- melt(AUC_1, id= c("Radiation"))
AUC_2 <- melt(AUC_2, id= c("Radiation"))
AUC_3 <- melt(AUC_3, id= c("Radiation"))
AUC_4 <- melt(AUC_4, id= c("Radiation"))
AUC_5 <- melt(AUC_5, id= c("Radiation"))
AUC_6 <- melt(AUC_6, id= c("Radiation"))
AUC_7 <- melt(AUC_7, id= c("Radiation"))
AUC_8 <- melt(AUC_8, id= c("Radiation"))
AUC_9 <- melt(AUC_9, id= c("Radiation"))
JHAUC <- rbind(AUC_1, AUC_2, AUC_3)
Hec8AUC <- rbind(AUC_4, AUC_5, AUC_6)
HecBAUC <- rbind(AUC_7, AUC_8, AUC_9)

###Aggregate to obtain the mean signal for each dose combination
JHMean <- dcast(JHAUC, Radiation ~ variable, mean)
Hec8Mean <- dcast(Hec8AUC, Radiation ~ variable, mean)
HecBMean <- dcast(HecBAUC, Radiation ~ variable, mean)

###Normalize all values to the max value for Fraction Surviving
JHMean <- melt(JHMean, id= c("Radiation"))
JHMean$value <- JHMean$value/max(JHMean$value)
Hec8Mean <- melt(Hec8Mean, id= c("Radiation"))
Hec8Mean$value <- Hec8Mean$value/max(Hec8Mean$value)
HecBMean <- melt(HecBMean, id= c("Radiation"))
HecBMean$value <- HecBMean$value/max(HecBMean$value)

###Create new Cell column with labels and rename remaining columns
JHMean$Cell <- "JHUEM2"
names(JHMean)[1:3] <- c("d2", "d1", "effect")
Hec8Mean$Cell <- "Hec108"
names(Hec8Mean)[1:3] <- c("d2", "d1", "effect")
HecBMean$Cell <- "Hec1B"
names(HecBMean)[1:3] <- c("d2", "d1", "effect")

###Combine all dataframes together
AllMean <- rbind(JHMean, Hec8Mean, HecBMean)

###Write all dataframes to excel file
Sheets <- list("All"= AllMean, "JHUEM2"= JHMean, "Hec108" = Hec8Mean, "Hec1B"= HecBMean)
write.xlsx(Sheets, "Synergy JHUEM2, Hec108, Hec1B -AMG232 + Rad, AP46.xlsx")


##Synergy analysis with BIGL
###Read in excel files containing AMG + Radiation responses for JHUEM2, Hec108, Hec1B
AMGSyn <- read_xlsx("Synergy JHUEM2, Hec108, Hec1B -AMG232 + Rad, AP46.xlsx", 1)

###Create separate AMG + Radiation tables for each cell line and rename columns for BIGL algorithms
JHAMG <- subset(AMGSyn, Cell== "JHUEM2", select = c("Cell", "d2", "d1", "effect"))
Hec8AMG <- subset(AMGSyn, Cell== "Hec108", select = c("Cell", "d2", "d1", "effect"))
HecBAMG <- subset(AMGSyn, Cell== "Hec1B", select = c("Cell", "d2", "d1", "effect"))

###Convert AMG and Radiation doses to numeric
JHAMG$d1 <- as.numeric(JHAMG$d1)
Hec8AMG$d1 <- as.numeric(Hec8AMG$d1)
HecBAMG$d1 <- as.numeric(HecBAMG$d1)
JHAMG$d2 <- as.numeric(JHAMG$d2)
Hec8AMG$d2 <- as.numeric(Hec8AMG$d2)
HecBAMG$d2 <- as.numeric(HecBAMG$d2)

###Fit marginal dose response curves to data using fitMarginal function in BIGL
Marg1Fit <- fitMarginals(JHAMG, method= "nlslm", names= c("AMG-232", "Radiation"))
Marg2Fit <- fitMarginals(Hec8AMG, method= "nlslm", names= c("AMG-232", "Radiation"))
Marg3Fit <- fitMarginals(HecBAMG, method= "nlslm", names= c("AMG-232", "Radiation"))

###Use fitSurface function in BIGL to fit predicted surface using Highest Single Agent (HSA) and Loewe models
JHAMGSurf <- fitSurface(JHAMG, Marg1Fit, null_model = "hsa", B.CP = 50, statistic = "both")
Hec8AMGSurf <- fitSurface(Hec8AMG, Marg2Fit, null_model = "hsa", B.CP = 50, statistic = "both")
HecBAMGSurf <- fitSurface(HecBAMG, Marg3Fit, null_model = "loewe", B.CP = 50, statistic = "both")

###Generate AMG isobolograms of predicted surfaces using isobologram function in BIGL and save as jpegs
jpeg("JHUEM2 + AMG isobologram, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(JHAMGSurf)
dev.off()

jpeg("Hec108 + AMG isobologram, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(Hec8AMGSurf)
dev.off()

jpeg("Hec1B + AMG isobologram, Loewe.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(HecBAMGSurf)
dev.off()

###Generate contour plots of AMG synergy scores using contour function in BIGL and save as jpegs
jpeg("JHUEM2 + AMG contour, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(JHAMGSurf, main= "Contour Plot for JHUEM2 MaxR")
dev.off()

jpeg("Hec108 + AMG contour, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(Hec8AMGSurf, main= "Contour Plot for Hec108 MaxR")
dev.off()

jpeg("Hec1B + AMG contour, Loewe.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(HecBAMGSurf, main= "Contour Plot for Hec1B MaxR")
dev.off()




#Analyzing Nutlin + radiation combination experiments

##Data prep for synergy analyses
###Read in Nutlin + Radiation data file
for (i in 1:length(sheet_names)){
  AUC[[i]] <-read_xlsx("JHUEM2, Hec108, Hec1B -Nutlin + Rad, AP47, 5-14-21.xlsx", sheet_names[i])
}
###Extract list of radiation doses and reverse order for use later
Rad <- as.character(rev(unique(AUC[[1]]$Radiation)))

###Loop through each dataframe, take the mean of 6 replicates, and restructure, making cell names into rownames 
for (i in 1:length(AUC)){
  AUC[[i]] <- melt(AUC[[i]], id = c("Dose", "Radiation"))
  AUC[[i]] <- dcast(AUC[[i]], Dose ~ Radiation, mean)
  row.names(AUC[[i]]) <- AUC[[i]]$Dose
  AUC[[i]] <- AUC[[i]][,-1]
}

###Loop through data tables and normalize each response value to Untreated sample
for (i in sheetnumber){
  AUC[[i]] <- AUC[[i]]/AUC[[i]]["0", "0"]
}

###For JHUEM2, Hec108 and Radiation doses of 2, 4, 8 Gy, divide fraction surviving by 2 because plated at twice the density
for (i in Rad){
  if (i == "2" | i == "4" | i == "8"){
    AUC[[1]][i] <- (AUC[[1]][i])/2
    AUC[[2]][i] <- (AUC[[2]][i])/2
    AUC[[3]][i] <- (AUC[[3]][i])/2
    AUC[[4]][i] <- (AUC[[4]][i])/2
    AUC[[5]][i] <- (AUC[[5]][i])/2
    AUC[[6]][i] <- (AUC[[6]][i])/2
  }
  else {
    AUC[[1]][i] <- AUC[[1]][i]
    AUC[[2]][i] <- AUC[[2]][i]
    AUC[[3]][i] <- AUC[[3]][i]
    AUC[[4]][i] <- AUC[[4]][i]
    AUC[[5]][i] <- AUC[[5]][i]
    AUC[[6]][i] <- AUC[[6]][i]
  }
}

###Transpose each table and extract to separate dataframe
AUC_1 <- as.data.frame(t(AUC[[1]]))
AUC_2 <- as.data.frame(t(AUC[[2]]))
AUC_3 <- as.data.frame(t(AUC[[3]]))
AUC_4 <- as.data.frame(t(AUC[[4]]))
AUC_5 <- as.data.frame(t(AUC[[5]]))
AUC_6 <- as.data.frame(t(AUC[[6]]))
AUC_7 <- as.data.frame(t(AUC[[7]]))
AUC_8 <- as.data.frame(t(AUC[[8]]))
AUC_9 <- as.data.frame(t(AUC[[9]]))

###Create new Radiation column for Excel labels
AUC_1$Radiation <- rownames(AUC_1)
AUC_2$Radiation <- rownames(AUC_2)
AUC_3$Radiation <- rownames(AUC_3)
AUC_4$Radiation <- rownames(AUC_4)
AUC_5$Radiation <- rownames(AUC_5)
AUC_6$Radiation <- rownames(AUC_6)
AUC_7$Radiation <- rownames(AUC_7)
AUC_8$Radiation <- rownames(AUC_8)
AUC_9$Radiation <- rownames(AUC_9)

###Pivot tables to long form and combine replicates for cell lines
AUC_1 <- melt(AUC_1, id= c("Radiation"))
AUC_2 <- melt(AUC_2, id= c("Radiation"))
AUC_3 <- melt(AUC_3, id= c("Radiation"))
AUC_4 <- melt(AUC_4, id= c("Radiation"))
AUC_5 <- melt(AUC_5, id= c("Radiation"))
AUC_6 <- melt(AUC_6, id= c("Radiation"))
AUC_7 <- melt(AUC_7, id= c("Radiation"))
AUC_8 <- melt(AUC_8, id= c("Radiation"))
AUC_9 <- melt(AUC_9, id= c("Radiation"))
JHAUC <- rbind(AUC_1, AUC_2, AUC_3)
Hec8AUC <- rbind(AUC_4, AUC_5, AUC_6)
HecBAUC <- rbind(AUC_7, AUC_8, AUC_9)

###Aggregate to obtain the mean signal for each dose combination
JHMean <- dcast(JHAUC, Radiation ~ variable, mean)
Hec8Mean <- dcast(Hec8AUC, Radiation ~ variable, mean)
HecBMean <- dcast(HecBAUC, Radiation ~ variable, mean)

###Normalize all values to the max value for Fraction Surviving
JHMean <- melt(JHMean, id= c("Radiation"))
JHMean$value <- JHMean$value/max(JHMean$value)
Hec8Mean <- melt(Hec8Mean, id= c("Radiation"))
Hec8Mean$value <- Hec8Mean$value/max(Hec8Mean$value)
HecBMean <- melt(HecBMean, id= c("Radiation"))
HecBMean$value <- HecBMean$value/max(HecBMean$value)

###Create new Cell column with labels and rename remaining columns
JHMean$Cell <- "JHUEM2"
names(JHMean)[1:3] <- c("d2", "d1", "effect")
Hec8Mean$Cell <- "Hec108"
names(Hec8Mean)[1:3] <- c("d2", "d1", "effect")
HecBMean$Cell <- "Hec1B"
names(HecBMean)[1:3] <- c("d2", "d1", "effect")

###Combine all dataframes together
AllMean <- rbind(JHMean, Hec8Mean, HecBMean)

###Write all dataframes to excel file
Sheets <- list("All"= AllMean, "JHUEM2"= JHMean, "Hec108" = Hec8Mean, "Hec1B"= HecBMean)
write.xlsx(Sheets, "Synergy JHUEM2, Hec108, Hec1B -Nutlin + Rad, AP47.xlsx")


##Synergy analysis with BIGL
###Read in excel files containing Nutlin + Radiation responses for JHUEM2, Hec108, Hec1B
NutlinSyn <- read_xlsx("Synergy JHUEM2, Hec108, Hec1B -Nutlin + Rad, AP46.xlsx", 1)

###Create separate Nutlin + Radiation tables for each cell line and rename columns for BIGL algorithms
JHNut <- subset(NutlinSyn, Cell= "JHUEM2", select = c("Cell", "d2", "d1", "effect"))
Hec8Nut <- subset(NutlinSyn, Cell= "Hec108", select = c("Cell", "d2", "d1", "effect"))
HecBNut <- subset(NutlinSyn, Cell= "Hec1B", select = c("Cell", "d2", "d1", "effect"))

###Convert Nutlin and Radiation doses to numeric
JHNut$d1 <- as.numeric(JHNutd1)
Hec8Nut$d1 <- as.numeric(Hec8Nut$d1)
HecBNut$d1 <- as.numeric(HecBNut$d1)
JHNut$d2 <- as.numeric(JHNut$d2)
Hec8Nut$d2 <- as.numeric(Hec8Nut$d2)
HecBNut$d2 <- as.numeric(HecBNut$d2)

###Fit marginal dose response curves to data using fitMarginal function in BIGL
Marg4Fit <- fitMarginals(JHNut, method= "nlslm", names= c("Nutlin", "Radiation"))
Marg5Fit <- fitMarginals(Hec8Nut, method= "nlslm", names= c("Nutlin", "Radiation"))
Marg6Fit <- fitMarginals(HecBNut, method= "nlslm", names= c("Nutlin", "Radiation"))

###Use fitSurface function in BIGL to fit predicted surface using Highest Single Agent (HSA) and Loewe models
JHNutSurf <- fitSurface(JHNut, Marg4Fit, null_model = "hsa", B.CP = 50, statistic = "both")
Hec8NutSurf <- fitSurface(Hec8Nut, Marg5Fit, null_model = "hsa", B.CP = 50, statistic = "both")
HecBNutSurf <- fitSurface(HecBNut, Marg6Fit, null_model = "loewe", B.CP = 50, statistic = "both")

###Generate Nutlin isobolograms of predicted surfaces using isobologram function in BIGL and save as jpegs
jpeg("JHUEM2 + Nutlin isobologram, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(JHNutSurf)
dev.off()

jpeg("Hec108 + Nutlin isobologram, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(Hec8NutSurf)
dev.off()

jpeg("Hec1B + Nutlin isobologram, Loewe.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(HecBNutSurf)
dev.off()

###Generate contour plots of Nutlin synergy scores using contour function in BIGL and save as jpegs
jpeg("JHUEM2 + Nutlin contour, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(JHNutSurf, main= "Contour Plot for JHUEM2 MaxR")
dev.off()

jpeg("Hec108 + Nutlin contour, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(Hec8NutSurf, main= "Contour Plot for Hec108 MaxR")
dev.off()

jpeg("Hec1B + Nutlin contour, Loewe.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(HecBNutSurf, main= "Contour Plot for Hec1B MaxR")
dev.off()

