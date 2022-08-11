#Load required packages
library(readxl)
library(tidyverse)
library(BIGL)

#Set working directory
setwd("/Users/pettya/Desktop/R Analysis")


#Synergy analysis with BIGL (Figure 5, Supplemental Figure 5)
##Read in excel files containing AMG, Nutlin + Radiation responses for JHUEM2, Hec108, Hec1B
AMGSyn <- read_xlsx("Combined AMG + Radiation, cell lines, 7-21-22.xlsx", 1)
NutlinSyn <- read_xlsx("Combined Nutlin + Radiation, cell lines, 7-21-22.xlsx", 1)

##Create separate AMG + Radiation tables for each cell line and rename columns for BIGL algorithms
JHAMG <- subset(AMGSyn, block_id == 1, select = c("block_id", "conc1", "conc2", "response"))
Hec8AMG <- subset(AMGSyn, block_id == 2, select = c("block_id", "conc1", "conc2", "response"))
HecBAMG <- subset(AMGSyn, block_id == 3, select = c("block_id", "conc1", "conc2", "response"))
names(JHAMG)[2:4] <- c("d1", "d2", "effect")
names(Hec8AMG)[2:4] <- c("d1", "d2", "effect")
names(HecBAMG)[2:4] <- c("d1", "d2", "effect")

##Convert AMG dose to numeric
JHAMG$d1 <- as.numeric(JHAMG$d1)
Hec8AMG$d1 <- as.numeric(Hec8AMG$d1)
HecBAMG$d1 <- as.numeric(HecBAMG$d1)

##Create separate Nutlin + Radiation tables for each cell line and rename columns for BIGL algorithms
JHNut <- subset(NutlinSyn, block_id == 1, select = c("block_id", "conc1", "conc2", "response"))
Hec8Nut <- subset(NutlinSyn, block_id == 2, select = c("block_id", "conc1", "conc2", "response"))
HecBNut <- subset(NutlinSyn, block_id == 3, select = c("block_id", "conc1", "conc2", "response"))
names(JHNut)[2:4] <- c("d1", "d2", "effect")
names(Hec8Nut)[2:4] <- c("d1", "d2", "effect")
names(HecBNut)[2:4] <- c("d1", "d2", "effect")

#Convert Nutlin dose to numeric
JHNut$d1 <- as.numeric(JHNut$d1)
Hec8Nut$d1 <- as.numeric(Hec8Nut$d1)
HecBNut$d1 <- as.numeric(HecBNut$d1)

##Fit marginal dose response curves to data using fitMarginals function in BIGL
Marg1Fit <- fitMarginals(JHAMG, method= "nlslm", names= c("AMG-232", "Radiation"))
Marg2Fit <- fitMarginals(Hec8AMG, method= "nlslm", names= c("AMG-232", "Radiation"))
Marg3Fit <- fitMarginals(HecBAMG, method= "nlslm", names= c("AMG-232", "Radiation"))
Marg4Fit <- fitMarginals(JHNut, method= "nlslm", names= c("Nutlin", "Radiation"))
Marg5Fit <- fitMarginals(Hec8Nut, method= "nlslm", names= c("Nutlin", "Radiation"))
Marg6Fit <- fitMarginals(HecBNut, method= "nlslm", names= c("Nutlin", "Radiation"))

##Use fitSurface function in BIGL to fit predicted surface using Highest Single Agent (HSA) and Loewe models
JHAMGSurf <- fitSurface(JHAMG, Marg1Fit, null_model = "hsa", B.CP = 50, statistic = "both")
Hec8AMGSurf <- fitSurface(Hec8AMG, Marg2Fit, null_model = "hsa", B.CP = 50, statistic = "both")
HecBAMGSurf <- fitSurface(HecBAMG, Marg3Fit, null_model = "loewe", B.CP = 50, statistic = "both")
JHNutSurf <- fitSurface(JHNut, Marg4Fit, null_model = "hsa", B.CP = 50, statistic = "both")
Hec8NutSurf <- fitSurface(Hec8Nut, Marg5Fit, null_model = "hsa", B.CP = 50, statistic = "both")
HecBNutSurf <- fitSurface(HecBNut, Marg6Fit, null_model = "loewe", B.CP = 50, statistic = "both")

##Generate isobolograms of predicted surfaces using isobologram function in BIGL and save as jpegs
jpeg("JHUEM2 + AMG isobologram, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(JHAMGSurf)
dev.off()

jpeg("Hec108 + AMG isobologram, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(Hec8AMGSurf)
dev.off()

jpeg("Hec1B + AMG isobologram, Loewe.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(HecBAMGSurf)
dev.off()

jpeg("JHUEM2 + Nutlin isobologram, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(JHNutSurf)
dev.off()

jpeg("Hec108 + Nutlin isobologram, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(Hec8NutSurf)
dev.off()

jpeg("Hec1B + Nutlin isobologram, Loewe.jpeg", width= 1200, height= 720, units= "px", res= 150)
isobologram(HecBNutSurf)
dev.off()


##Generate contour plots of synergy scores using contour function in BIGL and save as jpegs
jpeg("JHUEM2 + AMG contour, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(JHAMGSurf, main= "Contour Plot for JHUEM2 MaxR")
dev.off()

jpeg("Hec108 + AMG contour, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(Hec8AMGSurf, main= "Contour Plot for Hec108 MaxR")
dev.off()

jpeg("Hec1B + AMG contour, Loewe.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(HecBAMGSurf, main= "Contour Plot for Hec1B MaxR")
dev.off()

jpeg("JHUEM2 + Nutlin contour, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(JHNutSurf, main= "Contour Plot for JHUEM2 MaxR")
dev.off()

jpeg("Hec108 + Nutlin contour, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(Hec8NutSurf, main= "Contour Plot for Hec108 MaxR")
dev.off()

jpeg("Hec1B + Nutlin contour, Loewe.jpeg", width= 1200, height= 720, units= "px", res= 150)
contour(HecBNutSurf, main= "Contour Plot for Hec1B MaxR")
dev.off()


##Generate table of confidence intervals using plotConfInt function in BIGL and save as jpegs
jpeg("JHUEM2 + AMG confidence intervals, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
plotConfInt(JHAMGSurf, color= "effect-size")
dev.off()

jpeg("Hec108 + AMG confidence intervals, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
plotConfInt(Hec8AMGSurf, color= "maxR")
dev.off()

jpeg("Hec1B + AMG confidence intervals, Loewe.jpeg", width= 1200, height= 720, units= "px", res= 150)
plotConfInt(HecBAMGSurf, color= "maxR")
dev.off()

jpeg("JHUEM2 + Nutlin confidence intervals, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
plotConfInt(JHNutSurf, color= "effect-size")
dev.off()

jpeg("Hec108 + Nutlin confidence intervals, HSA.jpeg", width= 1200, height= 720, units= "px", res= 150)
plotConfInt(Hec8NutSurf, color= "maxR")
dev.off()

jpeg("Hec1B + Nutlin confidence intervals, Loewe.jpeg", width= 1200, height= 720, units= "px", res= 150)
plotConfInt(HecBNutSurf, color= "maxR")
dev.off()


