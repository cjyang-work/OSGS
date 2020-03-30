library(emmeans)
setwd("C:/Users/cyang/Desktop/revision")

#read in the raw phenotype.
pheno <- read.csv("Barley_NAM_phenotype_raw.csv", as.is=T)

str(pheno)
#'data.frame':   11164 obs. of  8 variables:
# $ Line     : chr  "HEB_01_001" "HEB_01_003" "HEB_01_004" "HEB_01_005" ...
# $ Location : chr  "Dundee" "Dundee" "Dundee" "Dundee" ...
# $ Treatment: chr  "N0" "N0" "N0" "N0" ...
# $ Year     : int  2014 2014 2014 2014 2014 2014 2014 2014 2014 2014 ...
# $ DTH      : int  69 78 69 69 74 76 74 69 75 64 ...
# $ YLD      : num  NA NA NA NA NA NA NA NA NA NA ...
# $ GW       : num  3.8 3.7 3.5 3.5 3.7 3.5 3.7 3.6 3.6 3.4 ...
# $ GL       : num  10.4 10.3 8.6 8.6 9.1 10.1 9.2 9.7 9.7 10.3 ...

#calculate the least square means for DTH.
DTH.lm <- lm(DTH~Line+Location+Treatment+Year, data=pheno)
DTH.lsm <- data.frame(emmeans(DTH.lm, ~Line))[,1:2]
DTH.lsm$Line <- as.character(DTH.lsm$Line)

#calculate the least square means for YLD.
YLD.lm <- lm(YLD~Line+Treatment+Year, data=pheno)
YLD.lsm <- data.frame(emmeans(YLD.lm, ~Line))[,1:2]
YLD.lsm$Line <- as.character(YLD.lsm$Line)

#read in the line information
Line.info <- read.csv("Barley_NAM_line_info.csv", as.is=T)

str(Line.info)
#'data.frame':   1420 obs. of  3 variables:
# $ Line  : chr  "HEB_01_001" "HEB_01_003" "HEB_01_004" "HEB_01_005" ...
# $ Family: int  1 1 1 1 1 1 1 1 1 1 ...
# $ Donor : chr  "HID003" "HID003" "HID003" "HID003" ...

#combine all the least square means.
LSM <- merge(Line.info[,1:2], DTH.lsm, by="Line", all=T); colnames(LSM)[3] <- "DTH"
LSM <- merge(LSM, YLD.lsm, by="Line", all=T); colnames(LSM)[4] <- "YLD"

str(LSM)
#'data.frame':   1420 obs. of  4 variables:
# $ Line  : chr  "HEB_01_001" "HEB_01_003" "HEB_01_004" "HEB_01_005" ...
# $ Family: int  1 1 1 1 1 1 1 1 1 1 ...
# $ DTH   : num  68.8 77 69.4 76.6 74.3 ...
# $ YLD   : num  229 369 268 260 276 ...

#export the least square means.
write.csv(LSM, "Barley_NAM_phenotype.csv", row.names=F, quote=F)



















