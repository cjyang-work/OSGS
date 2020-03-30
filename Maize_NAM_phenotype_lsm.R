library(emmeans)
setwd("C:/Users/cyang/Desktop/revision")

#read in the raw phenotype.
pheno.DTS <- read.csv("Maize_NAM_phenotype_raw_DTS.csv", as.is=T)
pheno.CL <- read.csv("Maize_NAM_phenotype_raw_CL.csv", as.is=T)

str(pheno.DTS)
#'data.frame':   38968 obs. of  8 variables:
# $ Line    : chr  "Z001E0003" "Z001E0009" "Z001E0012" "Z001E0013" ...
# $ Location: chr  "IL" "IL" "IL" "IL" ...
# $ Rep     : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Year    : int  2006 2006 2006 2006 2006 2006 2006 2006 2006 2006 ...
# $ Family  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Block   : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Entry   : int  3 9 12 13 14 15 17 19 35 36 ...
# $ Trait   : int  77 81 81 75 78 75 79 79 75 79 ...

str(pheno.CL)
#'data.frame':   51204 obs. of  7 variables:
# $ Line    : chr  "Z001E0005" "Z001E0015" "Z001E0015" "Z001E0039" ...
# $ Location: chr  "FL" "FL" "FL" "FL" ...
# $ Year    : int  2006 2006 2006 2006 2006 2006 2006 2006 2006 2006 ...
# $ Family  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Block   : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Entry   : int  5 15 15 39 40 40 50 52 52 57 ...
# $ Trait   : num  88.6 106 126 52.3 94.5 ...

#calculate the least square means for DTS.
DTS.lm <- lm(Trait~Line+Location+Location:Rep+Year+Block, data=pheno.DTS)
DTS.lsm <- data.frame(emmeans(DTS.lm, ~Line))[,1:2]
DTS.lsm$Line <- as.character(DTS.lsm$Line)

#calculate the least square means for CL.
CL.lm <- lm(Trait~Line+Location+Year+Block, data=pheno.CL)
CL.lsm <- data.frame(emmeans(CL.lm, ~Line))[,1:2]
CL.lsm$Line <- as.character(CL.lsm$Line)

#read in the line information
Line.info <- read.csv("Maize_NAM_line_info.csv", as.is=T)

str(Line.info)
#'data.frame':   4699 obs. of  2 variables:
# $ Line  : chr  "Z001E0001" "Z001E0002" "Z001E0003" "Z001E0004" ...
# $ Family: int  1 1 1 1 1 1 1 1 1 1 ...

#combine all the least square means.
LSM <- merge(Line.info[,1:2], DTS.lsm, by="Line", all.x=T, all.y=F); colnames(LSM)[3] <- "DTS"
LSM <- merge(LSM, CL.lsm, by="Line", all.x=T, all.y=F); colnames(LSM)[4] <- "CL"

str(LSM)
#'data.frame':   4699 obs. of  4 variables:
# $ Line  : chr  "Z001E0001" "Z001E0002" "Z001E0003" "Z001E0004" ...
# $ Family: int  1 1 1 1 1 1 1 1 1 1 ...
# $ DTS   : num  75.7 75.7 73.1 74 83.1 ...
# $ CL    : num  120 108 118 114 140 ...

#export the least square means.
write.csv(LSM, "Maize_NAM_phenotype.csv", row.names=F, quote=F)



















