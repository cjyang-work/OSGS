library(rrBLUP)
library(glmnet)
library(BGLR)
library(ggplot2)
library(reshape2)
library(AlphaSimR)
setwd("C:/Users/cyang/Desktop/OSGS")

#This file has 3 parts of scripts.
#1. Data preparation.
#2. Joint/independent prediction using rrBLUP/LASSO/BayesC.
#3. Simulating one generation of selection under OSGS vs GS.

################################
### Part 1. Data preparation ###
################################

#read in marker genotype.
geno.B <- read.csv("Barley_NAM_marker_genotype.csv", as.is=T, row.names=1)
geno.M <- read.csv("Maize_NAM_marker_genotype.csv", as.is=T, row.names=1)

#read in phenotype.
pheno.B <- read.csv("Barley_NAM_phenotype.csv", as.is=T)
pheno.M <- read.csv("Maize_NAM_phenotype.csv", as.is=T)

#remove lines with missing phenotype.
temp <- pheno.B[is.na(pheno.B$DTH) | is.na(pheno.B$YLD), 1]
geno.B <- geno.B[!(rownames(geno.B)%in%temp),]
pheno.B <- pheno.B[!(pheno.B[,1]%in%temp),]

temp <- pheno.M[is.na(pheno.M$DTS) | is.na(pheno.M$CL), 1]
geno.M <- geno.M[!(rownames(geno.M)%in%temp),]
pheno.M <- pheno.M[!(pheno.M[,1]%in%temp),]

#reorder lines in barley NAM phenotype to match with marker data.
pheno.B <- pheno.B[order(pheno.B$Family,pheno.B$Line),]
sum(pheno.B$Line==rownames(geno.B)) #1371

#read in marker information.
marker.B <- read.csv("Barley_NAM_marker_info.csv", as.is=T)
marker.M <- read.csv("Maize_NAM_marker_info.csv", as.is=T)

str(marker.B)
#'data.frame':	5704 obs. of  4 variables:
# $ NewMarkerName: chr  "M0001" "M0002" "M0003" "M0004" ...
# $ OldMarkerName: chr  "BOPA1_2511_533" "BOPA1_8670_388" "BOPA2_12_10420" "BOPA2_12_30653" ...
# $ Chr          : chr  "1H" "1H" "1H" "1H" ...
# $ PosGen       : num  0 0 0 0 0 0 0 0 0 0 ...

str(marker.M)
#'data.frame':   1106 obs. of  5 variables:
# $ NewMarkerName: chr  "M0001" "M0002" "M0003" "M0004" ...
# $ OldMarkerName: chr  "i0" "i1" "i2" "i3" ...
# $ Chr          : int  1 1 1 1 1 1 1 1 1 1 ...
# $ PosGen       : num  0 0.9 3.7 5.1 9.7 11.5 13.4 15.6 15.6 16.2 ...
# $ PanzeaMarker : chr  "PZA01271.1" "PZA03613.1" "PZA02129.1" "PZA02032.1" ...

#convert the marker genotype from data.frame to matrix.
geno.B <- as.matrix(geno.B)
geno.M <- as.matrix(geno.M)

#recode the marker genotype from 0/1/2 to -1/0/1
geno.B <- geno.B - 1
geno.M <- geno.M - 1

#split the marker genotype and phenotype by family.
family.B <- split(pheno.B[,1], pheno.B$Family)
family.M <- split(pheno.M[,1], pheno.M$Family)
geno.split.B <- lapply(1:length(family.B), FUN=function(x) geno.B[rownames(geno.B)%in%family.B[[x]],])
geno.split.M <- lapply(1:length(family.M), FUN=function(x) geno.M[rownames(geno.M)%in%family.M[[x]],])
pheno.split.B <- lapply(1:length(family.B), FUN=function(x) pheno.B[pheno.B[,1]%in%family.B[[x]],])
pheno.split.M <- lapply(1:length(family.M), FUN=function(x) pheno.M[pheno.M[,1]%in%family.M[[x]],])

#save all of the current data.
save.image("OSGS_dat1.RData")


##################################################################
### 2. Joint/independent prediction using rrBLUP/LASSO/BayesC. ###
##################################################################

###Joint analysis
#rrBLUP - Barley DTH
RR.DTH.1 <- mixed.solve(y=pheno.B$DTH,
                        Z=geno.B)
RR.mean.DTH.1 <- c(RR.DTH.1$beta)
RR.coef.DTH.1 <- c(RR.DTH.1$u)
RR.cfP.DTH.1 <- RR.coef.DTH.1
RR.cfS.DTH.1 <- RR.coef.DTH.1
RR.cfP.DTH.1[RR.cfP.DTH.1>0] <- 0 #high DTH is favorable in barley.
RR.cfS.DTH.1[RR.cfS.DTH.1<0] <- 0
RR.pred.DTH.1 <- data.frame(Observed=pheno.B$DTH-RR.mean.DTH.1,
                            Predicted.A=geno.B%*%RR.coef.DTH.1,
                            Predicted.P=geno.B%*%RR.cfP.DTH.1,
                            Predicted.S=geno.B%*%RR.cfS.DTH.1)

#rrBLUP - Barley YLD
RR.YLD.1 <- mixed.solve(y=pheno.B$YLD,
                        Z=geno.B)
RR.mean.YLD.1 <- c(RR.YLD.1$beta)
RR.coef.YLD.1 <- c(RR.YLD.1$u)
RR.cfP.YLD.1 <- RR.coef.YLD.1
RR.cfS.YLD.1 <- RR.coef.YLD.1
RR.cfP.YLD.1[RR.cfP.YLD.1>0] <- 0 #high YLD is favorable in barley.
RR.cfS.YLD.1[RR.cfS.YLD.1<0] <- 0
RR.pred.YLD.1 <- data.frame(Observed=pheno.B$YLD-RR.mean.YLD.1,
                            Predicted.A=geno.B%*%RR.coef.YLD.1,
                            Predicted.P=geno.B%*%RR.cfP.YLD.1,
                            Predicted.S=geno.B%*%RR.cfS.YLD.1)

#rrBLUP - Maize DTS
RR.DTS.1 <- mixed.solve(y=pheno.M$DTS,
                        Z=geno.M)
RR.mean.DTS.1 <- c(RR.DTS.1$beta)
RR.coef.DTS.1 <- c(RR.DTS.1$u)
RR.cfP.DTS.1 <- RR.coef.DTS.1
RR.cfS.DTS.1 <- RR.coef.DTS.1
RR.cfP.DTS.1[RR.cfP.DTS.1<0] <- 0 #low DTS is favorable in maize.
RR.cfS.DTS.1[RR.cfS.DTS.1>0] <- 0
RR.pred.DTS.1 <- data.frame(Observed=pheno.M$DTS-RR.mean.DTS.1,
                            Predicted.A=geno.M%*%RR.coef.DTS.1,
                            Predicted.P=geno.M%*%RR.cfP.DTS.1,
                            Predicted.S=geno.M%*%RR.cfS.DTS.1)

#rrBLUP - Maize CL
RR.CL.1 <- mixed.solve(y=pheno.M$CL,
                        Z=geno.M)
RR.mean.CL.1 <- c(RR.CL.1$beta)
RR.coef.CL.1 <- c(RR.CL.1$u)
RR.cfP.CL.1 <- RR.coef.CL.1
RR.cfS.CL.1 <- RR.coef.CL.1
RR.cfP.CL.1[RR.cfP.CL.1>0] <- 0 #high CL is favorable in maize.
RR.cfS.CL.1[RR.cfS.CL.1<0] <- 0
RR.pred.CL.1 <- data.frame(Observed=pheno.M$CL-RR.mean.CL.1,
                            Predicted.A=geno.M%*%RR.coef.CL.1,
                            Predicted.P=geno.M%*%RR.cfP.CL.1,
                            Predicted.S=geno.M%*%RR.cfS.CL.1)

#LASSO - Barley DTH
LASSO.DTH.1 <- cv.glmnet(x=geno.B,
                         y=pheno.B$DTH,
                         family="gaussian",
                         nfolds=5,
                         alpha=1,
                         standardize=F)
LASSO.mean.DTH.1 <- c(coef(LASSO.DTH.1, s="lambda.min")[1])
LASSO.coef.DTH.1 <- c(coef(LASSO.DTH.1, s="lambda.min")[-1])
LASSO.cfP.DTH.1 <- LASSO.coef.DTH.1
LASSO.cfS.DTH.1 <- LASSO.coef.DTH.1
LASSO.cfP.DTH.1[LASSO.cfP.DTH.1>0] <- 0 #high DTH is favorable in barley.
LASSO.cfS.DTH.1[LASSO.cfS.DTH.1<0] <- 0
LASSO.pred.DTH.1 <- data.frame(Observed=pheno.B$DTH-LASSO.mean.DTH.1,
                               Predicted.A=geno.B%*%LASSO.coef.DTH.1,
                               Predicted.P=geno.B%*%LASSO.cfP.DTH.1,
                               Predicted.S=geno.B%*%LASSO.cfS.DTH.1)

#LASSO - Barley YLD
LASSO.YLD.1 <- cv.glmnet(x=geno.B,
                         y=pheno.B$YLD,
                         family="gaussian",
                         nfolds=5,
                         alpha=1,
                         standardize=F)
LASSO.mean.YLD.1 <- c(coef(LASSO.YLD.1, s="lambda.min")[1])
LASSO.coef.YLD.1 <- c(coef(LASSO.YLD.1, s="lambda.min")[-1])
LASSO.cfP.YLD.1 <- LASSO.coef.YLD.1
LASSO.cfS.YLD.1 <- LASSO.coef.YLD.1
LASSO.cfP.YLD.1[LASSO.cfP.YLD.1>0] <- 0 #high YLD is favorable in barley.
LASSO.cfS.YLD.1[LASSO.cfS.YLD.1<0] <- 0
LASSO.pred.YLD.1 <- data.frame(Observed=pheno.B$YLD-LASSO.mean.YLD.1,
                               Predicted.A=geno.B%*%LASSO.coef.YLD.1,
                               Predicted.P=geno.B%*%LASSO.cfP.YLD.1,
                               Predicted.S=geno.B%*%LASSO.cfS.YLD.1)

#LASSO - Maize DTS
LASSO.DTS.1 <- cv.glmnet(x=geno.M,
                         y=pheno.M$DTS,
                         family="gaussian",
                         nfolds=5,
                         alpha=1,
                         standardize=F)
LASSO.mean.DTS.1 <- c(coef(LASSO.DTS.1, s="lambda.min")[1])
LASSO.coef.DTS.1 <- c(coef(LASSO.DTS.1, s="lambda.min")[-1])
LASSO.cfP.DTS.1 <- LASSO.coef.DTS.1
LASSO.cfS.DTS.1 <- LASSO.coef.DTS.1
LASSO.cfP.DTS.1[LASSO.cfP.DTS.1<0] <- 0 #low DTS is favorable in maize.
LASSO.cfS.DTS.1[LASSO.cfS.DTS.1>0] <- 0
LASSO.pred.DTS.1 <- data.frame(Observed=pheno.M$DTS-LASSO.mean.DTS.1,
                               Predicted.A=geno.M%*%LASSO.coef.DTS.1,
                               Predicted.P=geno.M%*%LASSO.cfP.DTS.1,
                               Predicted.S=geno.M%*%LASSO.cfS.DTS.1)

#LASSO - Maize CL
LASSO.CL.1 <- cv.glmnet(x=geno.M,
                         y=pheno.M$CL,
                         family="gaussian",
                         nfolds=5,
                         alpha=1,
                         standardize=F)
LASSO.mean.CL.1 <- c(coef(LASSO.CL.1, s="lambda.min")[1])
LASSO.coef.CL.1 <- c(coef(LASSO.CL.1, s="lambda.min")[-1])
LASSO.cfP.CL.1 <- LASSO.coef.CL.1
LASSO.cfS.CL.1 <- LASSO.coef.CL.1
LASSO.cfP.CL.1[LASSO.cfP.CL.1>0] <- 0 #high CL is favorable in maize.
LASSO.cfS.CL.1[LASSO.cfS.CL.1<0] <- 0
LASSO.pred.CL.1 <- data.frame(Observed=pheno.M$CL-LASSO.mean.CL.1,
                               Predicted.A=geno.M%*%LASSO.coef.CL.1,
                               Predicted.P=geno.M%*%LASSO.cfP.CL.1,
                               Predicted.S=geno.M%*%LASSO.cfS.CL.1)

#Phenotype is mean-centered and marker genotype is kept the same (no scaling).
#BayesCpi - Barley DTH
BayesCpi.DTH.1 <- BGLR(y=pheno.B$DTH-mean(pheno.B$DTH),
                       ETA=list(list(X=geno.B, model='BayesC', saveEffects=T)),
                       nIter=60000,
                       burnIn=10000,
                       thin=50,
                       verbose=F,
                       saveAt="output/DTH_1_")
BayesCpi.coef.DTH.1 <- readBinMat("output/DTH_1_ETA_1_b.bin", byrow=T)
BayesCpi.cfP.DTH.1 <- BayesCpi.coef.DTH.1
BayesCpi.cfS.DTH.1 <- BayesCpi.coef.DTH.1
BayesCpi.cfP.DTH.1[,BayesCpi.DTH.1$ETA[[1]]$b>0] <- 0 #high DTH is favorable in barley.
BayesCpi.cfS.DTH.1[,BayesCpi.DTH.1$ETA[[1]]$b<0] <- 0
BayesCpi.pred.DTH.1 <- data.frame(Observed=pheno.B$DTH-mean(pheno.B$DTH)-BayesCpi.DTH.1$mu,
                                  Predicted.A=rowSums(sapply(1:nrow(BayesCpi.coef.DTH.1), FUN=function(x) geno.B%*%BayesCpi.coef.DTH.1[x,]))/nrow(BayesCpi.coef.DTH.1),
                                  Predicted.P=rowSums(sapply(1:nrow(BayesCpi.cfP.DTH.1), FUN=function(x) geno.B%*%BayesCpi.cfP.DTH.1[x,]))/nrow(BayesCpi.cfP.DTH.1),
                                  Predicted.S=rowSums(sapply(1:nrow(BayesCpi.cfS.DTH.1), FUN=function(x) geno.B%*%BayesCpi.cfS.DTH.1[x,]))/nrow(BayesCpi.cfS.DTH.1))

#BayesCpi - Barley YLD
BayesCpi.YLD.1 <- BGLR(y=pheno.B$YLD-mean(pheno.B$YLD),
                       ETA=list(list(X=geno.B, model='BayesC', saveEffects=T)),
                       nIter=60000,
                       burnIn=10000,
                       thin=50,
                       verbose=F,
                       saveAt="output/YLD_1_")
BayesCpi.coef.YLD.1 <- readBinMat("output/YLD_1_ETA_1_b.bin", byrow=T)
BayesCpi.cfP.YLD.1 <- BayesCpi.coef.YLD.1
BayesCpi.cfS.YLD.1 <- BayesCpi.coef.YLD.1
BayesCpi.cfP.YLD.1[,BayesCpi.YLD.1$ETA[[1]]$b>0] <- 0 #high YLD is favorable in barley.
BayesCpi.cfS.YLD.1[,BayesCpi.YLD.1$ETA[[1]]$b<0] <- 0
BayesCpi.pred.YLD.1 <- data.frame(Observed=pheno.B$YLD-mean(pheno.B$YLD)-BayesCpi.YLD.1$mu,
                                  Predicted.A=rowSums(sapply(1:nrow(BayesCpi.coef.YLD.1), FUN=function(x) geno.B%*%BayesCpi.coef.YLD.1[x,]))/nrow(BayesCpi.coef.YLD.1),
                                  Predicted.P=rowSums(sapply(1:nrow(BayesCpi.cfP.YLD.1), FUN=function(x) geno.B%*%BayesCpi.cfP.YLD.1[x,]))/nrow(BayesCpi.cfP.YLD.1),
                                  Predicted.S=rowSums(sapply(1:nrow(BayesCpi.cfS.YLD.1), FUN=function(x) geno.B%*%BayesCpi.cfS.YLD.1[x,]))/nrow(BayesCpi.cfS.YLD.1))

#BayesCpi - Maize DTS
BayesCpi.DTS.1 <- BGLR(y=pheno.M$DTS-mean(pheno.M$DTS),
                       ETA=list(list(X=geno.M, model='BayesC', saveEffects=T)),
                       nIter=60000,
                       burnIn=10000,
                       thin=50,
                       verbose=F,
                       saveAt="output/DTS_1_")
BayesCpi.coef.DTS.1 <- readBinMat("output/DTS_1_ETA_1_b.bin", byrow=T)
BayesCpi.cfP.DTS.1 <- BayesCpi.coef.DTS.1
BayesCpi.cfS.DTS.1 <- BayesCpi.coef.DTS.1
BayesCpi.cfP.DTS.1[,BayesCpi.DTS.1$ETA[[1]]$b<0] <- 0 #low DTS is favorable in maize.
BayesCpi.cfS.DTS.1[,BayesCpi.DTS.1$ETA[[1]]$b>0] <- 0
BayesCpi.pred.DTS.1 <- data.frame(Observed=pheno.M$DTS-mean(pheno.M$DTS)-BayesCpi.DTS.1$mu,
                                  Predicted.A=rowSums(sapply(1:nrow(BayesCpi.coef.DTS.1), FUN=function(x) geno.M%*%BayesCpi.coef.DTS.1[x,]))/nrow(BayesCpi.coef.DTS.1),
                                  Predicted.P=rowSums(sapply(1:nrow(BayesCpi.cfP.DTS.1), FUN=function(x) geno.M%*%BayesCpi.cfP.DTS.1[x,]))/nrow(BayesCpi.cfP.DTS.1),
                                  Predicted.S=rowSums(sapply(1:nrow(BayesCpi.cfS.DTS.1), FUN=function(x) geno.M%*%BayesCpi.cfS.DTS.1[x,]))/nrow(BayesCpi.cfS.DTS.1))

#BayesCpi - Maize CL
BayesCpi.CL.1 <- BGLR(y=pheno.M$CL-mean(pheno.M$CL),
                       ETA=list(list(X=geno.M, model='BayesC', saveEffects=T)),
                       nIter=60000,
                       burnIn=10000,
                       thin=50,
                       verbose=F,
                       saveAt="output/CL_1_")
BayesCpi.coef.CL.1 <- readBinMat("output/CL_1_ETA_1_b.bin", byrow=T)
BayesCpi.cfP.CL.1 <- BayesCpi.coef.CL.1
BayesCpi.cfS.CL.1 <- BayesCpi.coef.CL.1
BayesCpi.cfP.CL.1[,BayesCpi.CL.1$ETA[[1]]$b>0] <- 0 #high CL is favorable in maize.
BayesCpi.cfS.CL.1[,BayesCpi.CL.1$ETA[[1]]$b<0] <- 0
BayesCpi.pred.CL.1 <- data.frame(Observed=pheno.M$CL-mean(pheno.M$CL)-BayesCpi.CL.1$mu,
                                  Predicted.A=rowSums(sapply(1:nrow(BayesCpi.coef.CL.1), FUN=function(x) geno.M%*%BayesCpi.coef.CL.1[x,]))/nrow(BayesCpi.coef.CL.1),
                                  Predicted.P=rowSums(sapply(1:nrow(BayesCpi.cfP.CL.1), FUN=function(x) geno.M%*%BayesCpi.cfP.CL.1[x,]))/nrow(BayesCpi.cfP.CL.1),
                                  Predicted.S=rowSums(sapply(1:nrow(BayesCpi.cfS.CL.1), FUN=function(x) geno.M%*%BayesCpi.cfS.CL.1[x,]))/nrow(BayesCpi.cfS.CL.1))


###Independent analysis
#rrBLUP - Barley DTH
RR.DTH.4 <- list()
RR.mean.DTH.4 <- list()
RR.coef.DTH.4 <- list()
RR.cfP.DTH.4 <- list()
RR.cfS.DTH.4 <- list()
RR.pred.DTH.4 <- list()

for(i in 1:25){
  RR.DTH.4[[i]] <- mixed.solve(y=pheno.split.B[[i]]$DTH,
                               Z=geno.split.B[[i]])
  RR.mean.DTH.4[[i]] <- c(RR.DTH.4[[i]]$beta)
  RR.coef.DTH.4[[i]] <- c(RR.DTH.4[[i]]$u)
  RR.cfP.DTH.4[[i]] <- RR.coef.DTH.4[[i]]
  RR.cfS.DTH.4[[i]] <- RR.coef.DTH.4[[i]]
  RR.cfP.DTH.4[[i]][RR.cfP.DTH.4[[i]]>0] <- 0 #high DTH is favorable in barley.
  RR.cfS.DTH.4[[i]][RR.cfS.DTH.4[[i]]<0] <- 0
  RR.pred.DTH.4[[i]] <- data.frame(Observed=pheno.split.B[[i]]$DTH-RR.mean.DTH.4[[i]],
                                   Predicted.A=geno.split.B[[i]]%*%RR.coef.DTH.4[[i]],
                                   Predicted.P=geno.split.B[[i]]%*%RR.cfP.DTH.4[[i]],
                                   Predicted.S=geno.split.B[[i]]%*%RR.cfS.DTH.4[[i]])
  cat("..",i)
}

#rrBLUP - Barley YLD
RR.YLD.4 <- list()
RR.mean.YLD.4 <- list()
RR.coef.YLD.4 <- list()
RR.cfP.YLD.4 <- list()
RR.cfS.YLD.4 <- list()
RR.pred.YLD.4 <- list()

for(i in 1:25){
  RR.YLD.4[[i]] <- mixed.solve(y=pheno.split.B[[i]]$YLD,
                               Z=geno.split.B[[i]])
  RR.mean.YLD.4[[i]] <- c(RR.YLD.4[[i]]$beta)
  RR.coef.YLD.4[[i]] <- c(RR.YLD.4[[i]]$u)
  RR.cfP.YLD.4[[i]] <- RR.coef.YLD.4[[i]]
  RR.cfS.YLD.4[[i]] <- RR.coef.YLD.4[[i]]
  RR.cfP.YLD.4[[i]][RR.cfP.YLD.4[[i]]>0] <- 0 #high YLD is favorable in barley.
  RR.cfS.YLD.4[[i]][RR.cfS.YLD.4[[i]]<0] <- 0
  RR.pred.YLD.4[[i]] <- data.frame(Observed=pheno.split.B[[i]]$YLD-RR.mean.YLD.4[[i]],
                                   Predicted.A=geno.split.B[[i]]%*%RR.coef.YLD.4[[i]],
                                   Predicted.P=geno.split.B[[i]]%*%RR.cfP.YLD.4[[i]],
                                   Predicted.S=geno.split.B[[i]]%*%RR.cfS.YLD.4[[i]])
  cat("..",i)
}

#rrBLUP - Maize DTS
RR.DTS.4 <- list()
RR.mean.DTS.4 <- list()
RR.coef.DTS.4 <- list()
RR.cfP.DTS.4 <- list()
RR.cfS.DTS.4 <- list()
RR.pred.DTS.4 <- list()

for(i in 1:25){
  RR.DTS.4[[i]] <- mixed.solve(y=pheno.split.M[[i]]$DTS,
                               Z=geno.split.M[[i]])
  RR.mean.DTS.4[[i]] <- c(RR.DTS.4[[i]]$beta)
  RR.coef.DTS.4[[i]] <- c(RR.DTS.4[[i]]$u)
  RR.cfP.DTS.4[[i]] <- RR.coef.DTS.4[[i]]
  RR.cfS.DTS.4[[i]] <- RR.coef.DTS.4[[i]]
  RR.cfP.DTS.4[[i]][RR.cfP.DTS.4[[i]]<0] <- 0 #low DTS is favorable in maize.
  RR.cfS.DTS.4[[i]][RR.cfS.DTS.4[[i]]>0] <- 0
  RR.pred.DTS.4[[i]] <- data.frame(Observed=pheno.split.M[[i]]$DTS-RR.mean.DTS.4[[i]],
                                   Predicted.A=geno.split.M[[i]]%*%RR.coef.DTS.4[[i]],
                                   Predicted.P=geno.split.M[[i]]%*%RR.cfP.DTS.4[[i]],
                                   Predicted.S=geno.split.M[[i]]%*%RR.cfS.DTS.4[[i]])
  cat("..",i)
}

#rrBLUP - Maize CL
RR.CL.4 <- list()
RR.mean.CL.4 <- list()
RR.coef.CL.4 <- list()
RR.cfP.CL.4 <- list()
RR.cfS.CL.4 <- list()
RR.pred.CL.4 <- list()

for(i in 1:25){
  RR.CL.4[[i]] <- mixed.solve(y=pheno.split.M[[i]]$CL,
                              Z=geno.split.M[[i]])
  RR.mean.CL.4[[i]] <- c(RR.CL.4[[i]]$beta)
  RR.coef.CL.4[[i]] <- c(RR.CL.4[[i]]$u)
  RR.cfP.CL.4[[i]] <- RR.coef.CL.4[[i]]
  RR.cfS.CL.4[[i]] <- RR.coef.CL.4[[i]]
  RR.cfP.CL.4[[i]][RR.cfP.CL.4[[i]]>0] <- 0 #high CL is favorable in maize.
  RR.cfS.CL.4[[i]][RR.cfS.CL.4[[i]]<0] <- 0
  RR.pred.CL.4[[i]] <- data.frame(Observed=pheno.split.M[[i]]$CL-RR.mean.CL.4[[i]],
                                  Predicted.A=geno.split.M[[i]]%*%RR.coef.CL.4[[i]],
                                  Predicted.P=geno.split.M[[i]]%*%RR.cfP.CL.4[[i]],
                                  Predicted.S=geno.split.M[[i]]%*%RR.cfS.CL.4[[i]])
  cat("..",i)
}

#LASSO - Barley DTH
LASSO.DTH.4 <- list()
LASSO.mean.DTH.4 <- list()
LASSO.coef.DTH.4 <- list()
LASSO.cfP.DTH.4 <- list()
LASSO.cfS.DTH.4 <- list()
LASSO.pred.DTH.4 <- list()

for(i in 1:25){
  LASSO.DTH.4[[i]] <- cv.glmnet(x=geno.split.B[[i]],
                                y=pheno.split.B[[i]]$DTH,
                                family="gaussian",
                                nfolds=5,
                                alpha=1,
                                standardize=F)
  LASSO.mean.DTH.4[[i]] <- c(coef(LASSO.DTH.4[[i]], s="lambda.min")[1])
  LASSO.coef.DTH.4[[i]] <- c(coef(LASSO.DTH.4[[i]], s="lambda.min")[-1])
  LASSO.cfP.DTH.4[[i]] <- LASSO.coef.DTH.4[[i]]
  LASSO.cfS.DTH.4[[i]] <- LASSO.coef.DTH.4[[i]]
  LASSO.cfP.DTH.4[[i]][LASSO.cfP.DTH.4[[i]]>0] <- 0 #high DTH is favorable in barley.
  LASSO.cfS.DTH.4[[i]][LASSO.cfS.DTH.4[[i]]<0] <- 0
  LASSO.pred.DTH.4[[i]] <- data.frame(Observed=pheno.split.B[[i]]$DTH-LASSO.mean.DTH.4[[i]],
                                      Predicted.A=geno.split.B[[i]]%*%LASSO.coef.DTH.4[[i]],
                                      Predicted.P=geno.split.B[[i]]%*%LASSO.cfP.DTH.4[[i]],
                                      Predicted.S=geno.split.B[[i]]%*%LASSO.cfS.DTH.4[[i]])
  cat("..",i)
}

#LASSO - Barley YLD
LASSO.YLD.4 <- list()
LASSO.mean.YLD.4 <- list()
LASSO.coef.YLD.4 <- list()
LASSO.cfP.YLD.4 <- list()
LASSO.cfS.YLD.4 <- list()
LASSO.pred.YLD.4 <- list()

for(i in 1:25){
  LASSO.YLD.4[[i]] <- cv.glmnet(x=geno.split.B[[i]],
                                y=pheno.split.B[[i]]$YLD,
                                family="gaussian",
                                nfolds=5,
                                alpha=1,
                                standardize=F)
  LASSO.mean.YLD.4[[i]] <- c(coef(LASSO.YLD.4[[i]], s="lambda.min")[1])
  LASSO.coef.YLD.4[[i]] <- c(coef(LASSO.YLD.4[[i]], s="lambda.min")[-1])
  LASSO.cfP.YLD.4[[i]] <- LASSO.coef.YLD.4[[i]]
  LASSO.cfS.YLD.4[[i]] <- LASSO.coef.YLD.4[[i]]
  LASSO.cfP.YLD.4[[i]][LASSO.cfP.YLD.4[[i]]>0] <- 0 #high YLD is favorable in barley.
  LASSO.cfS.YLD.4[[i]][LASSO.cfS.YLD.4[[i]]<0] <- 0
  LASSO.pred.YLD.4[[i]] <- data.frame(Observed=pheno.split.B[[i]]$YLD-LASSO.mean.YLD.4[[i]],
                                      Predicted.A=geno.split.B[[i]]%*%LASSO.coef.YLD.4[[i]],
                                      Predicted.P=geno.split.B[[i]]%*%LASSO.cfP.YLD.4[[i]],
                                      Predicted.S=geno.split.B[[i]]%*%LASSO.cfS.YLD.4[[i]])
  cat("..",i)
}

#LASSO - Maize DTS
LASSO.DTS.4 <- list()
LASSO.mean.DTS.4 <- list()
LASSO.coef.DTS.4 <- list()
LASSO.cfP.DTS.4 <- list()
LASSO.cfS.DTS.4 <- list()
LASSO.pred.DTS.4 <- list()

for(i in 1:25){
  LASSO.DTS.4[[i]] <- cv.glmnet(x=geno.split.M[[i]],
                                y=pheno.split.M[[i]]$DTS,
                                family="gaussian",
                                nfolds=5,
                                alpha=1,
                                standardize=F)
  LASSO.mean.DTS.4[[i]] <- c(coef(LASSO.DTS.4[[i]], s="lambda.min")[1])
  LASSO.coef.DTS.4[[i]] <- c(coef(LASSO.DTS.4[[i]], s="lambda.min")[-1])
  LASSO.cfP.DTS.4[[i]] <- LASSO.coef.DTS.4[[i]]
  LASSO.cfS.DTS.4[[i]] <- LASSO.coef.DTS.4[[i]]
  LASSO.cfP.DTS.4[[i]][LASSO.cfP.DTS.4[[i]]<0] <- 0 #low DTS is favorable in maize.
  LASSO.cfS.DTS.4[[i]][LASSO.cfS.DTS.4[[i]]>0] <- 0
  LASSO.pred.DTS.4[[i]] <- data.frame(Observed=pheno.split.M[[i]]$DTS-LASSO.mean.DTS.4[[i]],
                                      Predicted.A=geno.split.M[[i]]%*%LASSO.coef.DTS.4[[i]],
                                      Predicted.P=geno.split.M[[i]]%*%LASSO.cfP.DTS.4[[i]],
                                      Predicted.S=geno.split.M[[i]]%*%LASSO.cfS.DTS.4[[i]])
  cat("..",i)
}

#LASSO - Maize CL
LASSO.CL.4 <- list()
LASSO.mean.CL.4 <- list()
LASSO.coef.CL.4 <- list()
LASSO.cfP.CL.4 <- list()
LASSO.cfS.CL.4 <- list()
LASSO.pred.CL.4 <- list()

for(i in 1:25){
  LASSO.CL.4[[i]] <- cv.glmnet(x=geno.split.M[[i]],
                               y=pheno.split.M[[i]]$CL,
                               family="gaussian",
                               nfolds=5,
                               alpha=1,
                               standardize=F)
  LASSO.mean.CL.4[[i]] <- c(coef(LASSO.CL.4[[i]], s="lambda.min")[1])
  LASSO.coef.CL.4[[i]] <- c(coef(LASSO.CL.4[[i]], s="lambda.min")[-1])
  LASSO.cfP.CL.4[[i]] <- LASSO.coef.CL.4[[i]]
  LASSO.cfS.CL.4[[i]] <- LASSO.coef.CL.4[[i]]
  LASSO.cfP.CL.4[[i]][LASSO.cfP.CL.4[[i]]>0] <- 0 #high CL is favorable in maize.
  LASSO.cfS.CL.4[[i]][LASSO.cfS.CL.4[[i]]<0] <- 0
  LASSO.pred.CL.4[[i]] <- data.frame(Observed=pheno.split.M[[i]]$CL-LASSO.mean.CL.4[[i]],
                                     Predicted.A=geno.split.M[[i]]%*%LASSO.coef.CL.4[[i]],
                                     Predicted.P=geno.split.M[[i]]%*%LASSO.cfP.CL.4[[i]],
                                     Predicted.S=geno.split.M[[i]]%*%LASSO.cfS.CL.4[[i]])
  cat("..",i)
}

#BayesCpi - Barley DTH
BayesCpi.DTH.4 <- list()
BayesCpi.coef.DTH.4 <- list()
BayesCpi.cfP.DTH.4 <- list()
BayesCpi.cfS.DTH.4 <- list()
BayesCpi.pred.DTH.4 <- list()

for(i in 1:25){
  BayesCpi.DTH.4[[i]] <- BGLR(y=pheno.split.B[[i]]$DTH-mean(pheno.split.B[[i]]$DTH),
                              ETA=list(list(X=geno.split.B[[i]], model='BayesC', saveEffects=T)),
                              nIter=60000,
                              burnIn=10000,
                              thin=50,
                              verbose=F,
                              saveAt=paste("output/DTH_4_",i,"_", sep=""))
  BayesCpi.coef.DTH.4[[i]] <- readBinMat(paste("output/DTH_4_",i,"_ETA_1_b.bin",sep=""), byrow=T)
  BayesCpi.cfP.DTH.4[[i]] <- BayesCpi.coef.DTH.4[[i]]
  BayesCpi.cfS.DTH.4[[i]] <- BayesCpi.coef.DTH.4[[i]]
  BayesCpi.cfP.DTH.4[[i]][,BayesCpi.DTH.4[[i]]$ETA[[1]]$b>0] <- 0 #high DTH is favorable in barley.
  BayesCpi.cfS.DTH.4[[i]][,BayesCpi.DTH.4[[i]]$ETA[[1]]$b<0] <- 0
  BayesCpi.pred.DTH.4[[i]] <- data.frame(Observed=pheno.split.B[[i]]$DTH-mean(pheno.split.B[[i]]$DTH)-BayesCpi.DTH.4[[i]]$mu,
                                         Predicted.A=rowSums(sapply(1:nrow(BayesCpi.coef.DTH.4[[i]]), FUN=function(x) geno.split.B[[i]]%*%BayesCpi.coef.DTH.4[[i]][x,]))/nrow(BayesCpi.coef.DTH.4[[i]]),
                                         Predicted.P=rowSums(sapply(1:nrow(BayesCpi.cfP.DTH.4[[i]]), FUN=function(x) geno.split.B[[i]]%*%BayesCpi.cfP.DTH.4[[i]][x,]))/nrow(BayesCpi.cfP.DTH.4[[i]]),
                                         Predicted.S=rowSums(sapply(1:nrow(BayesCpi.cfS.DTH.4[[i]]), FUN=function(x) geno.split.B[[i]]%*%BayesCpi.cfS.DTH.4[[i]][x,]))/nrow(BayesCpi.cfS.DTH.4[[i]]))
  cat("..",i)
}

#BayesCpi - Barley YLD
BayesCpi.YLD.4 <- list()
BayesCpi.coef.YLD.4 <- list()
BayesCpi.cfP.YLD.4 <- list()
BayesCpi.cfS.YLD.4 <- list()
BayesCpi.pred.YLD.4 <- list()

for(i in 1:25){
  BayesCpi.YLD.4[[i]] <- BGLR(y=pheno.split.B[[i]]$YLD-mean(pheno.split.B[[i]]$YLD),
                              ETA=list(list(X=geno.split.B[[i]], model='BayesC', saveEffects=T)),
                              nIter=60000,
                              burnIn=10000,
                              thin=50,
                              verbose=F,
                              saveAt=paste("output/YLD_4_",i,"_", sep=""))
  BayesCpi.coef.YLD.4[[i]] <- readBinMat(paste("output/YLD_4_",i,"_ETA_1_b.bin",sep=""), byrow=T)
  BayesCpi.cfP.YLD.4[[i]] <- BayesCpi.coef.YLD.4[[i]]
  BayesCpi.cfS.YLD.4[[i]] <- BayesCpi.coef.YLD.4[[i]]
  BayesCpi.cfP.YLD.4[[i]][,BayesCpi.YLD.4[[i]]$ETA[[1]]$b>0] <- 0 #high YLD is favorable in barley.
  BayesCpi.cfS.YLD.4[[i]][,BayesCpi.YLD.4[[i]]$ETA[[1]]$b<0] <- 0
  BayesCpi.pred.YLD.4[[i]] <- data.frame(Observed=pheno.split.B[[i]]$YLD-mean(pheno.split.B[[i]]$YLD)-BayesCpi.YLD.4[[i]]$mu,
                                         Predicted.A=rowSums(sapply(1:nrow(BayesCpi.coef.YLD.4[[i]]), FUN=function(x) geno.split.B[[i]]%*%BayesCpi.coef.YLD.4[[i]][x,]))/nrow(BayesCpi.coef.YLD.4[[i]]),
                                         Predicted.P=rowSums(sapply(1:nrow(BayesCpi.cfP.YLD.4[[i]]), FUN=function(x) geno.split.B[[i]]%*%BayesCpi.cfP.YLD.4[[i]][x,]))/nrow(BayesCpi.cfP.YLD.4[[i]]),
                                         Predicted.S=rowSums(sapply(1:nrow(BayesCpi.cfS.YLD.4[[i]]), FUN=function(x) geno.split.B[[i]]%*%BayesCpi.cfS.YLD.4[[i]][x,]))/nrow(BayesCpi.cfS.YLD.4[[i]]))
  cat("..",i)
}

#BayesCpi - Maize DTS
BayesCpi.DTS.4 <- list()
BayesCpi.coef.DTS.4 <- list()
BayesCpi.cfP.DTS.4 <- list()
BayesCpi.cfS.DTS.4 <- list()
BayesCpi.pred.DTS.4 <- list()

for(i in 1:25){
  BayesCpi.DTS.4[[i]] <- BGLR(y=pheno.split.M[[i]]$DTS-mean(pheno.split.M[[i]]$DTS),
                              ETA=list(list(X=geno.split.M[[i]], model='BayesC', saveEffects=T)),
                              nIter=60000,
                              burnIn=10000,
                              thin=50,
                              verbose=F,
                              saveAt=paste("output/DTS_4_",i,"_", sep=""))
  BayesCpi.coef.DTS.4[[i]] <- readBinMat(paste("output/DTS_4_",i,"_ETA_1_b.bin",sep=""), byrow=T)
  BayesCpi.cfP.DTS.4[[i]] <- BayesCpi.coef.DTS.4[[i]]
  BayesCpi.cfS.DTS.4[[i]] <- BayesCpi.coef.DTS.4[[i]]
  BayesCpi.cfP.DTS.4[[i]][,BayesCpi.DTS.4[[i]]$ETA[[1]]$b<0] <- 0 #low DTS is favorable in maize.
  BayesCpi.cfS.DTS.4[[i]][,BayesCpi.DTS.4[[i]]$ETA[[1]]$b>0] <- 0
  BayesCpi.pred.DTS.4[[i]] <- data.frame(Observed=pheno.split.M[[i]]$DTS-mean(pheno.split.M[[i]]$DTS)-BayesCpi.DTS.4[[i]]$mu,
                                         Predicted.A=rowSums(sapply(1:nrow(BayesCpi.coef.DTS.4[[i]]), FUN=function(x) geno.split.M[[i]]%*%BayesCpi.coef.DTS.4[[i]][x,]))/nrow(BayesCpi.coef.DTS.4[[i]]),
                                         Predicted.P=rowSums(sapply(1:nrow(BayesCpi.cfP.DTS.4[[i]]), FUN=function(x) geno.split.M[[i]]%*%BayesCpi.cfP.DTS.4[[i]][x,]))/nrow(BayesCpi.cfP.DTS.4[[i]]),
                                         Predicted.S=rowSums(sapply(1:nrow(BayesCpi.cfS.DTS.4[[i]]), FUN=function(x) geno.split.M[[i]]%*%BayesCpi.cfS.DTS.4[[i]][x,]))/nrow(BayesCpi.cfS.DTS.4[[i]]))
  cat("..",i)
}

#BayesCpi - Maize CL
BayesCpi.CL.4 <- list()
BayesCpi.coef.CL.4 <- list()
BayesCpi.cfP.CL.4 <- list()
BayesCpi.cfS.CL.4 <- list()
BayesCpi.pred.CL.4 <- list()

for(i in 1:25){
  BayesCpi.CL.4[[i]] <- BGLR(y=pheno.split.M[[i]]$CL-mean(pheno.split.M[[i]]$CL),
                             ETA=list(list(X=geno.split.M[[i]], model='BayesC', saveEffects=T)),
                             nIter=60000,
                             burnIn=10000,
                             thin=50,
                             verbose=F,
                             saveAt=paste("output/CL_4_",i,"_", sep=""))
  BayesCpi.coef.CL.4[[i]] <- readBinMat(paste("output/CL_4_",i,"_ETA_1_b.bin",sep=""), byrow=T)
  BayesCpi.cfP.CL.4[[i]] <- BayesCpi.coef.CL.4[[i]]
  BayesCpi.cfS.CL.4[[i]] <- BayesCpi.coef.CL.4[[i]]
  BayesCpi.cfP.CL.4[[i]][,BayesCpi.CL.4[[i]]$ETA[[1]]$b>0] <- 0 #high CL is favorable in maize.
  BayesCpi.cfS.CL.4[[i]][,BayesCpi.CL.4[[i]]$ETA[[1]]$b<0] <- 0
  BayesCpi.pred.CL.4[[i]] <- data.frame(Observed=pheno.split.M[[i]]$CL-mean(pheno.split.M[[i]]$CL)-BayesCpi.CL.4[[i]]$mu,
                                        Predicted.A=rowSums(sapply(1:nrow(BayesCpi.coef.CL.4[[i]]), FUN=function(x) geno.split.M[[i]]%*%BayesCpi.coef.CL.4[[i]][x,]))/nrow(BayesCpi.coef.CL.4[[i]]),
                                        Predicted.P=rowSums(sapply(1:nrow(BayesCpi.cfP.CL.4[[i]]), FUN=function(x) geno.split.M[[i]]%*%BayesCpi.cfP.CL.4[[i]][x,]))/nrow(BayesCpi.cfP.CL.4[[i]]),
                                        Predicted.S=rowSums(sapply(1:nrow(BayesCpi.cfS.CL.4[[i]]), FUN=function(x) geno.split.M[[i]]%*%BayesCpi.cfS.CL.4[[i]][x,]))/nrow(BayesCpi.cfS.CL.4[[i]]))
  cat("..",i)
}


#function to plot predicted vs observed trait value.
plot.ObsPred <- function(dat, outFile){
  temp.cor <- cor(dat)
  colnames(dat)[2] <- paste("A (", format(round(temp.cor[2,1],2), nsmall=2), ")", sep="")
  colnames(dat)[3] <- paste("P (", format(round(temp.cor[3,1],2), nsmall=2), ")", sep="")
  colnames(dat)[4] <- paste("S (", format(round(temp.cor[4,1],2), nsmall=2), ")", sep="")
  
  dat <- melt(dat, id.vars="Observed")
  colnames(dat) <- c("Observed", "Type", "Predicted")
  
  temp.xy <- with(dat,range(c(range(Observed),range(Predicted))))

  dat.plot <- ggplot(data=dat) +
    geom_point(aes(x=Observed, y=Predicted, color=Type), size=0.5) +
    geom_abline(slope=1, intercept=0) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    theme(axis.line=element_line(color="#000000")) +
    coord_fixed(ratio=1, xlim=temp.xy, ylim=temp.xy) +
    scale_color_manual(values=c("#000000","#FF0000","#00FF00")) +
    theme(legend.key=element_rect(colour=NA, fill=NA)) +
    theme(legend.background=element_rect(colour=NA, fill=NA)) +
    theme(legend.text=element_text(size=8), legend.title=element_text(size=8))
  
  ggsave(filename=paste(outFile,".svg",sep=""),
         plot=dat.plot,
         width=3,
         height=3,
         units="in")
  
  ggsave(filename=paste(outFile,".png",sep=""),
         plot=dat.plot,
         width=3,
         height=3,
         units="in",
         dpi=300)
}

#function to plot Primary vs Secondary predicted trait value.
plot.PredPS <- function(dat, outFile){
  dat <- dat[,3:4]
  colnames(dat) <- c("Primary", "Secondary")
  
  temp.cor <- cor(dat)
  temp.xy <- with(dat,range(c(range(Primary),range(Secondary))))
  temp.label <- paste("r = ", format(round(temp.cor[2,1],2), nsmall=2),sep="")
  
  dat.plot <- ggplot(data=dat) +
    geom_point(aes(x=Secondary, y=Primary), size=0.5, color="#0000FF") +
    geom_abline(slope=1, intercept=0) +
    annotate("text", x=temp.xy[2]*0.6, y=temp.xy[1]*0.9, label=temp.label) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    theme(axis.line=element_line(color="#000000")) +
    coord_fixed(ratio=1, xlim=temp.xy, ylim=temp.xy) +
    guides(fill=F)
  
  ggsave(filename=paste(outFile,".svg",sep=""),
         plot=dat.plot,
         width=2,
         height=2,
         units="in")
  
  ggsave(filename=paste(outFile,".png",sep=""),
         plot=dat.plot,
         width=2,
         height=2,
         units="in",
         dpi=300)
}

#function to plot marker effects.
plot.mEff <- function(dat, outFile, is.LASSO=F){
  if(!is.null(dim(dat))){dat <- colSums(dat)/ncol(dat)}
  if(is.LASSO){dat <- dat[dat!=0]}
  
  temp.label <- round(100*sum(dat<0)/(sum(dat>0)+sum(dat<0)))
  temp.label <- paste(temp.label, ":", 100-temp.label, sep="")
  
  dat.plot <- ggplot() +
    geom_histogram(aes(x=dat), fill="#000000", alpha=0.7, bins=15) +
    geom_vline(xintercept=0, color="#800000") +
    annotate("text", x=max(dat)*0.7, y=Inf, vjust=1, label=temp.label) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    theme(axis.line=element_line(color="#000000")) +
    xlab("marker effects") +
    guides(fill=F)
  
  ggsave(filename=paste(outFile,".svg",sep=""),
         plot=dat.plot,
         width=2,
         height=2,
         units="in")
  
  ggsave(filename=paste(outFile,".png",sep=""),
         plot=dat.plot,
         width=2,
         height=2,
         units="in",
         dpi=300)
}

#plot joint analyses.
plot.ObsPred(RR.pred.DTH.1, "output2/plotA/RR_DTH_1")
plot.ObsPred(RR.pred.YLD.1, "output2/plotA/RR_YLD_1")
plot.ObsPred(RR.pred.DTS.1, "output2/plotA/RR_DTS_1")
plot.ObsPred(RR.pred.CL.1, "output2/plotA/RR_CL_1")
plot.ObsPred(LASSO.pred.DTH.1, "output2/plotA/LASSO_DTH_1")
plot.ObsPred(LASSO.pred.YLD.1, "output2/plotA/LASSO_YLD_1")
plot.ObsPred(LASSO.pred.DTS.1, "output2/plotA/LASSO_DTS_1")
plot.ObsPred(LASSO.pred.CL.1, "output2/plotA/LASSO_CL_1")
plot.ObsPred(BayesCpi.pred.DTH.1, "output2/plotA/BayesCpi_DTH_1")
plot.ObsPred(BayesCpi.pred.YLD.1, "output2/plotA/BayesCpi_YLD_1")
plot.ObsPred(BayesCpi.pred.DTS.1, "output2/plotA/BayesCpi_DTS_1")
plot.ObsPred(BayesCpi.pred.CL.1, "output2/plotA/BayesCpi_CL_1")

plot.PredPS(RR.pred.DTH.1, "output2/plotB/RR_DTH_1")
plot.PredPS(RR.pred.YLD.1, "output2/plotB/RR_YLD_1")
plot.PredPS(RR.pred.DTS.1, "output2/plotB/RR_DTS_1")
plot.PredPS(RR.pred.CL.1, "output2/plotB/RR_CL_1")
plot.PredPS(LASSO.pred.DTH.1, "output2/plotB/LASSO_DTH_1")
plot.PredPS(LASSO.pred.YLD.1, "output2/plotB/LASSO_YLD_1")
plot.PredPS(LASSO.pred.DTS.1, "output2/plotB/LASSO_DTS_1")
plot.PredPS(LASSO.pred.CL.1, "output2/plotB/LASSO_CL_1")
plot.PredPS(BayesCpi.pred.DTH.1, "output2/plotB/BayesCpi_DTH_1")
plot.PredPS(BayesCpi.pred.YLD.1, "output2/plotB/BayesCpi_YLD_1")
plot.PredPS(BayesCpi.pred.DTS.1, "output2/plotB/BayesCpi_DTS_1")
plot.PredPS(BayesCpi.pred.CL.1, "output2/plotB/BayesCpi_CL_1")

plot.mEff(RR.coef.DTH.1, "output2/plotC/RR_DTH_1", is.LASSO=F)
plot.mEff(RR.coef.YLD.1, "output2/plotC/RR_YLD_1", is.LASSO=F)
plot.mEff(RR.coef.DTS.1, "output2/plotC/RR_DTS_1", is.LASSO=F)
plot.mEff(RR.coef.CL.1, "output2/plotC/RR_CL_1", is.LASSO=F)
plot.mEff(LASSO.coef.DTH.1, "output2/plotC/LASSO_DTH_1", is.LASSO=T)
plot.mEff(LASSO.coef.YLD.1, "output2/plotC/LASSO_YLD_1", is.LASSO=T)
plot.mEff(LASSO.coef.DTS.1, "output2/plotC/LASSO_DTS_1", is.LASSO=T)
plot.mEff(LASSO.coef.CL.1, "output2/plotC/LASSO_CL_1", is.LASSO=T)
plot.mEff(BayesCpi.coef.DTH.1, "output2/plotC/BayesCpi_DTH_1", is.LASSO=F)
plot.mEff(BayesCpi.coef.YLD.1, "output2/plotC/BayesCpi_YLD_1", is.LASSO=F)
plot.mEff(BayesCpi.coef.DTS.1, "output2/plotC/BayesCpi_DTS_1", is.LASSO=F)
plot.mEff(BayesCpi.coef.CL.1, "output2/plotC/BayesCpi_CL_1", is.LASSO=F)

#plot independent analyses
for(i in 1:25){
  plot.ObsPred(RR.pred.DTH.4[[i]], paste("output2/plotA/RR_DTH_4_", i, sep=""))
  plot.ObsPred(RR.pred.YLD.4[[i]], paste("output2/plotA/RR_YLD_4_", i, sep=""))
  plot.ObsPred(RR.pred.DTS.4[[i]], paste("output2/plotA/RR_DTS_4_", i, sep=""))
  plot.ObsPred(RR.pred.CL.4[[i]], paste("output2/plotA/RR_CL_4_", i, sep=""))
  plot.ObsPred(LASSO.pred.DTH.4[[i]], paste("output2/plotA/LASSO_DTH_4_", i, sep=""))
  plot.ObsPred(LASSO.pred.YLD.4[[i]], paste("output2/plotA/LASSO_YLD_4_", i, sep=""))
  plot.ObsPred(LASSO.pred.DTS.4[[i]], paste("output2/plotA/LASSO_DTS_4_", i, sep=""))
  plot.ObsPred(LASSO.pred.CL.4[[i]], paste("output2/plotA/LASSO_CL_4_", i, sep=""))
  plot.ObsPred(BayesCpi.pred.DTH.4[[i]], paste("output2/plotA/BayesCpi_DTH_4_", i, sep=""))
  plot.ObsPred(BayesCpi.pred.YLD.4[[i]], paste("output2/plotA/BayesCpi_YLD_4_", i, sep=""))
  plot.ObsPred(BayesCpi.pred.DTS.4[[i]], paste("output2/plotA/BayesCpi_DTS_4_", i, sep=""))
  plot.ObsPred(BayesCpi.pred.CL.4[[i]], paste("output2/plotA/BayesCpi_CL_4_", i, sep=""))
  
  plot.PredPS(RR.pred.DTH.4[[i]], paste("output2/plotB/RR_DTH_4_", i, sep=""))
  plot.PredPS(RR.pred.YLD.4[[i]], paste("output2/plotB/RR_YLD_4_", i, sep=""))
  plot.PredPS(RR.pred.DTS.4[[i]], paste("output2/plotB/RR_DTS_4_", i, sep=""))
  plot.PredPS(RR.pred.CL.4[[i]], paste("output2/plotB/RR_CL_4_", i, sep=""))
  plot.PredPS(LASSO.pred.DTH.4[[i]], paste("output2/plotB/LASSO_DTH_4_", i, sep=""))
  plot.PredPS(LASSO.pred.YLD.4[[i]], paste("output2/plotB/LASSO_YLD_4_", i, sep=""))
  plot.PredPS(LASSO.pred.DTS.4[[i]], paste("output2/plotB/LASSO_DTS_4_", i, sep=""))
  plot.PredPS(LASSO.pred.CL.4[[i]], paste("output2/plotB/LASSO_CL_4_", i, sep=""))
  plot.PredPS(BayesCpi.pred.DTH.4[[i]], paste("output2/plotB/BayesCpi_DTH_4_", i, sep=""))
  plot.PredPS(BayesCpi.pred.YLD.4[[i]], paste("output2/plotB/BayesCpi_YLD_4_", i, sep=""))
  plot.PredPS(BayesCpi.pred.DTS.4[[i]], paste("output2/plotB/BayesCpi_DTS_4_", i, sep=""))
  plot.PredPS(BayesCpi.pred.CL.4[[i]], paste("output2/plotB/BayesCpi_CL_4_", i, sep=""))
  
  plot.mEff(RR.coef.DTH.4[[i]], paste("output2/plotC/RR_DTH_4_", i, sep=""), is.LASSO=F)
  plot.mEff(RR.coef.YLD.4[[i]], paste("output2/plotC/RR_YLD_4_", i, sep=""), is.LASSO=F)
  plot.mEff(RR.coef.DTS.4[[i]], paste("output2/plotC/RR_DTS_4_", i, sep=""), is.LASSO=F)
  plot.mEff(RR.coef.CL.4[[i]], paste("output2/plotC/RR_CL_4_", i, sep=""), is.LASSO=F)
  plot.mEff(LASSO.coef.DTH.4[[i]], paste("output2/plotC/LASSO_DTH_4_", i, sep=""), is.LASSO=T)
  plot.mEff(LASSO.coef.YLD.4[[i]], paste("output2/plotC/LASSO_YLD_4_", i, sep=""), is.LASSO=T)
  plot.mEff(LASSO.coef.DTS.4[[i]], paste("output2/plotC/LASSO_DTS_4_", i, sep=""), is.LASSO=T)
  plot.mEff(LASSO.coef.CL.4[[i]], paste("output2/plotC/LASSO_CL_4_", i, sep=""), is.LASSO=T)
  plot.mEff(BayesCpi.coef.DTH.4[[i]], paste("output2/plotC/BayesCpi_DTH_4_", i, sep=""), is.LASSO=F)
  plot.mEff(BayesCpi.coef.YLD.4[[i]], paste("output2/plotC/BayesCpi_YLD_4_", i, sep=""), is.LASSO=F)
  plot.mEff(BayesCpi.coef.DTS.4[[i]], paste("output2/plotC/BayesCpi_DTS_4_", i, sep=""), is.LASSO=F)
  plot.mEff(BayesCpi.coef.CL.4[[i]], paste("output2/plotC/BayesCpi_CL_4_", i, sep=""), is.LASSO=F)
}

temp.DTH <- temp.YLD <- temp.DTS <- temp.CL <- vector()
for(i in 1:25){
  temp.DTH <- rbind(temp.DTH, cor(RR.pred.DTH.4[[i]])[2:4,1])
  temp.YLD <- rbind(temp.YLD, cor(RR.pred.YLD.4[[i]])[2:4,1])
  temp.DTS <- rbind(temp.DTS, cor(RR.pred.DTS.4[[i]])[2:4,1])
  temp.CL <- rbind(temp.CL, cor(RR.pred.CL.4[[i]])[2:4,1])
}
temp.DTH <- data.frame(temp.DTH)
temp.YLD <- data.frame(temp.YLD)
temp.DTS <- data.frame(temp.DTS)
temp.CL <- data.frame(temp.CL)
temp.DTH$Group <- "DTH (Barley)"
temp.YLD$Group <- "YLD (Barley)"
temp.DTS$Group <- "DTS (Maize)"
temp.CL$Group <- "CL (Maize)"
cor.A <- rbind(temp.DTH[,c(1,4)],
               temp.YLD[,c(1,4)],
               temp.DTS[,c(1,4)],
               temp.CL[,c(1,4)])
cor.P <- rbind(temp.DTH[,c(2,4)],
               temp.YLD[,c(2,4)],
               temp.DTS[,c(2,4)],
               temp.CL[,c(2,4)])
cor.S <- rbind(temp.DTH[,c(3,4)],
               temp.YLD[,c(3,4)],
               temp.DTS[,c(3,4)],
               temp.CL[,c(3,4)])
colnames(cor.A)[1] <- colnames(cor.P)[1] <- colnames(cor.S)[1] <- "Correlation"
cor.A$Group <- factor(cor.A$Group, levels=c("DTH (Barley)","YLD (Barley)","DTS (Maize)","CL (Maize)"))
cor.P$Group <- factor(cor.P$Group, levels=c("DTH (Barley)","YLD (Barley)","DTS (Maize)","CL (Maize)"))
cor.S$Group <- factor(cor.S$Group, levels=c("DTH (Barley)","YLD (Barley)","DTS (Maize)","CL (Maize)"))

cor.A$Type <- "A"
cor.P$Type <- "P"
cor.S$Type <- "S"

cor.dat <- rbind(cor.A, cor.P, cor.S)
cor.dat$Type <- factor(cor.dat$Type, levels=c("A","P","S"))

cor.dat0 <- data.frame(Correlation=c(cor(RR.pred.DTH.1)[2:4,1],
                                     cor(RR.pred.YLD.1)[2:4,1],
                                     cor(RR.pred.DTS.1)[2:4,1],
                                     cor(RR.pred.CL.1)[2:4,1]),
                       Group=c(rep("DTH (Barley)",3),
                               rep("YLD (Barley)",3),
                               rep("DTS (Maize)",3),
                               rep("CL (Maize)",3)),
                       Type=rep(c("A","P","S"),4))

plot.APS <- ggplot() +
  geom_boxplot(data=cor.dat, aes(x=Type, y=Correlation), outlier.size=0.5) +
  facet_wrap(vars(Group), nrow=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(strip.background=element_blank()) +
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, color="#000000", fill=NA) +
  geom_point(data=cor.dat0, aes(x=Type, y=Correlation), color="#0000FF", shape=15, size=2)

ggsave(filename="cor_APS.svg",
       plot=plot.APS,
       width=6,
       height=2.25,
       units="in")
ggsave(filename="cor_APS.png",
       plot=plot.APS,
       width=6,
       height=2.25,
       units="in",
       dpi=600)


cor.all.APS <- replicate(12, vector())
for(i in 1:25){
  cor.all.APS[[1]] <- rbind(cor.all.APS[[1]],
                            cor(RR.pred.DTH.4[[i]])[1,2:4])
  cor.all.APS[[2]] <- rbind(cor.all.APS[[2]],
                            cor(LASSO.pred.DTH.4[[i]])[1,2:4])
  cor.all.APS[[3]] <- rbind(cor.all.APS[[3]],
                            cor(BayesCpi.pred.DTH.4[[i]])[1,2:4])
  cor.all.APS[[4]] <- rbind(cor.all.APS[[4]],
                            cor(RR.pred.YLD.4[[i]])[1,2:4])
  cor.all.APS[[5]] <- rbind(cor.all.APS[[5]],
                            cor(LASSO.pred.YLD.4[[i]])[1,2:4])
  cor.all.APS[[6]] <- rbind(cor.all.APS[[6]],
                            cor(BayesCpi.pred.YLD.4[[i]])[1,2:4])
  cor.all.APS[[7]] <- rbind(cor.all.APS[[7]],
                            cor(RR.pred.DTS.4[[i]])[1,2:4])
  cor.all.APS[[8]] <- rbind(cor.all.APS[[8]],
                            cor(LASSO.pred.DTS.4[[i]])[1,2:4])
  cor.all.APS[[9]] <- rbind(cor.all.APS[[9]],
                            cor(BayesCpi.pred.DTS.4[[i]])[1,2:4])
  cor.all.APS[[10]] <- rbind(cor.all.APS[[10]],
                            cor(RR.pred.CL.4[[i]])[1,2:4])
  cor.all.APS[[11]] <- rbind(cor.all.APS[[11]],
                            cor(LASSO.pred.CL.4[[i]])[1,2:4])
  cor.all.APS[[12]] <- rbind(cor.all.APS[[12]],
                            cor(BayesCpi.pred.CL.4[[i]])[1,2:4])
}

for(i in 1:12){
  cor.all.APS[[i]] <- cbind(cor.all.APS[[i]], Method=(i-1)%%3, Trait=(i-1)%/%3)
}

cor.all.APS <- do.call(rbind, cor.all.APS)
cor.all.APS <- data.frame(cor.all.APS)
colnames(cor.all.APS)[1:3] <- c("A","P","S")
cor.all.APS$Method <- as.factor(cor.all.APS$Method)
cor.all.APS$Method <- factor(cor.all.APS$Method,
                             levels=c(0,1,2),
                             labels=c("rrBLUP","LASSO","BayesCpi"))
cor.all.APS$Trait <- as.factor(cor.all.APS$Trait)
cor.all.APS$Trait <- factor(cor.all.APS$Trait,
                            levels=c(0,1,2,3),
                            labels=c("DTH","YLD","DTS","CL"))
write.csv(cor.all.APS, "cor_all_APS.csv", quote=F, row.names=F)
cor.all.APS <- melt(cor.all.APS, id.vars=c("Method", "Trait"))

cor.all.APS2 <- rbind(cor(RR.pred.DTH.1)[1,2:4],
                      cor(LASSO.pred.DTH.1)[1,2:4],
                      cor(BayesCpi.pred.DTH.1)[1,2:4],
                      cor(RR.pred.YLD.1)[1,2:4],
                      cor(LASSO.pred.YLD.1)[1,2:4],
                      cor(BayesCpi.pred.YLD.1)[1,2:4],
                      cor(RR.pred.DTS.1)[1,2:4],
                      cor(LASSO.pred.DTS.1)[1,2:4],
                      cor(BayesCpi.pred.DTS.1)[1,2:4],
                      cor(RR.pred.CL.1)[1,2:4],
                      cor(LASSO.pred.CL.1)[1,2:4],
                      cor(BayesCpi.pred.CL.1)[1,2:4])
colnames(cor.all.APS2) <- c("A","P","S")
cor.all.APS2 <- data.frame(cor.all.APS2,
                           Method=rep(c("rrBLUP","LASSO","BayesCpi"), 4),
                           Trait=c(rep("DTH",3),rep("YLD",3),rep("DTS",3),rep("CL",3)))
cor.all.APS2$Method <- as.factor(cor.all.APS2$Method)
cor.all.APS2$Method <- factor(cor.all.APS2$Method,
                              levels=c("rrBLUP","LASSO","BayesCpi"),
                              labels=c("rrBLUP","LASSO","BayesCpi"))
cor.all.APS2$Trait <- as.factor(cor.all.APS2$Trait)
cor.all.APS2$Trait <- factor(cor.all.APS2$Trait,
                             levels=c("DTH","YLD","DTS","CL"),
                             labels=c("DTH","YLD","DTS","CL"))
cor.all.APS2 <- melt(cor.all.APS2, id.vars=c("Method", "Trait"))

ggplot() +
  geom_boxplot(data=cor.all.APS, aes(x=variable, y=value), na.rm=T) +
  geom_point(data=cor.all.APS2, aes(x=variable, y=value), shape=15, size=1.5, color="#0000FF") +
  facet_grid(Trait~Method) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  scale_y_continuous(breaks=c(0,0.5,1.0), limits=c(-0.2,1.2)) +
  theme(axis.title.x=element_blank()) +
  ylab("Correlation")

ggsave(filename="cor_all_APS.svg",
       scale=1,
       width=4.5,
       height=6,
       units="in")
ggsave(filename="cor_all_APS.png",
       scale=1,
       width=4.5,
       height=6,
       units="in",
       dpi=600)

misc.summary <- replicate(12, vector())
for(i in 1:25){
  misc.summary[[1]] <- rbind(misc.summary[[1]],
                             c(nrow(RR.pred.DTH.4[[i]]),
                               sum(RR.coef.DTH.4[[i]]!=0)/length(RR.coef.DTH.4[[i]]),
                               sum(RR.coef.DTH.4[[i]]<0)/length(RR.coef.DTH.4[[i]]),
                               sum(RR.coef.DTH.4[[i]]>0)/length(RR.coef.DTH.4[[i]])))
  misc.summary[[2]] <- rbind(misc.summary[[2]],
                             c(nrow(LASSO.pred.DTH.4[[i]]),
                               sum(LASSO.coef.DTH.4[[i]]!=0)/length(LASSO.coef.DTH.4[[i]]),
                               sum(LASSO.coef.DTH.4[[i]]<0)/length(LASSO.coef.DTH.4[[i]]),
                               sum(LASSO.coef.DTH.4[[i]]>0)/length(LASSO.coef.DTH.4[[i]])))
  misc.summary[[3]] <- rbind(misc.summary[[3]],
                             c(nrow(BayesCpi.pred.DTH.4[[i]]),
                               sum(colSums(BayesCpi.coef.DTH.4[[i]])!=0)/ncol(BayesCpi.coef.DTH.4[[i]]),
                               sum(colSums(BayesCpi.coef.DTH.4[[i]])<0)/ncol(BayesCpi.coef.DTH.4[[i]]),
                               sum(colSums(BayesCpi.coef.DTH.4[[i]])>0)/ncol(BayesCpi.coef.DTH.4[[i]])))

  misc.summary[[4]] <- rbind(misc.summary[[4]],
                             c(nrow(RR.pred.YLD.4[[i]]),
                               sum(RR.coef.YLD.4[[i]]!=0)/length(RR.coef.YLD.4[[i]]),
                               sum(RR.coef.YLD.4[[i]]<0)/length(RR.coef.YLD.4[[i]]),
                               sum(RR.coef.YLD.4[[i]]>0)/length(RR.coef.YLD.4[[i]])))
  misc.summary[[5]] <- rbind(misc.summary[[5]],
                             c(nrow(LASSO.pred.YLD.4[[i]]),
                               sum(LASSO.coef.YLD.4[[i]]!=0)/length(LASSO.coef.YLD.4[[i]]),
                               sum(LASSO.coef.YLD.4[[i]]<0)/length(LASSO.coef.YLD.4[[i]]),
                               sum(LASSO.coef.YLD.4[[i]]>0)/length(LASSO.coef.YLD.4[[i]])))
  misc.summary[[6]] <- rbind(misc.summary[[6]],
                             c(nrow(BayesCpi.pred.YLD.4[[i]]),
                               sum(colSums(BayesCpi.coef.YLD.4[[i]])!=0)/ncol(BayesCpi.coef.YLD.4[[i]]),
                               sum(colSums(BayesCpi.coef.YLD.4[[i]])<0)/ncol(BayesCpi.coef.YLD.4[[i]]),
                               sum(colSums(BayesCpi.coef.YLD.4[[i]])>0)/ncol(BayesCpi.coef.YLD.4[[i]])))

  misc.summary[[7]] <- rbind(misc.summary[[7]],
                             c(nrow(RR.pred.DTS.4[[i]]),
                               sum(RR.coef.DTS.4[[i]]!=0)/length(RR.coef.DTS.4[[i]]),
                               sum(RR.coef.DTS.4[[i]]<0)/length(RR.coef.DTS.4[[i]]),
                               sum(RR.coef.DTS.4[[i]]>0)/length(RR.coef.DTS.4[[i]])))
  misc.summary[[8]] <- rbind(misc.summary[[8]],
                             c(nrow(LASSO.pred.DTS.4[[i]]),
                               sum(LASSO.coef.DTS.4[[i]]!=0)/length(LASSO.coef.DTS.4[[i]]),
                               sum(LASSO.coef.DTS.4[[i]]<0)/length(LASSO.coef.DTS.4[[i]]),
                               sum(LASSO.coef.DTS.4[[i]]>0)/length(LASSO.coef.DTS.4[[i]])))
  misc.summary[[9]] <- rbind(misc.summary[[9]],
                             c(nrow(BayesCpi.pred.DTS.4[[i]]),
                               sum(colSums(BayesCpi.coef.DTS.4[[i]])!=0)/ncol(BayesCpi.coef.DTS.4[[i]]),
                               sum(colSums(BayesCpi.coef.DTS.4[[i]])<0)/ncol(BayesCpi.coef.DTS.4[[i]]),
                               sum(colSums(BayesCpi.coef.DTS.4[[i]])>0)/ncol(BayesCpi.coef.DTS.4[[i]])))

  misc.summary[[10]] <- rbind(misc.summary[[10]],
                             c(nrow(RR.pred.CL.4[[i]]),
                               sum(RR.coef.CL.4[[i]]!=0)/length(RR.coef.CL.4[[i]]),
                               sum(RR.coef.CL.4[[i]]<0)/length(RR.coef.CL.4[[i]]),
                               sum(RR.coef.CL.4[[i]]>0)/length(RR.coef.CL.4[[i]])))
  misc.summary[[11]] <- rbind(misc.summary[[11]],
                             c(nrow(LASSO.pred.CL.4[[i]]),
                               sum(LASSO.coef.CL.4[[i]]!=0)/length(LASSO.coef.CL.4[[i]]),
                               sum(LASSO.coef.CL.4[[i]]<0)/length(LASSO.coef.CL.4[[i]]),
                               sum(LASSO.coef.CL.4[[i]]>0)/length(LASSO.coef.CL.4[[i]])))
  misc.summary[[12]] <- rbind(misc.summary[[12]],
                             c(nrow(BayesCpi.pred.CL.4[[i]]),
                               sum(colSums(BayesCpi.coef.CL.4[[i]])!=0)/ncol(BayesCpi.coef.CL.4[[i]]),
                               sum(colSums(BayesCpi.coef.CL.4[[i]])<0)/ncol(BayesCpi.coef.CL.4[[i]]),
                               sum(colSums(BayesCpi.coef.CL.4[[i]])>0)/ncol(BayesCpi.coef.CL.4[[i]])))
}

for(i in 1:12){
  misc.summary[[i]] <- cbind(misc.summary[[i]], Method=(i-1)%%3, Trait=(i-1)%/%3)
}

misc.summary <- do.call(rbind, misc.summary)
colnames(misc.summary)[1:4] <- c("n","A","P","S")
write.csv(misc.summary, "misc_summary.csv", quote=F, row.names=F)


#plot the comparison of predicted values from rrBLUP, LASSO, BayesCpi.
temp.DTH <- sapply(1:25, FUN=function(x) cor(RR.pred.DTH.4[[x]])[4,3])
temp.YLD <- sapply(1:25, FUN=function(x) cor(RR.pred.YLD.4[[x]])[4,3])
temp.DTS <- sapply(1:25, FUN=function(x) cor(RR.pred.DTS.4[[x]])[4,3])
temp.CL <- sapply(1:25, FUN=function(x) cor(RR.pred.CL.4[[x]])[4,3])

cor.PS <- data.frame(Correlation=c(temp.DTH,
                                   temp.YLD,
                                   temp.DTS,
                                   temp.CL),
                     Trait=c(rep("DTH",25),
                             rep("YLD",25),
                             rep("DTS",25),
                             rep("CL",25)))

cor.PS$Trait <- as.factor(cor.PS$Trait)
cor.PS$Trait <- factor(cor.PS$Trait, levels=c("DTH","YLD","DTS","CL"))

plot.PS <- ggplot() +
  geom_boxplot(data=cor.PS, aes(x=Trait, y=Correlation), outlier.size=0.5) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, color="#000000", fill=NA) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=0, vjust=0.5))

ggsave(filename="cor_PS.svg",
       plot=plot.PS,
       width=2,
       height=2,
       units="in")
ggsave(filename="cor_PS.png",
       plot=plot.PS,
       width=2,
       height=2,
       units="in",
       dpi=600)

cor.pred <- replicate(12, list(replicate(3, vector())))

for(i in 1:25){
  for(j in 1:3){
    cor.pred[[1]][[j]] <- c(cor.pred[[1]][[j]],
                            cor(RR.pred.DTH.4[[i]][,j+1], BayesCpi.pred.DTH.4[[i]][,j+1]))
    cor.pred[[2]][[j]] <- c(cor.pred[[2]][[j]],
                            cor(RR.pred.DTH.4[[i]][,j+1], LASSO.pred.DTH.4[[i]][,j+1]))
    cor.pred[[3]][[j]] <- c(cor.pred[[3]][[j]],
                            cor(LASSO.pred.DTH.4[[i]][,j+1], BayesCpi.pred.DTH.4[[i]][,j+1]))
    cor.pred[[4]][[j]] <- c(cor.pred[[4]][[j]],
                            cor(RR.pred.YLD.4[[i]][,j+1], BayesCpi.pred.YLD.4[[i]][,j+1]))
    cor.pred[[5]][[j]] <- c(cor.pred[[5]][[j]],
                            cor(RR.pred.YLD.4[[i]][,j+1], LASSO.pred.YLD.4[[i]][,j+1]))
    cor.pred[[6]][[j]] <- c(cor.pred[[6]][[j]],
                            cor(LASSO.pred.YLD.4[[i]][,j+1], BayesCpi.pred.YLD.4[[i]][,j+1]))
    cor.pred[[7]][[j]] <- c(cor.pred[[7]][[j]],
                            cor(RR.pred.DTS.4[[i]][,j+1], BayesCpi.pred.DTS.4[[i]][,j+1]))
    cor.pred[[8]][[j]] <- c(cor.pred[[8]][[j]],
                            cor(RR.pred.DTS.4[[i]][,j+1], LASSO.pred.DTS.4[[i]][,j+1]))
    cor.pred[[9]][[j]] <- c(cor.pred[[9]][[j]],
                            cor(LASSO.pred.DTS.4[[i]][,j+1], BayesCpi.pred.DTS.4[[i]][,j+1]))
    cor.pred[[10]][[j]] <- c(cor.pred[[10]][[j]],
                             cor(RR.pred.CL.4[[i]][,j+1], BayesCpi.pred.CL.4[[i]][,j+1]))
    cor.pred[[11]][[j]] <- c(cor.pred[[11]][[j]],
                             cor(RR.pred.CL.4[[i]][,j+1], LASSO.pred.CL.4[[i]][,j+1]))
    cor.pred[[12]][[j]] <- c(cor.pred[[12]][[j]],
                             cor(LASSO.pred.CL.4[[i]][,j+1], BayesCpi.pred.CL.4[[i]][,j+1]))
  }
}

for(i in 1:12){
  cor.pred[[i]] <- do.call(cbind, c(cor.pred[[i]], list(rep((i-1)%%3,25),rep((i-1)%/%3,25))))
  colnames(cor.pred[[i]]) <- c("A","P","S","Method","Trait")
}
cor.pred <- do.call(rbind, cor.pred)
cor.pred <- data.frame(cor.pred)

cor.pred$Method <- as.factor(cor.pred$Method)
cor.pred$Method <- factor(cor.pred$Method,
                          levels=c(0,1,2),
                          labels=c("rrBLUP - BayesCpi", "rrBLUP - LASSO", "LASSO - BayesCpi"))

cor.pred$Trait <- as.factor(cor.pred$Trait)
cor.pred$Trait <- factor(cor.pred$Trait,
                         levels=c(0,1,2,3),
                         labels=c("DTH", "YLD", "DTS", "CL"))

cor.pred <- melt(cor.pred, id.vars=c("Method","Trait"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  geom_boxplot(data=cor.pred, aes(x=variable, y=value)) +
  facet_grid(Trait~Method,) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.title.x=element_blank()) +
  ylab("Correlation")

ggsave(filename="Method_comparison.svg",
       scale=1,
       width=4.5,
       height=6,
       units="in")

ggsave(filename="Method_comparison.png",
       scale=1,
       width=4.5,
       height=6,
       units="in",
       dpi=600)

#Kolmogorov-Smirnov test for favorable primary-vs-secondary allele effects       
ks.DTH <- sapply(1:25, FUN=function(x) ks.test(abs(RR.coef.DTH.4[[x]][RR.coef.DTH.4[[x]]<0]), abs(RR.coef.DTH.4[[x]][RR.coef.DTH.4[[x]]>0]))$p.value)
ks.YLD <- sapply(1:25, FUN=function(x) ks.test(abs(RR.coef.YLD.4[[x]][RR.coef.YLD.4[[x]]<0]), abs(RR.coef.YLD.4[[x]][RR.coef.YLD.4[[x]]>0]))$p.value)
ks.DTS <- sapply(1:25, FUN=function(x) ks.test(abs(RR.coef.DTS.4[[x]][RR.coef.DTS.4[[x]]<0]), abs(RR.coef.DTS.4[[x]][RR.coef.DTS.4[[x]]>0]))$p.value)
ks.CL <- sapply(1:25, FUN=function(x) ks.test(abs(RR.coef.CL.4[[x]][RR.coef.CL.4[[x]]<0]), abs(RR.coef.CL.4[[x]][RR.coef.CL.4[[x]]>0]))$p.value)

ks.dat <- cbind(ks.DTH, ks.YLD, ks.DTS, ks.CL)
ks.dat[ks.dat==0] <- .Machine$double.neg.eps
ks.dat <- -log10(ks.dat)
colnames(ks.dat) <- c("DTH","YLD","DTS","CL")
ks.dat <- melt(ks.dat)
ks.dat <- ks.dat[,2:3]
colnames(ks.dat) <- c("Trait", "Value")

set.seed(32424)
ks.plot <- ggplot() +
  geom_jitter(data=ks.dat, aes(x=Trait, y=Value), width=0.3, color="#0000FF", size=2) +
  geom_vline(xintercept=c(1.5,2.5,3.5), color="#FFFFFF") +
  geom_hline(yintercept=-log10(0.05/25), color="#FF0000", alpha=0.5) +
  theme(axis.title.x=element_blank()) +
  ylab("-log10(p)") +
  theme(panel.grid=element_blank())

ggsave(filename="KS_test.svg",
       plot=ks.plot,
       scale=1.3,
       width=3,
       height=1.5,
       units="in")
ggsave(filename="KS_test.png",
       plot=ks.plot,
       scale=1.3,
       width=3,
       height=1.5,
       units="in",
       dpi=600)

#save all of the current data.
save.image("OSGS_dat2.RData")


##################################################################
### 3. Simulating one generation of selection under OSGS vs GS ###
##################################################################

#alphasim requires complete marker data, i.e. no +/-0.5
#we need to approximate the marker data.
#change -0.5 to -1.
#change 0.5 to 1.
#change single 0 to the prior marker.
geno.convert <- function(dat.geno, dat.marker){
  marker.by.chr <- lapply(unique(dat.marker$Chr), FUN=function(x) which(dat.marker$Chr==x))
  dat.geno[dat.geno==-0.5] <- -1
  dat.geno[dat.geno==0.5] <- 1
  geno.by.chr <- lapply(1:length(marker.by.chr), FUN=function(x) dat.geno[,marker.by.chr[[x]]])
  for(i in 1:length(geno.by.chr)){
    for(j in 1:nrow(geno.by.chr[[i]])){
      k <- 2
      while(k < ncol(geno.by.chr[[i]])){
        if(geno.by.chr[[i]][j,k]==0 & geno.by.chr[[i]][j,k]!=geno.by.chr[[i]][j,k-1] & geno.by.chr[[i]][j,k]!=geno.by.chr[[i]][j,k+1]) geno.by.chr[[i]][j,k] <- geno.by.chr[[i]][j,k-1]
        k <- k + 1
      }
    }
  }
  geno.by.chr.A <- geno.by.chr.B <- lapply(1:length(geno.by.chr), FUN=function(x) matrix(0, ncol=ncol(geno.by.chr[[x]]), nrow=nrow(geno.by.chr[[x]])))
  for(i in 1:length(geno.by.chr)){
    geno.by.chr.A[[i]][geno.by.chr[[i]]==1] <- 1
    geno.by.chr.B[[i]][geno.by.chr[[i]]==1] <- 1
    geno.by.chr.A[[i]][geno.by.chr[[i]]==0] <- 1
  }
  geno.by.chr <- lapply(1:length(geno.by.chr), FUN=function(x) rbind(geno.by.chr.A[[x]], geno.by.chr.B[[x]]))
  return(geno.by.chr)
}

#function given by R Chris Gaynor (AlphaSim author) because the current version is bugged.
newMapPop2 = function(genMap,haplotypes,inbred=FALSE,
                      ploidy=2L){
  stopifnot(length(genMap)==length(haplotypes))
  nRow = lapply(haplotypes,nrow)
  nRow = unlist(nRow)
  ploidy = as.integer(ploidy)
  if(length(nRow)>1L){
    if(any(nRow[1]!=nRow)){
      stop("Number of rows must be equal in haplotypes")
    }
    nRow = nRow[1]
  }
  if(inbred){
    nInd = nRow
  }else{
    if(nRow%%ploidy != 0L){
      stop("Number of haplotypes must be divisible by ploidy")
    }
    nInd = nRow/ploidy
  }
  nCol = lapply(haplotypes,ncol)
  nCol = unlist(nCol)
  segSites = lapply(genMap,length)
  segSites = unlist(segSites)
  if(!all(nCol == segSites)){
    stop("Number of segregating sites in haplotypes and genMap don't match")
  }
  output = vector("list",length(genMap))
  for(chr in 1:length(genMap)){
    geno = AlphaSimR:::packHaplo(as.matrix(haplotypes[[chr]]),
                                 ploidy=ploidy,inbred=inbred)
    output[[chr]] = new("MapPop",
                        nInd=as.integer(nInd),
                        nChr=1L,
                        ploidy=ploidy,
                        nLoci=as.integer(segSites[chr]),
                        geno=as.matrix(list(geno)),
                        genMap=as.matrix(genMap[chr]),
                        centromere=max(genMap[[chr]])/2)
  }
  output = do.call("cChr",output)
  return(output)
}

#function to pick top/bottom n-lines.
pick.line <- function(ebv, n, wt, top=T){
  ebv <- lapply(1:ncol(ebv), FUN=function(x) ebv[,x])
  for(i in 1:length(ebv)){
    names(ebv[[i]]) <- 1:length(ebv[[i]])
    ebv[[i]] <- if(top) rank(ebv[[i]], ties.method="random") else rank(-ebv[[i]], ties.method="random")
  }

  ebv <- if(wt<0) {
           sort(ebv[[2]], decreasing=T)
         } else {
           sort(wt*ebv[[3]]+(1-wt)*ebv[[4]], decreasing=T)
         }
  ebv <- sort(as.numeric(names(ebv)[1:n]))
  return(ebv)
}

#function to calculate percent of favorable primary/secondary and breeding values.
calc.fav <- function(geno, m.eff, wt, pop, top=T){
  fav <- which(m.eff!=0)
  
  if(top){
    favP <- which(m.eff<0)
    favS <- which(m.eff>0)
  } else {
    favP <- which(m.eff>0)
    favS <- which(m.eff<0)
  }
  
  out.favP <- (rowSums(geno[,favP]==-1) + 0.5*rowSums(geno[,favP]==0))/length(fav)
  out.favS <- (rowSums(geno[,favS]==1) + 0.5*rowSums(geno[,favS]==0))/length(fav)
  out.favA <- out.favP + out.favS
  
  if(wt<0) wt <- "GS"
  if(wt==9) wt <- "B"
  
  out.fav <- data.frame(favP=mean(out.favP),
                        favS=mean(out.favS),
                        favA=mean(out.favA),
                        ebv_mean=mean(geno%*%m.eff),
                        ebv_sd=sd(geno%*%m.eff),
                        wt=wt,
                        pop=pop)
  
  out.fav$wt <- as.character(out.fav$wt)
  
  return(out.fav)
}

#function to thin out linked markers.
thin.geno <- function(geno, GenPos){
  temp.GenPos <- lapply(1:length(GenPos), FUN=function(x) seq(0,max(GenPos[[x]])+0.01,0.01))
  thin.GenPos <- rep(list(vector()), length(temp.GenPos))
  for(i in 1:length(temp.GenPos)){
    for(j in 2:length(temp.GenPos[[i]])){
      thin.GenPos[[i]] <- c(thin.GenPos[[i]], which(GenPos[[i]] >= temp.GenPos[[i]][j-1] & GenPos[[i]] < temp.GenPos[[i]][j])[1])
    }
  }
  out.GenPos <- lapply(1:length(GenPos), FUN=function(x) temp.GenPos[[x]][-length(temp.GenPos[[x]])])
  out.GenPos <- lapply(1:length(GenPos), FUN=function(x) out.GenPos[[x]][!is.na(thin.GenPos[[x]])])
  out.marker.A <- lapply(1:length(GenPos), FUN=function(x) thin.GenPos[[x]][!is.na(thin.GenPos[[x]])])
  out.geno <- lapply(1:length(GenPos), FUN=function(x) geno[[x]][,out.marker.A[[x]]])
  
  temp.marker <- lapply(1:length(GenPos), FUN=function(x) c(out.marker.A[[x]], length(GenPos[[x]])+1))
  out.marker.B <- rep(list(vector()), length(temp.marker))
  for(i in 1:length(temp.marker)){
    for(j in 2:length(temp.marker[[i]])){
      out.marker.B[[i]] <- c(out.marker.B[[i]], temp.marker[[i]][j]-temp.marker[[i]][j-1])
    }
  }
  
  return(list(out.geno, out.GenPos, out.marker.A, out.marker.B))
}

#prepare barley NAM data for AlphaSim
n.ind.B <- sapply(1:25, FUN=function(x) nrow(geno.split.B[[x]]))
n.chr.B <- 7
gen.pos.B <- lapply(unique(marker.B$Chr), FUN=function(x) marker.B[marker.B$Chr==x,4]/100)
geno.founder.B <- lapply(1:25, FUN=function(x) geno.convert(geno.split.B[[x]], marker.B))
geno.founder.B <- lapply(1:25, FUN=function(x) thin.geno(geno.founder.B[[x]], gen.pos.B))

#prepare maize NAM data for AlphaSim
n.ind.M <- sapply(1:25, FUN=function(x) nrow(geno.split.M[[x]]))
n.chr.M <- 10
gen.pos.M <- lapply(unique(marker.M$Chr), FUN=function(x) marker.M[marker.M$Chr==x,4]/100)
geno.founder.M <- lapply(1:25, FUN=function(x) geno.convert(geno.split.M[[x]], marker.M))
geno.founder.M <- lapply(1:25, FUN=function(x) thin.geno(geno.founder.M[[x]], gen.pos.M))


#loop to iterate through 25 barley NAM families.
set.seed(93944)
for(i in 1:25){
  #create barley NAM base population in AlphaSim.
  founder.B <- newMapPop2(genMap=geno.founder.B[[i]][[2]],
                          haplotypes=geno.founder.B[[i]][[1]],
                          inbred=T,
                          ploidy=2L)
  SP <- SimParam$new(founder.B)
  founder.B <- newPop(founder.B, simParam=SP)
  base.B <- makeCross(pop=founder.B,
                      crossPlan=matrix(1:founder.B@nInd, nrow=founder.B@nInd/2, ncol=2, byrow=F),
                      simParam=SP)
  marker.rep <- unlist(geno.founder.B[[i]][[4]])
  
  #loop to iterate through different types of selection weights in OSGS (-0.1 is GS).
  for(j in seq(-0.1,1,0.1)){
    #identify top 4 lines.
    sel.line.DTH <- pick.line(ebv=RR.pred.DTH.4[[i]],
                              n=4,
                              wt=j,
                              top=T)
    sel.line.YLD <- pick.line(ebv=RR.pred.YLD.4[[i]],
                              n=4,
                              wt=j,
                              top=T)
    
    #cross the top 4 lines using all possible combinations.
    F1.DTH <- makeCross(pop=base.B,
                        crossPlan=t(combn(sel.line.DTH,2)),
                        simParam=SP)
    F1.YLD <- makeCross(pop=base.B,
                        crossPlan=t(combn(sel.line.YLD,2)),
                        simParam=SP)

    #create 10 double-haploids from each F1 (a total of 60 DHs).
    DH.DTH <- makeDH(pop=F1.DTH,
                     nDH=10,
                     simParam=SP)
    DH.YLD <- makeDH(pop=F1.YLD,
                     nDH=10,
                     simParam=SP)
    
    #extract the DH genotypes.
    geno.DTH <- pullSegSiteGeno(DH.DTH, chr=NULL, simParam=SP) - 1
    geno.YLD <- pullSegSiteGeno(DH.YLD, chr=NULL, simParam=SP) - 1
    
    #expand the DH genotypes to include markers that were previously thinned to prevent AlphaSim from crashing.
    geno.DTH <- do.call(cbind, sapply(1:length(marker.rep), FUN=function(x) replicate(marker.rep[x], geno.DTH[,x])))
    geno.YLD <- do.call(cbind, sapply(1:length(marker.rep), FUN=function(x) replicate(marker.rep[x], geno.YLD[,x])))

    #calculate the % favorable P/S/A alleles and EBV of the DH
    temp.DTH <- calc.fav(geno=geno.DTH, m.eff=RR.coef.DTH.4[[i]], wt=j, pop=i, top=T)
    temp.YLD <- calc.fav(geno=geno.YLD, m.eff=RR.coef.YLD.4[[i]], wt=j, pop=i, top=T)
  
    #combine the results for all families/selection weights.
    out.DTH <- if(!exists("out.DTH")) temp.DTH else rbind(out.DTH, temp.DTH)
    out.YLD <- if(!exists("out.YLD")) temp.YLD else rbind(out.YLD, temp.YLD)
  }
}

#loop to iterate through 25 maize NAM families.
set.seed(43234)
for(i in 1:25){
  #create maize NAM base population in AlphaSim.
  founder.M <- newMapPop2(genMap=geno.founder.M[[i]][[2]],
                          haplotypes=geno.founder.M[[i]][[1]],
                          inbred=T,
                          ploidy=2L)
  SP <- SimParam$new(founder.M)
  founder.M <- newPop(founder.M, simParam=SP)
  base.M <- makeCross(pop=founder.M,
                      crossPlan=matrix(1:founder.M@nInd, nrow=founder.M@nInd/2, ncol=2, byrow=F),
                      simParam=SP)
  marker.rep <- unlist(geno.founder.M[[i]][[4]])
  
  #loop to iterate through different types of selection weights in OSGS (-0.1 is GS).
  for(j in seq(-0.1,1,0.1)){
    #identify top 4 lines.
    sel.line.DTS <- pick.line(ebv=RR.pred.DTS.4[[i]],
                              n=4,
                              wt=j,
                              top=F)
    sel.line.CL <- pick.line(ebv=RR.pred.CL.4[[i]],
                             n=4,
                             wt=j,
                             top=T)
    
    #cross the top 4 lines using all possible combinations.
    F1.DTS <- makeCross(pop=base.M,
                        crossPlan=t(combn(sel.line.DTS,2)),
                        simParam=SP)
    F1.CL <- makeCross(pop=base.M,
                       crossPlan=t(combn(sel.line.CL,2)),
                       simParam=SP)
    
    #create 10 double-haploids from each F1 (a total of 60 DHs).
    DH.DTS <- makeDH(pop=F1.DTS,
                     nDH=10,
                     simParam=SP)
    DH.CL <- makeDH(pop=F1.CL,
                    nDH=10,
                    simParam=SP)
    
    #extract the DH genotypes.
    geno.DTS <- pullSegSiteGeno(DH.DTS, chr=NULL, simParam=SP) - 1
    geno.CL <- pullSegSiteGeno(DH.CL, chr=NULL, simParam=SP) - 1
    
    #expand the DH genotypes to include markers that were previously thinned to prevent AlphaSim from crashing.
    geno.DTS <- do.call(cbind, sapply(1:length(marker.rep), FUN=function(x) replicate(marker.rep[x], geno.DTS[,x])))
    geno.CL <- do.call(cbind, sapply(1:length(marker.rep), FUN=function(x) replicate(marker.rep[x], geno.CL[,x])))
    
    #calculate the % favorable P/S/A alleles and EBV of the DH
    temp.DTS <- calc.fav(geno=geno.DTS, m.eff=RR.coef.DTS.4[[i]], wt=j, pop=i, top=F)
    temp.CL <- calc.fav(geno=geno.CL, m.eff=RR.coef.CL.4[[i]], wt=j, pop=i, top=T)
    
    #combine the results for all families/selection weights.
    out.DTS <- if(!exists("out.DTS")) temp.DTS else rbind(out.DTS, temp.DTS)
    out.CL <- if(!exists("out.CL")) temp.CL else rbind(out.CL, temp.CL)
  }
}

#calculate the %favorable P/S/A alleles and EBV of the base population.
base.DTH <- lapply(1:25, FUN=function(x) calc.fav(geno=geno.split.B[[x]], m.eff=RR.coef.DTH.4[[x]], wt=9, pop=x, top=T))
base.YLD <- lapply(1:25, FUN=function(x) calc.fav(geno=geno.split.B[[x]], m.eff=RR.coef.YLD.4[[x]], wt=9, pop=x, top=T))
base.DTS <- lapply(1:25, FUN=function(x) calc.fav(geno=geno.split.M[[x]], m.eff=RR.coef.DTS.4[[x]], wt=9, pop=x, top=F))
base.CL <- lapply(1:25, FUN=function(x) calc.fav(geno=geno.split.M[[x]], m.eff=RR.coef.CL.4[[x]], wt=9, pop=x, top=T))

#combine all data.
base.DTH <- do.call(rbind, base.DTH)
base.YLD <- do.call(rbind, base.YLD)
base.DTS <- do.call(rbind, base.DTS)
base.CL <- do.call(rbind, base.CL)

out.DTH <- rbind(out.DTH, base.DTH)
out.YLD <- rbind(out.YLD, base.YLD)
out.DTS <- rbind(out.DTS, base.DTS)
out.CL <- rbind(out.CL, base.CL)

#standardize A/P/S and EBV by the base population and calculate the difference between OSGS and GS.
for(i in c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1")){
  temp.DTH <- data.frame(favP=(out.DTH[out.DTH$wt==i,1]-out.DTH[out.DTH$wt=="GS",1])/out.DTH[out.DTH$wt=="B",1],
                         favS=(out.DTH[out.DTH$wt==i,2]-out.DTH[out.DTH$wt=="GS",2])/out.DTH[out.DTH$wt=="B",2],
                         favA=(out.DTH[out.DTH$wt==i,3]-out.DTH[out.DTH$wt=="GS",3])/out.DTH[out.DTH$wt=="B",3],
                         ebv=(out.DTH[out.DTH$wt==i,4]-out.DTH[out.DTH$wt=="GS",4])/out.DTH[out.DTH$wt=="B",5],
                         wt=i)
  temp.YLD <- data.frame(favP=(out.YLD[out.YLD$wt==i,1]-out.YLD[out.YLD$wt=="GS",1])/out.YLD[out.YLD$wt=="B",1],
                         favS=(out.YLD[out.YLD$wt==i,2]-out.YLD[out.YLD$wt=="GS",2])/out.YLD[out.YLD$wt=="B",2],
                         favA=(out.YLD[out.YLD$wt==i,3]-out.YLD[out.YLD$wt=="GS",3])/out.YLD[out.YLD$wt=="B",3],
                         ebv=(out.YLD[out.YLD$wt==i,4]-out.YLD[out.YLD$wt=="GS",4])/out.YLD[out.YLD$wt=="B",5],
                         wt=i)
  temp.DTS <- data.frame(favP=(out.DTS[out.DTS$wt==i,1]-out.DTS[out.DTS$wt=="GS",1])/out.DTS[out.DTS$wt=="B",1],
                         favS=(out.DTS[out.DTS$wt==i,2]-out.DTS[out.DTS$wt=="GS",2])/out.DTS[out.DTS$wt=="B",2],
                         favA=(out.DTS[out.DTS$wt==i,3]-out.DTS[out.DTS$wt=="GS",3])/out.DTS[out.DTS$wt=="B",3],
                         ebv=(out.DTS[out.DTS$wt==i,4]-out.DTS[out.DTS$wt=="GS",4])/out.DTS[out.DTS$wt=="B",5],
                         wt=i)
  temp.CL <- data.frame(favP=(out.CL[out.CL$wt==i,1]-out.CL[out.CL$wt=="GS",1])/out.CL[out.CL$wt=="B",1],
                        favS=(out.CL[out.CL$wt==i,2]-out.CL[out.CL$wt=="GS",2])/out.CL[out.CL$wt=="B",2],
                        favA=(out.CL[out.CL$wt==i,3]-out.CL[out.CL$wt=="GS",3])/out.CL[out.CL$wt=="B",3],
                        ebv=(out.CL[out.CL$wt==i,4]-out.CL[out.CL$wt=="GS",4])/out.CL[out.CL$wt=="B",5],
                        wt=i)
  plot.DTH <- if(!exists("plot.DTH")) temp.DTH else rbind(plot.DTH, temp.DTH)
  plot.YLD <- if(!exists("plot.YLD")) temp.YLD else rbind(plot.YLD, temp.YLD)
  plot.DTS <- if(!exists("plot.DTS")) temp.DTS else rbind(plot.DTS, temp.DTS)
  plot.CL <- if(!exists("plot.CL")) temp.CL else rbind(plot.CL, temp.CL)
}

#prepare the data for plotting.
plot.DTH <- melt(plot.DTH, id.vars="wt")
plot.YLD <- melt(plot.YLD, id.vars="wt")
plot.DTS <- melt(plot.DTS, id.vars="wt")
plot.CL <- melt(plot.CL, id.vars="wt")

plot.DTH$wt <- as.character(plot.DTH$wt)
plot.YLD$wt <- as.character(plot.YLD$wt)
plot.DTS$wt <- as.character(plot.DTS$wt)
plot.CL$wt <- as.character(plot.CL$wt)

plot.DTH$variable <- as.character(plot.DTH$variable)
plot.YLD$variable <- as.character(plot.YLD$variable)
plot.DTS$variable <- as.character(plot.DTS$variable)
plot.CL$variable <- as.character(plot.CL$variable)

plot.DTH$fillSig <- "NS"
plot.YLD$fillSig <- "NS"
plot.DTS$fillSig <- "NS"
plot.CL$fillSig <- "NS"

#identify the outlier and set them to NA (otherwise the density plots are hard to visualize).
for(i in c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1")){
  for(j in c("favP","favS","favA","ebv")){
  
    temp <- plot.DTH$wt==i & plot.DTH$variable==j
    outlier <- boxplot(plot.DTH[temp,3], plot=FALSE)$out
    plot.DTH[temp,3][plot.DTH[temp,3]%in%outlier] <- NA
    
    temp <- plot.YLD$wt==i & plot.YLD$variable==j
    outlier <- boxplot(plot.YLD[temp,3], plot=FALSE)$out
    plot.YLD[temp,3][plot.YLD[temp,3]%in%outlier] <- NA
    
    temp <- plot.DTS$wt==i & plot.DTS$variable==j
    outlier <- boxplot(plot.DTS[temp,3], plot=FALSE)$out
    plot.DTS[temp,3][plot.DTS[temp,3]%in%outlier] <- NA
    
    temp <- plot.CL$wt==i & plot.CL$variable==j
    outlier <- boxplot(plot.CL[temp,3], plot=FALSE)$out
    plot.CL[temp,3][plot.CL[temp,3]%in%outlier] <- NA
  }
}

#summarize the outliers.
outlier <- list(DTH=data.frame(wt=seq(0,1,0.1)[(which(is.na(plot.DTH[,3]))-1)%%275%/%25 + 1],
                               variable=c("favP","favS","favA","ebv")[(which(is.na(plot.DTH[,3]))-1)%/%275+1],
                               pop=(which(is.na(plot.DTH[,3]))-1)%%100%%25 + 1),
                YLD=data.frame(wt=seq(0,1,0.1)[(which(is.na(plot.YLD[,3]))-1)%%275%/%25 + 1],
                               variable=c("favP","favS","favA","ebv")[(which(is.na(plot.YLD[,3]))-1)%/%275+1],
                               pop=(which(is.na(plot.YLD[,3]))-1)%%100%%25 + 1),
                DTS=data.frame(wt=seq(0,1,0.1)[(which(is.na(plot.DTS[,3]))-1)%%275%/%25 + 1],
                               variable=c("favP","favS","favA","ebv")[(which(is.na(plot.DTS[,3]))-1)%/%275+1],
                               pop=(which(is.na(plot.DTS[,3]))-1)%%100%%25 + 1),
                CL=data.frame(wt=seq(0,1,0.1)[(which(is.na(plot.CL[,3]))-1)%%275%/%25 + 1],
                              variable=c("favP","favS","favA","ebv")[(which(is.na(plot.CL[,3]))-1)%/%275+1],
                              pop=(which(is.na(plot.CL[,3]))-1)%%100%%25+ 1))

#remove rows that have been previously identified as outliers.
plot.DTH <- plot.DTH[!is.na(plot.DTH[,3]),]
plot.YLD <- plot.YLD[!is.na(plot.YLD[,3]),]
plot.DTS <- plot.DTS[!is.na(plot.DTS[,3]),]
plot.CL <- plot.CL[!is.na(plot.CL[,3]),]

#run t.test with Bonferonni correction and identify which is significantly different from 0.
for(i in c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1")){
  for(j in c("favP","favS","favA","ebv")){
    if(t.test(plot.DTH[plot.DTH$wt==i & plot.DTH$variable==j,3])$p.value < 0.05/11) plot.DTH[plot.DTH$wt==i & plot.DTH$variable==j,4] <- "S"
    if(t.test(plot.YLD[plot.YLD$wt==i & plot.YLD$variable==j,3])$p.value < 0.05/11) plot.YLD[plot.YLD$wt==i & plot.YLD$variable==j,4] <- "S"
    if(t.test(plot.DTS[plot.DTS$wt==i & plot.DTS$variable==j,3])$p.value < 0.05/11) plot.DTS[plot.DTS$wt==i & plot.DTS$variable==j,4] <- "S"
    if(t.test(plot.CL[plot.CL$wt==i & plot.CL$variable==j,3])$p.value < 0.05/11) plot.CL[plot.CL$wt==i & plot.CL$variable==j,4] <- "S"
  }
}

#create dummy datasets so that we can plot the Y-axis on a mirror scale (e.g. -0.5 to 0.5).
dummy.DTH <- data.frame(wt=NA,
                        variable=rep(c("favP","favS"),2),
                        value=c(-max(abs(plot.DTH[plot.DTH[,2]=="favP",3])),
                                -max(abs(plot.DTH[plot.DTH[,2]=="favS",3])),
                                max(abs(plot.DTH[plot.DTH[,2]=="favP",3])),
                                max(abs(plot.DTH[plot.DTH[,2]=="favS",3]))),
                        fillSig=NA)

dummy.YLD <- data.frame(wt=NA,
                        variable=rep(c("favP","favS"),2),
                        value=c(-max(abs(plot.YLD[plot.YLD[,2]=="favP",3])),
                                -max(abs(plot.YLD[plot.YLD[,2]=="favS",3])),
                                max(abs(plot.YLD[plot.YLD[,2]=="favP",3])),
                                max(abs(plot.YLD[plot.YLD[,2]=="favS",3]))),
                        fillSig=NA)

dummy.DTS <- data.frame(wt=NA,
                        variable=rep(c("favP","favS"),2),
                        value=c(-max(abs(plot.DTS[plot.DTS[,2]=="favP",3])),
                                -max(abs(plot.DTS[plot.DTS[,2]=="favS",3])),
                                max(abs(plot.DTS[plot.DTS[,2]=="favP",3])),
                                max(abs(plot.DTS[plot.DTS[,2]=="favS",3]))),
                        fillSig=NA)

dummy.CL <- data.frame(wt=NA,
                       variable=rep(c("favP","favS"),2),
                       value=c(-max(abs(plot.CL[plot.CL[,2]=="favP",3])),
                               -max(abs(plot.CL[plot.CL[,2]=="favS",3])),
                               max(abs(plot.CL[plot.CL[,2]=="favP",3])),
                               max(abs(plot.CL[plot.CL[,2]=="favS",3]))),
                       fillSig=NA)

#prepare the data for plotting, again.
plot.DTH$wt <- as.factor(plot.DTH$wt)
plot.DTH$wt <- factor(plot.DTH$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
plot.DTH <- plot.DTH[plot.DTH$variable%in%c("ebv","favP","favS"),]
plot.DTH$variable <- as.factor(plot.DTH$variable)
plot.DTH$variable <- factor(plot.DTH$variable,
                            levels=c("ebv","favP","favS"),
                            labels=c("ebv","favP","favS"))
plot.DTH$fillSig <- as.factor(plot.DTH$fillSig)
plot.DTH$fillSig <- factor(plot.DTH$fillSig,
                           levels=c("NS","S"),
                           labels=c("NS","S"))

plot.YLD$wt <- as.factor(plot.YLD$wt)
plot.YLD$wt <- factor(plot.YLD$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
plot.YLD <- plot.YLD[plot.YLD$variable%in%c("ebv","favP","favS"),]
plot.YLD$variable <- as.factor(plot.YLD$variable)
plot.YLD$variable <- factor(plot.YLD$variable,
                            levels=c("ebv","favP","favS"),
                            labels=c("ebv","favP","favS"))
plot.YLD$fillSig <- as.factor(plot.YLD$fillSig)
plot.YLD$fillSig <- factor(plot.YLD$fillSig,
                           levels=c("NS","S"),
                           labels=c("NS","S"))

plot.DTS$wt <- as.factor(plot.DTS$wt)
plot.DTS$wt <- factor(plot.DTS$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
plot.DTS <- plot.DTS[plot.DTS$variable%in%c("ebv","favP","favS"),]
plot.DTS$variable <- as.factor(plot.DTS$variable)
plot.DTS$variable <- factor(plot.DTS$variable,
                            levels=c("ebv","favP","favS"),
                            labels=c("ebv","favP","favS"))
plot.DTS$fillSig <- as.factor(plot.DTS$fillSig)
plot.DTS$fillSig <- factor(plot.DTS$fillSig,
                           levels=c("NS","S"),
                           labels=c("NS","S"))

plot.CL$wt <- as.factor(plot.CL$wt)
plot.CL$wt <- factor(plot.CL$wt,
                     levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                     labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
plot.CL <- plot.CL[plot.CL$variable%in%c("ebv","favP","favS"),]
plot.CL$variable <- as.factor(plot.CL$variable)
plot.CL$variable <- factor(plot.CL$variable,
                           levels=c("ebv","favP","favS"),
                           labels=c("ebv","favP","favS"))
plot.CL$fillSig <- as.factor(plot.CL$fillSig)
plot.CL$fillSig <- factor(plot.CL$fillSig,
                          levels=c("NS","S"),
                          labels=c("NS","S"))

#create the plots.
plot.DTH <- ggplot() +
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#000000", linetype=3) +
  geom_violin(data=plot.DTH, aes(x=wt, y=value, fill=fillSig, color=fillSig)) +
  facet_wrap(vars(variable), ncol=1, scales="free_y", strip.position="left", labeller=labeller(variable=c(favP="%P+",favS="%S+",favA="%A+",ebv="EBV"))) +
  theme(strip.placement.y="outside", strip.background=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  geom_point(data=dummy.DTH, aes(x=wt, y=value), na.rm=T) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

plot.YLD <- ggplot() +
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#000000", linetype=3) +
  geom_violin(data=plot.YLD, aes(x=wt, y=value, fill=fillSig, color=fillSig)) +
  facet_wrap(vars(variable), ncol=1, scales="free_y", strip.position="left", labeller=labeller(variable=c(favP="%P+",favS="%S+",favA="%A+",ebv="EBV"))) +
  theme(strip.placement.y="outside", strip.background=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  geom_point(data=dummy.YLD, aes(x=wt, y=value), na.rm=T) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

plot.DTS <- ggplot() +
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#000000", linetype=3) +
  geom_violin(data=plot.DTS, aes(x=wt, y=value, fill=fillSig, color=fillSig)) +
  facet_wrap(vars(variable), ncol=1, scales="free_y", strip.position="left", labeller=labeller(variable=c(favP="%P+",favS="%S+",favA="%A+",ebv="EBV"))) +
  theme(strip.placement.y="outside", strip.background=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  geom_point(data=dummy.DTS, aes(x=wt, y=value), na.rm=T) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

plot.CL <- ggplot() +
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#000000", linetype=3) +
  geom_violin(data=plot.CL, aes(x=wt, y=value, fill=fillSig, color=fillSig)) +
  facet_wrap(vars(variable), ncol=1, scales="free_y", strip.position="left", labeller=labeller(variable=c(favP="%P+",favS="%S+",favA="%A+",ebv="EBV"))) +
  theme(strip.placement.y="outside", strip.background=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  geom_point(data=dummy.CL, aes(x=wt, y=value), na.rm=T) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))


#export the plots
ggsave(filename="selection_DTH.png",
       plot=plot.DTH,
       scale=2,
       width=1.5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="selection_YLD.png",
       plot=plot.YLD,
       scale=2,
       width=1.5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="selection_DTS.png",
       plot=plot.DTS,
       scale=2,
       width=1.5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="selection_CL.png",
       plot=plot.CL,
       scale=2,
       width=1.5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="selection_DTH.svg",
       plot=plot.DTH,
       scale=1.5,
       width=1.5,
       height=3,
       units="in")

ggsave(filename="selection_YLD.svg",
       plot=plot.YLD,
       scale=1.5,
       width=1.5,
       height=3,
       units="in")

ggsave(filename="selection_DTS.svg",
       plot=plot.DTS,
       scale=1.5,
       width=1.5,
       height=3,
       units="in")

ggsave(filename="selection_CL.svg",
       plot=plot.CL,
       scale=1.5,
       width=1.5,
       height=3,
       units="in")


#plot the marker effects along chr and position.
plot.coef <- data.frame(marker.B,
                        effect=RR.coef.YLD.4[[1]])

plot.coef <- ggplot() + 
  geom_point(data=plot.coef, aes(x=PosGen, y=effect, color=as.character(sign(effect))), size=0.1) +
  facet_wrap(vars(Chr), nrow=1, scales="free_x") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(legend.position="bottom") +
  scale_color_manual(name="Favorable allele", values=c("#FF0000", "#00FF00"), labels=c("P", "S")) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  xlab("Genetic Position (cM)") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000")

ggsave(filename="coef_YLD_1.svg",
       plot=plot.coef,
       scale=1.3,
       width=4,
       height=2,
       units="in")


save.image("OSGS_dat3.RData")

