setwd("C:/Users/cyang/Desktop/OSGS/real_data")

#read in the marker information.
marker <- read.csv("Barley_NAM_marker_info.csv", as.is=T)
marker <- split(marker[,1], marker[,3])

#read in the marker genotype.
geno <- read.csv("Barley_NAM_marker_genotype_preimputed.csv", row.names=1)
geno <- as.matrix(geno)

#separate the marker genotype by chromosomes
geno <- lapply(1:length(marker), FUN=function(x) geno[,colnames(geno)%in%marker[[x]]])

#impute the first and last missing sites with the same genotype as nearest marker.
for(i in 1:length(marker)){
  for(j in 1:nrow(geno[[i]])){
    if(is.na(geno[[i]][j,1])){
      m <- 1; while(is.na(geno[[i]][j,m])){m <- m+1}
      geno[[i]][j,1] <- geno[[i]][j,m]
    }
    if(is.na(geno[[i]][j,ncol(geno[[i]])])){
      m <- ncol(geno[[i]]); while(is.na(geno[[i]][j,m])){m <- m-1}
      geno[[i]][j,ncol(geno[[i]])] <- geno[[i]][j,m]
    }
  }
}

#impute the missing sites with the same algorithm as maize NAM.
#0/NA/0 -> 0/0/0
#2/NA/2 -> 2/2/2
#0/NA/2 -> 0/1/2
#2/NA/0 -> 2/1/0
#0/NA/1 -> 0/0.5/1
#1/NA/0 -> 1/0.5/0
#1/NA/2 -> 1/1.5/2
#2/NA/1 -> 2/1.5/1
for(i in 1:length(marker)){
  for(j in 1:nrow(geno[[i]])){
    for(k in 1:ncol(geno[[i]])){
      if(is.na(geno[[i]][j,k])){
        m <- k - 1
        n <- k; while(is.na(geno[[i]][j,n])){n <- n+1}
        if(geno[[i]][j,m]==geno[[i]][j,n]){
          geno[[i]][j,k:(n-1)] <- geno[[i]][j,m]
        } else if(geno[[i]][j,m]==0 & geno[[i]][j,n]==2){
          geno[[i]][j,k:(n-1)] <- 1
        } else if(geno[[i]][j,m]==2 & geno[[i]][j,n]==0){
          geno[[i]][j,k:(n-1)] <- 1
        } else if(geno[[i]][j,m]==0 & geno[[i]][j,n]==1){
          geno[[i]][j,k:(n-1)] <- 0.5
        } else if(geno[[i]][j,m]==1 & geno[[i]][j,n]==0){
          geno[[i]][j,k:(n-1)] <- 0.5
        } else if(geno[[i]][j,m]==1 & geno[[i]][j,n]==2){
          geno[[i]][j,k:(n-1)] <- 1.5
        } else if(geno[[i]][j,m]==2 & geno[[i]][j,n]==1){
          geno[[i]][j,k:(n-1)] <- 1.5
        }        
      }
    }
  }
}

#combine the marker genotype for all chromosomes.
geno <- do.call(cbind, geno)

#convert the marker genotype into data frame.
geno <- data.frame(Line=rownames(geno), geno)
geno$Line <- as.character(geno$Line)
rownames(geno) <- NULL

#export the marker genotype.
write.csv(geno, "Barley_NAM_marker_genotype.csv", row.names=F, quote=F)







