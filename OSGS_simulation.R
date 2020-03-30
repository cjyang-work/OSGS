library(AlphaSimR)
library(rrBLUP)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(stringr)
setwd("C:/Users/cjyan/OneDrive/Desktop/OSGS/output_sim")

sessionInfo()
#R version 3.6.3 (2020-02-29)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 17763)
#
#Matrix products: default
#
#locale:
#[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
#[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
#[5] LC_TIME=English_United States.1252    
#
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
#[1] gdtools_0.2.1    gridExtra_2.3    ggplot2_3.2.1    reshape2_1.4.3   rrBLUP_4.6.1    
#[6] AlphaSimR_0.11.1 R6_2.4.1        
#
#loaded via a namespace (and not attached):
#[1] Rcpp_1.0.3        rstudioapi_0.11   magrittr_1.5      munsell_0.5.0     colorspace_1.4-1 
#[6] rlang_0.4.4       stringr_1.4.0     plyr_1.8.5        tools_3.6.3       grid_3.6.3       
#[11] gtable_0.3.0      withr_2.1.2       systemfonts_0.1.1 digest_0.6.24     lazyeval_0.2.2   
#[16] tibble_2.1.3      lifecycle_0.1.0   crayon_1.3.4      farver_2.0.3      labeling_0.3     
#[21] stringi_1.4.6     compiler_3.6.3    pillar_1.4.3      scales_1.1.0      svglite_1.2.3    
#[26] pkgconfig_2.0.3

#simulate base population.
OSGS.SimPop <- function(n.chr=10, #number of chromosomes
                        n.marker=seq(1000,550,-50), #vector of number of marker per chromosome
                        gen.len=seq(2,1.1,-0.1), #0.2cM/marker
                        pop.type="F2", #accepted value = "F2","BC1B","RBC1B"
                        n.line=200, #number of lines to simulate
                        n.self=4){ #number of generations to self

    #calculate initial values.
  gen.map <- lapply(1:n.chr, FUN=function(x) seq(0, gen.len[x]-gen.len[x]/n.marker[x], gen.len[x]/n.marker[x]))
  
  #create 2 founders.
  founder <- newMapPop2(genMap=gen.map,
                        haplotypes=lapply(1:n.chr, FUN=function(x) matrix(rep(c(0,1), n.marker[x]), nrow=2, ncol=n.marker[x])),
                        inbred=T,
                        ploidy=2L)
  SP <- SimParam$new(founder)
  founder <- newPop(founder, simParam=SP)
  
  #make F1.
  F1 <- makeCross(pop=founder,
                  crossPlan=matrix(c(1,2), nrow=1, ncol=2),
                  simParam=SP)
  
  #make F2/BC1B/RBC1B.
  if(pop.type=="F2"){
    G2 <- self(pop=F1,
               nProgeny=n.line,
               simParam=SP)
  } else if(pop.type=="BC1"){
    G2 <- makeCross2(females=founder,
                     males=F1,
                     crossPlan=matrix(rep(c(1,1), n.line), nrow=n.line, ncol=2, byrow=T),
                     simParam=SP)
  } else if(pop.type=="RBC1"){
    G2 <- makeCross2(females=founder,
                     males=F1,
                     crossPlan=matrix(rep(c(2,1), n.line), nrow=n.line, ncol=2, byrow=T),
                     simParam=SP)
  }
  
  #self F2/BC1B/RBC1B for n.self generations.
  temp <- G2
  for(i in 1:n.self){
    temp <- self(pop=temp,
                 nProgeny=1,
                 simParam=SP)
    assign(paste("G",i+2,sep=""),temp)
  }
  
  return(list(temp,
              SP))
}

#simulate QTL markers and effects.
OSGS.SimQTL <- function(n.chr=10, #number of chromosomes simulated
                        n.marker=seq(1000,550,-50), #vector of number of markers per chromosome
                        p.qtl=0.01, #proportion of markers as QTL
                        p.P=c(0.5,0.55,0.6,0.7,0.8,0.9), #vector of proportions of P-QTL
                        meanG=0, #mean of genetic effect
                        varG=1){ #variance of genetic effect
                        
  #simulate QTL marker positions based on p.qtl.
  qtl.pos <- lapply(1:n.chr, FUN=function(x) sort(sample(1:n.marker[x], floor(p.qtl*n.marker[x]), replace=F)))
  
  #identify total number of QTL and favorable-P/S.
  n.qtl <- length(unlist(qtl.pos))
  n.P <- floor(p.P*n.qtl)
  n.S <- n.qtl - n.P
  
  #simulate QTL marker effects such that var(P)=var(S) and var(A)=1/n.qtl.
  qtl.eff <- replicate(length(p.P),list(replicate(n.chr,list())))
  
  for(i in 1:length(p.P)){
    temp <- c(-abs(rnorm(n.P[i], meanG, sqrt(varG/(n.P[i]+n.S[i]-2*(n.S[i]-n.P[i])^2/pi/(n.P[i]+n.S[i]))))),
              abs(rnorm(n.S[i], meanG, sqrt(varG/(n.P[i]+n.S[i]-2*(n.S[i]-n.P[i])^2/pi/(n.P[i]+n.S[i]))))))
    for(j in 1:n.chr){
      temp.qtl <- sample(1:length(temp), length(qtl.pos[[j]]), replace=F)
      qtl.eff[[i]][[j]] <- temp[temp.qtl]
      temp <- temp[-temp.qtl]
    }
  }
  
  #summarize the marker effects.
  m.eff <- lapply(1:n.chr, FUN=function(x) rep(0, n.marker[x]))
  m.eff <- replicate(length(p.P), list(m.eff))
  
  for(i in 1:length(p.P)){
    for(j in 1:n.chr){
      m.eff[[i]][[j]][qtl.pos[[j]]] <- qtl.eff[[i]][[j]]
    }
    m.eff[[i]] <- unlist(m.eff[[i]])
  }
  
  return(m.eff)
}

#simulate selection based on OSGS/GS.  
OSGS.SimSel <- function(S0, #base pop
                        m.eff, #list of simulated marker effects
                        n.line=200, #number of lines
                        meanE=0, #mean of residual effect
                        varE=1, #variance of residual effect
                        wt=seq(0,1,0.1), #vector of selection weights on P
                        n.sel=5, #number of lines to be selected
                        n.DH=20, #number of double-haploids per cross
                        p.P=c(0.5,0.55,0.6,0.7,0.8,0.9), #vector of proportion P-QTLs previously simulated
                        RS.gen=5, #number of generations of recurrent selection
                        SP=SP){
  
  #extract marker data for generation S0 (base population).
  geno.S0 <- pullSegSiteGeno(S0, chr=NULL, simParam=SP) - 1
  
  #calculate the true breeding values (BVs) and trait value for S0.
  tbv.S0 <- lapply(1:length(m.eff), FUN=function(x) c(geno.S0%*%m.eff[[x]]))
  pheno.S0 <- lapply(1:length(m.eff), FUN=function(x) tbv.S0[[x]]+rnorm(n.line,meanE,sqrt(varE)))
  
  #predict marker effects in S0.
  GP.S0 <- lapply(1:length(m.eff), FUN=function(x) mixed.solve(y=pheno.S0[[x]], Z=geno.S0))
  u.S0 <- lapply(1:length(m.eff), FUN=function(x) c(GP.S0[[x]]$u)) #marker effect
  beta.S0 <- lapply(1:length(m.eff), FUN=function(x) c(GP.S0[[x]]$beta)) #mean
  
  #isolate P and S favorable marker effects.
  u.P.S0 <- u.S.S0 <- u.S0
  for(i in 1:length(m.eff)){
    u.P.S0[[i]][u.P.S0[[i]]>0] <- 0
    u.S.S0[[i]][u.S.S0[[i]]<0] <- 0
  }
  
  #calculate the % favorable P/S/A alleles and EBV in generation S0.
  out.S0 <- vector()
  for(i in 1:length(m.eff)){
    out.S0 <- rbind(out.S0,
                    calc.fav(geno=geno.S0,
                             m.eff=m.eff[[i]],
                             wt=9,
                             pop=i,
                             top=T))
  }

  #predict BVs in S0 based on A/P/S alleles.
  ebvA.S0 <- lapply(1:length(m.eff), FUN=function(x) c(geno.S0%*%u.S0[[x]]))
  ebvP.S0 <- lapply(1:length(m.eff), FUN=function(x) c(geno.S0%*%u.P.S0[[x]]))
  ebvS.S0 <- lapply(1:length(m.eff), FUN=function(x) c(geno.S0%*%u.S.S0[[x]]))
  
  #identify top lines based on GS and wt.
  sel.S0 <- replicate(length(m.eff), list(replicate(length(wt)+1, vector())))
  for(i in 1:length(m.eff)){
    for(j in 1:length(wt)){
      temp <- rank(ebvP.S0[[i]], ties.method="random")*wt[j] + rank(ebvS.S0[[i]], ties.method="random")*(1-wt[j])
      names(temp) <- 1:length(temp)
      temp <- sort(-temp)
      sel.S0[[i]][[j]] <- sort(as.numeric(names(temp)[1:n.sel]))
    }
    temp <- rank(ebvA.S0[[i]], ties.method="random")
    names(temp) <- 1:length(temp)
    temp <- sort(-temp)
    sel.S0[[i]][[j+1]] <- sort(as.numeric(names(temp)[1:n.sel]))
  }
  
  #make crosses among the selected S0 lines and create DH S1 lines.
  S1 <- replicate(length(m.eff), list(replicate(length(wt)+1, vector())))
  for(i in 1:length(m.eff)){
    for(j in 1:(length(wt)+1)){
      S1[[i]][[j]] <- makeCross(pop=S0,
                                crossPlan=t(combn(sel.S0[[i]][[j]],2)),
                                simParam=SP)
      S1[[i]][[j]] <- makeDH(pop=S1[[i]][[j]],
                             nDH=n.DH,
                             simParam=SP)
    }
  }
      
  #extract marker data for generation S1.
  geno.S1 <- replicate(length(m.eff), list(replicate(length(wt)+1, vector())))
  for(i in 1:length(m.eff)){
    for(j in 1:(length(wt)+1)){
      geno.S1[[i]][[j]] <- pullSegSiteGeno(S1[[i]][[j]], chr=NULL, simParam=SP) - 1
    }
  }
  
  #calculate the % favorable P/S/A alleles and EBV in generation S1.
  out.S1 <- vector()
  for(i in 1:length(m.eff)){
    for(j in 1:(length(wt)+1)){
      out.S1 <- rbind(out.S1, 
                      calc.fav(geno=geno.S1[[i]][[j]],
                               m.eff=m.eff[[i]],
                               wt=c(wt,-0.1)[j],
                               pop=i,
                               top=T))
    }
  }

  #calculate some initial values (S1) for the recurrent selection loop.
  temp.ebvA <- replicate(length(m.eff), vector())
  temp.ebvP <- temp.ebvS <- replicate(length(m.eff), list(replicate(length(wt), vector())))
  for(i in 1:length(m.eff)){
    for(j in 1:length(wt)){
      temp.ebvP[[i]][[j]] <- c(geno.S1[[i]][[j]]%*%u.P.S0[[i]])
      temp.ebvS[[i]][[j]] <- c(geno.S1[[i]][[j]]%*%u.S.S0[[i]])
    }
    temp.ebvA[[i]] <- c(geno.S1[[i]][[j+1]]%*%u.S0[[i]])
  }
  temp.S <- S1

  #compile all the outputs for recurrent and single-generation selection.
  out.RS <- out.SS <- rbind(data.frame(out.S0,
                                       generation="S0"),
                            data.frame(out.S1,
                                       generation="S1"))
  
  #loop through the remaining generations of recurrent selection.
  for(k in 2:RS.gen){
    
    #identify top lines.
    temp.sel <- replicate(length(m.eff), list(replicate(length(wt)+1, vector())))
    for(i in 1:length(m.eff)){
      for(j in 1:length(wt)){
        temp <- rank(temp.ebvP[[i]][[j]], ties.method="random")*wt[j] + rank(temp.ebvS[[i]][[j]], ties.method="random")*(1-wt[j])
        names(temp) <- 1:length(temp)
        temp.sel[[i]][[j]] <- sort(as.numeric(names(sort(-temp))[1:n.sel]))
      }
      temp <- rank(temp.ebvA[[i]], ties.method="random")
      names(temp) <- 1:length(temp)
      temp.sel[[i]][[j+1]] <- sort(as.numeric(names(sort(-temp))[1:n.sel]))
    }

    #make crosses among the selected lines and create DH.
    for(i in 1:length(m.eff)){
      for(j in 1:(length(wt)+1)){
        temp.S[[i]][[j]] <- makeCross(pop=temp.S[[i]][[j]],
                                      crossPlan=t(combn(temp.sel[[i]][[j]],2)),
                                      simParam=SP)
        temp.S[[i]][[j]] <- makeDH(pop=temp.S[[i]][[j]],
                                   nDH=n.DH,
                                   simParam=SP)
      }
    }
    
    #extract marker data.
    temp.geno <- replicate(length(m.eff), list(replicate(length(wt)+1, vector())))
    for(i in 1:length(m.eff)){
      for(j in 1:(length(wt)+1)){
        temp.geno[[i]][[j]] <- pullSegSiteGeno(temp.S[[i]][[j]], chr=NULL, simParam=SP) - 1
      }
    }
    
    #calculate the % favorable P/S/A alleles and true BV.
    temp.out <- vector()
    for(i in 1:length(m.eff)){
      for(j in 1:(length(wt)+1)){
        temp.out <- rbind(temp.out, 
                          calc.fav(geno=temp.geno[[i]][[j]],
                                   m.eff=m.eff[[i]],
                                   wt=c(wt,-0.1)[j],
                                   pop=i,
                                   top=T))
      }
    }
    
    out.RS <- rbind(out.RS,
                    data.frame(temp.out,
                               generation=paste("S",k,sep="")))
    
    temp.ebvA <- replicate(length(m.eff), vector())
    temp.ebvP <- temp.ebvS <- replicate(length(m.eff), list(replicate(length(wt), vector())))
    for(i in 1:length(m.eff)){
      for(j in 1:length(wt)){
        temp.ebvP[[i]][[j]] <- c(temp.geno[[i]][[j]]%*%u.P.S0[[i]])
        temp.ebvS[[i]][[j]] <- c(temp.geno[[i]][[j]]%*%u.S.S0[[i]])
      }
      temp.ebvA[[i]] <- c(temp.geno[[i]][[j+1]]%*%u.S0[[i]])
    }
  }

  return(list(out.SS,
              out.RS))
}

#summarize the selection outputs into data for plotting.
OSGS.SimCompare <- function(out.SimSel, RS){ #Here, we only use one list at a time.
  
  dat.base <- out.SimSel[out.SimSel$wt=="B",]
  dat.sel <- out.SimSel[!(out.SimSel$wt=="B"),]
  
  for(i in unique(out.SimSel$sim)){
    for(j in unique(out.SimSel$pop)){
      for(k in 1:3){
        dat.sel[dat.sel$sim==i & dat.sel$pop==j,k] <- dat.sel[dat.sel$sim==i & dat.sel$pop==j,k]/dat.base[dat.base$sim==i & dat.base$pop==j,k]
      }
      dat.sel[dat.sel$sim==i & dat.sel$pop==j,4] <- (dat.sel[dat.sel$sim==i & dat.sel$pop==j,4]-dat.base[dat.base$sim==i & dat.base$pop==j,4])/dat.base[dat.base$sim==i & dat.base$pop==j,5]
    }
  }
  
  dat.sel <- dat.sel[,-5]
  
  if(RS){
    dat.OSGS <- melt(dat.sel, id.vars=c("wt","pop","generation","sim"))
    dat.OSGS$Sig <- "NS"
    
    for(j in unique(dat.OSGS$pop)){
      for(k in unique(dat.OSGS$generation)){
        for(m in unique(dat.OSGS$variable)){
          temp <- unique(dat.OSGS$wt)
          temp <- t.test(dat.OSGS[dat.OSGS$wt==temp[1] & dat.OSGS$pop==j & dat.OSGS$generation==k & dat.OSGS$variable==m,"value"],
                         dat.OSGS[dat.OSGS$wt==temp[2] & dat.OSGS$pop==j & dat.OSGS$generation==k & dat.OSGS$variable==m,"value"], paired=T)$p.value < 0.05
          if(temp) dat.OSGS[dat.OSGS$pop==j & dat.OSGS$generation==k & dat.OSGS$variable==m,"Sig"] <- "S"
        }
      }
    }
    
  } else {
    
    dat.GS <- dat.sel[dat.sel$wt=="GS",]
    dat.OSGS <- dat.sel[!(dat.sel$wt=="GS"),]
    
    for(i in unique(out.SimSel$sim)){
      for(j in unique(out.SimSel$pop)){
        for(k in 1:4){
          dat.OSGS[dat.OSGS$sim==i & dat.OSGS$pop==j,k] <- dat.OSGS[dat.OSGS$sim==i & dat.OSGS$pop==j,k] - dat.GS[dat.GS$sim==i & dat.GS$pop==j,k]
        }
      }
    }
    
    #need to add significance
    dat.OSGS <- melt(dat.OSGS, id.vars=c("wt","pop","generation","sim"))
    dat.OSGS$Sig <- "NS"
    
    for(i in unique(dat.OSGS$wt)){
      for(j in unique(dat.OSGS$pop)){
        for(k in unique(dat.OSGS$generation)){
          for(m in unique(dat.OSGS$variable)){
            temp <- t.test(dat.OSGS[dat.OSGS$wt==i & dat.OSGS$pop==j & dat.OSGS$generation==k & dat.OSGS$variable==m,"value"])$p.value < 0.05/length(unique(dat.OSGS$wt))
            if(temp) dat.OSGS[dat.OSGS$wt==i & dat.OSGS$pop==j & dat.OSGS$generation==k & dat.OSGS$variable==m,"Sig"] <- "S"
          }
        }
      }
    }
    
  }
  
  return(dat.OSGS)
  
}

#calculate proportion P/S and ebv.
calc.fav <- function(geno,
                     m.eff,
                     wt,
                     pop,
                     top=T){
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

#summarize the estimated P/S, and h2.
OSGS.SimSummary <- function(out.SimPop=F2,
                            out.SimQTL=QTL,
                            meanE=0,
                            varE=1){

  #extract the marker genotype.
  geno <- pullSegSiteGeno(out.SimPop[[1]], chr=NULL, simParam=out.SimPop[[2]]) - 1
  
  #calculate the phenotypic trait value.
  pheno <- lapply(1:length(out.SimQTL), FUN=function(x) c(geno%*%out.SimQTL[[x]]) + rnorm(200,meanE,varE))
  
  #estimate marker effects.
  u <- lapply(1:length(out.SimQTL), FUN=function(x) c(mixed.solve(y=pheno[[x]], Z=geno)$u))
  
  #summarize the true and estimated P/S.
  PS <- lapply(1:length(out.SimQTL), FUN=function(x) c(sum(u[[x]]<0)/sum(u[[x]]!=0),
                                                       sum(u[[x]]>0)/sum(u[[x]]!=0)))  

  #summarize the h2.
  h2 <- lapply(1:length(out.SimQTL), FUN=function(x) var(c(geno%*%out.SimQTL[[x]]))/var(pheno[[x]]))
  
  return(list(PS, h2))
}

#replacement function by R Chris Gaynor (AlphaSim author) because the current version is bugged.
newMapPop2 = function(genMap,
                      haplotypes,
                      inbred=FALSE,
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


###############
### 2cM/QTL ###
###############
F2.out <- list(vector(), vector())
BC1.out <- list(vector(), vector())
RBC1.out <- list(vector(), vector())
QTL.out <- replicate(6, vector())
PS.out <- replicate(6, vector())
h2.out <- replicate(6, vector())

set.seed(31423)
for(i in 1:100){
  #simulate F2,BC1,RBC1 populations.
  F2 <- OSGS.SimPop(n.chr=10,
                    n.marker=seq(1000,550,-50),
                    gen.len=seq(2,1.1,-0.1),
                    pop.type="F2",
                    n.line=200,
                    n.self=4)
  BC1 <- OSGS.SimPop(n.chr=10,
                     n.marker=seq(1000,550,-50),
                     gen.len=seq(2,1.1,-0.1),
                     pop.type="BC1",
                     n.line=200,
                     n.self=4)
  RBC1 <- OSGS.SimPop(n.chr=10,
                      n.marker=seq(1000,550,-50),
                      gen.len=seq(2,1.1,-0.1),
                      pop.type="RBC1",
                      n.line=200,
                      n.self=4)
  
  #simulate QTL.
  QTL <- OSGS.SimQTL(n.chr=10,
                     n.marker=seq(1000,550,-50),
                     p.qtl=0.1,
                     p.P=c(0.5,0.55,0.6,0.7,0.8,0.9),
                     meanG=0,
                     varG=1)
  
  #simulate selection.
  F2.Sel <- OSGS.SimSel(S0=F2[[1]],
                        m.eff=QTL,
                        n.line=200,
                        meanE=0,
                        varE=1,
                        wt=seq(0,1,0.1),
                        n.sel=5,
                        n.DH=200/choose(5,2),
                        p.P=c(0.5,0.55,0.6,0.7,0.8,0.9),
                        RS.gen=5,
                        SP=F2[[2]])
  BC1.Sel <- OSGS.SimSel(S0=BC1[[1]],
                         m.eff=QTL,
                         n.line=200,
                         meanE=0,
                         varE=1,
                         wt=seq(0,1,0.1),
                         n.sel=5,
                         n.DH=200/choose(5,2),
                         p.P=c(0.5,0.55,0.6,0.7,0.8,0.9),
                         RS.gen=5,
                         SP=BC1[[2]])
  RBC1.Sel <- OSGS.SimSel(S0=RBC1[[1]],
                          m.eff=QTL,
                          n.line=200,
                          meanE=0,
                          varE=1,
                          wt=seq(0,1,0.1),
                          n.sel=5,
                          n.DH=200/choose(5,2),
                          p.P=c(0.5,0.55,0.6,0.7,0.8,0.9),
                          RS.gen=5,
                          SP=RBC1[[2]])
  
  #summarize output from single-gen selection (j=1) and recurrent selection (j=2).
  for(j in 1:2){
    F2.out[[j]] <- rbind(F2.out[[j]],
                         data.frame(F2.Sel[[j]],
                                    sim=i))
    BC1.out[[j]] <- rbind(BC1.out[[j]],
                          data.frame(BC1.Sel[[j]],
                                     sim=i))
    RBC1.out[[j]] <- rbind(RBC1.out[[j]],
                           data.frame(RBC1.Sel[[j]],
                                      sim=i))
  }
  
  #summarize the total QTL effect, i.e. P-parent true breeding value.
  for(j in 1:6){
    QTL.out[[j]] <- c(QTL.out[[j]], -sum(QTL[[j]]))
  }
  
  #summarize the estimated proportion of P/S and h2.
  temp.F2 <- OSGS.SimSummary(out.SimPop=F2, out.SimQTL=QTL, meanE=0, varE=1)
  temp.BC1 <- OSGS.SimSummary(out.SimPop=BC1, out.SimQTL=QTL, meanE=0, varE=1)
  temp.RBC1 <- OSGS.SimSummary(out.SimPop=RBC1, out.SimQTL=QTL, meanE=0, varE=1)
  
  for(j in 1:6){
    PS.out[[j]] <- rbind(PS.out[[j]],
                         c(temp.F2[[1]][[j]], temp.BC1[[1]][[j]], temp.RBC1[[1]][[j]]))
    h2.out[[j]] <- rbind(h2.out[[j]],
                         c(temp.F2[[2]][[j]], temp.BC1[[2]][[j]], temp.RBC1[[2]][[j]]))
  }

  message(i, " done.")
  
}

#save the selection output.    
save(F2.out, BC1.out, RBC1.out, file="sim_20200317.RData")

#prepare plot data.
##################

#identify the P-parent breeding value baseline.
  check.F2 <- check.BC1 <- check.RBC1 <- replicate(6,vector())
  for(i in 1:6){
    check.F2[[i]] <- c(mean(F2.out[[1]][F2.out[[1]][,6]=="B" & F2.out[[1]]$pop==i,4]),
                       mean(F2.out[[1]][F2.out[[1]][,6]=="B" & F2.out[[1]]$pop==i,5]))
    check.BC1[[i]] <- c(mean(BC1.out[[1]][BC1.out[[1]][,6]=="B" & BC1.out[[1]]$pop==i,4]),
                        mean(BC1.out[[1]][BC1.out[[1]][,6]=="B" & BC1.out[[1]]$pop==i,5]))
    check.RBC1[[i]] <- c(mean(RBC1.out[[1]][RBC1.out[[1]][,6]=="B" & RBC1.out[[1]]$pop==i,4]),
                         mean(RBC1.out[[1]][RBC1.out[[1]][,6]=="B" & RBC1.out[[1]]$pop==i,5]))
  }
  check.QTL <- lapply(1:6, FUN=function(x) mean(QTL.out[[x]]))

sapply(1:6, FUN=function(x) round((check.QTL[[x]]-check.F2[[x]][1])/check.F2[[x]][2],1))
#0.0 2.0 3.5 4.6 5.3 5.3
sapply(1:6, FUN=function(x) round((check.QTL[[x]]-check.BC1[[x]][1])/check.BC1[[x]][2],1))
#0.0 1.2 2.1 2.6 3.1 3.1
sapply(1:6, FUN=function(x) round((check.QTL[[x]]-check.RBC1[[x]][1])/check.RBC1[[x]][2],1))
#0.0 3.7 6.1 8.2 9.1 9.3

#extract the selection output for plotting.
F2.plot <- list(OSGS.SimCompare(F2.out[[1]], RS=F),
                OSGS.SimCompare(F2.out[[2]][F2.out[[2]]$pop==3 & F2.out[[2]]$wt%in%c("B","0.5","GS"),], RS=T))
BC1.plot <- list(OSGS.SimCompare(BC1.out[[1]], RS=F),
                 OSGS.SimCompare(BC1.out[[2]][BC1.out[[2]]$pop==3 & BC1.out[[2]]$wt%in%c("B","0.5","GS"),], RS=T))
RBC1.plot <- list(OSGS.SimCompare(RBC1.out[[1]], RS=F),
                  OSGS.SimCompare(RBC1.out[[2]][RBC1.out[[2]]$pop==3 & RBC1.out[[2]]$wt%in%c("B","0.5","GS"),], RS=T))

##################

#plot selection in F2.
#####################

#Plot simulation results for F2, pop 1-4.
dat.plot <- F2.plot[[1]][F2.plot[[1]]$pop%in%c(1:4) & F2.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_F2_pop1to4.png",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_F2_pop1to4.svg",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for F2, pop 1-6.
dat.plot <- F2.plot[[1]][F2.plot[[1]]$pop%in%c(1:6) & F2.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30","5"="80:20","6"="90:10"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_F2_pop1to6.png",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_F2_pop1to6.svg",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for F2, pop 3 with recurrent selection.
dat.plot <- F2.plot[[2]][F2.plot[[2]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0.5","GS"),
                      labels=c("0.5","GS"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))
dummy.plot <- data.frame(variable=c("ebv","favP","favS"),
                         value=c(3.5,1,1))

out.plot <- ggplot() +
  geom_hline(data=dummy.plot, aes(yintercept=value), color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=generation, y=value, fill=wt, color=wt)) +
  facet_wrap(vars(variable), scales="free_y", ncol=1, strip.position="left") +
  theme(strip.background=element_blank(), strip.placement="outside", axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  scale_color_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  xlab("Generation") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_F2_pop3_RS.png",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_F2_pop3_RS.svg",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in")
rm(dat.plot, out.plot, dummy.plot)


#####################

#plot selection in BC1.
######################

#Plot simulation results for BC1, pop 1-4.
dat.plot <- BC1.plot[[1]][BC1.plot[[1]]$pop%in%c(1:4) & BC1.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_BC1_pop1to4.png",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_BC1_pop1to4.svg",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for BC1, pop 1-6.
dat.plot <- BC1.plot[[1]][BC1.plot[[1]]$pop%in%c(1:6) & BC1.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30","5"="80:20","6"="90:10"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_BC1_pop1to6.png",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_BC1_pop1to6.svg",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for BC1, pop 3 with recurrent selection.
dat.plot <- BC1.plot[[2]][BC1.plot[[2]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0.5","GS"),
                      labels=c("0.5","GS"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))
dummy.plot <- data.frame(variable=c("ebv","favP","favS"),
                         value=c(2.1,1,1))

out.plot <- ggplot() +
  geom_hline(data=dummy.plot, aes(yintercept=value), color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=generation, y=value, fill=wt, color=wt)) +
  facet_wrap(vars(variable), scales="free_y", ncol=1, strip.position="left") +
  theme(strip.background=element_blank(), strip.placement="outside", axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  scale_color_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  xlab("Generation") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_BC1_pop3_RS.png",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_BC1_pop3_RS.svg",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in")
rm(dat.plot, out.plot, dummy.plot)
######################

#plot selection in RBC1.
#######################

#Plot simulation results for RBC1, pop 1-4.
dat.plot <- RBC1.plot[[1]][RBC1.plot[[1]]$pop%in%c(1:4) & RBC1.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_RBC1_pop1to4.png",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_RBC1_pop1to4.svg",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for RBC1, pop 1-6.
dat.plot <- RBC1.plot[[1]][RBC1.plot[[1]]$pop%in%c(1:6) & RBC1.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30","5"="80:20","6"="90:10"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_RBC1_pop1to6.png",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_RBC1_pop1to6.svg",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for RBC1, pop 3 with recurrent selection.
dat.plot <- RBC1.plot[[2]][RBC1.plot[[2]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0.5","GS"),
                      labels=c("0.5","GS"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))
dummy.plot <- data.frame(variable=c("ebv","favP","favS"),
                         value=c(6.1,1,1))

out.plot <- ggplot() +
  geom_hline(data=dummy.plot, aes(yintercept=value), color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=generation, y=value, fill=wt, color=wt)) +
  facet_wrap(vars(variable), scales="free_y", ncol=1, strip.position="left") +
  theme(strip.background=element_blank(), strip.placement="outside", axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  scale_color_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  xlab("Generation") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_RBC1_pop3_RS.png",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_RBC1_pop3_RS.svg",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in")
rm(dat.plot, out.plot, dummy.plot)
#######################

#plot estimated PS.
##################

for(i in 1:6){
  colnames(PS.out[[i]]) <- c("F2_P","F2_S","BC1_P","BC1_S","rBC1_P","rBC1_S")
  PS.out[[i]] <- melt(data.frame(PS.out[[i]]))
  PS.out[[i]] <- data.frame(str_split_fixed(as.character(PS.out[[i]]$variable), "_", 2),
                            PS.out[[i]]$value,
                            c("50:50","55:45","60:40","70:30","80:20","90:10")[i])
  colnames(PS.out[[i]]) <- c("Pop","PS","Proportion","QTL")
}

plot.PS <- do.call(rbind, PS.out)
plot.PS$Pop <- factor(plot.PS$Pop, levels=c("F2","BC1","rBC1"))

ggplot() +
  geom_boxplot(data=plot.PS, aes(x=Pop, y=Proportion, fill=PS), width=0.4, position=position_dodge(0.5)) +
  facet_wrap(vars(QTL), nrow=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  scale_fill_manual(values=c("#FF0000", "#00FF00")) +
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  theme(legend.title=element_blank(), legend.key=element_blank())

ggsave(filename="QTL_Proportion_2cM.png",
       scale=1,
       width=10,
       height=4,
       units="in",
       dpi=600)
##################

#plot h2.
########

plot.h2 <- do.call(cbind, out.h2)
colnames(plot.h2) <- paste(rep(c("A","B","C"),6),c(rep("A",3),rep("B",3),rep("C",3),rep("D",3),rep("E",3),rep("F",3)),sep="")
plot.h2 <- melt(data.frame(plot.h2))
plot.h2 <- data.frame(str_split_fixed(as.character(plot.h2$variable), "", 2),
                      plot.h2$value)
colnames(plot.h2) <- c("pop","qtl","h2")
plot.h2$pop <- factor(plot.h2$pop,
                      labels=c("F2","BC1","rBC1"))

ggplot() +
  geom_boxplot(data=plot.h2, aes(x=qtl, y=h2, fill=pop), width=0.4, position=position_dodge(0.5)) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  scale_x_discrete(labels=c("50:50", "55:45", "60:40", "70:30", "80:20", "90:10")) +
  theme(legend.title=element_blank(), legend.key=element_blank())

ggsave(filename="QTL_h2_2cM.png",
       scale=1,
       width=6,
       height=2,
       units="in",
       dpi=600)
########



################
### 20cM/QTL ###
################
F2B.out <- list(vector(), vector())
BC1B.out <- list(vector(), vector())
RBC1B.out <- list(vector(), vector())
QTLB.out <- replicate(6, vector())
PSB.out <- replicate(6, vector())
h2B.out <- replicate(6, vector())

set.seed(54693)
for(i in 1:100){
  F2B <- OSGS.SimPop(n.chr=10,
                    n.marker=seq(1000,550,-50),
                    gen.len=seq(2,1.1,-0.1),
                    pop.type="F2",
                    n.line=200,
                    n.self=4)
  BC1B <- OSGS.SimPop(n.chr=10,
                     n.marker=seq(1000,550,-50),
                     gen.len=seq(2,1.1,-0.1),
                     pop.type="BC1",
                     n.line=200,
                     n.self=4)
  RBC1B <- OSGS.SimPop(n.chr=10,
                      n.marker=seq(1000,550,-50),
                      gen.len=seq(2,1.1,-0.1),
                      pop.type="RBC1",
                      n.line=200,
                      n.self=4)
  
  QTL <- OSGS.SimQTL(n.chr=10,
                     n.marker=seq(1000,550,-50),
                     p.qtl=0.01,
                     p.P=c(0.5,0.55,0.6,0.7,0.8,0.9),
                     meanG=0,
                     varG=1)
  
  F2B.Sel <- OSGS.SimSel(S0=F2B[[1]],
                        m.eff=QTL,
                        n.line=200,
                        meanE=0,
                        varE=1,
                        wt=seq(0,1,0.1),
                        n.sel=5,
                        n.DH=200/choose(5,2),
                        p.P=c(0.5,0.55,0.6,0.7,0.8,0.9),
                        RS.gen=5,
                        SP=F2B[[2]])
  BC1B.Sel <- OSGS.SimSel(S0=BC1B[[1]],
                         m.eff=QTL,
                         n.line=200,
                         meanE=0,
                         varE=1,
                         wt=seq(0,1,0.1),
                         n.sel=5,
                         n.DH=200/choose(5,2),
                         p.P=c(0.5,0.55,0.6,0.7,0.8,0.9),
                         RS.gen=5,
                         SP=BC1B[[2]])
  RBC1B.Sel <- OSGS.SimSel(S0=RBC1B[[1]],
                          m.eff=QTL,
                          n.line=200,
                          meanE=0,
                          varE=1,
                          wt=seq(0,1,0.1),
                          n.sel=5,
                          n.DH=200/choose(5,2),
                          p.P=c(0.5,0.55,0.6,0.7,0.8,0.9),
                          RS.gen=5,
                          SP=RBC1B[[2]])
  
  #summarize output from single-gen selection (j=1) and recurrent selection (j=2).
  for(j in 1:2){
    F2B.out[[j]] <- rbind(F2B.out[[j]],
                          data.frame(F2B.Sel[[j]],
                                     sim=i))
    BC1B.out[[j]] <- rbind(BC1B.out[[j]],
                           data.frame(BC1B.Sel[[j]],
                                      sim=i))
    RBC1B.out[[j]] <- rbind(RBC1B.out[[j]],
                            data.frame(RBC1B.Sel[[j]],
                                       sim=i))
  }
  
  #summarize the total QTL effect, i.e. P-parent true breeding value.
  for(j in 1:6){
    QTLB.out[[j]] <- c(QTLB.out[[j]], -sum(QTLB[[j]]))
  }
  
  #summarize the estimated proportion of P/S and h2.
  temp.F2 <- OSGS.SimSummary(out.SimPop=F2B, out.SimQTL=QTL, meanE=0, varE=1)
  temp.BC1 <- OSGS.SimSummary(out.SimPop=F2B, out.SimQTL=QTL, meanE=0, varE=1)
  temp.RBC1 <- OSGS.SimSummary(out.SimPop=F2B, out.SimQTL=QTL, meanE=0, varE=1)
  
  for(j in 1:6){
    PSB.out[[j]] <- rbind(PSB.out[[j]],
                          c(temp.F2[[1]][[j]], temp.BC1[[1]][[j]], temp.RBC1[[1]][[j]]))
    h2B.out[[j]] <- rbind(h2B.out[[j]],
                          c(temp.F2[[2]][[j]], temp.BC1[[2]][[j]], temp.RBC1[[2]][[j]]))
  }
  
  message(i, " done.")
  
}

save(F2B.out, BC1B.out, RBC1B.out, file="sim_20200318.RData")

#prepare plot data.
##################

#identify the P-parent breeding value baseline.
check.F2B <- check.BC1B <- check.RBC1B <- replicate(6,vector())
for(i in 1:6){
  check.F2B[[i]] <- c(mean(F2B.out[[1]][F2B.out[[1]][,6]=="B" & F2B.out[[1]]$pop==i,4]),
                      mean(F2B.out[[1]][F2B.out[[1]][,6]=="B" & F2B.out[[1]]$pop==i,5]))
  check.BC1B[[i]] <- c(mean(BC1B.out[[1]][BC1B.out[[1]][,6]=="B" & BC1B.out[[1]]$pop==i,4]),
                       mean(BC1B.out[[1]][BC1B.out[[1]][,6]=="B" & BC1B.out[[1]]$pop==i,5]))
  check.RBC1B[[i]] <- c(mean(RBC1B.out[[1]][RBC1B.out[[1]][,6]=="B" & RBC1B.out[[1]]$pop==i,4]),
                        mean(RBC1B.out[[1]][RBC1B.out[[1]][,6]=="B" & RBC1B.out[[1]]$pop==i,5]))
}
check.QTLB <- lapply(1:6, FUN=function(x) mean(QTLB.out[[x]]))

sapply(1:6, FUN=function(x) round((check.QTLB[[x]]-check.F2B[[x]][1])/check.F2B[[x]][2],1))
#0.1 0.5 1.4 2.5 3.4 4.1
sapply(1:6, FUN=function(x) round((check.QTLB[[x]]-check.BC1B[[x]][1])/check.BC1B[[x]][2],1))
#0.2 0.2 0.8 1.3 2.0 2.3
sapply(1:6, FUN=function(x) round((check.QTLB[[x]]-check.RBC1B[[x]][1])/check.RBC1B[[x]][2],1))
#0.0 1.0 2.3 4.1 5.8 6.9

#extract the selection output for plotting.
F2B.plot <- list(OSGS.SimCompare(F2B.out[[1]], RS=F),
                 OSGS.SimCompare(F2B.out[[2]][F2B.out[[2]]$pop==3 & F2B.out[[2]]$wt%in%c("B","0.5","GS"),], RS=T))
BC1B.plot <- list(OSGS.SimCompare(BC1B.out[[1]], RS=F),
                  OSGS.SimCompare(BC1B.out[[2]][BC1B.out[[2]]$pop==3 & BC1B.out[[2]]$wt%in%c("B","0.5","GS"),], RS=T))
RBC1B.plot <- list(OSGS.SimCompare(RBC1B.out[[1]], RS=F),
                   OSGS.SimCompare(RBC1B.out[[2]][RBC1B.out[[2]]$pop==3 & RBC1B.out[[2]]$wt%in%c("B","0.5","GS"),], RS=T))

##################


#plot selection in F2B.
######################
#Plot simulation results for F2B, pop 1-4.
dat.plot <- F2B.plot[[1]][F2B.plot[[1]]$pop%in%c(1:4) & F2B.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_F2B_pop1to4.png",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_F2B_pop1to4.svg",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for F2B, pop 1-6.
dat.plot <- F2B.plot[[1]][F2B.plot[[1]]$pop%in%c(1:6) & F2B.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30","5"="80:20","6"="90:10"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_F2B_pop1to6.png",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_F2B_pop1to6.svg",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for F2B, pop 3 with recurrent selection.
dat.plot <- F2B.plot[[2]][F2B.plot[[2]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0.5","GS"),
                      labels=c("0.5","GS"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))
dummy.plot <- data.frame(variable=c("ebv","favP","favS"),
                         value=c(1.4,1,1))

out.plot <- ggplot() +
  geom_hline(data=dummy.plot, aes(yintercept=value), color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=generation, y=value, fill=wt, color=wt)) +
  facet_wrap(vars(variable), scales="free_y", ncol=1, strip.position="left") +
  theme(strip.background=element_blank(), strip.placement="outside", axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  scale_color_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  xlab("Generation") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_F2B_pop3_RS.png",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_F2B_pop3_RS.svg",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in")
rm(dat.plot, out.plot, dummy.plot)
######################
  
#plot selection in BC1B.
#######################
#Plot simulation results for BC1B, pop 1-4.
dat.plot <- BC1B.plot[[1]][BC1B.plot[[1]]$pop%in%c(1:4) & BC1B.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_BC1B_pop1to4.png",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_BC1B_pop1to4.svg",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for BC1B, pop 1-6.
dat.plot <- BC1B.plot[[1]][BC1B.plot[[1]]$pop%in%c(1:6) & BC1B.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30","5"="80:20","6"="90:10"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_BC1B_pop1to6.png",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_BC1B_pop1to6.svg",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for BC1B, pop 3 with recurrent selection.
dat.plot <- BC1B.plot[[2]][BC1B.plot[[2]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0.5","GS"),
                      labels=c("0.5","GS"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))
dummy.plot <- data.frame(variable=c("ebv","favP","favS"),
                         value=c(0.8,1,1))

out.plot <- ggplot() +
  geom_hline(data=dummy.plot, aes(yintercept=value), color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=generation, y=value, fill=wt, color=wt)) +
  facet_wrap(vars(variable), scales="free_y", ncol=1, strip.position="left") +
  theme(strip.background=element_blank(), strip.placement="outside", axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  scale_color_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  xlab("Generation") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_BC1B_pop3_RS.png",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_BC1B_pop3_RS.svg",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in")
rm(dat.plot, out.plot, dummy.plot)
#######################

#plot selection in RBC1B.
########################
#Plot simulation results for RBC1B, pop 1-4.
dat.plot <- RBC1B.plot[[1]][RBC1B.plot[[1]]$pop%in%c(1:4) & RBC1B.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_RBC1B_pop1to4.png",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_RBC1B_pop1to4.svg",
       plot=out.plot,
       scale=1.5,
       width=5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for RBC1B, pop 1-6.
dat.plot <- RBC1B.plot[[1]][RBC1B.plot[[1]]$pop%in%c(1:6) & RBC1B.plot[[1]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                      labels=c("0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))

out.plot <- ggplot() +
  geom_hline(yintercept=0, color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=wt, y=value, fill=Sig, color=Sig)) +
  facet_grid(variable~pop, scales="free_y", labeller=labeller(pop=c("1"="50:50","2"="55:45","3"="60:40","4"="70:30","5"="80:20","6"="90:10"))) +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  scale_color_manual(values=c("#999999","#FF6060"), name="OSGS vs GS: ", labels=c("Not Sig. Diff.","Sig. Diff.")) +
  xlab("Selection Weight") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_RBC1B_pop1to6.png",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_RBC1B_pop1to6.svg",
       plot=out.plot,
       scale=1.5,
       width=7.5,
       height=3,
       units="in")
rm(dat.plot, out.plot)

#Plot simulation results for RBC1B, pop 3 with recurrent selection.
dat.plot <- RBC1B.plot[[2]][RBC1B.plot[[2]]$variable%in%c("favP","favS","ebv_mean"),]

dat.plot$wt <- as.factor(dat.plot$wt)
dat.plot$wt <- factor(dat.plot$wt,
                      levels=c("0.5","GS"),
                      labels=c("0.5","GS"))
dat.plot$Sig <- as.factor(dat.plot$Sig)
dat.plot$Sig <- factor(dat.plot$Sig,
                       levels=c("NS","S"),
                       labels=c("NS","S"))
dat.plot$variable <- factor(dat.plot$variable,
                            levels=c("ebv_mean","favP","favS"),
                            labels=c("ebv","favP","favS"))
dummy.plot <- data.frame(variable=c("ebv","favP","favS"),
                         value=c(2.3,1,1))

out.plot <- ggplot() +
  geom_hline(data=dummy.plot, aes(yintercept=value), color="#000000", linetype=3) +
  geom_violin(data=dat.plot, aes(x=generation, y=value, fill=wt, color=wt)) +
  facet_wrap(vars(variable), scales="free_y", ncol=1, strip.position="left") +
  theme(strip.background=element_blank(), strip.placement="outside", axis.title.y=element_blank()) +
  theme(panel.grid=element_blank(), panel.background=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  scale_color_manual(values=c("#5555FF","#FF9900"), name="Selection: ", labels=c("OSGS (0.5)","GS")) +
  xlab("Generation") +
  theme(legend.position="bottom", legend.title=element_text(size=8), legend.text=element_text(size=8))

ggsave(filename="OSGS_sim_RBC1B_pop3_RS.png",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in",
       dpi=600)

ggsave(filename="OSGS_sim_RBC1B_pop3_RS.svg",
       plot=out.plot,
       scale=1.5,
       width=2,
       height=3,
       units="in")
rm(dat.plot, out.plot, dummy.plot)
########################





#misc summary.
#############

#summarize BV/P/S in S0 and S5 for F2,BC1,rBC1 60:40 populations. (Table 2)
F2.out.temp <- F2.out[[2]][F2.out[[2]]$pop==3 & F2.out[[2]]$generation%in%c("S0","S5"),]
F2.out.summary <- list()
for(i in 1:4){
  F2.out.summary[[i]] <- sapply(unique(F2.out.temp$wt),
                                FUN=function(x) mean(F2.out.temp[F2.out.temp$wt==x,i]))
}
F2.out.summary <- do.call(cbind, F2.out.summary)

BC1.out.temp <- BC1.out[[2]][BC1.out[[2]]$pop==3 & BC1.out[[2]]$generation%in%c("S0","S5"),]
BC1.out.summary <- list()
for(i in 1:4){
  BC1.out.summary[[i]] <- sapply(unique(BC1.out.temp$wt),
                                FUN=function(x) mean(BC1.out.temp[BC1.out.temp$wt==x,i]))
}
BC1.out.summary <- do.call(cbind, BC1.out.summary)

RBC1.out.temp <- RBC1.out[[2]][RBC1.out[[2]]$pop==3 & RBC1.out[[2]]$generation%in%c("S0","S5"),]
RBC1.out.summary <- list()
for(i in 1:4){
  RBC1.out.summary[[i]] <- sapply(unique(RBC1.out.temp$wt),
                                FUN=function(x) mean(RBC1.out.temp[RBC1.out.temp$wt==x,i]))
}
RBC1.out.summary <- do.call(cbind, RBC1.out.summary)

write.csv(cbind(F2.out.summary, BC1.out.summary, RBC1.out.summary),
          "RS_out.csv",
          quote=F, row.names=T)

#plot the distribution of simulated QTL.
plot.qtl <- list()
for(i in c(5000,5500,6000,7000,8000,9000)){
  j <- 10000 - i
  plot.qtl <- c(plot.qtl, list(c(-abs(rnorm(i, 0, sqrt(1/(i+j-2*(j-i)^2/pi/(i+j))))),
                                 abs(rnorm(j, 0, sqrt(1/(i+j-2*(j-i)^2/pi/(i+j))))))))
}

plot.qtl <- do.call(cbind, plot.qtl)
colnames(plot.qtl) <- c("A","B","C","D","E","F")
plot.qtl <- melt(data.frame(plot.qtl))
colnames(plot.qtl) <- c("QTL_Proportion", "Effect")
plot.qtl$QTL_Proportion <- factor(plot.qtl$QTL_Proportion,
                                  labels=c("50:50", "55:45", "60:40", "70:30", "80:20", "90:10"))
ggplot() +
  geom_histogram(data=plot.qtl, aes(x=Effect), bins=15) +
  facet_wrap(vars(QTL_Proportion), nrow=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#000000") +
  theme(legend.title=element_blank(), legend.key=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_x_continuous(breaks=c(-0.05,0,0.05))

ggsave(filename="QTL_effects.png",
       scale=1,
       width=6,
       height=2,
       units="in",
       dpi=600)

#############
