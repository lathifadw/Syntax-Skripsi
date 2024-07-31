library(sf)
library(terra)
library(sp)
library(spData)
library(shapefiles)
library(maptools) 
library(spdep)
library(classInt) 
library(rgeos) 
library(raster)
library(ggplot2) 
library(ggsn) 
library(lmtest) 
library(spatialreg) 
library(gtools) 
library(splm) 
library(plm) 
library(tseries)
library(tidyr)
library(CARBayes) 
library(CARBayesST) 
library(mvtnorm)
library(CAMAN)
library(tmap) 
library(INLA) 
library(tidyverse) 
library(orcutt)
library(HoRM)
library(RColorBrewer)
library(brinla) 
library(lattice)
library(ctv)
library(splancs) 


#Input the map
Indonesia <- getData('GADM', country='IDN', level=3)
Bandung <- Indonesia[Indonesia$NAME_2 == "Kota Bandung",]
plot(Bandung)
Kecamatan <- as.data.frame(Bandung$NAME_3);Kecamatan

#Input Data
setwd("D:/LATHIFA/KULIAH/SKRIPSI/DATA SKRIPSI") 
DiareISPA <- read.csv("Data Skripsi.csv", sep=";", dec=".")
head(DiareISPA)
summary(DiareISPA)

## MENGHUBUNGKAN DATA DENGAN PETA ##
# - Buat ID untuk menghubungkan Peta dengan Data
ID <- c(1:30)
Bandung$id <- ID
X <- coordinates(Bandung)[,1]
Y <- coordinates(Bandung)[,2]
DiareISPA$id <- rep(ID, 84)
DiareISPA$X <- rep(X, 84)
DiareISPA$Y <- rep(Y, 84)
# - Menghubungkan Peta dengan Data
BandungData <- fortify(Bandung)
id <- as.numeric(unique(BandungData$id))
DiareISPA$id <-rep(id,84)
data.shp <- merge(BandungData, DiareISPA, by="id", all.X=TRUE)
dataMap <- data.shp[order(data.shp$order), ]

## MATRIKS PEMBOBOT SPASIAL ##
# Queen
W <- poly2nb(Bandung, row.names=ID, queen=TRUE) 
WB <- nb2mat(W, style='B', zero.policy = TRUE) 
Wls<- nb2listw(W, zero.policy = TRUE)
Ws <- as(as_dgRMatrix_listw(Wls), "CsparseMatrix")

plot(Bandung, col="pink",main="Queen")
plot(Wls,CoordK,add=T, col="red", lwd=1)

## PENGUJIAN AUTOKORELASI SPATIALTEMPORAL ##
Wsmatrix <- as.matrix(Ws)
# - Diare
KasusDiare <- DiareISPA$Y1
Y1 <- as.matrix(KasusDiare)
MoranST_Diare <- MoranST.MC(Y1,Wsmatrix,100)
# - ISPA
KasusISPA <- DiareISPA$Y2
Y2 <- as.matrix(KasusISPA)
MoranST_ISPA <- MoranST.MC(Y2,Wsmatrix,100)

## PEMILIHAN MODEL TERBAIK ##
# - Setup INLA Output
control <- list(predictor = list(compute = TRUE,link=1),
                results = list(return.marginals.random = TRUE),
                compute = list(return.marginals.predictor=TRUE, hyperpar=TRUE, return.marginals=TRUE, 
                               dic=TRUE, mlik = TRUE, cpo = TRUE, po = TRUE, waic=TRUE, graph=TRUE, 
                               openmp.strategy="huge"))

# - HC Prior
halfcauchy <- "expression: lambda = 25; precision = exp(log_precision); 
logdens = -1.5*log_precision-log(pi*lambda)-log(1+1/(precision*lambda^2));
log_jacobian = log_precision;
return(logdens+log_jacobian);"
hcprior = list(prec = list(prior = halfcauchy))
hcprior <- hcprior
# - Modeling
# Membentuk Struktur
Data <- read.csv("Data Skripsi.csv", sep=";", dec=".")
k <- 2
n <- dim(Data)[1]
Y <- matrix(NA, n, k)
Y[1:n, 1] <- Data$Y1
Y[1:n, 2] <- Data$Y2

share.dat <- list(Y=matrix(NA, nrow=n*2, ncol=2))
share.dat$Y[1:n, 1] <- Y[,1]
share.dat$Y[n+(1:n), 2] <- Y[,2]

share.dat$Observed <- c(Data$Y1,Data$Y2)
share.dat$E <- c(Data$EiD,Data$EiI)

share.dat$Region <-c(Data$Kecamatan,Data$Kecamatan)
share.dat$Disease <- c(rep(1,n),rep(2,n))
share.dat$Time <-c (Data$IT1,Data$IT1) 
share.dat$Spasial <-c (Data$ID1,Data$ID1)

## Tahun
share.dat$Tahun=c(Data$Tahun,Data$Tahun)
## Spasial
share.dat$sharedSDiare <- c(Data$ID1, rep(NA,n))
share.dat$sharedSISPA <- c(rep(NA,n),Data$ID1)

share.dat$spatialSDiare<- c(Data$ID1, rep(NA,n))
share.dat$spatialSISPA <- c(rep(NA,n),Data$ID1)

share.dat$randomSDiare <- c(Data$ID1, rep(NA,n))
share.dat$randomSISPA <- c(rep(NA,n),Data$ID1)

share.dat$interacSDiare <- c(Data$ID1, rep(NA,n))
share.dat$interacSISPA <- c(rep(NA,n),Data$ID1)
## Temporal
share.dat$sharedTDiare <- c(Data$IT1, rep(NA,n))
share.dat$sharedTISPA <- c(rep(NA,n),Data$IT1)

share.dat$temporalTDiare <- c(Data$IT1, rep(NA,n))
share.dat$temporalTISPA <- c(rep(NA,n),Data$IT1)

share.dat$randomTDiare <- c(Data$IT1, rep(NA,n))
share.dat$randomTISPA <- c(rep(NA,n),Data$IT1)

share.dat$interacTDiare <- c(Data$IT1, rep(NA,n))
share.dat$interacTISPA <- c(rep(NA,n),Data$IT1)
## --- ##
share.dat$interacDiare <- c(1:n, rep(NA,n))
share.dat$interacISPA <- c(rep(NA,n),1:n)
share.dat$alphaDiare <- rep(c(1,0), each=n)
share.dat$alphaISPA <- rep(c(0,1), each=n)
## PHBS
share.dat$DiarePHBS <- c(Data$X1,rep(0,n))
share.dat$ISPAPHBS <- c(rep(0,n),Data$X1)
## Pemberian ASI Eksklusif
share.dat$DiareASI <- c(Data$X2,rep(0,n))
share.dat$ISPAASI <- c(rep(0,n),Data$X2)

## Menggunakan Random Walk Orde 1 ##
## - Model with covariate Only + Shared Component
Model1 <- Y ~ 0 + alphaDiare + alphaISPA + DiarePHBS:factor(Time) + ISPAPHBS:factor(Time) +  
  DiareASI + ISPAASI +
  f(sharedTDiare, model="rw1", constr=T, scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(sharedTISPA, copy="sharedTDiare", 
    hyper=list(theta=list(fixed=FALSE, range=c(0,Inf))))+ 
  f(sharedSDiare, model="besagproper2", graph=Ws, constr=T, hyper=hcprior) +
  f(sharedSISPA, copy="sharedSDiare", 
    hyper=list(theta=list(fixed=FALSE, param=c(1,1), range=c(0,Inf))))
RModel1 <- inla(Model1, family=c("nbinomial", "nbinomial"), data=share.dat,
                control.compute = control$compute,control.predictor = control$predictor,
                control.inla = list(int.strategy = "eb", 
                                    strategy = "simplified.laplace"))

## - Model with covariate + Shared Component + Temporal Trend RW1
Model2  <- Y ~ 0 + alphaDiare + alphaISPA + DiarePHBS:factor(Time) + ISPAPHBS:factor(Time) +  
  DiareASI + ISPAASI +
  f(temporalTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(temporalTISPA, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(sharedTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(sharedTISPA, copy="sharedTDiare", 
    hyper=list(theta=list(fixed=FALSE, range=c(0,Inf))))+ 
  f(sharedSDiare, model="besagproper2", graph=Ws, constr=T, hyper=hcprior) +
  f(sharedSISPA, copy="sharedSDiare", 
    hyper=list(theta=list(fixed=FALSE, param=c(1,1), range=c(0,Inf))))
RModel2 <- inla(Model2, family=c("nbinomial", "nbinomial"), data=share.dat, E=E,
                control.compute = control$compute,control.predictor = control$predictor,
                control.inla = list(int.strategy = "eb",
                                    strategy = "simplified.laplace"))

## - Model with covariate + Shared Component + Spatial Dependency
Model3 <- Y ~ 0 + alphaDiare + alphaISPA + DiarePHBS:factor(Time) + ISPAPHBS:factor(Time) +  
  DiareASI + ISPAASI +
  f(spatialSDiare, model="besagproper2", graph=Ws, constr=T, 
    hyper=hcprior) +
  f(spatialSISPA, model="besagproper2", graph=Ws, constr=T, 
    hyper=hcprior) +
  f(sharedTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(sharedTISPA, copy="sharedTDiare", 
    hyper=list(theta=list(fixed=FALSE, range=c(0,Inf))))+ 
  f(sharedSDiare, model="besagproper2", graph=Ws, constr=T, hyper=hcprior) +
  f(sharedSISPA, copy="sharedSDiare", 
    hyper=list(theta=list(fixed=FALSE, param=c(1,1), range=c(0,Inf))))
RModel3 <- inla(Model3, family=c("nbinomial", "nbinomial"), data=share.dat, E=E,
                control.compute = control$compute,control.predictor = control$predictor,
                control.inla = list(int.strategy = "eb",
                                    strategy = "simplified.laplace"))

## - Model with covariate + Shared Component + Interaction Type IV
Model4 <- Y ~ 0 + alphaDiare + alphaISPA + DiarePHBS:factor(Time) + ISPAPHBS:factor(Time) +  
  DiareASI + ISPAASI +
  f(interacSDiare,model="besagproper2", graph=Ws, initial=1,constr=T, 
    hyper=hcprior, group=interacTDiare, control.group=list(model="rw1", 
                                                           hyper=hcprior)) + 
  f(interacSISPA,model="besagproper2", graph=Ws, initial=1,constr=T,
    hyper=hcprior, group=interacTISPA, control.group=list(model="rw1", 
                                                          hyper=hcprior)) +
  f(sharedTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(sharedTISPA, copy="sharedTDiare", 
    hyper=list(theta=list(fixed=FALSE, range=c(0,Inf))))+ 
  f(sharedSDiare, model="besagproper2", graph=Ws, constr=T, hyper=hcprior) +
  f(sharedSISPA, copy="sharedSDiare", 
    hyper=list(theta=list(fixed=FALSE, param=c(1,1), range=c(0,Inf))))
RModel4 <- inla(Model4, family=c("nbinomial", "nbinomial"), data=share.dat, E=E,
                control.compute = control$compute,control.predictor = control$predictor,
                control.inla = list(int.strategy = "eb", 
                                    strategy = "simplified.laplace"))

## - Model with covariate + Shared Component + Temporal Trend RW1 + Spatial Dependency
Model5 <- Y ~ 0 + alphaDiare + alphaISPA + DiarePHBS:factor(Time) + ISPAPHBS:factor(Time) +  
  DiareASI + ISPAASI +
  f(temporalTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(temporalTISPA, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) + 
  f(spatialSDiare, model="besagproper2", graph=Ws, constr=T, 
    hyper=hcprior) +
  f(spatialSISPA, model="besagproper2", graph=Ws, constr=T, 
    hyper=hcprior) +
  f(sharedTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(sharedTISPA, copy="sharedTDiare", 
    hyper=list(theta=list(fixed=FALSE, range=c(0,Inf))))+ 
  f(sharedSDiare, model="besagproper2", graph=Ws, constr=T, hyper=hcprior) +
  f(sharedSISPA, copy="sharedSDiare", 
    hyper=list(theta=list(fixed=FALSE, param=c(1,1), range=c(0,Inf))))
RModel5 <- inla(Model5, family=c("nbinomial", "nbinomial"), data=share.dat, E=E,
                control.compute = control$compute,control.predictor = control$predictor,
                control.inla = list(int.strategy = "eb", 
                                    strategy = "simplified.laplace"))

## - Model with covariate + Shared Component + Temporal Trend RW1 + Interaction Type IV
Model6 <- Y ~ 0 + alphaDiare + alphaISPA + DiarePHBS:factor(Time) + ISPAPHBS:factor(Time) +  
  DiareASI + ISPAASI +
  f(temporalTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(temporalTISPA, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(interacSDiare,model="besagproper2", graph=Ws, initial=1,constr=T, 
    hyper=hcprior, group=interacTDiare, control.group=list(model="rw1", 
                                                           hyper=hcprior)) + 
  f(interacSISPA,model="besagproper2", graph=Ws, initial=1,constr=T,
    hyper=hcprior, group=interacTISPA, control.group=list(model="rw1", 
                                                          hyper=hcprior)) +
  f(sharedTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(sharedTISPA, copy="sharedTDiare", 
    hyper=list(theta=list(fixed=FALSE, range=c(0,Inf))))+ 
  f(sharedSDiare, model="besagproper2", graph=Ws, constr=T, hyper=hcprior) +
  f(sharedSISPA, copy="sharedSDiare", 
    hyper=list(theta=list(fixed=FALSE, param=c(1,1), range=c(0,Inf))))
RModel6 <- inla(Model6, family=c("nbinomial", "nbinomial"), data=share.dat, E=E,
                control.compute = control$compute,control.predictor = control$predictor,
                control.inla = list(int.strategy = "eb", 
                                    strategy = "simplified.laplace"))

## - Model with covariate + Shared Component + Spatial Dependency + Interaction Type IV
Model7 <- Y ~ 0 + alphaDiare + alphaISPA + DiarePHBS:factor(Time) + ISPAPHBS:factor(Time) +  
  DiareASI + ISPAASI +
  f(spatialSDiare, model="besagproper2", graph=Ws, constr=T, 
    hyper=hcprior) +
  f(spatialSISPA, model="besagproper2", graph=Ws, constr=T, 
    hyper=hcprior) +
  f(interacSDiare,model="besagproper2", graph=Ws, initial=1,constr=T, 
    hyper=hcprior, group=interacTDiare, control.group=list(model="rw1", 
                                                           hyper=hcprior)) + 
  f(interacSISPA,model="besagproper2", graph=Ws, initial=1,constr=T,
    hyper=hcprior, group=interacTISPA, control.group=list(model="rw1", 
                                                          hyper=hcprior)) +
  f(sharedTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(sharedTISPA, copy="sharedTDiare", 
    hyper=list(theta=list(fixed=FALSE, range=c(0,Inf))))+ 
  f(sharedSDiare, model="besagproper2", graph=Ws, constr=T, hyper=hcprior) +
  f(sharedSISPA, copy="sharedSDiare", 
    hyper=list(theta=list(fixed=FALSE, param=c(1,1), range=c(0,Inf))))
RModel7 <- inla(Model7, family=c("nbinomial", "nbinomial"), data=share.dat, E=E,
                control.compute = control$compute,control.predictor = control$predictor,
                control.inla = list(int.strategy = "eb", 
                                    strategy = "simplified.laplace"))

## - Model with covariate + Shared Component + Temporal Trend RW1 + Spatial Dependency + Interaction Type IV
Model8 <- Y ~ 0 + alphaDiare + alphaISPA + DiarePHBS:factor(Time) + ISPAPHBS:factor(Time) +  
  DiareASI + ISPAASI +
  f(temporalTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(temporalTISPA, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) + 
  f(spatialSDiare, model="besagproper2", graph=Ws, constr=T, 
    hyper=hcprior) +
  f(spatialSISPA, model="besagproper2", graph=Ws, constr=T, 
    hyper=hcprior) +
  f(interacSDiare,model="besagproper2", graph=Ws, initial=1,constr=T, 
    hyper=hcprior, group=interacTDiare, control.group=list(model="rw1", 
                                                           hyper=hcprior)) + 
  f(interacSISPA,model="besagproper2", graph=Ws, initial=1,constr=T,
    hyper=hcprior, group=interacTISPA, control.group=list(model="rw1", 
                                                          hyper=hcprior)) +
  f(sharedTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(sharedTISPA, copy="sharedTDiare", 
    hyper=list(theta=list(fixed=FALSE, range=c(0,Inf))))+ 
  f(sharedSDiare, model="besagproper2", graph=Ws, constr=T, hyper=hcprior) +
  f(sharedSISPA, copy="sharedSDiare", 
    hyper=list(theta=list(fixed=FALSE, param=c(1,1), range=c(0,Inf))))
RModel8 <- inla(Model8, family=c("nbinomial", "nbinomial"), data=share.dat, E=E,
                control.compute = control$compute,control.predictor = control$predictor,
                control.inla = list(int.strategy = "eb", 
                                    strategy = "simplified.laplace"))

summary(RModel1)
summary(RModel2)
summary(RModel3)
summary(RModel4)
summary(RModel5)
summary(RModel6)
summary(RModel7)
summary(RModel8)

## - R Korelasi 
pred1=RModel1$summary.fitted.values$mean*share.dat$E
cor1=cor(pred1,share.dat$Observed); cor1

pred2=RModel2$summary.fitted.values$mean*share.dat$E
cor2=cor(pred2,share.dat$Observed); cor2

pred3=RModel3$summary.fitted.values$mean*share.dat$E
cor3=cor(pred3,share.dat$Observed); cor3

pred4=RModel4$summary.fitted.values$mean*share.dat$E
cor4=cor(pred4,share.dat$Observed); cor4

pred5=RModel5$summary.fitted.values$mean*share.dat$E
cor5=cor(pred5,share.dat$Observed); cor5

pred6=RModel6$summary.fitted.values$mean*share.dat$E
cor6=cor(pred6,share.dat$Observed); cor6

pred7=RModel7$summary.fitted.values$mean*share.dat$E
cor7=cor(pred7,share.dat$Observed); cor7

pred8=RModel8$summary.fitted.values$mean*share.dat$E
cor8=cor(pred8,share.dat$Observed); cor8

pred5a=RModel5a$summary.fitted.values$mean*share.dat$E
cor5a=cor(pred5a,share.dat$Observed); cor5a

## EVALUASI KOMPONEN DARI MODEL TERBAIK ##
plot(RModel1)
plot(RModel2)
plot(RModel3)
plot(RModel4)
plot(RModel5)
plot(RModel6)
plot(RModel7)
plot(RModel8)

bri.hyperpar.summary(RModel5)

## UJI SIGNIFIKANSI vARIABEL ## 
summary(RModel5)

## MODEL FIX ##
ModelFix <- Y ~ 0 + alphaDiare + alphaISPA + DiarePHBS:factor(Time) + ISPAPHBS:factor(Time) + 
  f(temporalTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(temporalTISPA, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) + 
  f(spatialSDiare, model="besagproper2", graph=Ws, constr=T, 
    hyper=hcprior) +
  f(spatialSISPA, model="besagproper2", graph=Ws, constr=T, 
    hyper=hcprior) +
  f(sharedTDiare, model="rw1", constr=T,scale.model=FALSE, 
    cyclic=TRUE, hyper=hcprior) +
  f(sharedTISPA, copy="sharedTDiare", 
    hyper=list(theta=list(fixed=FALSE, range=c(0,Inf))))+ 
  f(sharedSDiare, model="besagproper2", graph=Ws, constr=T, hyper=hcprior) +
  f(sharedSISPA, copy="sharedSDiare", 
    hyper=list(theta=list(fixed=FALSE, param=c(1,1), range=c(0,Inf))))
RModelFix <- inla(ModelFix, family=c("nbinomial", "nbinomial"), data=share.dat, E=E,
                  control.compute = control$compute,control.predictor = control$predictor,
                  control.inla = list(int.strategy = "eb", 
                                       strategy = "simplified.laplace"))
summary(RModelFix)
bri.hyperpar.summary(RModelFix)

## MENGHITUNG TAKSIRAN RISIKO RELATIF ##  
setwd("D:/LATHIFA/KULIAH/SKRIPSI/SIDANG SKRIPSI") 
RR <- RModelFix$summary.fitted.values$mean
write.csv(RR,"Taksiran RR Diare ISPA.csv")
RRsc <- RModelFix$summary.random
SCTemporal <- data.frame(RRsc$sharedTDiare,RRsc$sharedTISPA)
write.csv(SCTemporal,"Taksiran RR SC Temporal.csv")
SCSpasial <- data.frame(RRsc$sharedSDiare,RRsc$sharedSISPA)
write.csv(SCSpasial,"Taksiran RR SC Spasial.csv")

Temporal <- data.frame(RRsc$temporalTDiare,RRsc$temporalTISPA)
write.csv(Temporal,"Taksiran RR Temporal.csv")
Spasial <- data.frame(RRsc$spatialSDiare,RRsc$spatialSISPA)
write.csv(Spasial,"Taksiran RR Spasial.csv")


## PETA SHARED COMPONENT SPASIAL ##
ggplot(BandungMap, aes(fill = RR.SC, x = long, y = lat, group = group)) +
  geom_polygon(color = "black", linewidth = 0.01) +
  facet_grid(. ~ Penyakit) +
  theme_bw() +
  ylab("Latitude") +
  xlab("Longitude") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(name = NULL, low = "white", high = "red")

## IDENTIFIKASI HOTSPOT MENGGUNAKAN NILAI PELUANG RISIKO RELATIF ##
marg <- RModelFix$marginals.fitted.values[[1]]
Exceedance <- sapply(RModelFix$marginals.fitted.values, FUN=function(marg){
  1 - inla.pmarginal(q = 1, marginal = marg)})

# Diare
Cluster1 <- c()
for(i in 1:2520){
  if (Exceedance[i] > 0.95) {
    Cluster1[i] <- 1
  } else{
    Cluster1[i] <- 0
  }
}
CLUSTER1 <- as.vector(Cluster1)
Pred.HS1 <- sum(CLUSTER1)/2520*100

# ISPA
Cluster2 <- numeric(2520)  
for(i in 2521:5040){
  if (Exceedance[i] > 0.95) {
    Cluster2[i - 2520] <- 1  
  } else{
    Cluster2[i - 2520] <- 0  
  }
}
CLUSTER2 <- as.vector(Cluster2)
Pred.HS2 <- sum(CLUSTER2)/2520*100

## PEMETAAN EXCEEDANCE POSTERIOR PROB DENGAN THRESHOLD (1) DAN CUTOFF (0.95)
# Diare
Exceed1 <- c()
for(i in 1:2520){
  Exceed1[i] <- Exceedance[i] 
}
EXCEED1       <- as.data.frame(Exceed1)
EXCEED1$id    <- Data$ID1 
EXCEED1$Bulan <- Data$Bulan
EXCEED1$Tahun <- Data$Tahun

BandungData2 <- fortify(Bandung)
id <- as.numeric(unique(BandungData2$id))
EXCEED1$id <- rep(id,84)
data.shp2.1 <- merge(BandungData2, EXCEED1, by="id", all.X=TRUE)
dataMap2.1 <- data.shp2.1[order(data.shp2.1$order), ]

# - Tahun 2017
BDGMAP.2017 <- data.shp2.1[data.shp2.1$Tahun==2017,]

pretty_breaks17 <- c(0.95)

labels17 <- c()
brks17 <- c(-1, pretty_breaks17,9.5)
for(idx in 1:length(brks17)){
  labels17 <- c(labels17,round(brks17[idx + 1], 2)) }

labels17 <- labels17[1:length(labels17)-1]
BDGMAP.2017$brks17 <- cut(BDGMAP.2017$Exceed1,breaks = brks17,include.lowest = TRUE,labels = labels17)

brks_scale17 <- levels(BDGMAP.2017$brks17)
labels_scale17 <- rev(brks_scale17)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2017` = "2017")

BDGMAP.2017 <- na.omit(BDGMAP.2017) 
BDGMAP.2017$TahunBulan <- with(BDGMAP.2017, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP.2017, aes(fill = brks17, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2017)

# - Tahun 2018
BDGMAP.2018 <- data.shp2.1[data.shp2.1$Tahun==2018,]

pretty_breaks18 <- c(0.95)

labels18 <- c()
brks18 <- c(-1, pretty_breaks18,9.5)
for(idx in 1:length(brks18)){
  labels18 <- c(labels18,round(brks18[idx + 1], 2)) }

labels18 <- labels18[1:length(labels18)-1]
BDGMAP.2018$brks18 <- cut(BDGMAP.2018$Exceed1,breaks = brks18,include.lowest = TRUE,labels = labels18)

brks_scale18 <- levels(BDGMAP.2018$brks18)
labels_scale18 <- rev(brks_scale18)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2018` = "2018")

BDGMAP.2018 <- na.omit(BDGMAP.2018) 
BDGMAP.2018$TahunBulan <- with(BDGMAP.2018, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP.2018, aes(fill = brks18, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2018)

# - Tahun 2019
BDGMAP.2019 <- data.shp2.1[data.shp2.1$Tahun==2019,]

pretty_breaks19 <- c(0.95)

labels19 <- c()
brks19 <- c(-1, pretty_breaks19,9.5)
for(idx in 1:length(brks19)){
  labels19 <- c(labels19,round(brks19[idx + 1], 2)) }

labels19 <- labels19[1:length(labels19)-1]
BDGMAP.2019$brks19 <- cut(BDGMAP.2019$Exceed1,breaks = brks19,include.lowest = TRUE,labels = labels19)

brks_scale19 <- levels(BDGMAP.2019$brks19)
labels_scale19 <- rev(brks_scale19)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2019` = "2019")

BDGMAP.2019 <- na.omit(BDGMAP.2019) 
BDGMAP.2019$TahunBulan <- with(BDGMAP.2019, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP.2019, aes(fill = brks19, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2019)

# - Tahun 2020
BDGMAP.2020 <- data.shp2.1[data.shp2.1$Tahun==2020,]

pretty_breaks20 <- c(0.95)

labels20 <- c()
brks20 <- c(-1, pretty_breaks20,9.5)
for(idx in 1:length(brks20)){
  labels20 <- c(labels20,round(brks20[idx + 1], 2)) }

labels20 <- labels20[1:length(labels20)-1]
BDGMAP.2020$brks20 <- cut(BDGMAP.2020$Exceed1,breaks = brks20,include.lowest = TRUE,labels = labels20)

brks_scale20 <- levels(BDGMAP.2020$brks20)
labels_scale20 <- rev(brks_scale20)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2020` = "2020")

BDGMAP.2020 <- na.omit(BDGMAP.2020) 
BDGMAP.2020$TahunBulan <- with(BDGMAP.2020, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP.2020, aes(fill = brks20, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2020)

# - Tahun 2021
BDGMAP.2021 <- data.shp2.1[data.shp2.1$Tahun==2021,]

pretty_breaks21 <- c(0.95)

labels21 <- c()
brks21 <- c(-1, pretty_breaks21,9.5)
for(idx in 1:length(brks21)){
  labels21 <- c(labels21,round(brks21[idx + 1], 2)) }

labels21 <- labels21[1:length(labels21)-1]
BDGMAP.2021$brks21 <- cut(BDGMAP.2021$Exceed1,breaks = brks21,include.lowest = TRUE,labels = labels21)

brks_scale21 <- levels(BDGMAP.2021$brks21)
labels_scale21 <- rev(brks_scale21)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2021` = "2021")

BDGMAP.2021 <- na.omit(BDGMAP.2021) 
BDGMAP.2021$TahunBulan <- with(BDGMAP.2021, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP.2021, aes(fill = brks21, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2021)

# - Tahun 2022
BDGMAP.2022 <- data.shp2.1[data.shp2.1$Tahun==2022,]

pretty_breaks22 <- c(0.95)

labels22 <- c()
brks22 <- c(-1, pretty_breaks22,9.5)
for(idx in 1:length(brks22)){
  labels22 <- c(labels22,round(brks22[idx + 1], 2)) }

labels22 <- labels22[1:length(labels22)-1]
BDGMAP.2022$brks22 <- cut(BDGMAP.2022$Exceed1,breaks = brks22,include.lowest = TRUE,labels = labels22)

brks_scale22 <- levels(BDGMAP.2022$brks22)
labels_scale22 <- rev(brks_scale22)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2022` = "2022")

BDGMAP.2022 <- na.omit(BDGMAP.2022) 
BDGMAP.2022$TahunBulan <- with(BDGMAP.2022, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP.2022, aes(fill = brks22, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2022)

# - Tahun 2023
BDGMAP.2023 <- data.shp2.1[data.shp2.1$Tahun==2023,]

pretty_breaks23 <- c(0.95)

labels23 <- c()
brks23 <- c(-1, pretty_breaks23,9.5)
for(idx in 1:length(brks23)){
  labels23 <- c(labels23,round(brks23[idx + 1], 2)) }

labels23 <- labels23[1:length(labels23)-1]
BDGMAP.2023$brks23 <- cut(BDGMAP.2023$Exceed1,breaks = brks23,include.lowest = TRUE,labels = labels23)

brks_scale23 <- levels(BDGMAP.2023$brks23)
labels_scale23 <- rev(brks_scale23)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2023` = "2023")

BDGMAP.2023 <- na.omit(BDGMAP.2023) 
BDGMAP.2023$TahunBulan <- with(BDGMAP.2023, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP.2023, aes(fill = brks23, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2023)

# ISPA
Exceed2 <- c()
for(i in 2521:5040){
  Exceed2[i - 2520] <- Exceedance[i]
} 

EXCEED2       <- as.data.frame(Exceed2)
EXCEED2$id    <- Data$ID1 
EXCEED2$Bulan <- Data$Bulan
EXCEED2$Tahun <- Data$Tahun

BandungData2 <- fortify(Bandung)
id <- as.numeric(unique(BandungData2$id))
EXCEED2$id <- rep(id,84)
data.shp2.2 <- merge(BandungData2, EXCEED2, by="id", all.X=TRUE)
dataMap2.2 <- data.shp2.2[order(data.shp2.2$order), ]

# - Tahun 2017
BDGMAP2.2017<- data.shp2.2[data.shp2.2$Tahun==2017,]

pretty_breaks17.2 <- c(0.95)

labels17.2 <- c()
brks17.2 <- c(-1, pretty_breaks17.2,9.5)
for(idx in 1:length(brks17.2)){
  labels17.2 <- c(labels17.2,round(brks17.2[idx + 1], 2)) }

labels17.2 <- labels17.2[1:length(labels17.2)-1]
BDGMAP2.2017$brks17.2 <- cut(BDGMAP2.2017$Exceed2,breaks = brks17.2,include.lowest = TRUE,labels = labels17.2)

brks_scale17.2 <- levels(BDGMAP2.2017$brks17.2)
labels_scale17.2 <- rev(brks_scale17.2)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2017` = "2017")

BDGMAP2.2017 <- na.omit(BDGMAP2.2017) 
BDGMAP2.2017$TahunBulan <- with(BDGMAP2.2017, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP2.2017, aes(fill = brks17.2, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2017)

# - Tahun 2018
BDGMAP2.2018<- data.shp2.2[data.shp2.2$Tahun==2018,]

pretty_breaks18.2 <- c(0.95)

labels18.2 <- c()
brks18.2 <- c(-1, pretty_breaks18.2,9.5)
for(idx in 1:length(brks18.2)){
  labels18.2 <- c(labels18.2,round(brks18.2[idx + 1], 2)) }

labels18.2 <- labels18.2[1:length(labels18.2)-1]
BDGMAP2.2018$brks18.2 <- cut(BDGMAP2.2018$Exceed2,breaks = brks18.2,include.lowest = TRUE,labels = labels18.2)

brks_scale18.2 <- levels(BDGMAP2.2018$brks18.2)
labels_scale18.2 <- rev(brks_scale18.2)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2018` = "2018")

BDGMAP2.2018 <- na.omit(BDGMAP2.2018) 
BDGMAP2.2018$TahunBulan <- with(BDGMAP2.2018, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP2.2018, aes(fill = brks18.2, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2018)

# - Tahun 2019
BDGMAP2.2019<- data.shp2.2[data.shp2.2$Tahun==2019,]

pretty_breaks19.2 <- c(0.95)

labels19.2 <- c()
brks19.2 <- c(-1, pretty_breaks19.2,9.5)
for(idx in 1:length(brks19.2)){
  labels19.2 <- c(labels19.2,round(brks19.2[idx + 1], 2)) }

labels19.2 <- labels19.2[1:length(labels19.2)-1]
BDGMAP2.2019$brks19.2 <- cut(BDGMAP2.2019$Exceed2,breaks = brks19.2,include.lowest = TRUE,labels = labels19.2)

brks_scale19.2 <- levels(BDGMAP2.2019$brks19.2)
labels_scale19.2 <- rev(brks_scale19.2)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2019` = "2019")

BDGMAP2.2019 <- na.omit(BDGMAP2.2019) 
BDGMAP2.2019$TahunBulan <- with(BDGMAP2.2019, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP2.2019, aes(fill = brks19.2, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2019)

# - Tahun 2020
BDGMAP2.2020<- data.shp2.2[data.shp2.2$Tahun==2020,]

pretty_breaks20.2 <- c(0.95)

labels20.2 <- c()
brks20.2 <- c(-1, pretty_breaks20.2,9.5)
for(idx in 1:length(brks20.2)){
  labels20.2 <- c(labels20.2,round(brks20.2[idx + 1], 2)) }

labels20.2 <- labels20.2[1:length(labels20.2)-1]
BDGMAP2.2020$brks20.2 <- cut(BDGMAP2.2020$Exceed2,breaks = brks20.2,include.lowest = TRUE,labels = labels20.2)

brks_scale20.2 <- levels(BDGMAP2.2020$brks20.2)
labels_scale20.2 <- rev(brks_scale20.2)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2020` = "2020")

BDGMAP2.2020 <- na.omit(BDGMAP2.2020) 
BDGMAP2.2020$TahunBulan <- with(BDGMAP2.2020, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP2.2020, aes(fill = brks20.2, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2020)

# - Tahun 2021
BDGMAP2.2021<- data.shp2.2[data.shp2.2$Tahun==2021,]

pretty_breaks21.2 <- c(0.95)

labels21.2 <- c()
brks21.2 <- c(-1, pretty_breaks21.2,9.5)
for(idx in 1:length(brks21.2)){
  labels21.2 <- c(labels21.2,round(brks21.2[idx + 1], 2)) }

labels21.2 <- labels21.2[1:length(labels21.2)-1]
BDGMAP2.2021$brks21.2 <- cut(BDGMAP2.2021$Exceed2,breaks = brks21.2,include.lowest = TRUE,labels = labels21.2)

brks_scale21.2 <- levels(BDGMAP2.2021$brks21.2)
labels_scale21.2 <- rev(brks_scale21.2)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2021` = "2021")

BDGMAP2.2021 <- na.omit(BDGMAP2.2021) 
BDGMAP2.2021$TahunBulan <- with(BDGMAP2.2021, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP2.2021, aes(fill = brks21.2, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2021)

# - Tahun 2022
BDGMAP2.2022<- data.shp2.2[data.shp2.2$Tahun==2022,]

pretty_breaks22.2 <- c(0.95)

labels22.2 <- c()
brks22.2 <- c(-1, pretty_breaks22.2,9.5)
for(idx in 1:length(brks22.2)){
  labels22.2 <- c(labels22.2,round(brks22.2[idx + 1], 2)) }

labels22.2 <- labels22.2[1:length(labels22.2)-1]
BDGMAP2.2022$brks22.2 <- cut(BDGMAP2.2022$Exceed2,breaks = brks22.2,include.lowest = TRUE,labels = labels22.2)

brks_scale22.2 <- levels(BDGMAP2.2022$brks22.2)
labels_scale22.2 <- rev(brks_scale22.2)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2022` = "2022")

BDGMAP2.2022 <- na.omit(BDGMAP2.2022) 
BDGMAP2.2022$TahunBulan <- with(BDGMAP2.2022, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP2.2022, aes(fill = brks22.2, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2022)

# - Tahun 2023
BDGMAP2.2023<- data.shp2.2[data.shp2.2$Tahun==2023,]

pretty_breaks23.2 <- c(0.95)

labels23.2 <- c()
brks23.2 <- c(-1, pretty_breaks23.2,9.5)
for(idx in 1:length(brks23.2)){
  labels23.2 <- c(labels23.2,round(brks23.2[idx + 1], 2)) }

labels23.2 <- labels23.2[1:length(labels23.2)-1]
BDGMAP2.2023$brks23.2 <- cut(BDGMAP2.2023$Exceed2,breaks = brks23.2,include.lowest = TRUE,labels = labels23.2)

brks_scale23.2 <- levels(BDGMAP2.2023$brks23.2)
labels_scale23.2 <- rev(brks_scale23.2)

Bulan <- c(
  `1` = "Januari",`2` = "Februari",`3` = "Maret",`4` = "April",`5` = "Mei",`6` = "Juni",
  `7` = "Juli",`8` = "Agustus",`9` = "September",`10` = "Oktober",`11` = "November",`12` = "Desember")

Tahun <- c(`2023` = "2023")

BDGMAP2.2023 <- na.omit(BDGMAP2.2023) 
BDGMAP2.2023$TahunBulan <- with(BDGMAP2.2023, paste(Tahun, Bulan, sep = "-"))

ggplot() +
  geom_polygon(data = BDGMAP2.2023, aes(fill = brks23.2, x = long,  y = lat, group = group), color="black", linewidth=0.01) +
  facet_wrap(~ Bulan, ncol = 6, labeller = labeller(Bulan = as_labeller(Bulan))) +
  theme_bw() + 
  ylab("Latitude") + 
  xlab("Longitude") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("Bukan Hotspot", "Hotspot"), values =  c("green","red")) +
  labs(fill = expression(paste("Penyebaran Hotspot"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ggtitle(2023)
