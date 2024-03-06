#####################################################
################## Niche Overlap#####################
#####################################################
### Install and activate libraries.
library(ecospat)    ###Niche comparison in environmental space
library(raster)     ###Raster manipulation
library(SDMTools)   ###Niche-related analyses
library(dismo)      ###Distribution modeling
library(viridis)
library(ggplot2)
library(gridExtra)

### We have six genetic lineages to be compared.
sp1 <- read.csv("CRP/CRP_joint.csv", h=T)
sp2 <- read.csv("NCA/NCA_joint.csv", h=T)
sp3 <- read.csv("NChi/NChi_joint.csv", h=T)
sp4 <- read.csv("SMO/SMO_joint.csv", h=T)
sp5 <- read.csv("SMS/SMS_joint.csv", h=T)
sp6 <- read.csv("Tux/Tux_joint.csv", h=T)

### Load previously cropped layers that are relevant for comparison of ecological niches to generate the principal components.
### For 6 genetic lineages to be compared.
varclim1 <- list.files("E:/chlorospingus/CRP/M_variables/Set_1/", pattern=".asc", full.names=T)
varclim2 <- list.files("E:/chlorospingus/NCA/M_variables/Set_1/", pattern=".asc", full.names=T)
varclim3 <- list.files("E:/chlorospingus/NChi/M_variables/Set_1/", pattern=".asc", full.names=T)
varclim4 <- list.files("E:/chlorospingus/SMO/M_variables/Set_1/", pattern=".asc", full.names=T)
varclim5 <- list.files("E:/chlorospingus/SMS/M_variables/Set_1/", pattern=".asc", full.names=T)
varclim6 <- list.files("E:/chlorospingus/Tux/M_variables/Set_1/", pattern=".asc", full.names=T)

varclim1 <- stack(varclim1)
varclim2 <- stack(varclim2)
varclim3 <- stack(varclim3)
varclim4 <- stack(varclim4)
varclim5 <- stack(varclim5)
varclim6 <- stack(varclim6)

climpunto1 <- rasterToPoints(varclim1[[1]], fun=NULL, spatial=TRUE)
climpunto2 <- rasterToPoints(varclim2[[1]], fun=NULL, spatial=TRUE)
climpunto3 <- rasterToPoints(varclim3[[1]], fun=NULL, spatial=TRUE)
climpunto4 <- rasterToPoints(varclim4[[1]], fun=NULL, spatial=TRUE)
climpunto5 <- rasterToPoints(varclim5[[1]], fun=NULL, spatial=TRUE)
climpunto6 <- rasterToPoints(varclim6[[1]], fun=NULL, spatial=TRUE)

### Extraction of environmental data for each layer in capas of each coordinate in varclim.
clim1 <- extract(varclim1, climpunto1)
clim2 <- extract(varclim2, climpunto2)
clim3 <- extract(varclim3, climpunto3)
clim4 <- extract(varclim4, climpunto4)
clim5 <- extract(varclim5, climpunto5)
clim6 <- extract(varclim6, climpunto6)

### Format clim to be a normal table with two columns x and y.
clim1 <- data.frame(coordinates(climpunto1),clim1)
clim2 <- data.frame(coordinates(climpunto2),clim2)
clim3 <- data.frame(coordinates(climpunto3),clim3)
clim4 <- data.frame(coordinates(climpunto4),clim4)
clim5 <- data.frame(coordinates(climpunto5),clim5)
clim6 <- data.frame(coordinates(climpunto6),clim6)

##### IM-C: Quitamos filas con NAs.
clim1 <- na.omit(clim1)
clim2 <- na.omit(clim2)
clim3 <- na.omit(clim3)
clim4 <- na.omit(clim4)
clim5 <- na.omit(clim5)
clim6 <- na.omit(clim6)

### Add climatic variables to data.
occ_sp1 <- na.exclude(ecospat.sample.envar(dfsp=sp1,colspxy=2:3, colspkept=2:3,dfvar=clim1, colvarxy=1:2,colvar="all",resolution=0.1))
occ_sp2 <- na.exclude(ecospat.sample.envar(dfsp=sp2,colspxy=2:3, colspkept=2:3,dfvar=clim2, colvarxy=1:2,colvar="all",resolution=0.1))
occ_sp3 <- na.exclude(ecospat.sample.envar(dfsp=sp3,colspxy=2:3, colspkept=2:3,dfvar=clim3, colvarxy=1:2,colvar="all",resolution=0.1))
occ_sp4 <- na.exclude(ecospat.sample.envar(dfsp=sp4,colspxy=2:3, colspkept=2:3,dfvar=clim4, colvarxy=1:2,colvar="all",resolution=0.1))
occ_sp5 <- na.exclude(ecospat.sample.envar(dfsp=sp5,colspxy=2:3, colspkept=2:3,dfvar=clim5, colvarxy=1:2,colvar="all",resolution=0.1))
occ_sp6 <- na.exclude(ecospat.sample.envar(dfsp=sp6,colspxy=2:3, colspkept=2:3,dfvar=clim6, colvarxy=1:2,colvar="all",resolution=0.1))

#######################################################################
#########################   Environmental PCA   #######################
#######################################################################
### Data for PCA to include all sites for the study area and for species presences.
### First clim and then the species in order.
data <- rbind(clim1[,3:9],clim2[,3:9],clim3[,3:9],clim4[,3:9],clim5[,3:9],clim6[,3:9],occ_sp1[,3:9],occ_sp2[,3:9],occ_sp3[,3:9],occ_sp4[,3:9],occ_sp5[,3:9],occ_sp6[,3:9])
data <- na.omit(data)

##### IM-C: Hago un rbind.
clim123456 <- rbind(clim1,clim2,clim3,clim4,clim5,clim6)

##### IM-C: Version modificada para niceoverplot.
row.w.6.env<-1-(nrow(clim6)/nrow(clim123456))
row.w.5.env<-1-(nrow(clim5)/nrow(clim123456))
row.w.4.env<-1-(nrow(clim4)/nrow(clim123456))
row.w.3.env<-1-(nrow(clim3)/nrow(clim123456))
row.w.2.env<-1-(nrow(clim2)/nrow(clim123456))
row.w.1.env<-1-(nrow(clim1)/nrow(clim123456))

### Vector with 0 weight for occurrences and 1 for sites of the study area.
row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(row.w.3.env, nrow(clim3)),rep(row.w.4.env, nrow(clim4)),rep(row.w.5.env, nrow(clim5)),rep(row.w.6.env, nrow(clim6)),rep(0, nrow(occ_sp1)),rep(0, nrow(occ_sp2)),rep(0, nrow(occ_sp3)),rep(0, nrow(occ_sp4)),rep(0, nrow(occ_sp5)),rep(0, nrow(occ_sp6)))

### PCA is done with all data from the study area.
### Presences are not used for calibration of PCA but coordinates are calculated. NUMBER OF AXIS=2.
pca.cal <-dudi.pca(na.omit(data, row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2))

row.sp1<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2))
row.sp3<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3))
row.sp4<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4))
row.sp5<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4)+nrow(occ_sp5))
row.sp6<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4)+nrow(occ_sp5)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4)+nrow(occ_sp5)+nrow(occ_sp6))

row.clim1 <- 1:nrow(clim1)
row.clim2<-(nrow(clim2)+1):(nrow(clim1)+nrow(clim2))
row.clim3<-(nrow(clim3)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3))
row.clim4<-(nrow(clim4)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4))
row.clim5<-( nrow(clim5)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5))
row.clim6<-( nrow(clim6)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6))

row.clim123456<-1:(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6))
row.sp1<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2))
row.sp3<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3))
row.sp4<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4))
row.sp5<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4)+nrow(occ_sp5))
row.sp6<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4)+nrow(occ_sp5)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(clim6)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4)+nrow(occ_sp5)+nrow(occ_sp6))

### Having the PCA results, we need the first and second eigenvector values for the background and the occurrence records per species.
### Coordinates in each axis of PCA of all the study area and each species.
scores.clim123456 <- pca.cal$li[row.clim123456,] 
scores.clim1<- pca.cal$li[row.clim1,]
scores.clim2<- pca.cal$li[row.clim2,]
scores.clim3<- pca.cal$li[row.clim3,]
scores.clim4<- pca.cal$li[row.clim4,]
scores.clim5<- pca.cal$li[row.clim5,]
scores.clim6<- pca.cal$li[row.clim6,]

scores.sp1<- pca.cal$li[row.sp1,]
scores.sp2<- pca.cal$li[row.sp2,]
scores.sp3<- pca.cal$li[row.sp3,]
scores.sp4<- pca.cal$li[row.sp4,]
scores.sp5<- pca.cal$li[row.sp5,]
scores.sp6<- pca.cal$li[row.sp6,]

##### IM-C: Contribucion de cada variable a cada componente del PCA.
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)

### Environmental space.
### Resolution of the environmental space based on the PCA values calculated for the background and occurrence records.
R <- 100 ### This is the same as using R=100 in the lines below.
### Create density surfaces of occurrence in the environmental space (two axis), considering the observed occurrence density and availability of condition in the background.
z6_grid <- ecospat.grid.clim.dyn (scores.clim123456, scores.clim6, scores.sp6, R=100) 
z5_grid <- ecospat.grid.clim.dyn (scores.clim123456, scores.clim5, scores.sp5, R=100) 
z4_grid <- ecospat.grid.clim.dyn (scores.clim123456, scores.clim4, scores.sp4, R=100) 
z3_grid <- ecospat.grid.clim.dyn (scores.clim123456, scores.clim3, scores.sp3, R=100) 
z2_grid <- ecospat.grid.clim.dyn (scores.clim123456, scores.clim2, scores.sp2, R=100) 
z1_grid <- ecospat.grid.clim.dyn (scores.clim123456, scores.clim1, scores.sp1, R=100) 

### Metrics for niche overlap. Calculation of D metric and its significance through a similarity test. We can define number of iterations for the test as rep <- 100 (which is the same as rep=100 in the lines below).
### Generation of values of niche overlap.
ecospat.niche.overlap (z1=z1_grid, z2=z2_grid, cor=TRUE) #crp and nca
ecospat.niche.overlap (z1=z1_grid, z2=z3_grid, cor=TRUE) #crp and nchi  
ecospat.niche.overlap (z1=z1_grid, z2=z4_grid, cor=TRUE) #crp and smo  
ecospat.niche.overlap (z1=z1_grid, z2=z5_grid, cor=TRUE) #crp and sms
ecospat.niche.overlap (z1=z1_grid, z2=z6_grid, cor=TRUE) #crp and tux  
ecospat.niche.overlap (z1=z2_grid, z2=z3_grid, cor=TRUE) #nca and nchi  
ecospat.niche.overlap (z1=z2_grid, z2=z4_grid, cor=TRUE) #nca and smo  
ecospat.niche.overlap (z1=z2_grid, z2=z5_grid, cor=TRUE) #nca and sms
ecospat.niche.overlap (z1=z2_grid, z2=z6_grid, cor=TRUE) #nca and tux  
ecospat.niche.overlap (z1=z3_grid, z2=z4_grid, cor=TRUE) #nchi and smo  
ecospat.niche.overlap (z1=z3_grid, z2=z5_grid, cor=TRUE) #nchi and sms  
ecospat.niche.overlap (z1=z3_grid, z2=z6_grid, cor=TRUE) #nchi and tux  
ecospat.niche.overlap (z1=z4_grid, z2=z5_grid, cor=TRUE) #smo and sms
ecospat.niche.overlap (z1=z4_grid, z2=z6_grid, cor=TRUE) #smo and tux
ecospat.niche.overlap (z1=z5_grid, z2=z6_grid, cor=TRUE) #sms and tux

ecospat.niche.overlap (z1=z2_grid, z2=z1_grid, cor=TRUE) #nca and crp
ecospat.niche.overlap (z1=z3_grid, z2=z1_grid, cor=TRUE) #nchi and crp  
ecospat.niche.overlap (z1=z4_grid, z2=z1_grid, cor=TRUE) #smo and crp  
ecospat.niche.overlap (z1=z5_grid, z2=z1_grid, cor=TRUE) #sms and crp
ecospat.niche.overlap (z1=z6_grid, z2=z1_grid, cor=TRUE) #tux and crp  
ecospat.niche.overlap (z1=z3_grid, z2=z2_grid, cor=TRUE) #nchi and nca  
ecospat.niche.overlap (z1=z4_grid, z2=z2_grid, cor=TRUE) #smo and nca  
ecospat.niche.overlap (z1=z5_grid, z2=z2_grid, cor=TRUE) #sms and nca
ecospat.niche.overlap (z1=z6_grid, z2=z2_grid, cor=TRUE) #tux and nca  
ecospat.niche.overlap (z1=z4_grid, z2=z3_grid, cor=TRUE) #smo and nchi  
ecospat.niche.overlap (z1=z5_grid, z2=z3_grid, cor=TRUE) #sms and nchi  
ecospat.niche.overlap (z1=z6_grid, z2=z3_grid, cor=TRUE) #tux and nchi  
ecospat.niche.overlap (z1=z5_grid, z2=z4_grid, cor=TRUE) #sms and smo
ecospat.niche.overlap (z1=z6_grid, z2=z4_grid, cor=TRUE) #tux and smo
ecospat.niche.overlap (z1=z6_grid, z2=z5_grid, cor=TRUE) #tux and sms


### Occurrence density graphics in the environmental space.
windows() 
par(mfrow=c(2,3))
ecospat.plot.niche (z1_grid, title="CRP", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z2_grid, title="NCA", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z3_grid, title="NChi", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z4_grid, title="SMO", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z5_grid, title="SMS", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z6_grid, title="Tux", name.axis1="PC1", name.axis2="PC2", cor=F)

### Individual niche plots. Modifications of plot.niche functions from Broenniman et al. 2012 and Silva et al. 2014.
n.groups <- 6
g.names <- c("CRP", "NCA", "NChi", "SMO", "SMS", "Tux")
g.colors <- c("cyan", "orange", "yellow3", "blue", "red", "green")
z <- c(list(z1_grid, z2_grid, z3_grid,z4_grid,z5_grid,z6_grid))
plot.niche.mod <- function(z, name.axis1 = "PC1", name.axis2 = "PC2",
                           cor = F, corte,  contornar = TRUE, 
                           densidade = TRUE, quantis = 10, 
                           back = TRUE, x = "red", title = "",
                           i) {  
  cor1 <- function(cores.i, n) {
    al <- seq(0,1,(1/n))
    cores <- numeric(length(n))
    for(i in 1:n) {    
      corespar <- col2rgb(cores.i)/255
      cores[i] <- rgb(corespar[1, ], corespar[2, ],
                      corespar[3, ], alpha = al[i])
    }
    return(cores)
  }
  a1 <- colorRampPalette(c("transparent",cor1(x, quantis)), alpha = TRUE)  
  xlim <- c(min(sapply(z, function(x){min(x$x)})),
            max(sapply(z, function(x){max(x$x)})))
  ylim <- c(min(sapply(z, function(x){min(x$y)})),
            max(sapply(z, function(x){max(x$y)})))
  graphics::image(z[[1]]$x, z[[1]]$y, as.matrix(z[[1]]$z.uncor), col = "white", 
                  ylim = ylim, xlim = xlim,
                  zlim = c(0.000001, max(as.matrix(z[[1]]$z.uncor), na.rm = T)), 
                  xlab = "PC1", ylab = "PC2", cex.lab = 1.5,
                  cex.axis = 1.4)
  abline(h = 0, v = 0, lty = 2)
  if (back) {
    contour(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$Z),
            add = TRUE, levels = quantile(z[[i]]$Z[z[[i]]$Z > 0],
                                          c(0, 0.5)), drawlabels = FALSE,
            lty = c(1, 2), col = x, lwd = 1)
  } 
  if (densidade) {
    image(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$z.uncor), col = a1(100), add = TRUE)
  }
  if(contornar){
    contour(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$z.uncor), 
            add = TRUE, levels = quantile(z[[i]]$z.uncor[z[[i]]$z.uncor > 0],
                                          seq(0, 1, (1 / quantis))),
            drawlabels = FALSE, lty = c(rep(2,(quantis - 1)), 1), 
            col = cor1(x, quantis), lwd = c(rep(1, (quantis - 1)), 2))
  }
  title(title)
  box()
}

### Continuos line represent 100% of the available environmental background and the dashed lines represent 50% of the most common conditions.
for(i in 1:n.groups) {
  plot.niche.mod(z, name.axis1 = "PC1", name.axis2 = "PC2",
                 cor = F, corte,  contornar = FALSE, 
                 densidade = TRUE, quantis = 6, 
                 back = TRUE, x = g.colors[i], title = g.names[i], i)
}

### Multiple niche plots. Function modification to allow multiple regions/species plots.
windows() 
par(mfrow=c(1,1))
plot.niche.all <- function(z, n.groups, g.names,
                           contornar = TRUE, 
                           densidade = TRUE,
                           quantis = 10,
                           back = TRUE, title = "",
                           g.colors, n = 5,
                           cor1) { 
cor1 <- function(cores.i, n) {
    al <- seq(0,1,(1/n))
    cores <- numeric(length(n))
    for(i in 1:n) {    
      corespar <- col2rgb(cores.i)/255
      cores[i] <- rgb(corespar[1, ], corespar[2, ],
                      corespar[3, ], alpha = al[i])
    }
    return(cores)
  }
  a <- list() 
  for(i in 1:n.groups) {
    a[[i]] <- colorRampPalette(c("transparent", cor1(g.colors[i], n)),
                               alpha = TRUE)  
  }
  xlim <- c(min(sapply(z, function(x){min(x$x)})),
            max(sapply(z, function(x){max(x$x)})))
  ylim <- c(min(sapply(z, function(x){min(x$y)})),
            max(sapply(z, function(x){max(x$y)})))
  image(z[[1]]$x, z[[1]]$y, as.matrix(z[[1]]$z.uncor), col = "white", 
        ylim = ylim, xlim = xlim,
        zlim = c(0.000001, max(as.matrix(z[[1]]$Z), na.rm = T)), 
        xlab = "PC1", ylab = "PC2", cex.lab = 1.5,
        cex.axis = 1.4)
  abline(h = 0, v = 0, lty = 2)
  box()
  if (back) {
    for(i in 1:n.groups) {
      contour(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$Z), add = TRUE,
              levels = quantile(z[[i]]$Z[z[[i]]$Z > 0], c(0, 1)),
              drawlabels = FALSE,lty = c(1, 2),
              col = g.colors[i], lwd = 1)
    }
  }
  if (densidade) {
    for(i in 1:n.groups) {
      image(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$z.uncor),
            col = a[[i]](100), add = TRUE)
    }
  }
  if(contornar){
    for(i in 1:n.groups) {
      contour(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$z.uncor), add = TRUE,
              levels = quantile(z[[i]]$z.uncor[z[[i]]$z.uncor > 0],
                                seq(0, 1, (1/quantis)))[quantis],
              drawlabels = FALSE, lty = rev(c(rep(2,(quantis - 1)), 1)),
              col = rev(cor1(g.colors[i], quantis)),
              lwd = rev(c(rep(1, (quantis - 1)), 2)))
    }
  }
}

### Strong contours represent 5% highest values of density and the thin lines represent 100% of the background available region.
plot.niche.all(z, n.groups, g.names,
               contornar = TRUE, 
               densidade = TRUE,
               quantis = 10,
               back = FALSE, title = "",
               g.colors, n = 2,
               cor1)
plot.niche.all(z, n.groups, g.names,
               contornar = TRUE, 
               densidade = TRUE,
               quantis = 10,
               back = TRUE, title = "",
               g.colors, n = 2,
               cor1)


#######################################################################
######################   Niche Equivalency Test    ####################
#######################################################################
### Are two niches identical? 
### Null hypothesis is that niches are identical and you reject it if the observed value is lower than expected under a null model.
a.dyn1<-ecospat.niche.equivalency.test(z1=z1_grid, z2=z2_grid, rep=100)
ecospat.plot.overlap.test(a.dyn1,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn1,"I","Equivalency")
a.dyn2<-ecospat.niche.equivalency.test(z1=z1_grid, z2=z3_grid, rep=100)
a.dyn3<-ecospat.niche.equivalency.test(z1=z1_grid, z2=z4_grid, rep=100)
a.dyn4<-ecospat.niche.equivalency.test(z1=z1_grid, z2=z5_grid, rep=100)
a.dyn5<-ecospat.niche.equivalency.test(z1=z1_grid, z2=z6_grid, rep=100)
a.dyn6<-ecospat.niche.equivalency.test(z1=z2_grid, z2=z3_grid, rep=100)
a.dyn7<-ecospat.niche.equivalency.test(z1=z2_grid, z2=z4_grid, rep=100)
a.dyn8<-ecospat.niche.equivalency.test(z1=z2_grid, z2=z5_grid, rep=100)
a.dyn9<-ecospat.niche.equivalency.test(z1=z2_grid, z2=z6_grid, rep=100)
a.dyn10<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z4_grid, rep=100)
a.dyn11<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z5_grid, rep=100)
a.dyn12<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z6_grid, rep=100)
a.dyn13<-ecospat.niche.equivalency.test(z1=z4_grid, z2=z5_grid, rep=100)
a.dyn14<-ecospat.niche.equivalency.test(z1=z4_grid, z2=z6_grid, rep=100)
a.dyn15<-ecospat.niche.equivalency.test(z1=z5_grid, z2=z6_grid, rep=100)

a.dyn16<-ecospat.niche.equivalency.test(z1=z2_grid, z2=z1_grid, rep=100)
a.dyn17<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z1_grid, rep=100)
a.dyn18<-ecospat.niche.equivalency.test(z1=z4_grid, z2=z1_grid, rep=100)
a.dyn19<-ecospat.niche.equivalency.test(z1=z5_grid, z2=z1_grid, rep=100)
a.dyn20<-ecospat.niche.equivalency.test(z1=z6_grid, z2=z1_grid, rep=100)
a.dyn21<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z2_grid, rep=100)
a.dyn22<-ecospat.niche.equivalency.test(z1=z4_grid, z2=z2_grid, rep=100)
a.dyn23<-ecospat.niche.equivalency.test(z1=z5_grid, z2=z2_grid, rep=100)
a.dyn24<-ecospat.niche.equivalency.test(z1=z6_grid, z2=z2_grid, rep=100)
a.dyn25<-ecospat.niche.equivalency.test(z1=z4_grid, z2=z3_grid, rep=100)
a.dyn26<-ecospat.niche.equivalency.test(z1=z5_grid, z2=z3_grid, rep=100)
a.dyn27<-ecospat.niche.equivalency.test(z1=z6_grid, z2=z3_grid, rep=100)
a.dyn28<-ecospat.niche.equivalency.test(z1=z5_grid, z2=z4_grid, rep=100)
a.dyn29<-ecospat.niche.equivalency.test(z1=z6_grid, z2=z4_grid, rep=100)
a.dyn30<-ecospat.niche.equivalency.test(z1=z6_grid, z2=z5_grid, rep=100)
