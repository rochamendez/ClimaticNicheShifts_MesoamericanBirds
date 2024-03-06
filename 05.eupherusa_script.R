###########################################
########     Script for ENM        ########
###########################################
### In case some of the code doesn’t work, check the quotation marks, they may be off or not adequate.
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

### Read csv archive of distributional data.
### We have five genetic lineages to be compared.
sp1 <- read.csv("cyanophrys/Cyanoprhys_Limpio.csv", h=T)
sp2 <- read.csv("eximia/Eximia_Limpio.csv", h=T)
sp3 <- read.csv("nigriventris/nig_joint.csv", h=T)
sp4 <- read.csv("poliocerca/Poliocerca_Limpia.csv", h=T)
sp5 <- read.csv("ridgwayi/Thalurania_Limpio.csv", h=T)

### Load previously cropped layers that are relevant for comparison of ecological niches to generate the principal components.
### For 5genetic lineages to be compared.
varclim1 <- list.files("E:/cyanophrys/M_variables/Set_1/", pattern=".asc", full.names=T)
varclim2 <- list.files("E:/eximia/M_variables/Set_1/", pattern=".asc", full.names=T)
varclim3 <- list.files("E:/nigriventris/M_variables/Set_1/", pattern=".asc", full.names=T)
varclim4 <- list.files("E:/poliocerca/M_variables/Set_1/", pattern=".asc", full.names=T)
varclim5 <- list.files("E:/ridgwayi/M_variables/Set_1/", pattern=".asc", full.names=T)

varclim1 <- stack(varclim1)
varclim2 <- stack(varclim2)
varclim3 <- stack(varclim3)
varclim4 <- stack(varclim4)
varclim5 <- stack(varclim5)

### Create a table for each pixel (x, y coordinates).
climpunto1 <- rasterToPoints(varclim1[[1]], fun=NULL, spatial=TRUE)
climpunto2 <- rasterToPoints(varclim2[[1]], fun=NULL, spatial=TRUE)
climpunto3 <- rasterToPoints(varclim3[[1]], fun=NULL, spatial=TRUE)
climpunto4 <- rasterToPoints(varclim4[[1]], fun=NULL, spatial=TRUE)
climpunto5 <- rasterToPoints(varclim5[[1]], fun=NULL, spatial=TRUE)

### Extraction of environmental data for each layer in capas of each coordinate in varclim.
clim1 <- extract(varclim1, climpunto1)
clim2 <- extract(varclim2, climpunto2)
clim3 <- extract(varclim3, climpunto3)
clim4 <- extract(varclim4, climpunto4)
clim5 <- extract(varclim5, climpunto5)

### Format clim to be a normal table with two columns x and y.
clim1 <- data.frame(coordinates(climpunto1),clim1)
clim2 <- data.frame(coordinates(climpunto2),clim2)
clim3 <- data.frame(coordinates(climpunto3),clim3)
clim4 <- data.frame(coordinates(climpunto4),clim4)
clim5 <- data.frame(coordinates(climpunto5),clim5)

##### IM-C: Quitamos filas con NAs.
clim1 <- na.omit(clim1)
clim2 <- na.omit(clim2)
clim3 <- na.omit(clim3)
clim4 <- na.omit(clim4)
clim5 <- na.omit(clim5)

### Add climatic variables to data.
occ_sp1 <- na.exclude(ecospat.sample.envar(dfsp=sp1,colspxy=2:3, colspkept=2:3,dfvar=clim1, colvarxy=1:2,colvar="all",resolution=0.1))
occ_sp2 <- na.exclude(ecospat.sample.envar(dfsp=sp2,colspxy=2:3, colspkept=2:3,dfvar=clim2, colvarxy=1:2,colvar="all",resolution=0.1))
occ_sp3 <- na.exclude(ecospat.sample.envar(dfsp=sp3,colspxy=2:3, colspkept=2:3,dfvar=clim3, colvarxy=1:2,colvar="all",resolution=0.1))
occ_sp4 <- na.exclude(ecospat.sample.envar(dfsp=sp4,colspxy=2:3, colspkept=2:3,dfvar=clim4, colvarxy=1:2,colvar="all",resolution=0.1))
occ_sp5 <- na.exclude(ecospat.sample.envar(dfsp=sp5,colspxy=2:3, colspkept=2:3,dfvar=clim5, colvarxy=1:2,colvar="all",resolution=0.1))

#######################################################################
#########################   Environmental PCA   #######################
#######################################################################
### Data for PCA to include all sites for the study area and for species presences
### First clim and then the species in order
### For 5 genetic lineages to be compared. 
data <- rbind(clim1[,3:8],clim2[,3:8],clim3[,3:8],clim4[,3:8],clim5[,3:8],occ_sp1[,3:8],occ_sp2[,3:8],occ_sp3[,3:8],occ_sp4[,3:8],occ_sp5[,3:8])
data <- na.omit(data)

##### IM-C: Hago un rbind.
### For 5 genetic lineages to be compared
clim12345 <- rbind(clim1,clim2,clim3,clim4,clim5)

##### IM-C: Version modificada para niceoverplot.
row.w.1.env<-1-(nrow(clim1)/nrow(clim12345))  # prevalence of clim1
row.w.2.env<-1-(nrow(clim2)/nrow(clim12345))  # prevalence of clim2
row.w.3.env<-1-(nrow(clim3)/nrow(clim12345))  # prevalence of clim3
row.w.4.env<-1-(nrow(clim4)/nrow(clim12345))  # prevalence of clim4
row.w.5.env<-1-(nrow(clim5)/nrow(clim12345))  # prevalence of clim5

### Vector with 0 weight for occurrences and 1 for sites of the study area.
row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(row.w.3.env, nrow(clim3)),rep(row.w.4.env, nrow(clim4)),rep(row.w.5.env, nrow(clim5)),rep(0, nrow(occ_sp1)),rep(0, nrow(occ_sp2)),rep(0, nrow(occ_sp3)),rep(0, nrow(occ_sp4)),rep(0, nrow(occ_sp5)))

### PCA is done with all data from the study area.
### Presences are not used for calibration of PCA but coordinates are calculated. NUMBER OF AXIS=2.
pca.cal <-dudi.pca(na.omit(data, row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2))

### For 5 genetic lineages to be compared.
row.clim1 <- 1:nrow(clim1)
row.clim2<-(nrow(clim2)+1):(nrow(clim1)+nrow(clim2))
row.clim3<-(nrow(clim3)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3))
row.clim4<-(nrow(clim4)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4))
row.clim5<-(nrow(clim5)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5))
row.clim12345<-1:(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5))

### For 5 genetic lineages to be compared.
row.sp1<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(occ_sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(occ_sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(occ_sp1)+nrow(occ_sp2))
row.sp3<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(occ_sp1)+nrow(occ_sp2)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3))
row.sp4<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4))
row.sp5<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(clim4)+nrow(clim5)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3)+nrow(occ_sp4)+nrow(occ_sp5))

### Having the PCA results, we need the first and second eigenvector values for the background and the occurrence records per species.
### Coordinates in each axis of PCA of all the study area and each species.
scores.clim12345 <- pca.cal$li[row.clim12345,] 
scores.clim1<- pca.cal$li[row.clim1,]
scores.clim2<- pca.cal$li[row.clim2,]
scores.clim3<- pca.cal$li[row.clim3,]
scores.clim4<- pca.cal$li[row.clim4,]
scores.clim5<- pca.cal$li[row.clim5,]
scores.sp1<- pca.cal$li[row.sp1,]
scores.sp2<- pca.cal$li[row.sp2,]
scores.sp3<- pca.cal$li[row.sp3,]
scores.sp4<- pca.cal$li[row.sp4,]
scores.sp5<- pca.cal$li[row.sp5,]

##### IM-C: Contribucion de cada variable a cada componente del PCA.
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)

### Environmental space.
### Resolution of the environmental space based on the PCA values calculated for the background and occurrence records.
### Create density surfaces of occurrence in the environmental space (two axis), considering the observed occurrence density and availability of condition in the background.
z1_grid <-ecospat.grid.clim.dyn (scores.clim12345, scores.clim1, scores.sp1, R=100) 
z2_grid <- ecospat.grid.clim.dyn (scores.clim12345, scores.clim2, scores.sp2, R=100) 
z3_grid <- ecospat.grid.clim.dyn (scores.clim12345, scores.clim3, scores.sp3, R=100) 
z4_grid <- ecospat.grid.clim.dyn (scores.clim12345, scores.clim4, scores.sp4, R=100) 
z5_grid <- ecospat.grid.clim.dyn (scores.clim12345, scores.clim5, scores.sp5, R=100) 

### Metrics for niche overlap. Calculation of D metric and its significance through a similarity test. We can define number of iterations for the test as rep <- 100 (which is the same as rep=100 in the lines below).
### Generation of values of niche overlap.
ecospat.niche.overlap (z1=z1_grid, z2=z2_grid, cor=TRUE) #cyano and eximia
ecospat.niche.overlap (z1=z1_grid, z2=z3_grid, cor=TRUE) #cyano and nigriventris  
ecospat.niche.overlap (z1=z1_grid, z2=z4_grid, cor=TRUE) #cyano and poliocerca  
ecospat.niche.overlap (z1=z1_grid, z2=z5_grid, cor=TRUE) #cyano and ridg
ecospat.niche.overlap (z1=z2_grid, z2=z3_grid, cor=TRUE) #eximia and nigri  
ecospat.niche.overlap (z1=z2_grid, z2=z4_grid, cor=TRUE) #eximia and poliocerca  
ecospat.niche.overlap (z1=z2_grid, z2=z5_grid, cor=TRUE) #Eximia and ridg
ecospat.niche.overlap (z1=z3_grid, z2=z4_grid, cor=TRUE) #nigri and polio  
ecospat.niche.overlap (z1=z3_grid, z2=z5_grid, cor=TRUE) #nigri and ridg  
ecospat.niche.overlap (z1=z4_grid, z2=z5_grid, cor=TRUE) #polio and ridg  

ecospat.niche.overlap (z1=z2_grid, z2=z1_grid, cor=TRUE) #eximia and cyano
ecospat.niche.overlap (z1=z3_grid, z2=z1_grid, cor=TRUE) #nigriventris and cyano
ecospat.niche.overlap (z1=z4_grid, z2=z1_grid, cor=TRUE) #poliocerca and cyano
ecospat.niche.overlap (z1=z5_grid, z2=z1_grid, cor=TRUE) #ridg and cyano
ecospat.niche.overlap (z1=z3_grid, z2=z2_grid, cor=TRUE) #nigri and eximia  
ecospat.niche.overlap (z1=z4_grid, z2=z2_grid, cor=TRUE) #poliocerca and eximia  
ecospat.niche.overlap (z1=z5_grid, z2=z2_grid, cor=TRUE) #ridg and eximia
ecospat.niche.overlap (z1=z4_grid, z2=z3_grid, cor=TRUE) #polio and nigri  
ecospat.niche.overlap (z1=z5_grid, z2=z3_grid, cor=TRUE) #ridg and nigri  
ecospat.niche.overlap (z1=z5_grid, z2=z4_grid, cor=TRUE) #ridg and polio  

### Occurrence density graphics in the environmental space.
windows() 
par(mfrow=c(1,3))
ecospat.plot.niche (z1_grid, title="CRP", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z2_grid, title="NCA", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z3_grid, title="NChi", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z4_grid, title="SMO", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z5_grid, title="SMS", name.axis1="PC1", name.axis2="PC2", cor=F)

### Individual niche plots. Modifications of plot.niche functions from Broenniman et al. 2012 and Silva et al. 2014.
n.groups <- 5
g.names <- c("E. cyanophrys", "E. eximia", "E. nigriventris", “E. poliocerca”, “E. ridgwayi”)
g.colors <- c("cyan", "darkblue", "red", "green", "orange")
z <- c(list(z1_grid, z2_grid, z3_grid,z4_grid,z5_grid))

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
### Continuos line represent 100% of the available environmental background and the dashed lines represent 50% of the most common conditions
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
###################   Niche Equivalency Test    #######################
#######################################################################
### Are two niches identical? 
### Null hypothesis is that niches are identical and you reject it if the observed value is lower than expected under a null model.
a.dyn1<-ecospat.niche.equivalency.test(z1=z1_grid, z2=z2_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn1,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn1,"I","Equivalency")
a.dyn2<-ecospat.niche.equivalency.test(z1=z1_grid, z2=z3_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn2,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn2,"I","Equivalency")
a.dyn3<-ecospat.niche.equivalency.test(z1=z1_grid, z2=z4_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn3,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn3,"I","Equivalency")
a.dyn4<-ecospat.niche.equivalency.test(z1=z1_grid, z2=z5_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn4,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn4,"I","Equivalency")
a.dyn5<-ecospat.niche.equivalency.test(z1=z2_grid, z2=z3_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn5,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn5,"I","Equivalency")
a.dyn6<-ecospat.niche.equivalency.test(z1=z2_grid, z2=z4_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn6,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn6,"I","Equivalency")
a.dyn7<-ecospat.niche.equivalency.test(z1=z2_grid, z2=z5_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn7,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn7,"I","Equivalency")
a.dyn8<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z4_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn8,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn8,"I","Equivalency")
a.dyn9<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z5_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn9,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn9,"I","Equivalency")
a.dyn10<-ecospat.niche.equivalency.test(z1=z4_grid, z2=z5_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn10,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn10,"I","Equivalency")

a.dyn11<-ecospat.niche.equivalency.test(z1=z2_grid, z2=z1_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn1,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn1,"I","Equivalency")
a.dyn12<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z2_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn2,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn2,"I","Equivalency")
a.dyn13<-ecospat.niche.equivalency.test(z1=z4_grid, z2=z1_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn3,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn3,"I","Equivalency")
a.dyn14<-ecospat.niche.equivalency.test(z1=z5_grid, z2=z1_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn4,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn4,"I","Equivalency")
a.dyn15<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z2_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn5,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn5,"I","Equivalency")
a.dyn16<-ecospat.niche.equivalency.test(z1=z4_grid, z2=z2_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn6,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn6,"I","Equivalency")
a.dyn17<-ecospat.niche.equivalency.test(z1=z5_grid, z2=z2_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn7,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn7,"I","Equivalency")
a.dyn18<-ecospat.niche.equivalency.test(z1=z4_grid, z2=z3_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn8,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn8,"I","Equivalency")
a.dyn19<-ecospat.niche.equivalency.test(z1=z5_grid, z2=z3_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn9,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn9,"I","Equivalency")
a.dyn20<-ecospat.niche.equivalency.test(z1=z5_grid, z2=z4_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn10,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn10,"I","Equivalency")