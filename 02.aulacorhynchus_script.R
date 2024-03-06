####################################################
##################Niche Overlap#####################
####################################################

### Install and activate libraries.
library(ecospat)    ###Niche comparison in environmental space.
library(raster)     ###Raster manipulation.
library(SDMTools)   ###Niche-related analyses.
library(dismo)      ###Distribution modeling.
library(viridis)
library(ggplot2)
library(gridExtra)

### Read csv archive of distributional data.
### We have three genetic lineages to be compared.
sp1 <- read.csv("EMNCA/EMNCA_joint.csv", h=T)
sp2 <- read.csv("SCA/SCA_joint.csv", h=T)
sp3 <- read.csv("SMS/SMS_joint.csv", h=T)

### Load previously cropped layers that are relevant for comparison of ecological niches to generate the principal components.
varclim1 <- list.files("E:/aulacorhynchus/EMNCA/M_variables/Set_1/", pattern=".asc", full.names=T)
varclim2 <- list.files("E:/aulacorhynchus/SCA/M_variables/Set_1/", pattern=".asc", full.names=T)
varclim3 <- list.files("E:/aulacorhynchus/SMS/M_variables/Set_1/", pattern=".asc", full.names=T)

### Stack variables.
varclim1 <- stack(varclim1)
varclim2 <- stack(varclim2)
varclim3 <- stack(varclim3)

### Create a table for each pixel (x, y coordinates)
climpunto1 <- rasterToPoints(varclim1[[1]], fun=NULL, spatial=TRUE)
climpunto2 <- rasterToPoints(varclim2[[1]], fun=NULL, spatial=TRUE)
climpunto3 <- rasterToPoints(varclim3[[1]], fun=NULL, spatial=TRUE)

### Extraction of environmental data for each layer in capas of each coordinate in clipt.
clim1 <- extract(varclim1, climpunto1)
clim2 <- extract(varclim2, climpunto2)
clim3 <- extract(varclim3, climpunto3)

### Add climatic variables to data.
occ_sp1 <- na.exclude(ecospat.sample.envar(dfsp=sp1,colspxy=2:3, colspkept=2:3,dfvar=clim1, colvarxy=1:2,colvar="all",resolution=0.1))
occ_sp2 <- na.exclude(ecospat.sample.envar(dfsp=sp2,colspxy=2:3, colspkept=2:3,dfvar=clim2, colvarxy=1:2,colvar="all",resolution=0.1))
occ_sp3 <- na.exclude(ecospat.sample.envar(dfsp=sp3,colspxy=2:3, colspkept=2:3,dfvar=clim3, colvarxy=1:2,colvar="all",resolution=0.1))

#######################################################################
########################   Environmental PCA   ########################
#######################################################################
### Data for PCA to include all sites for the study area and for species presences.
### First clim and then the species in order.
### For our three genetic groups.
data <- rbind(clim1[,3:9],clim2[,3:9],clim3[,3:9],occ_sp1[,3:9],occ_sp2[,3:9],occ_sp3[,3:9])
data <- na.omit(data)

### Vector with 0 weight for occurrences and 1 for sites of the study area.
row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(row.w.3.env, nrow(clim3)),rep(0, nrow(occ_sp1)),rep(0, nrow(occ_sp2)),rep(0, nrow(occ_sp3)))

### PCA is done with all data from the study area.
### Presences are not used for calibration of PCA but coordinates are calculated. NUMBER OF AXIS=2.
pca.cal <-dudi.pca(na.omit(data, row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2))

##### IM-C: Tres áreas M distintas.
row.clim1 <- 1:nrow(clim1)
row.clim2<-(nrow(clim2)+1):(nrow(clim1)+nrow(clim2))
row.clim3<-(nrow(clim3)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3))

### For 3 genetic lineages to be compared.
row.clim123<-1:(nrow(clim1)+nrow(clim2)+nrow(clim3))
row.sp1<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(occ_sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(occ_sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(occ_sp1)+nrow(occ_sp2))
row.sp3<-(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(occ_sp1)+nrow(occ_sp2)+1):(nrow(clim1)+nrow(clim2)+nrow(clim3)+nrow(occ_sp1)+nrow(occ_sp2)+nrow(occ_sp3))

### Having the PCA results, we need the first and second eigenvector values for the background and the occurrence records per species
### Coordinates in each axis of PCA of all the study area and each species.
scores.clim123 <- pca.cal$li[row.clim123,] 
scores.clim1<- pca.cal$li[row.clim1,]
scores.clim2<- pca.cal$li[row.clim2,]
scores.clim3<- pca.cal$li[row.clim3,]
scores.sp1<- pca.cal$li[row.sp1,]
scores.sp2<- pca.cal$li[row.sp2,]
scores.sp3<- pca.cal$li[row.sp3,]

##### IM-C: Contribucion de cada variable a cada componente del PCA.
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)

### Environmental space.
### Resolution of the environmental space based on the PCA values calculated for the background and occurrence records
### Create density surfaces of occurrence in the environmental space (two axis), considering the observed occurrence density and availability of condition in the background
z1_grid <-ecospat.grid.clim.dyn (scores.clim123, scores.clim1, scores.sp1, R=100) 
z2_grid <- ecospat.grid.clim.dyn (scores.clim123, scores.clim2, scores.sp2, R=100) 
z3_grid <- ecospat.grid.clim.dyn (scores.clim123, scores.clim3, scores.sp3, R=100) 

### Metrics for niche overlap. Calculation of D metric and its significance through a similarity test.
### Generation of values of niche overlap.
ecospat.niche.overlap (z1=z1_grid, z2=z2_grid, cor=TRUE)  #EMNCA and SCA
ecospat.niche.overlap (z1=z1_grid, z2=z3_grid, cor=TRUE) #EMNCA and SMS  
ecospat.niche.overlap (z1=z2_grid, z2=z3_grid, cor=TRUE) #SCA and SMS

ecospat.niche.overlap (z1=z2_grid, z2=z1_grid, cor=TRUE)  #SCA and EMNCA
ecospat.niche.overlap (z1=z3_grid, z2=z1_grid, cor=TRUE) #SMS AND EMNCA  
ecospat.niche.overlap (z1=z3_grid, z2=z2_grid, cor=TRUE) #SMS AND SCA  

### Occurrence density graphics in the environmental space.
windows() 
par(mfrow=c(1,3))
ecospat.plot.niche (z1_grid, title="EMNCA", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z2_grid, title="SCA", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z3_grid, title="SMS", name.axis1="PC1", name.axis2="PC2", cor=F)

### Individual niche plots. Modifications of plot.niche functions from Broenniman et al. 2012 and Silva et al. 2014.
n.groups <- 3
g.names <- c("EMNCA", "SCA", "SMS")
g.colors <- c("darkblue", "cyan", "red")
z <- c(list(z1_grid, z2_grid, z3_grid))


### Plot E-space.
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

### Continuos line represents 100% of the available environmental background and the dashed lines represent 50% of the most common conditions.
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
a.dyn<-ecospat.niche.equivalency.test(z1=z1_grid , z2=z2_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn,"I","Equivalency")

a.dyn2<-ecospat.niche.equivalency.test(z1=z1_grid, z2=z3_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn2,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn2,"I","Equivalency")

a.dyn3<-ecospat.niche.equivalency.test(z1=z2_grid, z2=z3_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn3,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn3,"I","Equivalency")

a.dyn4<-ecospat.niche.equivalency.test(z1=z3_grid , z2=z1_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn,"I","Equivalency")

a.dyn5<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z1_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn2,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn2,"I","Equivalency")

a.dyn6<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z2_grid, rep=1000)
ecospat.plot.overlap.test(a.dyn3,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn3,"I","Equivalency")