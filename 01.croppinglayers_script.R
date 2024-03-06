##################################
####Script for layer cropping####
#################################
### In case some of the code doesn’t work, check the quotation marks, they may be off or not adequate

### Check working directory
### Install and activate packages
install.packages("raster")
install.packages("rgdal")
install.packages("sf")
install.packages("tidyverse")
install.packages("lattice")
install.packages("maptools")
install.packages("phyloclim")
install.packages("dismo")
install.packages("rJava")
install.packages("ENMeval")
install.packages("ecospat")
install.packages("grDevices")
library(sf)
library(raster)
library(rgdal)
library(tidyverse)
library(rgeos)
library(scales)
library(ecospat)
library(rJava)
library(dismo)
library(viridis)
library(gridExtra)

### Create a folder in which you can download and work your data

### Mesoamerican extension is : -118.32, -76.9, 7-0, 32.7189.
### The order of coordinates is: (West, East, South, North).
### Extent of SMS ecoregion is: -105.46, -95.082, 15.84, 21.01

### Extent of SMS for Aulacorhynchus is: -101.74, -95.62, 15.77, 18.45
### Extent of Eastern Mexico-NCA for Aulacorhynchus is: -100.81, -85.15, 12.66, 23.01
### Extent of SCA for Aulacorhynchus is: -84.99, -80.14, 7.38, 10.55

### Extent of TMBV for Cardellina is: -104.66, -96.70, 18.30, 20.83
### Extent of SMS for Cardellina is: -101.74, -95.08, 15.84, 19.13
### Extent of SMOCC for Cardellina is: -109.43, -104.31, 22.86, 30.88

### Extent of E. cyanophrys is: -97.55, -96.05, 15.84, 16.39
### Extent of E. eximia is: -97.19, -81.92, 8.51, 19.15
### Extent of E. poliocerca is: -100.66, -96.76, 16.08, 18.12
### Extent of E. ridgwayi is: -105.59, -103.24, 19.0024, 21.53
### Extent of E. nigriventris is: -84.99, -80.41, 8.34, 10.57

### Extent of Eastern Mexico for Chlorospingus is: -102.38, -95.62, 16.69, 25.68
### Extent of SCA for Chlorospingus is: -85.67, -81.90, 8.44, 11.04
### Extent of NCA for Chlorospingus is: -92.15, -84.69, 12.85, 16.07
### Extent of NChi for Chlorospingus is: -95.03, -91.25, 15.05, 17.45
### Extent of Tux for Chlorospingus is: -95.75, -94.64, 17.97, 18.69
### Extent of SMS for Chlorospingus is: -101.74, -95.62, 15.77, 18.45

###Extent for Aulacorhynchus lineages
sms_aula <- extent(-101.74, -95.62, 15.77, 18.45) 
emnca_aula <- extent (
sca_aula <- extent (

###Extent for Cardellina lineages
tmvb <- extent (
sms_carde <- extent (
smocc <- extent (

###Extent for Eupherusa lineages
cyano <- extent (
eximia <- extent (
polio <- extent (
rid <- extent (
nigri <- extent (

###Extent for Chlorospingus lineages
sca_chloro <- extent (
nca_chloro <- extent (
nchi <- extent (
tux <- extent (
sms_chloro <- extent (

### Here we only present the step-by-step for the SMS region for Aulacorhynchus. Steps should be repeated for each lineage extension and saved accordingly.  

### Create an object with a WGS84 projection
### To check the list of coordinate system references visit: http://spatialreference.org/ref/   https://epsg.io/
### WGS84 (code   EPSG:4326)
crs.wg <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

### List all tif archives found in a folder). 
bioclim.vars <- list.files("E:/GIS/Env_variables/Presente_2.1_30s", pattern=".tif$", full.names=T)

### Stack archives
bioclim.vars <- stack(bioclim.vars)

### Crop archives to the extension of the SMS region for Aulacorhynchus (same as Chlorospingus)
bioclim.vars.crop <- crop(bioclim.vars, sms_aula)

### Project to the defined extension. When projecting categorical data should use "ngb" (nearest neighbor) method, if data is continuous (like elevation) use "bilinear" method
bioclim.vars.WG <- projectRaster(bioclim.vars.crop, crs=crs.wg)

### Repeat for each lineage extension.
### One can check by plotting
### Repeat previous steps for LGM, Mid Holocene and LIG layers
LGMCC <- list.files("E:/GIS/Env_variables/LGM_CCSM4_2.5m", pattern=".tif$", full.names=T)
LGMCC <- stack(LGMCC)
LGMCC.crop <- crop(LGMCC, sms_aula)
LGMCC.WG <- projectRaster(LGMCC.crop, crs=crs.wg)
LGMMIR <- list.files("E:/GIS/Env_variables/LGM_MIROC_2.5m", pattern=".tif$", full.names=T)
LGMMIR <- stack(LGMMIR)
LGMMIR.crop <- crop(LGMMIR, sms_aula)
LGMMIR.WG <- projectRaster(LGMMIR.crop, crs=crs.wg)
LIG <- list.files("E:/GIS/Env_variables/LIG_30s", pattern=".bil$", full.names=T)
LIG <- stack(LIG)
LIG.crop <- crop(LIG, sms_aula)
LIG.WG <- projectRaster(LIG.crop, crs=crs.wg)
MHCC <- list.files("E:/GIS/Env_variables/MidHolocene_CCSM4_30s", pattern=".tif$", full.names=T)
MHCC <- stack(MHCC)
MHCC.crop <- crop(MHCC, sms_aula)
MHCC.WG <- projectRaster(MHCC.crop, crs=crs.wg)
MHMIR <- list.files("E:/GIS/Env_variables/MidHolocene_MIROC_30s", pattern=".tif$", full.names=T)
MHMIR <- stack(MHMIR)
MHMIR.crop <- crop(MHMIR, ext)
MHMIR.WG <- projectRaster(MHMIR.crop, crs=crs.wg)

### Check cropping by plotting
plot(bioclim.vars.WG$bio1)
plot(LGMCC.WG$bio1)
plot(LGMMIR.WG$bio1)
plot(LIG.WG$bio1)
plot(MHCC.WG$bio1)
plot(MHMIR.WG$bio1)

### Reescale LGM rasters, so cell size is the same
LGMCC.WG <- resample(LGMCC.WG, bioclim.vars.WG, method="bilinear")
LGMMIR.WG <-resample(LGMMIR.WG, bioclim.vars.WG, method="bilinear")

### Plot shapefiles. In this first example, I’m going to model Aulacorhynchus SMS lineage 
AulaSMS <- readOGR(dsn="D:/ms/nichos_SMS", layer="Pol_AulaSMS")
AulaSMS2 <- spTransform(AulaSMS, CRS=crs.wg) 
bioclim.vars.AulaSMS <- mask(bioclim.vars.WG, AulaSMS2)

### Layers with correlation of 0.8 for Aulacorhynchus according to "SDMToolbox" analysis
capas <- c("bio1","bio2","bio3","bio4","bio12","bio14","bio15")

###Here follows the variables for each taxon used in the study###
### Layers with correlation of 0.8 for Chlorospingus according to "SDMToolbox" analysis
#capas <- c("bio1", "bio2", "bio3", "bio12", "bio15", "bio18", "bio19")
### Layers with correlation of 0.8 for Cardellina according to "SDMToolbox" analysis
#capas <- c("bio1", "bio2", "bio3", "bio12", "bio14", "bio15", "bio18", "bio19")
### Layers with correlation of 0.8 for Eupherusa according to "SDMToolbox" analysis
#capas <- c("bio1", "bio2", "bio3", "bio12", "bio15", "bio18")
 
bioclim.vars.AulaSMS <- subset(bioclim.vars.AulaSMS, capas)
LGMCC.SMS <- mask(LGMCC.WG, AulaSMS2)
LGMMIR.SMS <- mask(LGMMIR.WG, AulaSMS2)
LIG.SMS <- mask(LIG.WG, AulaSMS2)
MHCC.SMS <- mask(MHCC.WG, AulaSMS2)
MHMIR.SMS <- mask(MHMIR.WG, AulaSMS2)

### Check by plotting
plot(bioclim.vars.AulaSMS$bio1)
plot(LGMCC.SMS$bio1)

### Save rasters of study area
writeRaster(bioclim.vars.AulaSMS, filename=paste0("aulacorhynchus/sms/M_variables/Set_1/", names(bioclim.vars.AulaSMS)), format="ascii", bylayer=T, overwrite=T)
writeRaster(LGMCC.SMS, filename=paste0("aulacorhynchus/sms/G_variables/Set_1/", names(LGMCC.SMS)), format="ascii", bylayer=T, overwrite=T)
writeRaster(LGMMIR.SMS, filename=paste0("aulacorhynchus/sms/G_variables/Set_2/", names(LGMMIR.SMS)), format="ascii", bylayer=T, overwrite=T)
writeRaster(LIG.SMS, filename=paste0("aulacorhynchus/sms/G_variables/Set_3/", names(LIG.SMS)), format="ascii", bylayer=T, overwrite=T)
writeRaster(MHCC.SMS, filename=paste0("aulacorhynchus/sms/G_variables/Set_4/", names(MHCC.SMS)), format="ascii", bylayer=T, overwrite=T)
writeRaster(MHMIR.SMS, filename=paste0("aulacorhynchus/sms/G_variables/Set_5/", names(MHMIR.SMS)), format="ascii", bylayer=T, overwrite=T)

###Repeat for all lineages!###
