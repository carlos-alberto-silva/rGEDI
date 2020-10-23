![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig1.png)<br/>

[![CRAN](https://www.r-pkg.org/badges/version/rGEDI)](https://cran.r-project.org/package=rGEDI)
![Github](https://img.shields.io/badge/Github-0.1.9-green.svg)
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rGEDI)
[![Build Status](https://travis-ci.com/carlos-alberto-silva/rGEDI.svg?token=Jqizwyc6gBxNafNccTdU&branch=master)](https://travis-ci.com/carlos-alberto-silva/rGEDI)

**rGEDI: An R Package for NASA's Global Ecosystem Dynamics Investigation (GEDI) Data Visualizing and Processing.**

Authors: Carlos Alberto Silva, Caio Hamamura, Ruben Valbuena, Steven Hancock, Adrian Cardil, Eben N. Broadbent, Danilo R. A. de Almeida, Celso H. L. Silva Junior and Carine Klauberg  

The rGEDI package provides functions for i) downloading, ii) visualizing, iii) clipping, iv) gridding, iv) simulating and v) exporting GEDI data.

# Getting Started

## Installation
```r
#The CRAN version:
install.packages("rGEDI")

#The development version:
library(devtools)
devtools::install_git("https://github.com/carlos-alberto-silva/rGEDI", dependencies = TRUE)

# loading rGEDI package
library(rGEDI)

```    

## Find GEDI data within your study area (GEDI finder tool)
```r
# Study area boundary box coordinates
ul_lat<- -44.0654
lr_lat<- -44.17246
ul_lon<- -13.76913
lr_lon<- -13.67646

# Specifying the date range
daterange=c("2019-07-01","2020-05-22")

# Get path to GEDI data
gLevel1B<-gedifinder(product="GEDI01_B",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel2A<-gedifinder(product="GEDI02_A",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel2B<-gedifinder(product="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
```
## Downloading GEDI data
```r
# Set output dir for downloading the files
outdir=getwd()

# Downloading GEDI data
gediDownload(filepath=gLevel1B,outdir=outdir)
gediDownload(filepath=gLevel2A,outdir=outdir)
gediDownload(filepath=gLevel2B,outdir=outdir)

#** Herein, we are using only a GEDI sample dataset for this tutorial.
# downloading zip file
download.file("https://github.com/carlos-alberto-silva/rGEDI/releases/download/datasets/examples.zip",destfile=file.path(outdir, "examples.zip"))

# unzip file 
unzip(file.path(outdir,"examples.zip"))

```

## Reading GEDI data
```r
# Reading GEDI data
gedilevel1b<-readLevel1B(level1Bpath = paste0(outdir,"\\GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.h5"))
gedilevel2a<-readLevel2A(level2Apath = paste0(outdir,"\\GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.h5"))
gedilevel2b<-readLevel2B(level2Bpath = paste0(outdir,"\\GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.h5"))
```

## Get GEDI Pulse Geolocation (GEDI Level1B)
```r
level1bGeo<-getLevel1BGeo(level1b=gedilevel1b,select=c("elevation_bin0"))
head(level1bGeo)

##           shot_number latitude_bin0 latitude_lastbin longitude_bin0 longitude_lastbin elevation_bin0
##  1: 19640002800109382     -13.75903        -13.75901      -44.17219         -44.17219       784.8348
##  2: 19640003000109383     -13.75862        -13.75859      -44.17188         -44.17188       799.0491
##  3: 19640003200109384     -13.75821        -13.75818      -44.17156         -44.17156       814.4647
##  4: 19640003400109385     -13.75780        -13.75777      -44.17124         -44.17124       820.1437
##  5: 19640003600109386     -13.75738        -13.75736      -44.17093         -44.17093       821.7012
##  6: 19640003800109387     -13.75697        -13.75695      -44.17061         -44.17061       823.2526

# Converting shot_number as "integer64" to "character"
level1bGeo$shot_number<-paste0(level1bGeo$shot_number)

# Converting level1bGeo as data.table to SpatialPointsDataFrame
library(sp)
level1bGeo_spdf<-SpatialPointsDataFrame(cbind(level1bGeo$longitude_bin0, level1bGeo$latitude_bin0),
                                        data=level1bGeo)

# Exporting level1bGeo as ESRI Shapefile
raster::shapefile(level1bGeo_spdf,paste0(outdir,"\\GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub"))
```
<img align="right" src="https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig2.PNG"  width="400">

```

library(leaflet)
library(leafsync)

leaflet() %>%
  addCircleMarkers(level1bGeo$longitude_bin0,
                   level1bGeo$latitude_bin0,
                   radius = 1,
                   opacity = 1,
                   color = "red")  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addLegend(colors = "red", labels= "Samples",title ="GEDI Level1B")


```

## Get GEDI Full-waveform (GEDI Level1B)
```r
# Extracting GEDI full-waveform for a giving shotnumber
wf <- getLevel1BWF(gedilevel1b, shot_number="19640521100108408")

par(mfrow = c(2,1), mar=c(4,4,1,1), cex.axis = 1.5)

plot(wf, relative=FALSE, polygon=TRUE, type="l", lwd=2, col="forestgreen",
     xlab="Waveform Amplitude", ylab="Elevation (m)")
grid()
plot(wf, relative=TRUE, polygon=FALSE, type="l", lwd=2, col="forestgreen",
     xlab="Waveform Amplitude (%)", ylab="Elevation (m)")
grid()
```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig3.png)

## Get GEDI Elevation and Height Metrics (GEDI Level2A)
```r
# Get GEDI Elevation and Height Metrics
level2AM<-getLevel2AM(gedilevel2a)
head(level2AM[,c("beam","shot_number","elev_highestreturn","elev_lowestmode","rh100")])

##          beam       shot_number elev_highestreturn elev_lowestmode rh100
##  1: BEAM0000 19640002800109382           740.7499        736.3301  4.41
##  2: BEAM0000 19640003000109383           756.0878        746.7614  9.32
##  3: BEAM0000 19640003200109384           770.3423        763.1509  7.19
##  4: BEAM0000 19640003400109385           775.9838        770.6652  5.31
##  5: BEAM0000 19640003600109386           777.8409        773.0841  4.75
##  6: BEAM0000 19640003800109387           778.7181        773.6990  5.01

# Converting shot_number as "integer64" to "character"
level2AM$shot_number<-paste0(level2AM$shot_number)

# Converting Elevation and Height Metrics as data.table to SpatialPointsDataFrame
level2AM_spdf<-SpatialPointsDataFrame(cbind(level2AM$lon_lowestmode,level2AM$lat_lowestmode),
                                        data=level2AM)

# Exporting Elevation and Height Metrics as ESRI Shapefile
raster::shapefile(level2AM_spdf,paste0(outdir,"\\GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub"))
```

## Plot waveform with RH metrics
```r
shot_number = "19640521100108408"

png("fig8.png", width = 8, height = 6, units = 'in', res = 300)
plotWFMetrics(gedilevel1b, gedilevel2a, shot_number, rh=c(25, 50, 75, 90))
dev.off()
```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig8.png)
## Get GEDI Vegetation Biophysical Variables (GEDI Level2B)
```r
level2BVPM<-getLevel2BVPM(gedilevel2b)
head(level2BVPM[,c("beam","shot_number","pai","fhd_normal","omega","pgap_theta","cover")])

##          beam       shot_number         pai fhd_normal omega pgap_theta       cover
##   1: BEAM0000 19640002800109382 0.007661204  0.6365142     1  0.9961758 0.003823273
##   2: BEAM0000 19640003000109383 0.086218357  2.2644432     1  0.9577964 0.042192958
##   3: BEAM0000 19640003200109384 0.299524575  1.8881851     1  0.8608801 0.139084846
##   4: BEAM0000 19640003400109385 0.079557180  1.6625489     1  0.9609926 0.038997617
##   5: BEAM0000 19640003600109386 0.018724868  1.5836401     1  0.9906789 0.009318732
##   6: BEAM0000 19640003800109387 0.017654873  1.2458609     1  0.9912092 0.008788579

# Converting shot_number as "integer64" to "character"
level2BVPM$shot_number<-paste0(level2BVPM$shot_number)

# Converting GEDI Vegetation Profile Biophysical Variables as data.table to SpatialPointsDataFrame
level2BVPM_spdf<-SpatialPointsDataFrame(cbind(level2BVPM$longitude_lastbin,level2BVPM$latitude_lastbin),data=level2BVPM)

# Exporting GEDI Vegetation Profile Biophysical Variables as ESRI Shapefile
raster::shapefile(level2BVPM_spdf,paste0(outdir,"\\GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub_VPM"))

```

## Get Plant Area Index (PAI) and Plant Area Volume Density (PAVD) Profiles (GEDI Level2B)
```r
level2BPAIProfile<-getLevel2BPAIProfile(gedilevel2b)
head(level2BPAIProfile[,c("beam","shot_number","pai_z0_5m","pai_z5_10m")])

##          beam       shot_number   pai_z0_5m   pai_z5_10m
##   1: BEAM0000 19640002800109382 0.007661204 0.0000000000
##   2: BEAM0000 19640003000109383 0.086218357 0.0581122264
##   3: BEAM0000 19640003200109384 0.299524575 0.0497199222
##   4: BEAM0000 19640003400109385 0.079557180 0.0004457365
##   5: BEAM0000 19640003600109386 0.018724868 0.0000000000
##   6: BEAM0000 19640003800109387 0.017654873 0.0000000000

level2BPAVDProfile<-getLevel2BPAVDProfile(gedilevel2b)
head(level2BPAVDProfile[,c("beam","shot_number","pavd_z0_5m","pavd_z5_10m")])

##          beam       shot_number  pavd_z0_5m  pavd_z5_10m
##   1: BEAM0000 19640002800109382 0.001532241 0.0007661204
##   2: BEAM0000 19640003000109383 0.005621226 0.0086218351
##   3: BEAM0000 19640003200109384 0.049960934 0.0299524590
##   4: BEAM0000 19640003400109385 0.015822290 0.0079557188
##   5: BEAM0000 19640003600109386 0.003744974 0.0018724868
##   6: BEAM0000 19640003800109387 0.003530974 0.0017654872

# Converting shot_number as "integer64" to "character"
level2BPAIProfile$shot_number<-paste0(level2BPAIProfile$shot_number)
level2BPAVDProfile$shot_number<-paste0(level2BPAVDProfile$shot_number)

# Converting PAI and PAVD Profiles as data.table to SpatialPointsDataFrame
level2BPAIProfile_spdf<-SpatialPointsDataFrame(cbind(level2BPAIProfile$lon_lowestmode,level2BPAIProfile$lat_lowestmode),
                                        data=level2BPAIProfile)
level2BPAVDProfile_spdf<-SpatialPointsDataFrame(cbind(level2BPAVDProfile$lon_lowestmode,level2BPAVDProfile$lat_lowestmode),
                                               data=level2BPAVDProfile)

# Exporting PAI and PAVD Profiles as ESRI Shapefile
raster::shapefile(level2BPAIProfile_spdf,paste0(outdir,"\\GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub_PAIProfile"))
raster::shapefile(level2BPAVDProfile_spdf,paste0(outdir,"\\GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub_PAVDProfile"))

```

## Plot Plant Area Index (PAI) and Plant Area Volume Density (PAVD) Profiles 
```r
#specify GEDI beam
beam="BEAM0101"

# Plot Level2B PAI Profile
gPAIprofile<-plotPAIProfile(level2BPAIProfile, beam=beam, elev=TRUE)

# Plot Level2B PAVD Profile
gPAVDprofile<-plotPAVDProfile(level2BPAVDProfile, beam=beam, elev=TRUE)

```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig9.png)


## Clip GEDI data (h5 files; gedi.level1b, gedi.level2a and gedi.level2b objects)
```r
## Clip GEDI data by coordinates
# Study area boundary box
xmin = -44.15036
xmax = -44.10066
ymin = -13.75831
ymax = -13.71244

## clipping GEDI data within boundary box
level1b_clip_bb <- clipLevel1B(gedilevel1b, xmin, xmax, ymin, ymax,output=paste0(outdir,"//level1b_clip_bb.h5"))
level2a_clip_bb <- clipLevel2A(gedilevel2a, xmin, xmax, ymin, ymax, output=paste0(outdir,"//level2a_clip_bb.h5"))
level2b_clip_bb <- clipLevel2B(gedilevel2b, xmin, xmax, ymin, ymax,output=paste0(outdir,"//level2b_clip_bb.h5"))

## Clipping GEDI data by geometry
# specify the path to shapefile for the study area
polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package="rGEDI")

# Reading shapefile as SpatialPolygonsDataFrame object
library(rgdal)
polygon_spdf<-readOGR(polygon_filepath)
head(polygon_spdf@data) # column id name "id"
split_by="id"

# Clipping h5 files
level1b_clip_gb <- clipLevel1BGeometry(gedilevel1b,polygon_spdf,output=paste0(outdir,"//level1b_clip_gb.h5"), split_by=split_by)
level2a_clip_gb <- clipLevel2AGeometry(gedilevel2a,polygon_spdf,output=paste0(outdir,"//level2a_clip_gb.h5"), split_by=split_by)
level2b_clip_gb <- clipLevel2BGeometry(gedilevel2b,polygon_spdf,output=paste0(outdir,"//level2b_clip_gb.h5"), split_by=split_by)
```
## Clip GEDI data (data.table objects)
```r
## Clipping GEDI data within boundary box
level1bGeo_clip_bb <-clipLevel1BGeo(level1bGeo, xmin, xmax, ymin, ymax)
level2AM_clip_bb <- clipLevel2AM(level2AM, xmin, xmax, ymin, ymax)
level2BVPM_clip_bb <- clipLevel2BVPM(level2BVPM, xmin, xmax, ymin, ymax)
level1BPAIProfile_clip_bb <- clipLevel2BPAIProfile(level2BPAIProfile, xmin, xmax, ymin, ymax)
level2BPAVDProfile_clip_bb <- clipLevel2BPAVDProfile(level2BPAVDProfile, xmin, xmax, ymin, ymax)

## Clipping GEDI data by geometry
level1bGeo_clip_gb <- clipLevel1BGeoGeometry(level1bGeo,polygon_spdf, split_by=split_by)
level2AM_clip_gb <- clipLevel2AMGeometry(level2AM,polygon_spdf, split_by=split_by)
level2BVPM_clip_gb <- clipLevel2BVPMGeometry(level2BVPM,polygon_spdf, split_by=split_by)
level1BPAIProfile_clip_gb <- clipLevel2BPAIProfileGeometry(level2BPAIProfile,polygon_spdf, split_by=split_by)
level2BPAVDProfile_clip_gb <- clipLevel2BPAVDProfileGeometry(level2BPAVDProfile,polygon_spdf, split_by=split_by)


## View GEDI clipped data by bbox
m1<-leaflet() %>%
  addCircleMarkers(level2AM$lon_lowestmode,
                   level2AM$lat_lowestmode,
                   radius = 1,
                   opacity = 1,
                   color = "red")  %>%
  addCircleMarkers(level2AM_clip_bb$lon_lowestmode,
                   level2AM_clip_bb$lat_lowestmode,
                   radius = 1,
                   opacity = 1,
                   color = "green")  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery)  %>%
  addLegend(colors = c("red","green"), labels= c("All samples","Clip bbox"),title ="GEDI Level2A") 

## View GEDI clipped data by geometry
# color palette
pal <- colorFactor(
  palette = c('blue', 'green', 'purple', 'orange',"white","black","gray","yellow"),
  domain = level2AM_clip_gb$poly_id
)

m2<-leaflet() %>%
  addCircleMarkers(level2AM$lon_lowestmode,
                   level2AM$lat_lowestmode,
                   radius = 1,
                   opacity = 1,
                   color = "red")  %>%
  addCircleMarkers(level2AM_clip_gb$lon_lowestmode,
                   level2AM_clip_gb$lat_lowestmode,
                   radius = 1,
                   opacity = 1,
                   color = pal(level2AM_clip_gb$poly_id))  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addPolygons(data=polygon_spdf,weight=1,col = 'white',
              opacity = 1, fillOpacity = 0) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addLegend(pal = pal, values = level2AM_clip_gb$poly_id,title ="Poly IDs" ) 

sync(m1, m2)
```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig4.png)

## Compute descriptive statistics of GEDI Level2A and Level2B data
```r
# Define your own function
mySetOfMetrics = function(x)
{
metrics = list(
    min =min(x), # Min of x
    max = max(x), # Max of x
    mean = mean(x), # Mean of x
    sd = sd(x)# Sd of x
  )
  return(metrics)
}

# Computing the maximum of RH100 stratified by polygon
rh100max_st<-polyStatsLevel2AM(level2AM_clip_gb,func=max(rh100), id="poly_id")
head(rh100max_st)

##    poly_id   max
## 1:       2 12.81
## 2:       1 12.62
## 3:       5  9.96
## 4:       6  8.98
## 5:       4 10.33
## 6:       8  8.72

# Computing a serie statistics for GEDI metrics stratified by polygon
rh100metrics_st<-polyStatsLevel2AM(level2AM_clip_gb,func=mySetOfMetrics(rh100),
id="poly_id")
head(rh100metrics_st)

##    poly_id  min   max     mean       sd
## 1:       2 4.08 12.81 5.508639 1.452143
## 2:       1 3.78 12.62 5.514930 1.745507
## 3:       5 4.12  9.96 5.100122 1.195272
## 4:       6 4.64  8.98 5.595294 1.024171
## 5:       4 4.38 10.33 7.909500 1.757200
## 6:       8 4.45  8.72 6.136471 1.097468

# Computing the max of the Total Plant Area Index
pai_max<-polyStatsLevel2BVPM(level2BVPM_clip_gb,func=max(pai), id=NULL)
pai_max

##          max
#   1: 1.224658

# Computing a serie of statistics of Canopy Cover stratified by polygon
cover_metrics_st<-polyStatsLevel2BVPM(level2BVPM_clip_gb,func=mySetOfMetrics(cover),
id="poly_id")
head(cover_metrics_st)

##     poly_id          min       max       mean         sd
##  1:       2 0.0010017310 0.3479594 0.05156159 0.05817241
##  2:       1 0.0003717059 0.3812594 0.04829096 0.06346548
##  3:       5 0.0020242794 0.4262614 0.03577852 0.06407325
##  4:       6 0.0028748326 0.2392146 0.03094646 0.05577988
##  5:       4 0.0022404396 0.3501986 0.11343149 0.09354305
##  6:       8 0.0050588539 0.1457105 0.04784596 0.04427151
```

## Compute Grids with descriptive statistics of GEDI-derived Elevation and Height Metrics (Level2A)

<img align="right" src="https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig5.png" width="300">

```r
# Computing a serie of statistics of GEDI RH100 metric
rh100metrics<-gridStatsLevel2AM(level2AM = level2AM, func=mySetOfMetrics(rh100), res=0.005)

# View maps
library(rasterVis)
library(viridis)

rh100maps<-levelplot(rh100metrics,
                     layout=c(1, 4),
                     margin=FALSE,
                     xlab = "Longitude (degree)", ylab = "Latitude (degree)",
                     colorkey=list(
                       space='right',
                       labels=list(at=seq(0, 18, 2), font=4),
                       axis.line=list(col='black'),
                       width=1),
                     par.settings=list(
                       strip.border=list(col='gray'),
                       strip.background=list(col='gray'),
                       axis.line=list(col='gray')
                     ),
                     scales=list(draw=TRUE),
                     col.regions=viridis,
                     at=seq(0, 18, len=101),
                     names.attr=c("rh100 min","rh100 max","rh100 mean", "rh100 sd"))

# Exporting maps 
png("fig6.png", width = 6, height = 8, units = 'in', res = 300)
rh100maps
dev.off()



```

## Compute Grids with descriptive statistics of GEDI-derived Canopy Cover and Vertical Profile Metrics (Level2B)

<img align="right" src="https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig6.png" width="300">

```r
# Computing a serie of statistics of Total Plant Area Index
level2BVPM$pai[level2BVPM$pai==-9999]<-NA # assing NA to -9999
pai_metrics<-gridStatsLevel2BVPM(level2BVPM = level2BVPM, func=mySetOfMetrics(pai), res=0.005)

# View maps
pai_maps<-levelplot(pai_metrics,
                    layout=c(1, 4),
                    margin=FALSE,
                    xlab = "Longitude (degree)", ylab = "Latitude (degree)",
                    colorkey=list(
                      space='right',
                      labels=list(at=seq(0, 1.5, 0.2), font=4),
                      axis.line=list(col='black'),
                      width=1),
                    par.settings=list(
                      strip.border=list(col='gray'),
                      strip.background=list(col='gray'),
                      axis.line=list(col='gray')
                    ),
                    scales=list(draw=TRUE),
                    col.regions=viridis,
                    at=seq(0, 1.5, len=101),
                    names.attr=c("PAI min","PAI max","PAI mean", "PAI sd"))

# Exporting maps 
png("fig6.png", width = 6, height = 8, units = 'in', res = 300)
pai_maps
dev.off()




```

## Simulating GEDI full-waveform data from Airborne Laser Scanning (ALS) 3-D point cloud and extracting canopy derived metrics
```r
# Specifying the path to ALS data
lasfile_amazon <- file.path(outdir, "Amazon.las")
lasfile_savanna <- file.path(outdir, "Savanna.las")

# Reading and plot ALS file
library(lidR)
library(plot3D)
las_amazon<-readLAS(lasfile_amazon)
las_savanna<-readLAS(lasfile_savanna)

# Extracting plot center geolocations
xcenter_amazon = mean(las_amazon@bbox[1,])
ycenter_amazon = mean(las_amazon@bbox[2,])
xcenter_savanna = mean(las_savanna@bbox[1,])
ycenter_savanna = mean(las_savanna@bbox[2,])

# Simulating GEDI full-waveform
wf_amazon<-gediWFSimulator(input=lasfile_amazon,output=paste0(getwd(),"//gediWF_amazon_simulation.h5"),coords = c(xcenter_amazon, ycenter_amazon))
wf_savanna<-gediWFSimulator(input=lasfile_savanna,output=paste0(getwd(),"//gediWF_savanna_simulation.h5"),coords = c(xcenter_savanna, ycenter_savanna))

# Plotting ALS and GEDI simulated full-waveform
png("gediWf.png", width = 8, height = 6, units = 'in', res = 300)

par(mfrow=c(2,2), mar=c(4,4,0,0), oma=c(0,0,1,1),cex.axis = 1.2)
scatter3D(las_amazon@data$X,las_amazon@data$Y,las_amazon@data$Z,pch = 16,colkey = FALSE, main="",
          cex = 0.5,bty = "u",col.panel ="gray90",phi = 30,alpha=1,theta=45,
          col.grid = "gray50", xlab="UTM Easting (m)", ylab="UTM Northing (m)", zlab="Elevation (m)")

plot(wf_amazon, relative=TRUE, polygon=TRUE, type="l", lwd=2, col="forestgreen",
     xlab="", ylab="Elevation (m)", ylim=c(90,140))
grid()
scatter3D(las_savanna@data$X,las_savanna@data$Y,las_savanna@data$Z,pch = 16,colkey = FALSE, main="",
          cex = 0.5,bty = "u",col.panel ="gray90",phi = 30,alpha=1,theta=45,
          col.grid = "gray50", xlab="UTM Easting (m)", ylab="UTM Northing (m)", zlab="Elevation (m)")

plot(wf_savanna, relative=TRUE, polygon=TRUE, type="l", lwd=2, col="green",
xlab="Waveform Amplitude (%)", ylab="Elevation (m)", ylim=c(815,835))
grid()
dev.off()
```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig7.png)

## Extracting GEDI full-waveform derived metrics without adding noise to the full-waveform
```
wf_amazon_metrics<-gediWFMetrics(input=wf_amazon,
                                outRoot=file.path(getwd(), "amazon"))
wf_savanna_metrics<-gediWFMetrics(input=wf_savanna,
                                outRoot=file.path(getwd(), "savanna"))

metrics<-rbind(wf_amazon_metrics,wf_savanna_metrics)
rownames(metrics)<-c("Amazon","Savanna")
head(metrics[,1:8])

#         #wave ID true ground true top ground slope ALS cover gHeight maxGround inflGround
# Amazon         0      -1e+06   133.29       -1e+06        -1   94.97     99.84      95.19
# Savanna        0      -1e+06   831.51       -1e+06        -1  822.18    822.21     822.25
```
## Extracting GEDI full-waveform derived metrics after adding noise to the full-waveform
```
wf_amazon_metrics_noise<-gediWFMetrics(input=wf_amazon,
                         outRoot=file.path(getwd(), "amazon"),
                         linkNoise= c(3.0103,0.95),
                         maxDN= 4096,
                         sWidth= 0.5,
                         varScale= 3)

wf_savanna_metrics_noise<-gediWFMetrics(
                        input=wf_savanna,
                        outRoot=file.path(getwd(), "savanna"),
                        linkNoise= c(3.0103,0.95),
                        maxDN= 4096,
                        sWidth= 0.5,
                        varScale= 3)

metrics_noise<-rbind(wf_amazon_metrics_noise,wf_savanna_metrics_noise)
rownames(metrics_noise)<-c("Amazon","Savanna")
head(metrics_noise[,1:8])

#         #wave ID true ground true top ground slope ALS cover gHeight maxGround inflGround
# Amazon         0      -1e+06   133.29       -1e+06        -1   99.17     99.99      95.39
# Savanna        0      -1e+06   831.36       -1e+06        -1  822.15    822.21     822.18

```

## Always close gedi objects, so HDF5 files won't be blocked!
```{r cleanup, echo=TRUE, results="hide", error=TRUE}
close(wf_amazon)
close(wf_savanna)
close(gedilevel1b)
close(gedilevel2a)
close(gedilevel2b)
```


# References
Dubayah, R., Blair, J.B., Goetz, S., Fatoyinbo, L., Hansen, M., Healey, S., Hofton, M., Hurtt, G.,         Kellner, J., Luthcke, S., & Armston, J. (2020) The Global Ecosystem Dynamics Investigation:         High-resolution laser ranging of the Earth’s forests and topography. Science of Remote             Sensing, p.100002. https://doi.org/10.1016/j.srs.2020.100002

Hancock, S., Armston, J., Hofton, M., Sun, X., Tang, H., Duncanson, L.I., Kellner,
       J.R. and Dubayah, R., 2019. The GEDI simulator: A large-footprint waveform lidar simulator
       for calibration and validation of spaceborne missions. Earth and Space Science.
       https://doi.org/10.1029/2018EA000506

Silva, C. A.; Saatchi, S.; Alonso, M. G. ; Labriere, N. ; Klauberg, C. ; Ferraz, A. ; Meyer, V. ;        Jeffery, K. J. ; Abernethy, K. ; White, L. ; Zhao, K. ; Lewis, S. L. ; Hudak, A. T. (2018)         Comparison of Small- and Large-Footprint Lidar Characterization of Tropical Forest                 Aboveground Structure and Biomass: A Case Study from Central Gabon. IEEE Journal of Selected       Topics in Applied Earth Observations and Remote Sensing, p. 1-15.
      https://ieeexplore.ieee.org/document/8331845

GEDI webpage. Accessed on February 15 2020 https://gedi.umd.edu/   
GEDI01_Bv001. Accessed on February 15 2020 https://lpdaac.usgs.gov/products/gedi01_bv001/   
GEDI02_Av001. Accessed on February 15 2020 https://lpdaac.usgs.gov/products/gedi02_av001/  
GEDI02_Bv001. Accessed on February 15 2020 https://lpdaac.usgs.gov/products/gedi02_bv001/  
GEDI Finder. Accessed on February 15 2020 https://lpdaacsvc.cr.usgs.gov/services/gedifinder

# Acknowledgements
The University of Maryland and NASA's Goddard Space Flight Center for developing GEDI mission.

We gratefully acknowledge funding from NASA’s Carbon Monitoring Systems, grant NNH15ZDA001N-CMS. Project entitled "Future Mission Fusion for High Biomass Forest Carbon Accounting" led by Dr. Laura Duncanson (lduncans@umd.edu, University of Maryland) and Dr. Lola Fatoyinbo (lola.fatoyinbo@nasa.gov, NASA's Goddard Space Flight Center).

The Brazilian National Council for Scientific and Technological Development (CNPq) for funding the project entitled "Mapping fuel load and simulation of fire behaviour and spread in the Cerrado biome using modeling and remote sensing technologies" and leaded by Prof. Dr. Carine Klauberg (carine_klauberg@hotmail.com) and Dr. Carlos Alberto Silva
(carlos_engflorestal@outlook.com).

# Reporting Issues 
Please report any issue reagardling the rGEDI package herein https://groups.yahoo.com/neo/groups/rGEDI

# Citing rGEDI
Silva,C.A; Hamamura,C.; Valbuena, R.; Hancock,S.; Cardil,A.; Broadbent, E. N.; Almeida,D.R.A.; Silva Junior, C.H.L; Klauberg, C. rGEDI: NASA's Global Ecosystem Dynamics Investigation (GEDI) Data Visualization and Processing.
version 0.1.9, accessed on October. 22 2020, available at: <https://CRAN.R-project.org/package=rGEDI>

# Disclaimer
**rGEDI package has not been developted by the GEDI team. It comes with no guarantee, expressed or implied, and the authors hold no responsibility for its use or reliability of its outputs.**

