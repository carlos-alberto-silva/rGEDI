![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig1.png)<br/>

![Github](https://img.shields.io/badge/CRAN-0.0.1-green.svg)
![Github](https://img.shields.io/badge/Github-0.0.1-green.svg)
[![Rdoc](http://www.rdocumentation.org/badges/version/rGEDI)](http://www.rdocumentation.org/packages/ForestGapR)
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
![R_Forge](https://img.shields.io/badge/R_Forge-0.0.1-green.svg) 
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rGEDI)

**rGEDI: An R Package for NASA's Global Ecosystem Dynamics Investigation (GEDI) data 
downloading, visualizing and processing.**

Authors: Carlos A. Silva, Caio Hamamura, Ruben Valbuena, Steve Hancock, Adrian Cardil, Eben Broadbent, Danilo R. A. de Almeida, Celso H. L. Silva Junior and Carine Klauberg  

The rGEDI package provides functions for i) downloading, ii) visualizing, iii) clipping, iv) exporting, iv) Gridding and v) Simulating GEDI data.

# Installation
```r
#The CRAN version:
install.packages("rGEDI")

#The development version:
library(devtools)
devtools::install_github("carlos-alberto-silva/rGEDI")
```    

# Getting Started

## GEDIfinder tool
```r
# Find GEDI data within your study area
# Study area boundary box
xmin<- -44.17246
ymin<- -44.0654
xmax<- -13.76913
ymax<- -13.67646

# Get path to GEDI data
gLevel1B<-gediFinder(level="GEDI01_B",xleft, xright, ybottom, ytop)
gLevel2A<-gediFinder(level="GEDI02_A",xleft, xright, ybottom, ytop)
gLevel2B<-gediFinder(level="GEDI02_B",xleft, xright, ybottom, ytop)
```
## Downloading GEDI data
```r
# Set output dir for downloading the files
outdir=tempdir()

# Downloading GEDI data
LPDAACDataPool(filepath=gLevel1B,outdir)
LPDAACDataPool(filepath=gLevel2A,outdir)
LPDAACDataPool(filepath=gLevel2B,outdir)
```

## Reading GEDI data
```r
# As it takes time for downloading GEDI data, herein we will be using only samples 
# specify the path to GEDI samples
GEDI01_B_urlfile="https://github.com/carlos-alberto-silva/rGEDI/blob/master/inst/extdata/GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.h5"
GEDI02_A_urlfile="https://github.com/carlos-alberto-silva/rGEDI/blob/master/inst/extdata/GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.h5"
GEDI02_B_urlfile="https://github.com/carlos-alberto-silva/rGEDI/blob/master/inst/extdata/GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.h5"

# Download GEDI samples for the example. The files will be downloaded in the current working directory [see getwd()]
download.file(GEDI01_B_urlfile, destfile = paste0(getwd(),"//",basename(GEDI01_B_urlfile)))
download.file(GEDI02_A_urlfile, destfile = paste0(getwd(),"//",basename(GEDI02_A_urlfile)))
download.file(GEDI02_B_urlfile, destfile = paste0(getwd(),"//",basename(GEDI02_B_urlfile)))

# Reading GEDI data
gedilevel1b<-readLevel1B(level1bpath = paste0(getwd(),"//",basename(GEDI01_B_urlfile)))
gedilevel2a<-readLevel1B(level2apath = paste0(getwd(),"//",basename(GEDI02_A_urlfile)))
gedilevel2b<-readLevel1B(level2bpath = paste0(getwd(),"//",basename(GEDI02_B_urlfile)))
```

## Get GEDI Pulse Full-Waveform Geolocation (GEDI Level1B)
```r
level1BGeo<-getLevel1BGeo(level1b,select=NULL)
head(level1BGeo)

##           shot_number latitude_bin0 latitude_lastbin longitude_bin0 longitude_lastbin
##  1: 19640002800109382     -13.75903        -13.75901      -44.17219         -44.17219
##  2: 19640003000109383     -13.75862        -13.75859      -44.17188         -44.17188
##  3: 19640003200109384     -13.75821        -13.75818      -44.17156         -44.17156
##  4: 19640003400109385     -13.75780        -13.75777      -44.17124         -44.17124
##  5: 19640003600109386     -13.75738        -13.75736      -44.17093         -44.17093
##  6: 19640003800109387     -13.75697        -13.75695      -44.17061         -44.17061
    
library(leaflet)
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
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig2.PNG)

## Get GEDI Pulse Full-Waveform (GEDI Level1B)
```r
# Extracting GEDI full-waveform for a giving shotnumber
wf <- getLevel1BWF(level1b, shot_number="19640521100108408")

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
level2AM<-getLevel2AM(level2a)
head(level2AM[,c("beam","shot_number","rh100")]) 

##         beam       shot_number rh100
##  1: BEAM0000 19640002800109382  4.41
##  2: BEAM0000 19640003000109383  9.32
##  3: BEAM0000 19640003200109384  7.19
##  4: BEAM0000 19640003400109385  5.31
##  5: BEAM0000 19640003600109386  4.75
##  6: BEAM0000 19640003800109387  5.01
```

## Get GEDI Plant Area Index (PAI) Profile (GEDI Level2B)
```r
level2BPAIProfile<-getLevel2BPAIProfile(level2b)
head(level2BPAIProfile[,c("beam","shot_number","pai_z0_5m","pai_z5_10m")])

##          beam       shot_number   pai_z0_5m   pai_z5_10m
##   1: BEAM0000 19640002800109382 0.007661204 0.0000000000
##   2: BEAM0000 19640003000109383 0.086218357 0.0581122264
##   3: BEAM0000 19640003200109384 0.299524575 0.0497199222
##   4: BEAM0000 19640003400109385 0.079557180 0.0004457365
##   5: BEAM0000 19640003600109386 0.018724868 0.0000000000
##   6: BEAM0000 19640003800109387 0.017654873 0.0000000000
```

## Get GEDI Plant Area Volume Density (PAVD) Index (GEDI Level2B)
```r
level2BPAVDProfile<-getLevel2BPAVDProfile(level2b)
head(level2BPAVDProfile[,c("beam","shot_number","pavd_z0_5m","pavd_z5_10m")])

##          beam       shot_number  pavd_z0_5m  pavd_z5_10m
##   1: BEAM0000 19640002800109382 0.001532241 0.0007661204
##   2: BEAM0000 19640003000109383 0.005621226 0.0086218351
##   3: BEAM0000 19640003200109384 0.049960934 0.0299524590
##   4: BEAM0000 19640003400109385 0.015822290 0.0079557188
##   5: BEAM0000 19640003600109386 0.003744974 0.0018724868
##   6: BEAM0000 19640003800109387 0.003530974 0.0017654872
```

## Get GEDI Vegetation Profile Biophysical Variables (GEDI Level2B)
```r
level2BVPM<-getLevel2BVPM(level2b)
head(level2BVPM[,c("beam","shot_number","pai","fhd_normal","omega","pgap_theta","cover")])

##          beam       shot_number         pai fhd_normal omega pgap_theta       cover
##   1: BEAM0000 19640002800109382 0.007661204  0.6365142     1  0.9961758 0.003823273
##   2: BEAM0000 19640003000109383 0.086218357  2.2644432     1  0.9577964 0.042192958
##   3: BEAM0000 19640003200109384 0.299524575  1.8881851     1  0.8608801 0.139084846
##   4: BEAM0000 19640003400109385 0.079557180  1.6625489     1  0.9609926 0.038997617
##   5: BEAM0000 19640003600109386 0.018724868  1.5836401     1  0.9906789 0.009318732
##   6: BEAM0000 19640003800109387 0.017654873  1.2458609     1  0.9912092 0.008788579

```
# Clip GEDI data (h5 files; gedi.level1b, gedi.level2a and gedi.level2b objects)
```r
## Clip GEDI data by coordinates
# Study area boundary box
xmin = -44.15036
xmax = -44.10066
ymin = -13.75831
ymax = -13.71244

## clipping GEDI data within boundary box
level1b_clip_bb <- clipLevel1BGeometry(level1b, xmin, xmax, ymin, ymax,output=paste0(getwd(),"//level1b_clip_gb.h5"))
level2a_clip_bb <- clipLevel2AGeometry(level2a, xmin, xmax, ymin, ytop,output=paste0(getwd(),"//level2a_clip_gb.h5"))
level2b_clip_bb <- clipLevel2BGeometry(level2b, xmin, xmax, ymin, ytop,output=paste0(getwd(),"//level2b_clip_gb.h5"))

## Clip GEDI data by geometry
# specify the path to shapefile for the study area
polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package="rGEDI")

# Reading shapefile as SpatialPolygonsDataFrame object
library(rgdal)
polygon_spdf<-readOGR(polygons_filepath)

# Clip h5 files
level1b_clip_gb <- clipLevel1BGeometry(level1b,polygon_spdf,output=paste0(getwd(),"//level1b_clip_gb.h5"), split_by="id")
level2a_clip_gb <- clipLevel2AGeometry(level2a,polygon_spdf,output=paste0(getwd(),"//level2a_clip_gb.h5"), split_by="id")
level2b_clip_gb <- clipLevel2BGeometry(level2b,polygon_spdf,output=paste0(getwd(),"//level2b_clip_gb.h5"), split_by="id")

```
# Clip GEDI data (data.table objects)
```r
## clipping GEDI data within boundary box
level1BGeo_clip_bb <-clipLevel1BGeo(level1BGeo, xmin, xmax, ymin, ymax)
level2AM_clip_bb <- clipLevel2AM(level2AM,xleft, xmin, xmax, ymin, ymax)
level2BVPM_clip_bb <- clipLevel2BVPM(level2BVPM, xmin, xmax, ymin, ymax)
level1BPAIProfile_clip_bb <- clipLevel1BPAIProfile(level1BPAIProfile, xmin, xmax, ymin, ymax)
level2BPAVDProfile_clip_bb <- clipLevel2BPAVDProfile(level2BPAVDProfile, xmin, xmax, ymin, ymax)

## Clip GEDI data by geometry
level1BGeo_clip_gb <- clipLevel1BGeo(level1BGeo,polygon_spdf, split_by="id")
level2AM_clip_gb <- clipLevel2AM(level2AM,polygon_spdf, split_by="id")
level2BVPM_clip_gb <- clipLevel2BVPM(level2BVPM,polygon_spdf, split_by="id")
level1BPAIProfile_clip_gb <- clipLevel1BPAIProfile(level1BPAIProfile,polygon_spdf, split_by="id")
level2BPAVDProfile_clip_gb <- clipLevel2BPAVDProfile(level2BPAVDProfile,polygon_spdf, split_by="id")

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
RH100max_st<-polyStatsLevel2AM(level2AM_clip_gb,func=max(RH100), id="poly_id")
head(RH100max_st)

# Computing a serie statistics for GEDI metrics stratified by polygon
RH100metrics_st<-polyStatsLevel2AM(level2AM_clip,func=mySetOfMetrics(RH100),
                      id="poly_id")
head(RH100metrics_st)

# Computing the max of the Total Plant Area Index 
pai_max<-polyStatsLevel2BVPM(level2BVPM_clip,func=max(pai), id=NULL)
pai_max

# Computing the serie of statistics of Foliage Clumping Index stratified by polygon
omega_metrics_st<-polyStatsLevel2BVPM(level2BVPM_clip,func=mySetOfMetrics(omega),
                      id="poly_id")
head(omega_metrics_st)
```
## Compute Grids with descriptive statistics of GEDI-derived Elevation and Height Metrics (Level2A)
```r
# Computing the serie of statistics of GEDI RH100 metric
#'RH100metrics<-gridStatsLevel2AM(level2AM = level2AM, func=mySetOfMetrics(RH100), res=0.0005)

# View maps
library(rasterVis)
library(viridis)
library(gridExtra)

rh100maps<-levelplot(RH100metrics,
                     layout=c(4, 1),
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
rh100maps
```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig5.png)

## Compute Grids with descriptive statistics of GEDI-derived Canopy Cover and Vertical Profile Metrics (Level2B)
```r
# Computing the max of the Total Plant Area Index only
level2BVPM$pai[level2BVPM$pai==-9999]<-NA # assing NA to -9999
pai_metrics<-gridStatsLevel2BVPM(level2BVPM = level2BVPM, func=mySetOfMetrics(pai), res=0.0005)

# View maps
pai_maps<-levelplot(pai_metrics,
                    layout=c(4, 1),
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

```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig6.png)


## Simulating GEDI full-waveform data from Airborne Laser Scanning (ALS) 3-D point cloud and extracting canopy derived metrics
```r
# specify the path to ALS data
lasfile_amazon <- system.file("extdata", "Amazon.las", package="rGEDI")
lasfile_cerrado <- system.file("extdata", "Cerrado.las", package="rGEDI")

# Reading and plot ALS file
library(lidR)
require(plot3D)
las_amazon<-readLAS(lasfile_amazon)
las_cerrado<-readLAS(lasfile_cerrado)

# Extracting plot center geolocations
xcenter_amazon = mean(las_amazon@bbox[1,])
ycenter_amazon = mean(las_amazon@bbox[2,])
xcenter_cerrado = mean(las_cerrado@bbox[1,])
ycenter_cerrado = mean(las_cerrado@bbox[2,])

# Simulating GEDI full-waveform
wf_amazon<-gediWFSimulator(input=lasfile_amazon,output=paste0(getwd(),"//gediWF_amazon_simulation.h5"),coords = c(xcenter_amazon, ycenter_amazon))
wf_cerrado<-gediWFSimulator(input=lasfile_cerrado,output=paste0(getwd(),"//gediWF_cerrado_simulation.h5"),coords = c(xcenter_cerrado, ycenter_cerrado))

# Plot ALS and GEDI simulated full-waveform
png("gediWf.png", width = 8, height = 6, units = 'in', res = 300)

par(mfrow=c(2,2), mar=c(4,4,0,0), oma=c(0,0,1,1),cex.axis = 1.2)
scatter3D(las_amazon@data$X,las_amazon@data$Y,las_amazon@data$Z,pch = 16,colkey = FALSE, main="",
          cex = 0.5,bty = "u",col.panel ="gray90",phi = 30,alpha=1,theta=45,
          col.grid = "gray50", xlab="UTM Easting (m)", ylab="UTM Northing (m)", zlab="Elevation (m)")

plot(wf_amazon, relative=TRUE, polygon=TRUE, type="l", lwd=2, col="forestgreen",
     xlab="", ylab="Elevation (m)", ylim=c(90,140))
grid()
scatter3D(las_cerrado@data$X,las_cerrado@data$Y,las_cerrado@data$Z,pch = 16,colkey = FALSE, main="",
          cex = 0.5,bty = "u",col.panel ="gray90",phi = 30,alpha=1,theta=45,
          col.grid = "gray50", xlab="UTM Easting (m)", ylab="UTM Northing (m)", zlab="Elevation (m)")

plot(wf_cerrado, relative=TRUE, polygon=TRUE, type="l", lwd=2, col="green",
xlab="Waveform Amplitude (%)", ylab="Elevation (m)", ylim=c(815,835))
grid()
dev.off()

# Extracting GEDI feull-waveform derived metrics
wf_amazon_metrics<-gediWFMetrics(input=wf_amazon@h5$filename,outRoot=getwd())
wf_cerrado_metrics<-gediWFMetrics(input=wf_cerrado@h5$filename,outRoot=getwd())

metrics<-rbind(wf_amazon_metrics,wf_cerrado_metrics)
rownames(metrics)<-c("Amazon","Cerrado")
head(metrics[,1:8])
```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig7.png)


# References
Dubayah, R., Blair, J.B., Goetz, S., Fatoyinbo, L., Hansen, M., Healey, S., Hofton, M., Hurtt, G.,         Kellner, J., Luthcke, S., & Armston, J. (2020) The Global Ecosystem Dynamics Investigation:         High-resolution laser ranging of the Earth’s forests and topography. Science of Remote             Sensing, p.100002.

Hancock, S., Armston, J., Hofton, M., Sun, X., Tang, H., Duncanson, L.I., Kellner,
       J.R. and Dubayah, R., 2019. The GEDI simulator: A large‐footprint waveform lidar simulator
       for calibration and validation of spaceborne missions. Earth and Space Science.
       https://doi.org/10.1029/2018EA000506

Silva, C. A.; Saatchi, S.; Alonso, M. G. ; Labriere, N. ; Klauberg, C. ; Ferraz, A. ; Meyer, V. ;        Jeffery, K. J. ; Abernethy, K. ; White, L. ; Zhao, K. ; Lewis, S. L. ; Hudak, A. T. (2018)         Comparison of Small- and Large-Footprint Lidar Characterization of Tropical Forest                 Aboveground Structure and Biomass: A Case Study from Central Gabon. IEEE Journal of Selected       Topics in Applied Earth Observations and Remote Sensing, p. 1-15.
      https://ieeexplore.ieee.org/document/8331845

GEDI webpage. Accessed on February 15 2020 https://gedi.umd.edu/   
GEDI01_Bv001. Accessed on February 15 2020 https://lpdaac.usgs.gov/products/gedi01_bv001/    GEDI02_Av001. Accessed on February 15 2020 https://lpdaac.usgs.gov/products/gedi02_av001/  
GEDI02_Bv001. Accessed on February 15 2020 https://lpdaac.usgs.gov/products/gedi02_bv001/  
GEDI Finder. Accessed on February 15 2020 https://lpdaacsvc.cr.usgs.gov/services/gedifinder

# Acknowledgements
University of Maryland and NASA Goddard Space Flight Center for developing GEDI mission

Brazilian National Council for Scientific and Technological Development (CNPq) for funding the project entitled "Mapping fuel load and simulation of fire behaviour and spread in the Cerrado biome using modeling and remote sensing technologies" and" leaded by Prof. Dr. Carine Klauberg and
Dr. Carlos Alberto Silva.

**rGEDI package has been not developted by the GEDI team. The authors assumes no responsibility whatsoever for its use, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.**
