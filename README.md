![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig1.png)<br/>

![Github](https://img.shields.io/badge/CRAN-0.0.1-green.svg)
![Github](https://img.shields.io/badge/Github-0.0.1-green.svg)
[![Rdoc](http://www.rdocumentation.org/badges/version/rGEDI)](http://www.rdocumentation.org/packages/ForestGapR)
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
![R_Forge](https://img.shields.io/badge/R_Forge-0.0.1-green.svg) 
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rGEDI)

rGEDI: An R Package for NASA's Global Ecosystem Dynamics Investigation (GEDI) data 
downloading, visualizing and processing. 

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
xmin<--45.62009
ymin<--15.71073
xmax<--44.51999
ymax<--14.46284

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
# specify the path to GEDI data
level1bpath <- system.file("extdata", "GEDIexample_level01B.h5", package="rGEDI")
level2apath <- system.file("extdata", "GEDIexample_level02A.h5", package="rGEDI")
level2bpath <- system.file("extdata", "GEDIexample_level02B.h5", package="rGEDI")

# Reading GEDI data
gedilevel1b<-readLevel1B(level1bpath)
gedilevel2a<-readLevel1B(level2apath)
gedilevel2b<-readLevel1B(level2bpath)
```

## Get GEDI Pulse Full-Waveform Geolocation (GEDI Level1B)
```r
level1bGeo<-getLevel1BGeo(gedilevel1b,select=c("latitude_bin0", "longitude_bin0","shot_number"))
head(level1bGeo)

    ##   id  min   max   mean        sd
    ##   1 2.06 65.40 32.48820  9.996999
    ##   3 2.47 57.26 37.95028 12.054305
    ##   2 6.92 59.78 37.23889  5.176369
    
library(leaflet)
leaflet() %>%
  addCircleMarkers(llevel1bGeo@data$longitude_bin0,
                   level1bGeo@data$latitude_bin0,
                   radius = 1,
                   opacity = 1,
                   color = "green")  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery)
```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig2.png)

## Get GEDI Pulse Full-Waveform (GEDI Level1B)
```r
# Extracting GEDI full-waveform for a giving shotnumber
wf <- getLevel1BWF(gedilevel1b, shot_number="19850022900500000")

# Plot full-waveform
par(mfrow = c(1,2), cex.axis = 1.5)

plot(wf, relative=FALSE, polygon=TRUE, type="l", lwd=2, col="forestgreen",
xlab="Waveform Amplitude", ylab="Elevation (m)")

plot(wf, relative=TRUE, polygon=FALSE, type="l", lwd=2, col="forestgreen",
xlab="Waveform Amplitude (%)", ylab="Elevation (m)")
```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig3.png)


## Get GEDI Elevation and Height Metrics (GEDI Level2A)
```r
# Get GEDI Elevation and Height Metrics
level2AM<-getLevel2AM(level2a)
head(level2AM)

    ##   id  min   max   mean        sd
    ##   1 2.06 65.40 32.48820  9.996999
    ##   3 2.47 57.26 37.95028 12.054305
    ##   2 6.92 59.78 37.23889  5.176369
```

## Get GEDI Plant Area Index (PAI) Profile (GEDI Level2B)
```r
level2BPAIProfile<-getLevel2BPAIProfile(level2b)
head(level2BPAIProfile)

    ##   id  min   max   mean        sd
    ##   1 2.06 65.40 32.48820  9.996999
    ##   3 2.47 57.26 37.95028 12.054305
    ##   2 6.92 59.78 37.23889  5.176369
```
## Get GEDI Plant Area Volume Density (PAVD) Index (GEDI Level2B)
```r
level2BPAVDProfile<-getLevel2BPAVDProfile(level2b)
head(level2BPAVDProfile)

    ##   id  min   max   mean        sd
    ##   1 2.06 65.40 32.48820  9.996999
    ##   3 2.47 57.26 37.95028 12.054305
    ##   2 6.92 59.78 37.23889  5.176369
```
## Get GEDI Vegetation Profile Biophysical Variables (GEDI Level2B)
```r
level2BVPM<-getLevel2BVPM(level2b)
head(level2BVPM)
```
# Clip GEDI data 
```r
## Clip GEDI data by coordinates
# Study area boundary box
xmin<--45.62009
ymin<--15.71073
xmax<--44.51999
ymax<--14.46284

## clipping GEDI data within boundary box

# Clip h5 files
level1b_clip_bb <- clipLevel1B(level1b,xleft, xright, ybottom, ytop)
level2a_clip_bb <- clipLevel2A(level2a,xleft, xright, ybottom, ytop)
level2b_clip_bb <- clipLevel2B(level2b,xleft, xright, ybottom, ytop)

# Clip gedi.level1b, gedi.level2a and gedi.level2b objects
level1BGeo_clip_bb <-clipLevel1BGeo(level1BGeo,xleft, xright, ybottom, ytop)
level2AM_clip_bb <- clipLevel2AM(level2AM,xleft, xright, ybottom, ytop)
level2BVPM_clip_bb <- clipLevel2BVPM(level2BVPM,xleft, xright, ybottom, ytop)
level1BPAIProfile_clip_bb <- clipLevel1BPAIProfile(level1BPAIProfile,xleft, xright, ybottom, ytop)
level2BPAVDProfile_clip_bb <- clipLevel2BPAVDProfile(level2BPAVDProfile,xleft, xright, ybottom, ytop)

## Clip GEDI data by geometry
# specify the path to shapefile for the study area
polygon_filepath <- system.file("extdata", "shp_np.shp", package="rGEDI")

# Reading shapefile as SpatialPolygonsDataFrame object
library(rgdal)
polygon_spdf<-readOGR(polygons_filepath)

# Clip h5 files
level1b_clip_gb <- clipLevel1BGeometry(level1b,polygon_spdf)
level2a_clip_gb <- clipLevel2AGeometry(level2a,polygon_spdf)
level2b_clip_gb <- clipLevel2BGeometry(level2b,polygon_spdf)

# Clip gedi.level1b, gedi.level2a and gedi.level2b objects
level1BGeo_clip_gb <- clipLevel1BGeo(level1BGeo,polygon_spdf)
level2AM_clip_gb <- clipLevel2AM(level2AM,polygon_spdf)
level2BVPM_clip_gb <- clipLevel2BVPM(level2BVPM,polygon_spdf)
level1BPAIProfile_clip_gb <- clipLevel1BPAIProfile(level1BPAIProfile,polygon_spdf)
level2BPAVDProfile_clip_gb <- clipLevel2BPAVDProfile(level2BPAVDProfile,polygon_spdf)

## visualizing clipped GEDI data
library(sp)
library(raster)
library(mapview)
library(leafsync)

# Setting spatial coordinates to create a spatial object
coordinates(level2AM_clip_bb) <- ~x+y
coordinates(level2AM_clip_gb) <- ~x+y
coordinates(level2BVPM_clip_bb) <- ~x+y
coordinates(level2BVPM_clip_gb) <- ~x+y

# Setting projection attributes to spatial object
proj4string(level2AM_clip_bb) <- CRS("+init=epsg:4674")
proj4string(level2AM_clip_gb) <- CRS("+init=epsg:4674")
proj4string(level2BVPM_clip_bb) <- CRS("+init=epsg:4674")
proj4string(level2BVPM_clip_gb) <- CRS("+init=epsg:4674")

# view GEDI-derived Canopy Height (RH100) and Plant Area Index
m1 <- mapview(level2AM_clip_bb, zcol = "RH100", burst = TRUE)
m2 <- mapview(level2AM_clip_bb, zcol = "RH100")
m3 <- mapview(level2BVPM_clip_bb, zcol = "pai", map.types = "Esri.WorldImagery")
m4 <- mapview(level2BVPM_clip_gb, zcol = "pai")

sync(m1, m2, m3, m4) # 4 panels synchronised
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
RH100max_st<-polyStatsLevel2AM(level2AM_clip,func=max(RH100), id="id")
head(RH100max_st)

# Computing a serie statistics for GEDI metrics stratified by polygon
RH100metrics_st<-polyStatsLevel2AM(level2AM_clip,func=mySetOfMetrics(RH100),
                      id=level2AM_clip_gb@data$id)
head(RH100metrics_st)

# Computing the max of the Total Plant Area Index 
pai_max<-polyStatsLevel2BVPM(level2BVPM_clip,func=max(pai), id=NULL)
pai_max

# Computing the serie of statistics of Foliage Clumping Index stratified by polygon
omega_metrics_st<-polyStatsLevel2BVPM(level2BVPM_clip,func=mySetOfMetrics(omega),
                      id=level2BM_clip_gb@data$id)
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
          layout=c(2, 2),
          margin=FALSE,
          colorkey=list(
            space='right',
            labels=list(at=seq(0, 50, 5), font=4),
            axis.line=list(col='black'),
            width=1),
          par.settings=list(
            strip.border=list(col='transparent'),
            strip.background=list(col='transparent'),
            axis.line=list(col='transparent')
          ),
          scales=list(draw=TRUE),
          col.regions=viridis,
          at=seq(0, 50, len=101),
          names.attr=c("min of RH100","max of RH100","mean of RH100", "sd of RH100"))

rh100maps
```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig5.png)


## Simulating GEDI full-waveform data from Airborne Laser Scanning (ALS) 3-D point cloud
and extracting canopy derived metrics
```r
# specify the path to ALS data
LASfile <- system.file("extdata", "LASexample1.las", package="rGEDI")

# Simulate GEDI full-waveform
wf<-gediWFSimulator(input=LASfile,output="gediSimulation.h5")

# Adding noise to GEDI full-waveform
wfn<-gediWFNoise(input=wf,output="gediSimulation_noise.h5")

# Extracting GEDI feull-waveform derived metrics
wfmetrics<-gediWFMetrics(input=wfn,outRoot=getwd())
head(wfmetrics)
```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig6.png)

### References
Dubayah, R., Blair, J.B., Goetz, S., Fatoyinbo, L., Hansen, M., Healey, S., Hofton, M., Hurtt, G.,         Kellner, J., Luthcke, S., & Armston, J. (2020) The Global Ecosystem Dynamics Investigation:         High-resolution laser ranging of the Earth’s forests and topography. Science of Remote             Sensing, p.100002.

Hancock, S., Armston, J., Hofton, M., Sun, X., Tang, H., Duncanson, L.I., Kellner,
       J.R. and Dubayah, R., 2019. The GEDI simulator: A large‐footprint waveform lidar simulator
       for calibration and validation of spaceborne missions. Earth and Space Science.
       https://doi.org/10.1029/2018EA000506


Silva, C. A.; Saatchi, S.; Alonso, M. G. ; Labriere, N. ; Klauberg, C. ; Ferraz, A. ; Meyer, V. ;        Jeffery, K. J. ; Abernethy, K. ; White, L. ; Zhao, K. ; Lewis, S. L. ; Hudak, A. T. (2018)         Comparison of Small- and Large-Footprint Lidar Characterization of Tropical Forest                 Aboveground Structure and Biomass: A Case Study from Central Gabon. IEEE Journal of Selected       Topics in Applied Earth Observations and Remote Sensing, p. 1-15.
      https://ieeexplore.ieee.org/document/8331845

GEDI webpage. Accessed on February 15 2020 <https://gedi.umd.edu/>
GEDI01_Bv001. Accessed on February 15 2020 <https://lpdaac.usgs.gov/products/gedi01_bv001/>
GEDI02_Av001. Accessed on February 15 2020 <https://lpdaac.usgs.gov/products/gedi02_av001/>
GEDI02_Bv001. Accessed on February 15 2020 <https://lpdaac.usgs.gov/products/gedi02_bv001/>

# Acknowledgements
University of Maryland and NASA Goddard Space Flight Center for developing GEDI mission.
Brazilian National Council for Scientific and Technological Development (CNPq) for funding the project entitled "Mapping fuel load and simulation of fire behaviour and spread in the Cerrado biome using modeling and remote sensing technologies" and" leaded by Prof. Dr. Carine Klauberg and
Dr. Carlos Alberto Silva.

rGEDI package has been not developted by the GEDI team. The authors assumes no responsibility whatsoever for its use, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.
