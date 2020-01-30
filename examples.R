#require(h5)
require(hdf5r)
#install.packages("hdf5r")
#h5_data <- hdf5r::H5File$new(level1bpath, mode = 'r')
#h5_data[["BEAM0000/beam"]][]
#Level1b<- h5::h5file(level1bpath, 'a')


#### function readLevel1b
level1bpath<-"C:\\Users\\carlo\\Documents\\GEDI_package\\GEDI01_B_2019109163004_O01985_T02206_02_003_01.h5"

level1b<-readLevel1b(level1bpath)


### plot waveform
windows()
par(mfrow=c(1,2))
par(cex.axis=1.5)
plot(level1b,shot_number="19850022900500000",relative=TRUE,polygon=TRUE,type="l", lwd=2, col="forestgreen")
grid()
plot(level1b,shot_number="19850022900500000",relative=FALSE,polygon=TRUE,type="l", lwd=2, col="forestgreen")
grid()

### level1b to dt
level1b_df<-level1b2dt(level1b,select="all")
head(level1b_df)

### plot Level1b as dt
require(rgdal)
worldshp<-readOGR("C:\\Users\\carlo\\Downloads\\countries_shp\\countries.shp")

windows()
plot(worldshp)
points(level1b_df$longitude_bin0,level1b_df$latitude_bin0, pch=".", col="green")

x<-as.numeric(spdfs@dt$longitude_bin0)
y<-as.numeric(spdfs@dt$latitude_bin0)
# Load the library
library(leaflet)

leaflet() %>%
  addCircleMarkers(x, y,
                   radius = 3,
                   opacity = 100,
                   color = "white")  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery)
