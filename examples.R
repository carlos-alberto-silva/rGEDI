#require(h5)
require(hdf5r)
#install.packages("hdf5r")
#h5_data <- hdf5r::H5File$new(level1bpath, mode = 'r')
#h5_data[["BEAM0000/beam"]][]
#Level1b<- h5::h5file(level1bpath, 'a')


#### function readLevel1b
level1bpath<-"C:\\Users\\carlo\\Documents\\GEDI_package\\GEDI01_B_2019109163004_O01985_T02206_02_003_01.h5"
level1b<-readLevel1B(level1bpath)


### plot waveform
windows()
par(mfrow=c(1,2))
par(cex.axis=1.5)
plot(level1b,shot_number="19850022900500000",relative=TRUE,polygon=TRUE,type="l", lwd=2, col="forestgreen")
grid()
plot(level1b,shot_number="19850022900500000",relative=FALSE,polygon=TRUE,type="l", lwd=2, col="forestgreen")
grid()

### level1b to dt
level1bdt2<-level1B2dt(level1b,select=c("latitude_bin0","latitude_lastbin","longitude_bin0","longitude_lastbin"))
#require(rgdal)
#writeOGR(points,"C:\\Users\\carlo\\Downloads\\countries_shp","gedi1b", drive="ESRI Shapefile")
### plot Level1b as dt
#require(rgdal)
#worldshp<-readOGR("C:\\Users\\carlo\\Downloads\\countries_shp\\countries.shp")
#windows()
#plot(worldshp)
#points(level1b_df$longitude_bin0,level1b_df$latitude_bin0, pch=".", col="green")

### clip level1bdt
polygon_spdf<-readOGR("C:\\Users\\carlo\\Downloads\\countries_shp\\study_area.shp")

# Rectangle area for cliping
xmin = -116.4683
xmax = -116.3177
ymin = 46.75208
ymax = 46.84229

ext<-extent(c(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax))
clipLevel1Bdtest<-clipLevel1Bdt(level1bdt2,extent=ext)
clipLevel1BdtGeometrytest<-clipLevel1BdtGeometry(level1bdt2,polygon_spdf = polygon_spdf)

# Load the library
library(leaflet)
leaflet() %>%
  addCircleMarkers(clipLevel1BdtGeometrytest@dt$longitude_bin0,
                   clipLevel1BdtGeometrytest@dt$latitude_bin0,
                   radius = 1,
                   opacity = 1,
                   color = "red")  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addPolygons(data=polygon_spdf,weight=1,col = 'white',
              opacity = 1, fillOpacity = 0) %>%
  addProviderTiles(providers$Esri.WorldImagery)
