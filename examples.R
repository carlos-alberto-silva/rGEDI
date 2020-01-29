#require(h5)
require(hdf5r)
#install.packages("hdf5r")
#h5_data <- hdf5r::H5File$new(level1bpath, mode = 'r')
#h5_data[["BEAM0000/beam"]][]
#Level1b<- h5::h5file(level1bpath, 'a')

level1bpath<-"C:\\Users\\carlo\\Documents\\GEDI_package\\GEDI01_B_2019109163004_O01985_T02206_02_003_01.h5"

GEDI01_B<-readLevel1b(level1bpath)


par(mfrow=c(1,2))
par(cex.axis=1.5)
plot(GEDI01_B,shot_number="19850022900500000",relative=TRUE,polygon=TRUE,type="l", lwd=2, col="red")
grid()
plot(GEDI01_B,shot_number="19850022900500000",relative=FALSE,polygon=TRUE,type="l", lwd=2, col="forestgreen")
grid()

spdfs<-level1bSPDF(GEDI01_B@h5)

head(spdfs@data)

plot(spdfs)
level1b<-GEDI01_B
y<-kkk@level1b.dt[,1]
x<-kkk@level1b.dt[,2]

# Load the library
library(leaflet)

leaflet() %>%
  addCircleMarkers(x, y,
                   radius = 3,
                   opacity = 100,
                   color = "white")  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery)

plot(kkk@level1b.dt[,1],kkk@level1b.dt[,2])

level1b<- h5::h5file(level1bpath, 'a')
class(level1b)

yBEAM0000<-level1b["BEAM0000"]["geolocation"]["latitude_bin0"][]
xBEAM0000<-level1b["BEAM0000"]["geolocation"]["longitude_bin0"][]
