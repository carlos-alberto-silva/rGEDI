#require(h5)
library(hdf5r)
library(raster)
library(sp)
library(rgdal)


#install.packages("hdf5r")
#h5_data <- hdf5r::H5File$new(level1bpath, mode = 'r')
#h5_data[["BEAM0000/beam"]][]
#Level1b<- h5::h5file(level1bpath, 'a')


#### function readLevel1b
level1bpath<-"E:\\GEDI01_B_2019109163004_O01985_T02206_02_003_01.h5"
level1b<-readLevel1B(level1bpath)

### plot waveform
x<-getLevel1BWF(level1b,shot_number="19850022900500000")
windows()
par(mfrow=c(1,2))
par(cex.axis=1.5)
plot(x,relative=FALSE,polygon=TRUE,type="l", lwd=2, col="forestgreen", xlab="", ylab="Elevation (m)")

### level1b to dt
level1bGeo<-getLevel1BGeo(level1b,select=c("latitude_bin0","latitude_lastbin","longitude_bin0","longitude_lastbin","shot_number"))
#require(rgdal)
head(level1bGeo)

#write.table(level1bdt2@dt,"C:\\Users\\carlo\\Downloads\\countries_shp\\gedi1b_2.txt", sep=" ")
#writeOGR(points,"C:\\Users\\carlo\\Downloads\\countries_shp","gedi1b_2", drive="ESRI Shapefile")
### plot Level1b as dt
#require(rgdal)
worldshp<-rgdal::readOGR("C:\\Users\\carlo\\Downloads\\countries_shp\\countries.shp")
#windows()
#plot(worldshp)
#points(level1b_df$longitude_bin0,level1b_df$latitude_bin0, pch=".", col="green")

### clip level1bdt
require(rgdal)
polygon_spdf<-raster::shapefile("C:\\Users\\carlo\\Downloads\\countries_shp\\study_area.shp")

require(raster)
# Rectangle area for cliping
xleft = -116.4683
xright = -116.4583
ybottom = 46.75208
ytop = 46.84229

# clip by extent boundary box
level1b_clip<-clipLevel1Bh5(level1b,xleft, xright, ybottom, ytop)

# clip by geometry
level1b_clip<-clipLevel1Bh5Geometry(level1b,polygon_spdf = polygon_spdf)
geo<-getLevel1BGeo(clipLevel1BGeometrytest)

windows()
# Load the library
library(leaflet)
leaflet() %>%
  addCircleMarkers(geo$longitude_bin0,
                   geo$latitude_bin0,
                   radius = 1,
                   opacity = 1,
                   color = "red")  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addPolygons(data=polygon_spdf,weight=1,col = 'white',
              opacity = 1, fillOpacity = 0) %>%
  addProviderTiles(providers$Esri.WorldImagery)

###########################
# read level2A

level2a<-readLevel2A("C:\\Users\\carlo\\Documents\\GEDI_package\\GEDI02_A_2019108002011_O01959_T03909_02_001_01.h5")
level2a2<-readLevel2A("E:\\GEDI02_A_2019108015252_O01960_T03910_02_001_01.h5")

rhmetrics1<-getLevel2AM(level2a)
rhmetrics2<-getLevel2AM(level2a2)

rhmetrics<-data.table::rbindlist(list(rhmetrics1,rhmetrics2))
maps<-level2AGridStats(x=rhmetrics,func=mySetOfMetrics(rh100),res = 0.5)

windows()
plot(maps)

#spdfg<-sp::SpatialPointsDataFrame(rhmetrics[,4:3],data=rhmetrics[,4:3])
#sp::plot(spdfg)

###########################
# read level2B

#level2b<-readLevel2B("C:\\Users\\carlo\\Documents\\GEDI_package\\GEDI02_B_2019108045815_O01962_T01066_02_001_01.h5")

level2b1<-readLevel2B("E:\\GEDI02_B_2019108002011_O01959_T03909_02_001_01.h5")
level2b2<-readLevel2B("E:\\GEDI02_B_2019108015252_O01960_T03910_02_001_01.h5")
level2b3<-readLevel2B("E:\\GEDI02_B_2019108032534_O01961_T03911_02_001_01.h5")
level2b4<-readLevel2B("E:\\GEDI02_B_2019108045815_O01962_T01066_02_001_01.h5")
level2b5<-readLevel2B("E:\\GEDI02_B_2019108063056_O01963_T01067_02_001_01.h5")
level2b6<-readLevel2B("E:\\GEDI02_B_2019108080338_O01964_T05337_02_001_01.h5")
level2b7<-readLevel2B("E:\\GEDI02_B_2019108110900_O01966_T02493_02_001_01.h5")

head(vpm_metrics1)
vpm_metrics1<-getLevel2BVPM(level2b1)
vpm_metrics2<-getLevel2BVPM(level2b2)
vpm_metrics3<-getLevel2BVPM(level2b3)
vpm_metrics4<-getLevel2BVPM(level2b4)
vpm_metrics5<-getLevel2BVPM(level2b5)
vpm_metrics6<-getLevel2BVPM(level2b6)
vpm_metrics7<-getLevel2BVPM(level2b7)

x<-data.table::rbindlist(list(vpm_metrics1,vpm_metrics2,vpm_metrics3,
                              vpm_metrics4,vpm_metrics5,vpm_metrics6,vpm_metrics6))

ids<-1:nrow(x)
x<-x[ids[!x[,"pai"]==-9999],]

maps<-level2BVPMGridStats(x=x,func=mean(pai), res = 0.5)

plot(maps)
