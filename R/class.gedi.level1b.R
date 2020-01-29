#setClass("gedi.level1b", representation(h5="H5File",level1b.spdf='SpatialPointsDataFrame'))
require(h5)

setClass("gedi.level1b", representation(h5="H5File"))
setClass("gedi.level1bSPDF", slots=c(spdf="SpatialPointsDataFrame"))


setMethod("plot", signature("gedi.level1b", y = "missing"), function(x,shot_number,relative=TRUE,polygon=FALSE,...) {
    level1b<-x@H5File
    groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                       list.groups(level1b, recursive = F)), value = T)
    k<-"BEAM1011"
    for ( k in groups_id){
      gid<-max(level1b[[paste0(k,"/shot_number")]][]==shot_number)
      if (gid==1) {i=k}
    }

  shot_number_i<-level1b[[paste0(i,"/shot_number")]][]
  shot_number_id<-which(shot_number_i[]==shot_number)
  elevation_bin0<-level1b[[paste0(i,"/geolocation/elevation_bin0")]][]
  elevation_lastbin<-level1b[[paste0(i,"/geolocation/elevation_lastbin")]][]
  rx_sample_count<-level1b[[paste0(i,"/rx_sample_count")]][]
  rx_sample_start_index<-level1b[[paste0(i,"/rx_sample_start_index")]][]
  rxwaveform_i<-level1b[[paste0(i,"/rxwaveform")]][rx_sample_start_index[shot_number_id]:(rx_sample_start_index[shot_number_id]+rx_sample_count[shot_number_id]-1)]
  rxwaveform_inorm<-(rxwaveform_i-min(rxwaveform_i))/(max(rxwaveform_i)-min(rxwaveform_i))*100
  elevation_bin0_i<-elevation_bin0[shot_number_id]
  elevation_lastbin_i<-elevation_lastbin[shot_number_id]
  z=rev(seq(elevation_lastbin_i,elevation_bin0_i,(elevation_bin0_i-elevation_lastbin_i)/rx_sample_count[shot_number_id]))[-1]

  if (relative==TRUE){x=rxwaveform_inorm } else{
      x=rxwaveform_i
    }
  if (polygon==TRUE){

    xstart<-x[which(z==min(z, na.rm=T))]
    xend<-x[which(z==max(z, na.rm=T))]

    xl<-c(min(x),min(x),xstart,rev(x),xend,min(x))
    yl<-c(max(z, na.rm=T),min(z, na.rm=T),min(z, na.rm=T),rev(z),max(z, na.rm=T),max(z, na.rm=T))

    plot(xl,yl,...)
    polygon(xl,yl,...)
  } else {
    plot(x=x,y=z,...)
   }
  }
)


setMethod("plot", signature("gedi.level1bSPDF", y = "missing"), function(x,...) {
    spdf_xy<-x@spdf
    leaflet(spdf_xy) %>%
    #addCircleMarkers(spdf_xy, y,
    #                 radius = 3,
    #                 opacity = 100,
    #                 color = "white")  %>%
    addScaleBar(options = list(imperial = FALSE)) %>%
    addProviderTiles(providers$Esri.WorldImagery)
})
