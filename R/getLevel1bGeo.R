#'Read LVIS Level1b data
#'
#'@description This function reads LVIS level1b data
#'
#'
#'@param level1b H5File object from level1b
#'@param select character vector of desired columns to extract from GEDI hdf5 format, default c("latitude_bin0", "latitude_lastbin", "longitude_bin0", "longitude_lastbin", "shot_number")
#'
#'@return gedi.level1b class; S4 object of class gedi.level1b;
#'@author Carlos Alberto Silva. This function calls \emph{h5file} function from h5 package (Author: Mario Annau)
#'@seealso \code{\link[h5]{h5file}} in the \emph{h5} package.
#'
#'@import hdf5r
#'@import utils
#'@import sp
#'@importFrom hdf5r H5File
#'@importFrom data.table data.table
#'@importFrom sp SpatialPointsDataFrame
#'@export
getxGeo<-function(x,select=c("latitude_bin0", "latitude_lastbin", "longitude_bin0", "longitude_lastbin", "shot_number")) {
  x<-x@h5

  datasets<-hdf5r::list.datasets(x, recursive = T)
  datasets_names<-basename(datasets)

  selected<-datasets_names %in% select

  for ( i in select){
    if  ( i =="shot_number"){
      assign(i,bit64::as.integer64(NaN))
    } else {
      assign(i,numeric())
    }
  }

  dtse2<-datasets[selected][!grepl("geolocation/shot_number",datasets[selected])]

  # Set progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(dtse2), style = 3)
  i.s=0


  # i = "BEAM0000/shot_number"
  for ( i in dtse2){
    #print(i)
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    name_i<-basename(i)
    if ( name_i =="shot_number"){

      assign(name_i, bit64::c.integer64(get(name_i),x[[i]][]))

    } else {

      assign(name_i, c(get(name_i), x[[i]][]))

    }

  }

  if ( "shot_number" %in% select){
    x.dt<-data.table::data.table(get("shot_number"))[-1,]
    colnames(x.dt)<-"shot_number"
    select2<-select[!select[]=="shot_number"]

  } else{
    x.dt<-data.table::data.table(get(select[1]))
    colnames(x.dt)<-select[1]
    select2<-select[-1]
  }

  for ( i in select2){
    x.dt[,i]<-get(i)
  }

  #x.spdf<-sp::SpatialPointsDataFrame(cbind(as.numeric(x.dt$longitude_bin0),as.numeric(x.dt$latitude_bin0)),data=x.dt)
  #x.dt<- methods::new("gedi.x.dt", dt = x.dt)
  close(pb)
  return(x.dt)
}

