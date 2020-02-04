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
getLevel1BGeo<-function(level1b,select=c("latitude_bin0", "latitude_lastbin", "longitude_bin0", "longitude_lastbin", "shot_number")) {
  level1b<-level1b@h5

  datasets<-hdf5r::list.datasets(level1b, recursive = T)
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

      assign(name_i, bit64::c.integer64(get(name_i),level1b[[i]][]))

    } else {

      assign(name_i, c(get(name_i), level1b[[i]][]))

    }

  }

  if ( "shot_number" %in% select){
    level1b.dt<-data.table::data.table(get("shot_number"))[-1,]
    colnames(level1b.dt)<-"shot_number"
    select2<-select[!select[]=="shot_number"]

  } else{
    level1b.dt<-data.table::data.table(get(select[1]))
    colnames(level1b.dt)<-select[1]
    select2<-select[-1]
  }

  for ( i in select2){
    level1b.dt[,i]<-get(i)
  }

  #level1b.spdf<-sp::SpatialPointsDataFrame(cbind(as.numeric(level1b.dt$longitude_bin0),as.numeric(level1b.dt$latitude_bin0)),data=level1b.dt)
  #level1b.dt<- methods::new("gedi.level1b.dt", dt = level1b.dt)
  close(pb)
  return(level1b.dt)
}

