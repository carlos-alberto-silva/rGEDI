#'Read LVIS level2a data
#'
#'@description This function reads LVIS level2a data
#'
#'
#'@param level2a H5File object from level2a
#'@return H5File; S4 object of class H5File;
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
level2A2dt<-function(level2a,select=c("latitude_bin0","latitude_lastbin","shot_number")) {

  level2a<-level2a@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(level2a, recursive = F)), value = T)
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)

  datasets<-hdf5r::list.datasets(level2a, recursive = T)
  datasets_names<-basename(datasets)

  selected<-datasets_names %in% select

  for ( i in select){
    if  ( i =="shot_number"){
      assign(i,bit64::as.integer64(NaN))
    } else {
      assign(i,numeric())
    }
  }

  #i="BEAM0000/geolocation/shot_number"
  i.s=0
  #i="BEAM0010"
  dtse2<-datasets[selected][!grepl("geolocation/shot_number",datasets[selected])]


  # i = "BEAM0000/shot_number"
  for ( i in dtse2){
    #print(i)
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    name_i<-basename(i)
    if ( name_i =="shot_number"){

      assign(name_i, bit64::c.integer64(get(name_i),level2a[[i]][]))

    } else {

      assign(name_i, c(get(name_i), level2a[[i]][]))

    }

  }

  if ( "shot_number" %in% select){
    level2a.dt<-data.table::data.table(get("shot_number"))[-1,]
    colnames(level2a.dt)<-"shot_number"
    select2<-select[!select[]=="shot_number"]

  } else{
    level2a.dt<-data.table::data.table(get(select[1]))
    colnames(level2a.dt)<-select[1]
    select2<-select[-1]
  }

  for ( i in select2){
    level2a.dt[,i]<-get(i)
  }

  format(level2a.dt$shot_number[1], scientific = F)
  format(shot_number[2], scientific = F)

  head(level2a.dt)
  length(latitude_bin0)
  length(latitude_bin0)
  length(latitude_lastbin)
  nrow(shot_number)


  format(shot_number[2], scientific = F)


  #level2a.spdf<-sp::SpatialPointsDataFrame(cbind(as.numeric(level2a.dt$longitude_bin0),as.numeric(level2a.dt$latitude_bin0)),data=level2a.dt)
  level2a.dt<- methods::new("gedi.level2a.dt", dt = level2a.dt)
  close(pb)
  return(level2a.dt@dt)
}

