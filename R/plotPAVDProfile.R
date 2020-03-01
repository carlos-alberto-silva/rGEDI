#'Plot Plant Area Volume Density Profile
#'
#'@description This functions plots Plant Area Volume Density profile (GEDI level2B)
#'
#'@param level2BPAVDProfile A GEDI Level2B object (output of \code{\link[rGEDI:getLevel2BPAVDProfile]{getLevel2BPAVDProfile}} function). A S4 object of class "data.table".
#'@param beam select GEDI beam. Default is "BEAM0101". See details section.
#'@param elev if TRUE, elevation will be used for plotting the PAVD profile. Otherwise,
#'height will be used instead.
#'
#'@return Returns a ggplot object. See \code{\link[ggplot2:ggplot]{ggplot}}
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
#'
#'@details list of GEDI beams. See output of \code{\link[rGEDI:getLevel2BPAVDProfile]{getLevel2BPAVDProfile}} function.
#'\itemize{
#'\item \emph{BEAM0000}
#'\item \emph{BEAM0001}
#'\item \emph{BEAM0010}
#'\item \emph{BEAM0011}
#'\item \emph{BEAM0101}
#'\item \emph{BEAM0110}
#'\item \emph{BEAM1000}
#'\item \emph{BEAM1011}
#'}
#'
#'@examples
#'\dontrun{
#'# specify the path to download GEDI example dataset
#'outdir<-getwd()
#'
#'# downloading GEDI example dataset (zip file)
#'download.file(
#'              paste0(
#'                     "https://github.com/carlos-alberto-silva/rGEDI/",
#'                     "releases/download/examples/examples.zip"
#'              ),
#'              destfile=paste0(outdir,"/examples.zip"))
#'
#'# unzip the file
#'unzip(paste0(outdir,"\\examples.zip"))
#'
#'# specify the path to GEDI level2B data
#'level2bpath = paste0(outdir,"\\GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.h5")
#'
#'# Reading GEDI level1B file
#'level2b<-readLevel2b(gedilevel2b)
#'
#'# Get Plant Area Volume Density profile
#'level2BPAVDProfile<-getLevel2BPAVDProfile(level2b)
#'
#'# Plot Level2B PAVD Profile
#'gprofile<-plotPAVDProfile(level2BPAVDProfile, beam="BEAM0101", elev=TRUE)
#'
#'
#'}
#'@export
plotPAVDProfile<-function(level2BPAVDProfile, beam="BEAM0101", elev=TRUE){
  #require(ggplot2)

  rids<-1:nrow(level2BPAVDProfile)
  rids<-rids[level2BPAVDProfile$beam==beam]
  level2BPAVDProfile_sub<-level2BPAVDProfile[rids,]
  level2BPAVDProfile_sub$height_bin0[level2BPAVDProfile_sub$height_bin0<0]<-0

  n0<-nrow(level2BPAVDProfile_sub)
  dft<-data.table::melt(level2BPAVDProfile_sub[,c(2,6,8,9:38)], id.vars=c("shot_number","elev_lowestmode", "height_bin0"), variable.name="pavd", value.name="value")
  dft$rowids<-rep(1:n0,30)
  df <- as.data.frame(lapply(dft, rep, rep(5,nrow(dft))))
  n<-nrow(df)
  seqs<-seq(0,150,5)
  hids<-NULL
  for ( i in 1:30){
    hids<-c(hids,rep(seq(seqs[i]+1,seqs[i+1]),n0))
  }
  df$value[df$value<0]<-0
  df$hids<-hids
  df<-df[df$value>0,]
  dif<-(df$elev_lowestmode+df$height_bin0) - (df$hids+df$elev_lowestmode)
  df<-df[dif>0,]
  xp<-((df$rowids*60)-60)/1000

  if( elev==TRUE){
    yp<-round(df$elev_lowestmode+df$hids)
    xsl<-((1:nrow(level2BPAVDProfile_sub)*60)-60)/1000
    gg <- ggplot2::ggplot()+
      geom_tile(aes(x=xp, y=yp,fill= df$value))+
      geom_line(aes(x=xsl, y=round(level2BPAVDProfile_sub$height_bin0+level2BPAVDProfile_sub$elev_lowestmode)),color = "black") +
      geom_line(aes(x=xsl, y=round(level2BPAVDProfile_sub$elev_lowestmode)),color = "black") +
      scale_fill_gradientn(colours = brewer.pal(n = 8, name = "Greens"))+
      xlab("Distance Along Track (km)") + ylab("Elevation (m)") +
      labs(fill=expression(PAVD~(m^2/m^3)))+
      theme(panel.border = element_rect(colour = "gray70", fill=NA, size=0.2))
    print(gg)
  } else {
    yp<-round(df$hids)
    xsl<-((1:nrow(level2BPAVDProfile_sub)*60)-60)/1000

    #require(ggplot2)
    gg <- ggplot()+
      geom_tile(aes(x=xp, y=yp,fill= df$value))+
      geom_line(aes(x=xsl, y=round(level2BPAVDProfile_sub$height_bin0)),color = "black") +
      scale_fill_gradientn(colours = brewer.pal(n = 8, name = "Greens"))+
      xlab("Distance Along Track (km)") + ylab("Height (m)") +
      labs(fill=expression(PAVD~(m^2/m^3)))+
      print(gg)
  }
  return(gg)
}

