#'Plot Plant Area Volume Density Profile
#'
#'@description This functions plots Plant Area Volume Density profile (GEDI level2B)
#'
#'@param level2BPAVDProfile A GEDI Level2B object (output of \code{\link[rGEDI:getLevel2BPAVDProfile]{getLevel2BPAVDProfile}} function). A S4 object of class "data.table".
#'@param beam Select GEDI beam. Default is "BEAM0101". See details section.
#'@param elev If TRUE, elevation will be used for plotting the PAVD profile. Otherwise,
#'height will be used instead.
#'
#'@return Returns a ggplot object. See \code{\link[ggplot2:ggplot]{ggplot}} package.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
#'
#'@details list of GEDI beams. See the output of \code{\link[rGEDI:getLevel2BPAVDProfile]{getLevel2BPAVDProfile}} function.
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
#'# specify the path to GEDI level2B data (zip file)
#'level2B_fp_zip <- system.file("extdata",
#'                   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Bpath <- unzip(level2B_fp_zip,exdir = dirname(level2B_fp_zip))
#'
#'# Reading GEDI level2B data (h5 file)
#'level2b<-readLevel2B(level2Bpath=level2Bpath)
#'
#'# Get Plant Area Volume Density profile
#'level2BPAVDProfile<-getLevel2BPAVDProfile(level2b)
#'
#'# Plot Level2B PAVD Profile
#'gprofile<-plotPAVDProfile(level2BPAVDProfile, beam="BEAM0101", elev=TRUE)
#'
#'
#'close(level2b)
#'@export
plotPAVDProfile<-function(level2BPAVDProfile, beam="BEAM0101", elev=TRUE){
  #require(ggplot2)
  #require(RColorBrewer)

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
  #df<-df[df$value>0,]

  if( elev==TRUE){
    dif<-(df$elev_lowestmode+df$height_bin0) - (df$hids+df$elev_lowestmode)
    df<-df[dif>0,]
    xp<-((df$rowids*60)-60)/1000
    yp<-round(df$elev_lowestmode+df$hids)

    xsl<-unique(xp)
    yl1<-tapply(yp,df$rowids,max) +0.5
    yl2<-tapply(yp,df$rowids,min) -0.5


    #xsl<-((1:nrow(level2BPAVDProfile_sub)*60)-60)/1000
    #yl1<-round(level2BPAVDProfile_sub$height_bin0+level2BPAVDProfile_sub$elev_lowestmode)
    #yl2<-round(level2BPAVDProfile_sub$elev_lowestmode)

    gg <- ggplot2::ggplot()+
      geom_tile(aes(x=xp, y=yp,fill= df$value))+
      scale_fill_gradientn(colours = brewer.pal(n = 8, name = "Greens"))+
      xlab("Distance Along Track (km)") + ylab("Elevation (m)")+
      geom_line(mapping = aes(x = xsl,y=yl1, color = "Canopy \nTop Height (m)"))+#,size=1) +
      geom_line(mapping = aes(x = xsl,y=yl2, color = "Ground \nElevation (m)"))+#,size=1) +
      scale_color_manual(name="",values = c("forestgreen", "black"))+
      theme(panel.border = element_rect(colour = "gray70", fill=NA, size=0.2))+
      labs(fill=expression(PAVD~(m^2/m^3)))+
      theme(legend.key.height=unit(1, "cm"))

      print(gg)

  } else {

    dif<-df$height_bin0 - df$hids
    df<-df[dif>0,]
    xp<-((df$rowids*60)-60)/1000
    yp<-df$hids-0.5

    yl<-tapply(df$hids,df$rowids,max)#+0.5
    xl<-((unique(df$rowids)*60)-60)/1000

    #xl<-((1:nrow(level2BPAVDProfile_sub)*60)-60)/1000
    #yl<-round(level2BPAVDProfile_sub$height_bin0)

    gg <- ggplot()+
      geom_tile(aes(x=xp, y=yp,fill= df$value))+
      geom_line(mapping = aes(x = xl,y=yl, color = "Canopy \nTop Height (m)"))+#,size=1) +
      scale_fill_gradientn(colours = brewer.pal(n = 8, name = "Greens"))+
      xlab("Distance Along Track (km)") + ylab("Height (m)") +
      theme(panel.border = element_rect(colour = "gray70", fill=NA, size=0.2))+
      labs(fill=expression(PAVD~(m^2/m^3)))+
      theme(legend.key.height=unit(1, "cm"))+
      scale_color_manual(name="",values = c("forestgreen", "black"))

      print(gg)
  }
  return(gg)
}

