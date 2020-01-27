#'Display LVIS Waveform
#'
#'@description This function plots LVIS level1b and level 2 data
#'
#'@usage plotWaveform(level1b,level2,shotnum,plot2=TRUE,...)
#'
#'@param level1b h5file; S4 object of class H5File
#'@param level2 dataframe containing LVIS level 2 data
#'@param shotnum LVIS shot number to display
#'@param plot2 if TRUE, plot both Level1b and Level2. If FALSE, plot only Level1b data
#'@param ... passing arguments on to the plot function
#'@return Returns 2-D scatterplot of the LVIS waveform
#'@author Carlos Alberto Silva.
#'@examples
#'
#'# LVIS level1b file path
#'level1b_filepath_zip <- system.file("extdata", "LVIS_Mondah_level1b.zip", package="rLVIS")
#'unzip(level1b_filepath_zip, exdir = tempdir())
#'level1b_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.h5")
#'
#'# Reading LVIS level1b file
#'level1b<-readLevel1b(level1bpath=level1b_filepath)
#'
#'# LVIS level2 file path
#'level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
#'unzip(level2_filepath_zip, exdir = tempdir())
#'level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level2.txt")
#'
#'# Reading LVIS level1b file
#'level2_spdf<-readLevel2(level2path=level2_filepath, spdf=TRUE, glatlon=TRUE)
#'
#'#'Plotting LVIS waveforms
#'plotWaveform(level1b=level1b,level2=level2_spdf,
#'              shotnum=10964985,plot2=TRUE,xlab="Relative amplitude (%)", ylab="Height (m)")
#'
#'@importFrom grDevices colorRamp colors rgb
#'@importFrom graphics abline grid legend par plot polygon
#'@importFrom utils read.table
#'@export
plotWaveform<-function(level1b,level2,shotnum=10964985,plot2=TRUE,...) {

  add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }

  par(cex.axis=2)
  all_shotnums<-level1b['SHOTNUMBER'][]
  wave_idx<-which(all_shotnums==shotnum)
  wave_idx[wave_idx==TRUE]
  waveform<-level1b['RXWAVE'][wave_idx][1:1024]

  waveform<-(waveform/max(waveform))*100

  Z0<-level1b['Z0'][wave_idx]
  Z1023<-level1b['Z1023'][wave_idx]
  zstretch=rev(seq(Z1023,Z0,(Z0-Z1023)/1023))
  mytxtpathsub<-subset(level2@data,level2@data$SHOTNUMBER==shotnum)

  ZG<-mytxtpathsub[,6]
  ZT<-mytxtpathsub[,9]

  zstretch<-zstretch - ZG
  ZT<-ZT-ZG
  ZG<-0

  RH10<-mytxtpathsub[,10] #+ ZG
  RH25<-mytxtpathsub[,13] #+ ZG
  RH50<-mytxtpathsub[,18] #+ ZG
  RH75<-mytxtpathsub[,23] #+ ZG
  RH98<-mytxtpathsub[,30] #+ ZG
  RH100<-mytxtpathsub[,32] #+ ZG
  zmin<-ZG-(ZT-RH10)/6
  zmax=ZT + (ZT - RH10) /18
  x=zstretch>=zmin
  y=zstretch<=zmax
  z<-x==y
  waveform_crop=NULL
  zstrech_crop=NULL

  for ( i in 1:length(z)){
    if ( z[i]==TRUE) {
      waveform_crop[i]<-level1b['RXWAVE'][wave_idx][i]
      zstrech_crop[i]<-zstretch[i]
    }
  }

  #windows()
  if(plot2==TRUE) {
    par(mfrow=c(1,2))
    par(cex.axis=1.5)
    plot(waveform,zstretch, type="l", lwd=2, col="forestgreen")
    grid()
    polygon(c(waveform,min(zstretch),min(zstretch)),c(zstretch,min(zstretch),max(zstretch)),col="forestgreen")
    abline(h=ZG, lwd=2, col="blue")
    abline(h=RH10, lwd=2, col="green")
    abline(h=RH25, lwd=2, col="yellow")
    abline(h=RH50, lwd=2, col="gray")
    abline(h=RH75, lwd=2, col="orange")
    abline(h=RH98, lwd=2, lty=2,col="black")
    abline(h=RH100, lwd=2, col="red")

  }

  #browser()
  waveform_crop<-(as.numeric(waveform_crop)/as.numeric(summary(waveform_crop))[6])*100
  waveform_crop<- waveform_crop-as.numeric(summary(waveform_crop))[1]

  #polygon(c(0,as.numeric(summary(waveform_crop))[1],waveform_crop,as.numeric(summary(waveform_crop))[1],0,0),
  #        c(as.numeric(summary(zstrech_crop))[1],as.numeric(summary(zstrech_crop))[1],zstrech_crop,as.numeric(summary(zstrech_crop))[6],as.numeric(summary(zstrech_crop))[6],
  #          as.numeric(summary(zstrech_crop))[1]),col="forestgreen")

  xstart<-waveform_crop[which(zstrech_crop==min(zstrech_crop, na.rm=T))]
  xend<-waveform_crop[which(zstrech_crop==max(zstrech_crop, na.rm=T))]

  xl<-c(0,0,xstart,rev(waveform_crop),xend,0)
  yl<-c(max(zstrech_crop, na.rm=T),min(zstrech_crop, na.rm=T),min(zstrech_crop, na.rm=T),rev(zstrech_crop),max(zstrech_crop, na.rm=T),max(zstrech_crop, na.rm=T))

  par(cex.axis=1.5)
  plot(xl,yl, lwd=2, col="forestgreen", type="l",...)#,xlim=c(0,100))# ylim=c(min(zstretch),max(zstretch)))
  grid()
  #points(xl,yl)
  polygon(xl,yl,col="forestgreen")

  abline(h=ZG, lwd=2, col="blue")
  abline(h=RH10, lwd=2, col="green")
  abline(h=RH25, lwd=2, col="yellow")
  abline(h=RH50, lwd=2, col="gray")
  abline(h=RH75, lwd=2, col="orange")
  abline(h=RH98, lwd=2, lty=2,col="black")
  abline(h=RH100, lwd=2, col="red")
  par(xpd=TRUE)
  add_legend("topright",horiz=T,legend=c("ZG","RH10","RH25","RH50","RH75","RH98","RH100"),
             col=c("blue","green","yellow","gray","orange","black","red"),lwd=2,lty=c(rep(1,5),2,1),bty="n",cex=0.8)
  par(xpd=FALSE)
  #return(data.frame(SHOTNUMBER=rep(shotnum,1024),waveform,zstretch,min(yl, na.rm=T),max(yl, na.rm=T)))
}
