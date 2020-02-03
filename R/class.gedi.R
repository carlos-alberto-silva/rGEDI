#' Class for GEDI Level1B derived waveform
#'
#' @slot dt Object of class data.table
#'
#' @export
setClass(
  Class="gedi.waveform",
  slots=list(df = "data.frame")
)


setMethod("plot", signature("gedi.waveform", y = "missing"), function(x,relative=TRUE,polygon=FALSE,...) {
  waveform<-x@df
  z<-waveform$elevation
  if (relative==TRUE){
    wf<-(waveform$rxwaveform-min(waveform$rxwaveform))/
      (max(waveform$rxwaveform)-min(waveform$rxwaveform))*100
    } else{
      wf=waveform$rxwaveform
  }

    if (polygon==TRUE){
      wfstart<-wf[which(z==min(z, na.rm=T))]
      wfend<-wf[which(z==max(z, na.rm=T))]
      wfl<-c(min(wf),min(wf),wfstart,rev(wf),wfend,min(wf))
      zl<-c(max(z, na.rm=T),min(z, na.rm=T),min(z, na.rm=T),rev(z),max(z, na.rm=T),max(z, na.rm=T))
      plot(wfl,zl,...)
      suppressWarnings(polygon(wfl,zl,...))
    } else {
      plot(x=wf,y=z,...)
    }
})

