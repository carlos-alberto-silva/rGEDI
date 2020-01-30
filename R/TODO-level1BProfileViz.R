
beam=c("BEAM0000","BEAM0001","BEAM0010","BEAM0011","BEAM0101","BEAM0110",
       "BEAM1000","BEAM1011")

beam="BEAM0000"


myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
# Color by height
library("RColorBrewer")

windows()
level1BProfileViz<-function(level1b,extent,
                            beam=c("BEAM0000","BEAM0001","BEAM0010","BEAM0011","BEAM0101","BEAM0110",
                                    "BEAM1000","BEAM1011")){

  shot_numberx<-clipLevel1Bdtest@dt$shot_number[1]
  shot_numberx<-level1bx[["BEAM0000/shot_number"]][]
  lenShot<-length(shot_numberx)

  distx<-seq(0,lenShot*60,60)[-length(distx)]
  head(distx)

  xp<-NULL
  yp<-NULL
  zp<-NULL

  for ( i in 1:length(shot_numberx)){
    print(i)

    xyext_i<-PlotWaveform(level1b,shot_number=shot_numberx[i],relative=TRUE,polygon=FALSE,return=TRUE,plotWave = FALSE)
    mz<-mean(xyext_i[,2])
    xyext_i<-subset(xyext_i,xyext_i[,2]>mz-30 & xyext_i[,2]<mz+30)

    zp<-c(zp,xyext_i[,2])
    xp<-c(xp,rep(distx[i],nrow(xyext_i)))
    yp<-c(yp,xyext_i[,1])
  }

  hist(yp)
  head(yp)
  col <- myColorRamp(brewer.pal(5, "Greens"),yp)

  plot(xp,zp,col=col, xlim=c(0,1000))

  length(distx)
  length(shot_numberx)

  }
