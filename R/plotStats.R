#'Scatterplot of a 1:1 comparison
#'
#'@description Compute and plot stastical parameters of 1:1 comparison
#'
#'@aliases plotStats
#'
#'@param x h5file; S4 object of class H5File
#'@param y dataframe containing LVIS level 2 data
#'@param colours h5file; S4 object of class H5File
#'@param n size
#'@param xlim dataframe containing LVIS level 2 data
#'@param ylim h5file; S4 object of class H5File
#'@param colfill dataframe containing LVIS level 2 data
#'@param xlab h5file; S4 object of class H5File
#'@param ylab dataframe containing LVIS level 2 data
#'@param title h5file; S4 object of class H5File
#'@param legend.position dataframe containing LVIS level 2 data
#'@param stats.position dataframe containing LVIS level 2 data
#'@param base_size dataframe containing LVIS level 2 data
#'@param stat.size dataframe containing LVIS level 2 data
#'@param legend.size dataframe containing LVIS level 2 data
#'@param fit.line.col dataframe containing LVIS level 2 data
#'@param x_axis dataframe containing LVIS level 2 data
#'@param y_axis dataframe containing LVIS level 2 data
#'
#'@return Returns A list containing statistical parameters and the Scatterplot of the 1:1 comparison
#'
#'@author Carlos Alberto Silva.
#'
#'@examples
#'# Importing libraries
#'library(raster)
#'library(rasterVis)
#'library(viridis)
#'library(gridExtra)
#'
#'# Importing dataset
#'sf_agb1ha_path <- system.file("extdata", "sf_agb_1ha.tif", package="rLVIS")
#'lf_agb1ha_path <- system.file("extdata", "lf_agb_1ha.tif", package="rLVIS")
#'
#'sf_agb<-raster(sf_agb1ha_path)
#'lf_agb<-raster(lf_agb1ha_path)
#'
#'# Ploting AGB maps
#'s <- stack(sf_agb,lf_agb)
#'agb.maps<-levelplot(s,
#'          layout=c(1, 2),
#'          margin=FALSE,
#'          colorkey=list(
#'            space='right',
#'            labels=list(at=seq(0, 500, 50), font=4),
#'            axis.line=list(col='black'),
#'            width=1
#'          ),
#'          par.settings=list(
#'            strip.border=list(col='transparent'),
#'            strip.background=list(col='transparent'),
#'            axis.line=list(col='transparent')
#'          ),
#'          scales=list(draw=TRUE),
#'          col.regions=viridis,
#'          at=seq(0, 500, len=101),
#'          names.attr=c("SF_AGB","LF_AGB"))
#'
#'# Ploting Stats
#'colours<-viridis(10)
#'base_size=15
#'legend.position= c(0.85, 0.3)
#'stat.size=5
#'stats.position=c(100,400,50)
#'base_size=base_size
#'xlim=c(0,500)
#'legend.size=c(8,15,10,2)
#'ylim=c(0,500)
#'ylab="LF-derived AGB (Mg/ha)"
#'xlab="SF-derived AGB (Mg/ha)"
#'fit.line.col=c("black","gray")
#'title="SF vs LF lidar"
#'x_axis=TRUE
#'y_axis=TRUE
#'
#'\dontrun{
#' x11()
#'}
#'agb.comp<-plotStats(y=getValues(lf_agb),
#'                  x=getValues(sf_agb),
#'                  colours=colours,
#'                  legend.position= legend.position,
#'                  stat.size=stat.size,
#'                  stats.position=stats.position,
#'                  base_size=base_size,
#'                  xlim=xlim,
#'                  legend.size=legend.size,
#'                  ylim=ylim,
#'                  ylab=ylab,
#'                  xlab=xlab,
#'                  fit.line.col=fit.line.col,
#'                  title=title,
#'                  x_axis=x_axis,
#'                  y_axis=y_axis)
#'
#'# Combining plots
#'grid.arrange(agb.maps,agb.comp$plotg, nrow = 1)
#'
#'@import ggplot2
#'@importFrom stats na.omit lm pf cor
#'@export
plotStats<-function(x,y,n,colours=c("blue","green","yellow","red"),xlim=NULL,ylim=NULL,colfill="gray",xlab=NULL,ylab=NULL,title="GEDI vs ALS",legend.position = c(0.8, 0.2),stats.position=c(50,400,5),
                    base_size = 16,stat.size=3.5, legend.size=c(8,8,3,1),fit.line.col=c("black","red"),x_axis=TRUE, y_axis=TRUE){

  density = NA

  xy<-na.omit(cbind(x,y))
  x<-xy[,1]
  y<-xy[,2]

  get_density <- function(x, y,h=200, n = 500) {
    dens <- MASS::kde2d(x = x, y = y,h=h, n = n)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }

  theme_set(theme_bw(base_size = base_size))
  dat<-na.omit(data.frame(cbind(y,x)))
  dat$labels<-rep(title, nrow(dat))
  dat$density <- get_density(dat$x, dat$y)
  #mid<-mean(dat$density)
  g<-ggplot(dat) + geom_point(aes(x, y, color = density))
  g<-g+ scale_color_gradientn(colours=colours) + theme(legend.position = legend.position)
  g<-g+ facet_wrap(~ labels) + theme(legend.text=element_text(size=legend.size[1]),
                                     legend.title=element_text(size=legend.size[2])) +
    guides(colour=guide_colourbar(barwidth=legend.size[4], barheight=legend.size[3]))+
    theme(strip.background = element_rect(fill=colfill),strip.text = element_text(size=15)) +
    labs(x =xlab,y=ylab) + geom_smooth(data=dat, formula=y~x,method="lm", color=fit.line.col[1],aes(x, y)) +
    geom_abline(slope=1, intercept=0, color=fit.line.col[2], size = 1)
  #+ theme_bw(base_size = base_size)

  if (!is.null(xlim)) { g<-g + xlim(xlim) }
  if (!is.null(ylim)) { g<-g + ylim(ylim) }

  rmse <- sqrt( sum(( y - x )^2)/length(x) ) # Root mean square error
  rmseR <- round(100 * sqrt( sum(( y - x )^2)/length(x) ) / mean( x ),2)
  bias <- mean( y - x ) # bias
  biasR <-round( 100 * mean( y - x ) / mean(x),2)

  m<-lm(y~x)
  ms<-summary(m)

  l <- list(r2 = format(ms$r.squared, digits = 2),
            pval = format(pf(ms$fstatistic[1],ms$fstatistic[2],ms$fstatistic[3],lower.tail=FALSE), digits=3)
  )

  eq <- substitute(italic(r)^2 == r2,l)

  eqstr <- as.character(as.expression(eq))
  g<-g + annotate("text", x=stats.position[1], y=stats.position[2], label= eqstr,parse=TRUE,size = stat.size) + #fontface = "bold"
    annotate("text", x = stats.position[1], y=stats.position[2]+stats.position[3], label = paste("RMSE(%)=",rmseR),size = stat.size)+
    annotate("text", x = stats.position[1], y=stats.position[2]+stats.position[3]*2, label = paste("Bias(%)=",biasR),size =stat.size)

  if( x_axis==FALSE & y_axis==F ) {
  g<-g+theme(
        axis.text.x=element_blank(),
        axis.text.y=element_blank()
        )
  }

  if( x_axis==FALSE & y_axis==TRUE ) {
    g<-g+theme(
      axis.text.x=element_blank()
    )
  }


  if( x_axis==TRUE & y_axis==FALSE ) {
    g<-g+theme(
      axis.text.y=element_blank()
    )
  }

  StatInfo<-data.frame( Stat=c("rmse","rmseR","bias","biasR","r","r2"),
                        Values=c(rmse,rmseR,bias,biasR,cor(x,y),ms$r.squared))
  print(g)
  return(list(stats=StatInfo,plotg=g))
}


