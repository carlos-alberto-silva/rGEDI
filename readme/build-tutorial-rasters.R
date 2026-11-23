# Offline demonstrations for native GEDI Level 3 and Level 4B products.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(rGEDI)
  library(terra); library(sf); library(data.table)
})
dir.create("readme/results",recursive=TRUE,showWarnings=FALSE)
emit <- function(x,n) writeLines(capture.output(print(x)),file.path("readme/results",n))
level3 <- readLevel3(system.file("extdata","GEDI03_Brasilia.tif",package="rGEDI"))
level4b <- readLevel4B(system.file("extdata","GEDI04_B_Brasilia.tif",package="rGEDI"))
grid_aoi <- st_read(system.file("extdata","brasilia.shp",package="rGEDI"),quiet=TRUE)
native_aoi <- st_transform(grid_aoi,crs(level3))
bounds <- as.vector(ext(vect(native_aoi)))
level3_bbox <- clipLevel3(level3,bounds)
level4b_bbox <- clipLevel4B(level4b,bounds)
level3_aoi <- clipLevel3(level3,vect(native_aoi))
level4b_aoi <- clipLevel4B(level4b,vect(native_aoi))
level3_values <- getLevel3(level3,geometry=vect(native_aoi))
level4b_values <- getLevel4B(level4b,geometry=vect(native_aoi))
emit(head(level3_values[is.finite(rh100_mean)],3),"level3-values.txt")
emit(head(level4b_values[is.finite(MU)],3),"level4b-values.txt")
level4b_stats <- polyStats(level4b,vect(native_aoi),fun=mean,na.rm=TRUE)
emit(level4b_stats,"level4b-statistics.txt")
emit(data.table(product=c("Level 3","Level 4B"),rows=c(nrow(level3),nrow(level4b)),
  columns=c(ncol(level3),ncol(level4b)),layers=c(nlyr(level3),nlyr(level4b)),
  resolution_m=c(res(level3)[1],res(level4b)[1])),"native-raster-info.txt")
emit(data.table(output=c("L3 bbox","L3 polygon","L4B bbox","L4B polygon"),
  cells=vapply(list(level3_bbox,level3_aoi,level4b_bbox,level4b_aoi),ncell,numeric(1)),
  valid=vapply(list(level3_bbox,level3_aoi,level4b_bbox,level4b_aoi),
    function(x) sum(is.finite(values(x[[1]]))),numeric(1))),"native-clips.txt")
outline <- vect(st_transform(grid_aoi,32723))
draw <- function(full,clipped,reader,name,unit,palette) {
  full <- if (reader=="3") projectLevel3(full,"EPSG:32723",method="near") else projectLevel4B(full,"EPSG:32723",method="near")
  clipped <- if (reader=="3") projectLevel3(clipped,"EPSG:32723",method="near") else projectLevel4B(clipped,"EPSG:32723",method="near")
  zlim <- range(values(full),na.rm=TRUE)
  old <- par(mfrow=c(1,2),mar=c(3,3,4,5),cex=.95);on.exit(par(old))
  for (i in 1:2) {
    obj <- if(i==1) full else clipped
    title <- if(i==1) paste0("A  Native Level ",reader," metric") else "B  Clipped to Brasilia boundary"
    plotfun <- if(reader=="3") plotLevel3 else plotLevel4B
    plotfun(obj,main=title,col=viridisLite::viridis(100,option=palette),range=zlim,
      xlim=c(xmin(full),xmax(full)),ylim=c(ymin(full),ymax(full)),axes=FALSE,
      plg=list(title=unit))
    plot(outline,add=TRUE,border="#D97936",lwd=1.5)
    viewport <- par("usr"); graphics::clip(viewport[1],viewport[2],viewport[3],viewport[4])
    u<-as.vector(ext(full));x<-u[1]+.07*diff(u[1:2]);y<-u[3]+.07*diff(u[3:4])
    segments(x,y,x+20000,y,lwd=3)
    text(x+10000,y+.035*diff(u[3:4]),"20 km",cex=.85)
    nx<-u[2]-.08*diff(u[1:2]);ny<-u[4]-.06*diff(u[3:4])
    arrows(nx,ny-.1*diff(u[3:4]),nx,ny,length=.09,lwd=1.5)
    text(nx,ny+.03*diff(u[3:4]),"N",cex=1)
    mtext("Orange outline: IBGE Brasilia boundary",side=1,line=1,cex=.8)
  }
}
png("readme/fig-level3-products.png",width=2100,height=800,res=190)
draw(level3[["rh100_mean"]],level3_aoi[["rh100_mean"]],"3","rh100_mean","RH100 (m)","viridis");dev.off()
png("readme/fig-level4b-products.png",width=2100,height=800,res=190)
draw(level4b[["MU"]],level4b_aoi[["MU"]],"4B","MU","AGBD (Mg/ha)","magma");dev.off()
