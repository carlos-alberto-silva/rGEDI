# Offline demonstrations for native GEDI Level 3 and Level 4B products.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(rGEDI)
  library(terra); library(sf); library(data.table)
})
dir.create("readme/results",recursive=TRUE,showWarnings=FALSE)
emit <- function(x,n) writeLines(capture.output(print(x)),file.path("readme/results",n))
level3 <- readLevel3(system.file("extdata","GEDI03_Sequoia.tif",package="rGEDI"))
level4b <- readLevel4B(system.file("extdata","GEDI04_B_Sequoia.tif",package="rGEDI"))
grid_aoi <- st_read(system.file("extdata","sequoia.shp",package="rGEDI"),quiet=TRUE)
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
# Figure display is projected; all extraction above uses the native grids.
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
outline <- st_transform(grid_aoi,32611)
draw <- function(full,clipped,product,unit,palette) {
  projector <- if(product=="3") projectLevel3 else projectLevel4B
  full <- projector(full,"EPSG:32611",method="near")
  clipped <- projector(clipped,"EPSG:32611",method="near")
  e <- as.vector(ext(full));dx<-e[2]-e[1];dy<-e[4]-e[3]
  limits <- range(values(full),na.rm=TRUE)
  panel <- function(r,title) {
    d <- as.data.frame(r,xy=TRUE,na.rm=TRUE);names(d)[3]<-"value"
    ggplot(d,aes(x,y,fill=value))+geom_tile(width=res(r)[1],height=res(r)[2])+
      geom_sf(data=outline,inherit.aes=FALSE,fill=NA,colour="white",linewidth=1.1)+
      geom_sf(data=outline,inherit.aes=FALSE,fill=NA,colour="#50575A",linewidth=.45)+
      scale_fill_viridis_c(option=palette,limits=limits,name=unit,
        guide=guide_colourbar(barwidth=grid::unit(7,"cm"),barheight=grid::unit(.35,"cm")))+
      coord_sf(crs=32611,xlim=e[1:2],ylim=c(e[3]-.12*dy,e[4]),datum=NA,expand=FALSE)+
      annotate("segment",x=e[1]+.08*dx,xend=e[1]+.08*dx+20000,y=e[3]-.07*dy,yend=e[3]-.07*dy,linewidth=.8)+
      annotate("text",x=e[1]+.08*dx+10000,y=e[3]-.035*dy,label="20 km",size=3.5)+
      annotate("rect",xmin=e[2]-.13*dx,xmax=e[2]-.03*dx,ymin=e[4]-.20*dy,ymax=e[4]-.015*dy,fill="white",colour="#E1E4E0",linewidth=.2)+
      annotate("segment",x=e[2]-.08*dx,xend=e[2]-.08*dx,y=e[4]-.18*dy,yend=e[4]-.07*dy,
        arrow=grid::arrow(type="closed",length=grid::unit(2,"mm")))+
      annotate("text",x=e[2]-.08*dx,y=e[4]-.035*dy,label="N",size=4)+
      labs(title=title,subtitle="Sequoia National Forest")+
      theme_void(base_size=13)+theme(plot.title=element_text(face="bold"),legend.position="bottom",
        panel.background=element_rect(fill="#F2F3F1",colour=NA),plot.background=element_rect(fill="white",colour=NA))
  }
  panel(full,paste0("A  Native Level ",product))+panel(clipped,"B  Clipped to the forest")+
    plot_layout(guides="collect")+plot_annotation(caption="Grey outline: USFS administrative boundary. Light grey: no displayed estimate. Native nominal 1 km cells; nearest-neighbour projection for display.") &
    theme(legend.position="bottom",plot.caption=element_text(size=9,hjust=0))
}
for (product in c("3","4B")) {
  p <- if(product=="3") draw(level3[["rh100_mean"]],level3_aoi[["rh100_mean"]],product,"RH100 (m)","viridis") else
    draw(level4b[["MU"]],level4b_aoi[["MU"]],product,"AGBD (Mg/ha)","magma")
  stem <- paste0("readme/fig-level",tolower(product),"-products")
  ggsave(paste0(stem,".png"),p,width=11,height=7,dpi=250,bg="white")
  ggsave(paste0(stem,".pdf"),p,width=11,height=7,bg="white")
}
