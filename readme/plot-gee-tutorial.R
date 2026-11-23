# Draw the package's own cloud tutorial outputs; no manuscript inputs.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload",quietly=TRUE)) pkgload::load_all(".",quiet=TRUE) else library(rGEDI)
  library(sf);library(terra);library(data.table);library(ggplot2);library(patchwork)
})
root <- Sys.getenv("RGEDI_EXAMPLE_CACHE",file.path(tempdir(),"rGEDI-package-examples"))
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face="bold")))
save_plot <- function(p,name,w,h) {
  ggsave(paste0("readme/",name,".png"),p,width=w,height=h,dpi=300,bg="white")
  ggsave(paste0("readme/",name,".pdf"),p,width=w,height=h,bg="white")
}
# Global native Level 4B, displayed on Equal Earth without spatial analysis.
r <- readLevel4B(file.path(root,"global_l4b_complete.tif"));r[r<0] <- NA
rr <- projectLevel4B(r,"EPSG:8857",method="near")
d <- as.data.table(as.data.frame(rr,xy=TRUE,na.rm=TRUE));setnames(d,3,"agbd")
world <- st_as_sf(rnaturalearthdata::countries110) |> st_make_valid()
world <- st_transform(world, 8857)
reference <- st_transform(st_sfc(st_linestring(matrix(c(0,0,44.915764,0),ncol=2,byrow=TRUE)),crs=4326),8857)
scale_length <- diff(st_coordinates(reference)[,1])
global <- ggplot()+geom_sf(data=world,fill="#E6E8E3",colour="#C6CCC6",linewidth=.15)+
  geom_tile(data=d,aes(x,y,fill=agbd),width=res(rr)[1],height=res(rr)[2])+
  scale_fill_viridis_c(option="magma",name="AGBD (Mg/ha)",limits=c(0,400),oob=scales::squish,
    breaks=c(0,100,200,300,400),labels=c("0","100","200","300",">=400"))+
  coord_sf(crs=8857,xlim=c(-17500000,17500000),ylim=c(-8100000,8100000),datum=NA,expand=FALSE)+
  annotate("segment",x=-15500000,xend=-15500000+scale_length,y=-6500000,yend=-6500000,linewidth=1)+
  annotate("text",x=-15500000+scale_length/2,y=-7050000,label="5,000 km at equator",size=3.3)+
  annotate("segment",x=0,xend=0,y=6400000,yend=7400000,arrow=grid::arrow(type="closed",length=grid::unit(2,"mm")))+
  annotate("text",x=0,y=7750000,label="N",size=3.5)+
  labs(title="GEDI Level 4B biomass through Google Earth Engine",
    subtitle="Native MU estimates; 20 km nearest-neighbour display preview",
    caption="Equal Earth projection. Grey: no displayed estimate. This preview is not a new biomass product.")+
  theme_void(base_size=12)+theme(legend.position="bottom",plot.title=element_text(face="bold"),
    plot.subtitle=element_text(size=11),plot.caption=element_text(hjust=0,size=9),legend.key.width=grid::unit(1.5,"cm"))
save_plot(global,"fig-gee-global",11.5,6.5)

pred <- fread("readme/results/shenandoah_predictions.csv")
pred[,stage:=factor(stage,levels=c("Training","Validation"))]
summary <- pred[,.(n=.N,R2=1-sum((predicted-agbd)^2)/sum((agbd-mean(agbd))^2),
  RMSE=sqrt(mean((predicted-agbd)^2)),bias=mean(predicted-agbd)),by=stage]
summary[,label:=sprintf("n = %d\nR\u00b2 = %.2f\nRMSE = %.1f Mg/ha\nBias = %.1f Mg/ha",n,R2,RMSE,bias)]
maximum <- ceiling(max(pred$agbd,pred$predicted)/100)*100
p <- ggplot(pred,aes(agbd,predicted))+geom_abline(slope=1,intercept=0,linetype=2,colour="#68777C")+
  geom_point(colour="#267F82",size=1.15,alpha=.55)+geom_smooth(method="lm",formula=y~x,se=FALSE,colour="#A42A73",linewidth=.7)+
  geom_text(data=summary,aes(x=.04*maximum,y=.96*maximum,label=label),inherit.aes=FALSE,hjust=0,vjust=1,size=3.5)+
  facet_wrap(~stage)+coord_equal(xlim=c(0,maximum),ylim=c(0,maximum),expand=FALSE)+
  labs(x="Reference GEDI AGBD (Mg/ha)",y="Predicted AGBD (Mg/ha)",title="Biomass model training and spatial validation",
    caption="Dashed: 1:1. Magenta: fitted relationship. Validation uses held-out 2 km blocks across all sampled years.")+
  theme(strip.text=element_text(face="bold"),plot.caption=element_text(hjust=0,size=9))
save_plot(p,"fig-gee-model-validation",10.5,5.5)
output <- summary[,.(stage,n,R2=round(R2,3),RMSE=round(RMSE,2),bias=round(bias,2))]
writeLines(capture.output(print(output)),"readme/results/shenandoah-evaluation.txt")
importance <- fread("readme/results/shenandoah_importance.csv")
rfe <- fread("readme/results/shenandoah_rfe.csv")
top <- importance[order(-importance)][1:15]
top[,parameter:=factor(parameter,levels=rev(parameter))]
chosen <- sum(importance$selected)
a <- ggplot(top,aes(importance,parameter,fill=selected))+geom_col(width=.7)+
  scale_fill_manual(values=c("FALSE"="#C9D0D0","TRUE"="#705196"),labels=c("Other","Selected"),name=NULL)+
  labs(title="A  Predictor importance",x="Relative importance",y=NULL)+theme(legend.position="bottom")
b <- ggplot(rfe,aes(nvariables,oob_rmse))+geom_line(colour="#267F82",linewidth=.8)+
  geom_vline(xintercept=chosen,colour="#A42A73",linetype=2)+
  labs(title="B  Recursive feature elimination",x="Predictors retained",y="Training OOB RMSE (Mg/ha)")
save_plot(a+b+plot_annotation(caption="Selection uses training blocks only. The final GEE forest uses the selected predictors."),"fig-gee-variable-selection",10.5,4.5)

area <- st_read(system.file("extdata","shenandoah.shp",package="rGEDI"),quiet=TRUE) |> st_transform(32617)
footprints <- fread("readme/results/shenandoah_footprints.csv")
points <- st_as_sf(footprints,coords=c("lon_lowestmode","lat_lowestmode"),crs=4326) |> st_transform(32617)
biomass <- rast(file.path(root,"shenandoah_agbd_2021.tif"))
df <- as.data.frame(biomass,xy=TRUE,na.rm=TRUE);names(df)[3]<-"agbd"
dem <- rast(file.path(root,"shenandoah_dem.tif"))
terrain <- terra::terrain(dem,v=c("slope","aspect"),unit="radians")
hill <- terra::shade(terrain[[1]],terrain[[2]],angle=40,direction=315)
background <- as.data.frame(hill,xy=TRUE,na.rm=TRUE);names(background)[3]<-"shade"
bb <- st_bbox(area);dx<-bb[3]-bb[1];dy<-bb[4]-bb[2]
decorate <- function(p) p+geom_sf(data=area,fill=NA,colour="#263A39",linewidth=.45)+
  coord_sf(xlim=c(bb[1]-.04*dx,bb[3]+.04*dx),ylim=c(bb[2]-.04*dy,bb[4]+.04*dy),datum=NA,expand=FALSE)+
  annotate("segment",x=bb[3]-.35*dx,xend=bb[3]-.35*dx+10000,y=bb[2]+.05*dy,yend=bb[2]+.05*dy,linewidth=.8)+
  annotate("label",x=bb[3]-.35*dx+5000,y=bb[2]+.08*dy,label="10 km",size=3,linewidth=0)+
  annotate("segment",x=bb[1]+.1*dx,xend=bb[1]+.1*dx,y=bb[4]-.15*dy,yend=bb[4]-.05*dy,
           arrow=grid::arrow(type="closed",length=grid::unit(2,"mm")))+
  annotate("text",x=bb[1]+.1*dx,y=bb[4]-.02*dy,label="N",size=3.5)+theme_void(base_size=12)+
  theme(legend.position="bottom",legend.key.width=grid::unit(1.1,"cm"),plot.title=element_text(face="bold"),plot.subtitle=element_text(size=10))
a <- ggplot()+geom_tile(data=background,aes(x,y,fill=shade),width=res(hill)[1],height=res(hill)[2])+scale_fill_gradient(low="#6A746E",high="#F0F2ED",guide="none")+
  geom_sf(data=points,aes(colour=partition),size=.6,alpha=.7)+
  scale_colour_manual(values=c("Training"="#228A87","Validation"="#CC763A"),name=NULL)+
  labs(title="A  Shenandoah National Park",subtitle="NASADEM hillshade and GEDI model samples")
b <- ggplot(df)+geom_tile(aes(x,y,fill=agbd),width=res(biomass)[1],height=res(biomass)[2])+
  scale_fill_viridis_c(option="magma",name="AGBD (Mg/ha)",limits=c(0,500),oob=scales::squish,
    breaks=c(0,100,200,300,400,500),labels=c("0","100","200","300","400",">=500"))+
  labs(title="B  Wall-to-wall prediction",subtitle="2021 ancillary predictors; 100 m display grid")
maps <- decorate(a)+decorate(b)+plot_annotation(caption="NPS park boundary simplified by 20 m. Final map refitted on all 2,348 GEDI samples. Predictions are calibrated to GEDI AGBD, not field biomass.")
save_plot(maps,"fig-gee-shenandoah-map",10.5,8)
raw <- readRDS(file.path(root,"shenandoah_raw_2019.rds"))
writeLines(capture.output(print(head(raw[,.(beam,shot_number,agbd,agbd_se,sensitivity)],3))),"readme/results/level4a-footprints.txt")
