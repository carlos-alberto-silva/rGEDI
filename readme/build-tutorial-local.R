# All displayed outputs use measured package data and public rGEDI functions.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(rGEDI)
  library(sf); library(terra); library(data.table); library(ggplot2); library(patchwork)
})
dir.create("readme/results", showWarnings = FALSE, recursive = TRUE)
work <- tempfile("rgedi-readme-");dir.create(work)
emit <- function(x, name) writeLines(capture.output(print(x)), file.path("readme/results", name))
open_zip <- function(pattern, reader) reader(unzip(system.file("extdata", pattern, package = "rGEDI"), exdir = work))
a <- open_zip("GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip", readLevel2A)
b <- open_zip("GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip", readLevel2B)
l1 <- open_zip("GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip", readLevel1B)
level1b_geo <- getLevel1BGeo(l1)
level2a_metrics <- getLevel2AM(a); level2b_vpm <- getLevel2BVPM(b)
emit(head(level1b_geo[, .(shot_number, longitude_bin0, latitude_bin0)], 3), "level1b-geolocation.txt")
emit(head(level2a_metrics[, .(beam,shot_number,rh50,rh98,rh100)],3),"level2a-metrics.txt")
emit(head(level2b_vpm[, .(beam,shot_number,pai,cover,fhd_normal)],3),"level2b-metrics.txt")
polygons <- st_read(system.file("extdata","stands_cerrado.shp",package="rGEDI"),quiet=TRUE)
polygons$unit_id <- seq_len(nrow(polygons))
level2b_good <- level2b_vpm[l2b_quality_flag == 1 & sensitivity >= .9]
footprints <- clip(level2b_good, polygons, split_by = "unit_id")
structure_metrics <- function(pai, cover) list(n=length(pai),mean_pai=mean(pai),mean_cover=mean(cover))
polygon_summary <- polyStatsLevel2BVPM(footprints, ~structure_metrics(pai,cover), id="poly_id")
emit(polygon_summary,"polygon-summary.txt")

# The real ALS examples originally distributed in the rGEDI datasets release.
waves <- list(); metrics <- list(); clouds <- list()
for (site in c("Amazon","Savanna")) {
  cloud <- lidR::readLAS(system.file("extdata",paste0(site,".laz"),package="rGEDI"))
  xyz <- as.data.table(cloud@data)
  centre <- c(mean(range(xyz$X)),mean(range(xyz$Y)))
  clouds[[site]] <- xyz[,.(X=X-centre[1],Y=Y-centre[2],Z=Z-min(Z))][X^2+Y^2<=12.5^2]
  for (noise in c(0,.001)) {
    obj <- gediWFSimulator(cloud,output=file.path(work,paste0(site,"_",noise,".h5")),coords=centre,noise=noise,seed=11)
    m <- gediWFMetrics(obj,rhRes=1,ground_method="lowest_peak")
    wf <- getLevel1BWF(obj,m$shot_number[1])@dt
    wf[,height:=elevation-m$gHeight[1]]
    tag <- paste(site,noise)
    waves[[tag]] <- wf
    metrics[[tag]] <- m[,.(site=site,noise=noise,cover,rh50,rh90,rh98,rh100,waveEnergy)]
    close(obj)
  }
}
ans <- rbindlist(metrics)
emit(ans,"simulation-metrics.txt");fwrite(ans,"readme/results/simulation-metrics.csv")
draw_als <- function() {
  old <- par(mfrow=c(2,2),mar=c(4,4,3,1),oma=c(0,0,0,1));on.exit(par(old))
  for (site in c("Amazon","Savanna")) {
    d <- clouds[[site]][seq(1,.N,by=3)]
    plot3D::scatter3D(d$X,d$Y,d$Z,colvar=d$Z,col=viridisLite::viridis(100),
      pch=16,cex=.45,clab="Height (m)",ticktype="detailed",colkey=list(length=.5,width=.5,cex.axis=.7,cex.clab=.75),theta=45,phi=25,bty="u",col.panel="#F6F7F5",
      col.grid="#DDE2DE",xlab="x (m)",ylab="y (m)",zlab="",
      xlim=c(-13,13),ylim=c(-13,13),zlim=c(0,45),clim=c(0,45),main=paste(site,"ALS: 25 m footprint"))
    clean <- waves[[paste(site,0)]];noisy <- waves[[paste(site,.001)]]
    plot(clean$rxwaveform/max(clean$rxwaveform)*100,clean$height,type="l",lwd=2,
      col="#267F82",xlim=c(0,105),ylim=c(-5,45),xlab="Peak return (%)",ylab="Height above ground (m)",
      main=paste(site,"simulated waveform"))
    lines(noisy$rxwaveform/max(noisy$rxwaveform)*100,noisy$height,col="#CA7137",lwd=1)
    abline(h=0,col="#65756B",lty=2)
    rh98 <- metrics[[paste(site,0)]]$rh98[1]
    abline(h=rh98,col="#A42A73",lty=2)
    text(103,rh98+1.5,sprintf("RH98 %.1f m",rh98),adj=1,col="#A42A73",cex=.85)
    legend("topright",c("Clean","Noise = 0.001"),col=c("#267F82","#CA7137"),lty=1,lwd=2,bty="n",cex=.85)
  }
}
png("readme/fig-package-als.png",width=2100,height=1700,res=200);draw_als();dev.off()
pdf("readme/fig-package-als.pdf",width=10.5,height=8.5);draw_als();dev.off()
close(a);close(b);close(l1)
unlink(work,recursive=TRUE)
