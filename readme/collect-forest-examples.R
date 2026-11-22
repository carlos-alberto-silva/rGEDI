# Access and process GEDI through rGEDI. Keep new caches on D:.
library(rGEDI)
library(data.table)
library(sf)
library(terra)
root <- Sys.getenv('RGEDI_FIGURE_DATA',unset='D:/rGEDI_manuscript/Manuscript_v07/forest_examples_v15')
dir.create(root,recursive=TRUE,showWarnings=FALSE)
ee <- ee_initialize(project=Sys.getenv('EE_PROJECT'),authenticate=FALSE,quiet=TRUE)

# Join the actual L2A metric to the L2B profile shots.
cache <- file.path(root,'appalachian_rh.csv')
if (!file.exists(cache)) {
  refs <- gedifinder('GEDI02_A',orbit=2629,version='003',return='table')
  refs <- refs[as.integer(granule_part)==3]
  stopifnot(nrow(refs)==1)
  obj <- readLevel2A(refs$url)
  metrics <- getLevel2AM(obj,beams='BEAM0101',rh=c(0,25,50,75,90,98,100))
  close(obj)
  metrics[, shot_number:=as.character(shot_number)]
  fwrite(metrics,cache)
  cat('Saved matched Appalachian L2A RH metrics.\n')
}

sites <- list(
 list(name='Sierra Nevada',label='Montane conifer forest',box=c(-120.48,-120.28,39.30,39.50),forest=c(42L,43L),project='CA_SierraNevada_6_2022',als_year=2022L),
 list(name='Palm Coast',label='Southeastern evergreen',box=c(-81.35,-81.18,29.44,29.61),forest=42L,project='FL_2017PalmCoast_1_C22',als_year=2017L))
traits <- ee$Image('USGS/NLCD_RELEASES/2019_REL/NLCD/2019')$select('landcover')$addBands(ee$Terrain$slope(ee$Image('USGS/SRTMGL1_003'))$rename('slope'))
for (site in sites) {
  folder <- file.path(root,gsub(' ','_',tolower(site$name)))
  dir.create(folder,showWarnings=FALSE)
  if (file.exists(file.path(folder,'shots.csv'))) next
  b <- site$box
  ref <- gedifinder('GEDI02_A',ul_lat=b[4],ul_lon=b[1],lr_lat=b[3],lr_lon=b[2],daterange=c('2019-05-01','2019-10-01'),source='gee')
  candidates <- getLevel2AM(readLevel2A(ref),rh=c(0,25,50,75,90,98,100),quality=TRUE,max_features=30000L)
  candidates <- candidates[quality_flag==1 & degrade_flag==0 & sensitivity>=.95 & solar_elevation<0 & beam %in% c(5,6,8,11) & between(rh98,12,58)]
  stopifnot(nrow(candidates)>0)
  if (nrow(candidates)>160) candidates <- candidates[round(seq(1,.N,length.out=160))]
  candidates[, shot_number:=as.character(shot_number)]
  context <- extractEE(traits,candidates[,.(shot_number,lon_lowestmode,lat_lowestmode)],scale=30)
  candidates <- merge(candidates,context[,.(shot_number,landcover,slope)],by='shot_number')
  candidates <- candidates[landcover %in% site$forest & slope<=20]
  stopifnot(nrow(candidates)>0)
  candidates[, priority:=abs(rh98-median(rh98))]
  setorder(candidates,priority)
  fwrite(candidates,file.path(folder,'candidates.csv'))
  chosen <- candidates[1]
  shot <- chosen$shot_number
  orbit <- as.integer(substr(shot,1,nchar(shot)-13))
  part <- as.integer(substr(shot,nchar(shot)-10,nchar(shot)-8))
  beam <- paste0('BEAM',paste(rev(as.integer(intToBits(as.integer(chosen$beam)))[1:4]),collapse=''))
  refs <- gedifinder('GEDI02_A',orbit=orbit,version='003',return='table')[as.integer(granule_part)==part]
  stopifnot(nrow(refs)==1)
  obj <- readLevel2A(refs$url)
  metrics <- getLevel2AM(obj,beams=beam,rh=c(0,25,50,75,90,98,100))
  close(obj)
  metrics[, shot_number:=as.character(shot_number)]
  actual <- metrics[shot_number==shot]
  stopifnot(nrow(actual)==1)
  refs1 <- gedifinder('GEDI01_B',orbit=orbit,version='003',return='table')[as.integer(granule_part)==part]
  stopifnot(nrow(refs1)==1)
  obj <- readLevel1B(refs1$url)
  waveform <- getLevel1BWF(obj,shot)
  g <- obj@h5[[beam]]
  index <- which(as.character(g[['shot_number']][])==shot)
  noise <- g[['noise_mean_corrected']][index]
  close(obj)
  signal <- pmax(as.numeric(waveform@dt$rxwaveform)-as.numeric(noise),0)
  stopifnot(length(signal)>100,all(is.finite(signal)),max(signal)>0)
  bins <- data.table(shot=shot,landscape=site$name,height_m=waveform@dt$elevation-actual$elev_lowestmode,amplitude_pct=100*signal/max(signal))
  date <- as.Date(strptime(substr(refs1$granule_id,10,16),'%Y%j'))
  record <- data.table(shot=shot,beam=beam,granule=sub('[.]h5$','',refs$granule_id),date=as.character(date),
    lon=actual$lon_lowestmode,lat=actual$lat_lowestmode,ground=actual$elev_lowestmode,
    rh0=actual$rh0,rh25=actual$rh25,rh50=actual$rh50,rh75=actual$rh75,rh90=actual$rh90,rh98=actual$rh98,rh100=actual$rh100,
    landscape=site$name,nlcd=chosen$landcover,slope=chosen$slope,l1b_granule=sub('[.]h5$','',refs1$granule_id),
    vegetation=site$label,als_year=site$als_year,project=site$project)
  fwrite(record,file.path(folder,'shots.csv'));fwrite(bins,file.path(folder,'bins.csv'))
  cat('Saved real GEDI waveform and L2A metrics for',site$name,'shot',shot,'\n')
}
