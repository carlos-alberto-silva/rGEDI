library(rGEDI)
library(data.table)
library(sf)
root <- 'D:/rGEDI_manuscript/Manuscript_v07/forest_examples_v15'
park <- st_read('readme/data/grsm_boundary.geojson',quiet=TRUE)
raw <- fread(file.path(root,'appalachian_rh.csv'),colClasses=c(shot_number='character'))
raw <- raw[is.finite(lon_lowestmode)&is.finite(lat_lowestmode)]
single <- if (file.exists('readme/data/appalachian_single_granule.csv')) fread('readme/data/appalachian_single_granule.csv',colClasses=c(shot_number='character')) else rGEDI::clip(raw,park)
if (length(unique(single$beam)) != 8L) {
  cat("Reading all eight tracks through rGEDI getTrack()...\n")
  refs <- gedifinder("GEDI02_A",orbit=2629,version="003",return="table")[as.integer(granule_part)==3]
  obj <- readLevel2A(refs$url)
  tracks <- getTrack(obj)
  close(obj)
  raw <- as.data.table(tracks)
  setnames(raw,c("longitude","latitude"),c("lon_lowestmode","lat_lowestmode"))
  raw[,shot_number:=as.character(shot_number)]
  single <- if (file.exists('readme/data/appalachian_single_granule.csv')) fread('readme/data/appalachian_single_granule.csv',colClasses=c(shot_number='character')) else rGEDI::clip(raw,park)
}
stopifnot(length(unique(single$beam))==8)
fwrite(single,'readme/data/appalachian_single_granule.csv')
cat('Single granule contains all eight beams in the park:',nrow(single),'footprints\n')
overview_file <- file.path(root,'appalachian_multigranule_leafon2019.csv')
if (!file.exists(overview_file)) {
  ee_initialize(project=Sys.getenv('EE_PROJECT'),authenticate=FALSE,quiet=TRUE)
  b<-st_bbox(park)
  ref <- gedifinder('GEDI02_A',ul_lat=unname(b['ymax']),ul_lon=unname(b['xmin']),lr_lat=unname(b['ymin']),lr_lon=unname(b['xmax']),
                    daterange=c('2019-05-01','2019-10-01'),source='gee')
  overview <- getLevel2AM(readLevel2A(ref),rh=98,quality=TRUE,max_features=200000L,chunk_size=10000L)
  overview <- overview[quality_flag==1 & degrade_flag==0]
  overview[,shot_number:=as.character(shot_number)]
  kept <- rGEDI::clip(overview,park)
  overview[,inside:=shot_number %in% kept$shot_number]
  overview[,granule_key:=paste(substr(shot_number,1,nchar(shot_number)-13),substr(shot_number,nchar(shot_number)-10,nchar(shot_number)-8),sep='_')]
  fwrite(overview,overview_file)
  cat('Overview:',nrow(overview),'screened footprints from',uniqueN(overview$granule_key),'granules\n')
}
file.copy(overview_file,'readme/data/appalachian_multigranule.csv',overwrite=TRUE)
wavefile <- file.path(root,'appalachian_inset_waveform.csv')
if (!file.exists(wavefile)) {
  best <- fread('readme/data/appalachian_profile_segment.csv',colClasses=c(shot_number='character'))
  chosen <- best[ceiling(.N/2)]
  refs <- gedifinder('GEDI01_B',orbit=2629,version='003',return='table')[as.integer(granule_part)==3]
  obj <- readLevel1B(refs$url)
  wf <- getLevel1BWF(obj,chosen$shot_number)
  g<-obj@h5[[chosen$beam]]
  i<-which(as.character(g[['shot_number']][])==chosen$shot_number)
  noise<-as.numeric(g[['noise_mean_corrected']][i])
  signal<-pmax(as.numeric(wf@dt$rxwaveform)-noise,0)
  close(obj)
  stopifnot(length(signal)>100,max(signal)>0)
  out<-data.table(shot_number=chosen$shot_number,height=wf@dt$elevation-chosen$l2a_ground,amplitude=100*signal/max(signal),rh98=chosen$rh98,distance=chosen$distance)
  fwrite(out,wavefile)
  cat('Inset waveform verified for shot',chosen$shot_number,'with',nrow(out),'return bins\n')
}
file.copy(wavefile,'readme/data/appalachian_inset_waveform.csv',overwrite=TRUE)


