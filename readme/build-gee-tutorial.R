# Complete upscaling workflow, also printed step by step in README.md.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload",quietly=TRUE)) pkgload::load_all(".",quiet=TRUE) else library(rGEDI)
  library(sf);library(terra);library(data.table)
})
main <- function() {
  root <- Sys.getenv("RGEDI_EXAMPLE_CACHE",file.path(path.expand("~"),"rGEDI-cache"))
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  dir.create("readme/results",recursive=TRUE,showWarnings=FALSE)
  emit <- function(x,name) writeLines(trimws(capture.output(print(x)),which="right"),file.path("readme/results",name))
  if (!nzchar(Sys.getenv("EE_PROJECT"))) stop("Set EE_PROJECT to an authenticated project.")
  ee <- ee_initialize(project=Sys.getenv("EE_PROJECT"),authenticate=FALSE)
  ee$data$setDeadline(180000L)

  # 1. Retrieve GEDI observations and screen each acquisition year.
  aoi <- st_read(system.file("extdata","shenandoah.shp",package="rGEDI"),quiet=TRUE)
  box <- st_bbox(aoi)
  source_years <- c(2019L,2021L,2022L)
  yearly_gedi <- list();counts <- list()
  for (year in source_years) {
    raw_file <- file.path(root,paste0("shenandoah_raw_",year,".rds"))
    if (file.exists(raw_file)) raw <- readRDS(raw_file) else {
      address <- gedifinder("GEDI04_A",ul_lat=box[["ymax"]],ul_lon=box[["xmin"]],
        lr_lat=box[["ymin"]],lr_lon=box[["xmax"]],version="002.1",source="gee",
        daterange=c(sprintf("%s-05-01",year),sprintf("%s-10-01",year)))
      level4a <- readLevel4A(address)
      lazy <- clip(getLevel4A(level4a,quality=TRUE,collect=FALSE,max_granules=500L),aoi)
      raw <- getLevel4A(lazy,quality=TRUE,max_features=10000L,max_granules=500L,
        cols=c("beam","shot_number","agbd","agbd_se","sensitivity","solar_elevation","l4_quality_flag","degrade_flag"))
      saveRDS(raw,raw_file)
    }
    d <- clip(raw,aoi)
    d <- d[is.finite(agbd) & agbd>=0 & sensitivity>=.9 & sensitivity<=1 &
      solar_elevation<0 & as.character(beam) %in% c("5","6","8","11","BEAM0101","BEAM0110","BEAM1000","BEAM1011")]
    d[,`:=`(sample_id=seq_len(.N)+year*100000L,year=year)]
    yearly_gedi[[as.character(year)]] <- d
    counts[[as.character(year)]] <- data.table(year=year,quality_screened=nrow(raw),night_power=nrow(d))
  }
  cohort_counts <- rbindlist(counts)
  emit(cohort_counts,"upscaling-cohorts.txt")

  # 2. Extract annual ancillary predictors at the screened footprints.
  all <- list()
  for (year in source_years) {
    dest <- file.path(root,paste0("shenandoah_predictors_",year,".rds"))
    if (file.exists(dest)) d <- readRDS(dest) else {
      d <- yearly_gedi[[as.character(year)]]
      predictors <- ee_build_AlphaEarth_embedding_terrain_stack(aoi,year,year,add_lonlat=FALSE)
      bands <- c(sprintf("A%02d",0:63),"elevation","slope","aspect")
      coordinates <- d[,.(sample_id,lon_lowestmode,lat_lowestmode)]
      d <- extractEE(predictors$select(as.list(bands)),d,scale=30,chunk_size=1000L)
      d <- merge(d,coordinates,by="sample_id",all.x=TRUE,sort=FALSE)
      d <- d[complete.cases(d[,c("agbd",bands),with=FALSE])]
      d[,`:=`(aspect_sin=sin(aspect*pi/180),aspect_cos=cos(aspect*pi/180))]
      saveRDS(d,dest)
    }
    all[[as.character(year)]] <- d
  }
  predictor_names <- c(sprintf("A%02d",0:63),"elevation","slope","aspect_sin","aspect_cos")
  predictor_preview <- head(rbindlist(all)[,.(year,agbd,A00,A01,elevation,slope,aspect_sin,aspect_cos)],6)
  emit(predictor_preview,"upscaling-predictors.txt")
  stack <- ee_build_AlphaEarth_embedding_terrain_stack(aoi,2021,2021,add_lonlat=FALSE)
  aspect <- stack$select("aspect")$multiply(pi/180)
  stack <- stack$addBands(aspect$sin()$rename("aspect_sin"))$addBands(aspect$cos()$rename("aspect_cos"))
  band_names <- unlist(reticulate::py_to_r(stack$bandNames()$getInfo()))
  predictor_inventory <- data.table(group=c("AlphaEarth embeddings","Terrain predictors"),
    bands=c(sum(grepl("^A[0-9]{2}$",band_names)),sum(predictor_names %in% c("elevation","slope","aspect_sin","aspect_cos"))))
  emit(predictor_inventory,"upscaling-bands.txt")
  display_file <- file.path(root,"shenandoah_predictors_display_2021.tif")
  if (!file.exists(display_file)) map_download(stack$select(as.list(c("A00","A20","A40","elevation","slope"))),
    display_file,region=aoi,scale=100,crs="EPSG:32617")

  # 3. Thin observations, then partition complete spatial blocks.
  complete_counts <- rbindlist(lapply(all,function(d) data.table(year=d$year[1],complete_predictors=nrow(d))))
  set.seed(726)
  all <- lapply(all,function(d) rGEDI::sample(d,spacedSampling(size=min(1000L,nrow(d)),radius=60)))
  data <- rbindlist(all,fill=TRUE)
  xy <- st_coordinates(st_transform(st_as_sf(data,coords=c("lon_lowestmode","lat_lowestmode"),crs=4326),32617))
  data[,block:=paste(floor(xy[,1]/2000),floor(xy[,2]/2000),sep="_")]
  set.seed(726)
  blocks <- base::sample(unique(data$block),ceiling(.3*length(unique(data$block))))
  train <- data[!block %in% blocks];valid <- data[block %in% blocks]
  stopifnot(!length(intersect(train$block,valid$block)))
  retained_counts <- data[,.(retained=.N),by=year]
  sample_counts <- merge(merge(cohort_counts,complete_counts,by="year"),retained_counts,by="year")
  emit(sample_counts,"upscaling-samples.txt")
  partition_counts <- data[,.(footprints=.N,blocks=uniqueN(block)),by=.(partition=fifelse(block %in% blocks,"Validation","Training"))]
  emit(partition_counts,"upscaling-partitions.txt")
  fwrite(data[,.(sample_id,year,beam,shot_number,lon_lowestmode,lat_lowestmode,agbd,block,
    partition=fifelse(block %in% blocks,"Validation","Training"))],"readme/results/shenandoah_footprints.csv")

  # 4. Select variables using training blocks only.
  selection_file <- file.path(root,"shenandoah_selection.rds")
  previous_file <- file.path(root,"shenandoah_model_inputs.rds")
  reuse_selection <- FALSE
  if (file.exists(selection_file) && file.exists(previous_file)) {
    previous <- readRDS(previous_file)
    previous_train <- previous$training[!block %in% previous$heldout_blocks]
    fields <- c("agbd",predictor_names)
    reuse_selection <- identical(train$sample_id,previous_train$sample_id) &&
      isTRUE(all.equal(as.data.frame(train[,..fields]),as.data.frame(previous_train[,..fields])))
  }
  if (reuse_selection) selection <- readRDS(selection_file) else {
    selection <- varSel(train[,..predictor_names],train$agbd,method="rfe",threshold=0,
      seed=726,ntree=150,nodesize=5)
    saveRDS(selection,selection_file)
  }
  keep <- selection$selvars
  selected_variables <- as.data.table(selection$importance)[selected==TRUE,.(parameter,importance)][order(-importance)]
  emit(selected_variables,"upscaling-selected.txt")
  fwrite(selection$importance,"readme/results/shenandoah_importance.csv")
  fwrite(selection$test,"readme/results/shenandoah_rfe.csv")

  # 5. Fit and validate the local model that will be mapped.
  fit_rf <- fit_model(train[,..keep],train$agbd,method="randomForest",test="none",
    seed=726,ntree=200,maxnodes=128,nodesize=5)
  emit(fit_rf$model,"upscaling-model.txt")
  pred_train <- as.numeric(predict(fit_rf,train))
  pred_valid <- as.numeric(predict(fit_rf,valid))
  predictions <- rbindlist(list(data.table(sample_id=train$sample_id,agbd=train$agbd,predicted=pred_train,stage="Training"),
    data.table(sample_id=valid$sample_id,agbd=valid$agbd,predicted=pred_valid,stage="Validation")))
  fwrite(predictions,"readme/results/shenandoah_predictions.csv")
  train_metrics <- as.data.table(fit_metrics(train$agbd,pred_train,unit="Mg/ha"))
  valid_metrics <- as.data.table(fit_metrics(valid$agbd,pred_valid,unit="Mg/ha"))
  metric_table <- merge(train_metrics,valid_metrics,by=c("stat","unit"),suffixes=c("_training","_validation"),sort=FALSE)
  emit(metric_table,"upscaling-metrics.txt")
  fwrite(rbind(copy(train_metrics)[,stage:="Training"],copy(valid_metrics)[,stage:="Validation"]),"readme/results/shenandoah_accuracy.csv")
  saveRDS(fit_rf,file.path(root,"shenandoah_local_rf.rds"))
  saveRDS(list(training=data,selected=keep,heldout_blocks=blocks),file.path(root,"shenandoah_model_inputs.rds"))

  # 6. Map this same forest; map_create handles conversion internally.
  image <- map_create(model=fit_rf,stack=stack$select(as.list(keep)),aoi=aoi,name="AGBD")
  dest <- file.path(root,"shenandoah_local_rf_agbd_2021.tif")
  map_download(image,dest,region=aoi,scale=100,crs="EPSG:32617",overwrite=TRUE)
  mapped <- terra::rast(dest)
  v <- terra::values(mapped,mat=FALSE);v<-v[is.finite(v)]
  map_summary <- data.table(resolution_m=res(mapped)[1],mapped_cells=length(v),
    min_agbd=min(v),mean_agbd=mean(v),max_agbd=max(v))
  emit(map_summary,"upscaling-map.txt")
  fwrite(data.table(n=nrow(data),train=nrow(train),validation=nrow(valid),predictors=length(keep),trees=200L,
    map_year=2021L,map_resolution_m=100L,model="Local randomForest; no EE retraining"),"readme/results/shenandoah_run.csv")
  cat("LOCAL MODEL UPSCALING COMPLETE\n")
}
main()
