# Retrieve a named NPS boundary and year-matched land cover for real GEDI shots.
library(rGEDI)
library(sf)
library(data.table)
raw <- 'readme/data/grsm_boundary.esri.json'
url <- paste0('https://services1.arcgis.com/fBc8EJBxQRMcHlei/ArcGIS/rest/services/',
              'GRSM_PARK_BOUNDARY/FeatureServer/1/query?',
              'where=1%3D1&outFields=*&outSR=4326&returnGeometry=true&f=json')
if (!file.exists(raw)) download.file(url,raw,mode='wb',quiet=TRUE)
boundary <- st_read(raw,quiet=TRUE)
boundary <- st_sf(unit_id='Great Smoky Mountains National Park',
                  geometry=st_union(st_make_valid(st_transform(boundary,32617))))
# Simplification is for a compact display boundary, in metres.
boundary <- st_simplify(boundary,dTolerance=10,preserveTopology=TRUE)
st_write(st_transform(boundary,4326),'readme/data/grsm_boundary.geojson',delete_dsn=TRUE,quiet=TRUE)
obj <- readLevel2B('readme/data/appalachian_level2b.h5')
d <- getLevel2BPAVDProfile(obj)
g <- obj@h5[[unique(d$beam)]]
d[, `:=`(cover=g[['cover']][],pai=g[['pai']][],degrade=g[['degrade_flag']][],sensitivity=g[['sensitivity']][])]
close(obj)
good <- d[l2b_quality_flag==1 & algorithmrun_flag==1 & degrade==0 & between(sensitivity,.9,1) & between(cover,0,1) & pai>=0]
ee <- ee_initialize(project=Sys.getenv('EE_PROJECT'),authenticate=FALSE,quiet=TRUE)
nlcd <- ee$Image('USGS/NLCD_RELEASES/2019_REL/NLCD/2019')$select('landcover')
values <- extractEE(nlcd,good[,.(shot_number=as.character(shot_number),lon_lowestmode,lat_lowestmode)],scale=30)
stopifnot(nrow(values)==nrow(good),!anyDuplicated(values$shot_number))
fwrite(values[,.(shot_number,landcover)],'readme/data/appalachian_nlcd2019.csv')
cat('Saved NPS boundary and NLCD land cover for',nrow(values),'screened GEDI shots.\n')
