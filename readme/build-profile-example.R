# Run from the repository root. Process measured GEDI data with rGEDI.
library(rGEDI)
library(data.table)
library(sf)
library(ggplot2)
library(patchwork)
obj <- readLevel2B('readme/data/appalachian_level2b.h5')
d <- getLevel2BPAVDProfile(obj)
g <- obj@h5[[unique(d$beam)]]
stopifnot(identical(as.character(d$shot_number),as.character(g[['shot_number']][])))
d[, `:=`(cover=g[['cover']][],pai=g[['pai']][],degrade=g[['degrade_flag']][],sensitivity=g[['sensitivity']][])]
close(obj)
setorder(d,delta_time)
xy <- st_coordinates(st_transform(st_as_sf(d,coords=c('lon_lowestmode','lat_lowestmode'),crs=4326),32617))
d[, distance := c(0,cumsum(sqrt(rowSums((xy[-1,]-xy[-nrow(xy),])^2))))/1000]
good <- d[l2b_quality_flag==1 & algorithmrun_flag==1 & degrade==0 & between(sensitivity,.9,1) & between(cover,0,1) & pai>=0]
park <- st_read('readme/data/grsm_boundary.geojson',quiet=TRUE)
inside <- copy(rGEDI::clip(good,park,split_by='unit_id'))
landcover <- fread('readme/data/appalachian_nlcd2019.csv',colClasses=c(shot_number='character'))
inside[, shot_number := as.character(shot_number)]
inside <- merge(inside,landcover,by='shot_number',all.x=TRUE,sort=FALSE)
l2a <- fread('readme/data/appalachian_l2a.csv',colClasses=c(shot_number='character'))
stopifnot(all(inside$shot_number %in% l2a[l2a_quality_flag_rel3==1]$shot_number))
inside <- merge(inside,l2a[,.(shot_number,rh98,l2a_ground=elev_lowestmode)],by='shot_number',all.x=TRUE,sort=FALSE)

inside[, forest_type := fcase(landcover==41,'Deciduous',landcover==42,'Evergreen',landcover==43,'Mixed',default='Other cover')]
inside[, return_height := elev_highestreturn-elev_lowestmode]
setorder(inside,distance)
inside[, run := cumsum(c(TRUE,diff(distance)>.105))]
runs <- inside[,.(n=.N,lo=min(distance),hi=max(distance)),by=run][order(-n)]
best <- inside[run==runs$run[1]]
metrics <- function(pai,cover,rh98) list(n=length(pai),mean_pai=mean(pai),mean_cover_percent=100*mean(cover),mean_rh98=mean(rh98))
summary <- polyStatsLevel2BVPM(inside,~metrics(pai,cover,rh98),id='poly_id')
forest_summary <- polyStatsLevel2BVPM(inside,~metrics(pai,cover,rh98),id='forest_type')
stopifnot(sum(summary$n)==nrow(inside),sum(forest_summary$n)==nrow(inside),nrow(best)>=10)
fwrite(summary,'readme/profile-polygon-summary.csv')
fwrite(forest_summary,'readme/profile-forest-type-summary.csv')
fwrite(inside,'readme/data/appalachian_profiles_screened.csv')
fwrite(best,'readme/data/appalachian_profile_segment.csv')
d[, status := fifelse(as.character(shot_number) %in% best$shot_number,'Profile segment',fifelse(as.character(shot_number) %in% inside$shot_number,'Screened inside park',fifelse(shot_number %in% good$shot_number,'Screened outside park','Excluded by screen')))]
pts <- st_transform(st_as_sf(d,coords=c('lon_lowestmode','lat_lowestmode'),crs=4326),32617)
park <- st_transform(park,32617)
b <- st_bbox(park)
pad <- 4000
map <- ggplot()+geom_sf(data=park,fill='#EEF2E8',colour='#506A59',linewidth=.65)+
  geom_sf(data=pts,aes(colour=status),size=1.15)+
  scale_colour_manual(values=c('Profile segment'='#CA5C2C','Screened inside park'='#245E86','Screened outside park'='#92A4AB','Excluded by screen'='#D5DADB'),name=NULL)+
  coord_sf(xlim=c(b['xmin']-pad,b['xmax']+pad),ylim=c(b['ymin']-pad,b['ymax']+pad),datum=NA,expand=FALSE)+
  annotate('segment',x=b['xmin']+2000,xend=b['xmin']+22000,y=b['ymin']-1500,yend=b['ymin']-1500,linewidth=1.2)+
  annotate('text',x=b['xmin']+12000,y=b['ymin']+1000,label='20 km',size=3.5)+
  annotate('segment',x=b['xmax']-1500,xend=b['xmax']-1500,y=b['ymax']-7000,yend=b['ymax']-2000,arrow=grid::arrow(length=grid::unit(2.5,'mm')))+
  annotate('text',x=b['xmax']-1500,y=b['ymax']+300,label='N',size=4)+
  labs(title='A  Clip to the national park boundary',subtitle='Great Smoky Mountains National Park')+
  guides(colour=guide_legend(ncol=2,override.aes=list(size=2.5)))+
  theme_void(base_size=11)+theme(legend.position='bottom',legend.text=element_text(size=9),plot.title=element_text(face='bold',size=12))
cols <- grep('^pavd_z[0-9]',names(best),value=TRUE)
m <- melt(best[,c('distance','return_height','elev_lowestmode',cols),with=FALSE],id.vars=c('distance','return_height','elev_lowestmode'),variable.name='bin',value.name='density')
m[, height := as.numeric(sub('^pavd_z([0-9]+)_.*','\\1',bin))+2.5]
m <- m[height<60 & is.finite(density) & density>=0]
native <- dcast(m,height~distance,value.var='density')
heights <- native$height
distances <- as.numeric(names(native)[-1])
values <- as.matrix(native[,-1])
newx <- distances
newz <- seq(0,60,by=.25)
horizontal <- t(apply(values,1,function(v) approx(distances,v,newx,rule=1)$y))
smooth <- apply(horizontal,2,function(v) approx(heights,v,newz,rule=2)$y)
display <- data.table(distance=rep(newx,each=length(newz)),height=rep(newz,times=length(newx)),density=as.vector(smooth))
display[, canopy := approx(best$distance,best$return_height,distance)$y]
display[, ground := approx(best$distance,best$elev_lowestmode,distance)$y]
display <- display[height<=canopy]
edges <- c(newx[1]-(newx[2]-newx[1])/2,(newx[-1]+newx[-length(newx)])/2,
           tail(newx,1)+(tail(newx,1)-newx[length(newx)-1])/2)
display[, left := edges[match(distance,newx)]]
display[, right := edges[match(distance,newx)+1L]]
limit <- unname(quantile(m[density>0]$density,.98))
theme_set(theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold',size=12),legend.position='bottom'))
profile <- ggplot(display,aes(distance,ground+height))+geom_rect(aes(xmin=left,xmax=right,ymin=ground+height-.125,ymax=ground+pmin(height+.125,canopy),fill=density))+
  geom_line(data=best,aes(y=l2a_ground+rh98,colour='L2A RH98 elevation'),linewidth=.7)+
  geom_line(data=best,aes(y=elev_lowestmode,colour='L2B ground'),linewidth=.8)+
  scale_colour_manual(values=c('L2A RH98 elevation'='#6F8D71','L2B ground'='#253F38'),name=NULL)+
  scale_fill_gradientn(colours=c('#F7FBF6','#D8EED2','#8BCA7F','#3A9652','#005A30'),limits=c(0,limit),oob=scales::squish,breaks=c(0,.1,.2),name=expression(PAVD~(m^2/m^3)))+
  coord_cartesian(expand=FALSE)+
  labs(title=sprintf('C  Terrain-following PAVD profile (%d footprints)',nrow(best)),subtitle=sprintf('%.2f km within the park, 31 May 2019',diff(range(best$distance))),x='Distance along track (km)',y='Ellipsoidal elevation (m)')+
  theme(panel.grid=element_blank(),legend.key.width=grid::unit(1.1,'cm'))
profile_metrics <- function(x) list(mean=mean(x),q25=unname(quantile(x,.25)),q75=unname(quantile(x,.75)))
forest_profiles <- melt(inside[forest_type %in% c('Deciduous','Evergreen','Mixed'),
                       c('forest_type',cols),with=FALSE],id.vars='forest_type',
                       variable.name='bin',value.name='density')
forest_profiles[, height := as.numeric(sub('^pavd_z([0-9]+)_.*','\\1',bin))+2.5]
forest_profiles <- forest_profiles[height<60 & is.finite(density) & density>=0]
forest_profiles[, profile_group := paste(forest_type,height,sep='|')]
vertical <- polyStatsLevel2BVPM(forest_profiles,~profile_metrics(density),id='profile_group')
vertical[, forest_type := sub('[|].*','',profile_group)]
vertical[, height := as.numeric(sub('.*[|]','',profile_group))]
fwrite(vertical,'readme/profile-by-forest-type.csv')
palette <- c(Deciduous='#247C78',Evergreen='#3D6796',Mixed='#9B5B8A','Other cover'='#B38C50')
forest_counts <- inside[forest_type %in% c('Deciduous','Evergreen','Mixed'),.N,by=forest_type]
facet_labels <- setNames(sprintf('%s\n(n = %d)',forest_counts$forest_type,forest_counts$N),forest_counts$forest_type)
profile_max <- max(vertical$q75)*1.04
vert <- ggplot(vertical,aes(mean,height,colour=forest_type,fill=forest_type))+
  geom_ribbon(aes(xmin=q25,xmax=q75),orientation='y',alpha=.20,colour=NA)+
  geom_path(linewidth=.9)+facet_wrap(~forest_type,nrow=1,labeller=as_labeller(facet_labels))+
  scale_y_continuous(limits=c(0,60),breaks=c(0,20,40,60))+
  scale_x_continuous(limits=c(0,profile_max),breaks=c(0,.1,.2))+
  scale_colour_manual(values=palette,guide='none')+scale_fill_manual(values=palette,guide='none')+
  labs(title='D  L2B PAVD profiles by forest type',subtitle='Mean and interquartile range of native bins',
       x=expression(PAVD~(m^2/m^3)),y='Height above ground (m)')+
  theme(panel.grid.major.x=element_blank(),strip.text=element_text(size=10),axis.text=element_text(size=9))
counts <- inside[,.N,by=forest_type]
labels <- setNames(sprintf('%s\n(n = %d)',counts$forest_type,counts$N),counts$forest_type)
boxdata <- melt(inside[forest_type!='Other cover',.(forest_type,`Canopy cover (%)`=100*cover,`PAI (m2/m2)`=pai,`RH98 (m)`=rh98)],id.vars='forest_type')
box <- ggplot(boxdata,aes(forest_type,value,fill=forest_type))+geom_boxplot(width=.6,outlier.shape=NA,alpha=.8)+
  geom_jitter(width=.12,size=.9,alpha=.5)+facet_wrap(~variable,scales='free_y',nrow=1)+scale_fill_manual(values=palette,guide='none')+
  scale_x_discrete(labels=labels)+labs(title='E  L2A and L2B metrics by forest type',subtitle='NLCD forests (n = 72, 18 and 57, respectively)',x=NULL,y=NULL)+
  theme(axis.text.x=element_text(size=9,angle=45,hjust=1),strip.text=element_text(size=9),panel.grid.major.x=element_blank())


# A shows several granules. B shows the eight tracks of one native granule.
overview <- fread('readme/data/appalachian_multigranule.csv',colClasses=c(shot_number='character'))
overview[, selection := fifelse(inside,'Clipped to park AOI','Outside park AOI')]
overview_pts <- st_transform(st_as_sf(overview,coords=c('lon_lowestmode','lat_lowestmode'),crs=4326),32617)
single <- fread('readme/data/appalachian_single_granule.csv',colClasses=c(shot_number='character'))
single_pts <- st_transform(st_as_sf(single,coords=c('lon_lowestmode','lat_lowestmode'),crs=4326),32617)
stopifnot(length(unique(single$beam))==8L)
v <- st_coordinates(st_transform(st_as_sf(best,coords=c('lon_lowestmode','lat_lowestmode'),crs=4326),32617))
direction <- (v[nrow(v),]-v[1,]);direction<-direction/sqrt(sum(direction^2))
mid <- colMeans(v)
allxy <- st_coordinates(single_pts)
along <- as.vector(sweep(allxy,2,mid)%*%direction)
detail_pts <- single_pts[abs(along)<=1600,]
stopifnot(length(unique(detail_pts$beam))==8L)
zb <- st_bbox(detail_pts)+c(-250,-250,250,250)
zoom_area <- st_as_sfc(zb)
detail_pts$selection <- ifelse(detail_pts$beam==unique(best$beam),'Profile beam','Other ground tracks')
footprints <- st_buffer(detail_pts,12.5)
wave <- fread('readme/data/appalachian_inset_waveform.csv',colClasses=c(shot_number='character'))
chosen <- best[shot_number==wave$shot_number[1]]
stopifnot(nrow(chosen)==1L,nrow(wave)>100)
chosen_point <- st_transform(st_as_sf(chosen,coords=c('lon_lowestmode','lat_lowestmode'),crs=4326),32617)
map <- ggplot()+geom_sf(data=park,fill='#F0F3E9',colour='#536E5B',linewidth=.65)+
  geom_sf(data=overview_pts,aes(colour=selection),size=.25,alpha=.65)+
  geom_sf(data=st_sf(geometry=zoom_area),fill=NA,colour='#C65D2B',linewidth=.75)+
  scale_colour_manual(values=c('Clipped to park AOI'='#2A7087','Outside park AOI'='#CCD4D6'),name=NULL)+
  coord_sf(xlim=c(b['xmin']-4000,b['xmax']+4000),ylim=c(b['ymin']-4000,b['ymax']+4000),datum=NA,expand=FALSE)+
  annotate('segment',x=b['xmin']+2000,xend=b['xmin']+22000,y=b['ymin']-1500,yend=b['ymin']-1500,linewidth=1.1)+
  annotate('text',x=b['xmin']+12000,y=b['ymin']+1200,label='20 km',size=3.6)+
  annotate('segment',x=b['xmax']-1500,xend=b['xmax']-1500,y=b['ymax']-7000,yend=b['ymax']-2000,arrow=grid::arrow(length=grid::unit(2.5,'mm')))+
  annotate('text',x=b['xmax']-1500,y=b['ymax']+300,label='N',size=4)+
  labs(title='A  Multiple granules clipped to the park AOI',subtitle=sprintf('May to September 2019: %d granules; box locates panel B',uniqueN(overview$granule_key)))+
  theme_void(base_size=11)+theme(legend.position='bottom',plot.title=element_text(face='bold',size=12))
xy <- as.data.table(st_drop_geometry(detail_pts))
coords <- st_coordinates(detail_pts);xy[,`:=`(x=coords[,1],y=coords[,2])]
beam_labels <- xy[order(delta_time),.SD[.N],by=beam]
beam_labels[, label_x := pmin(pmax(x,zb['xmin']+800),zb['xmax']-800)]
beam_labels[, label_y := pmin(pmax(y,zb['ymin']+300),zb['ymax']-300)]
zoom <- ggplot()+geom_sf(data=park,fill='#F7F9F3',colour='#C3CDBE',linewidth=.5)+
  geom_sf(data=footprints,aes(fill=selection),colour=NA)+
  geom_segment(data=beam_labels,aes(x,y,xend=label_x,yend=label_y),colour='#7B8B8A',linewidth=.25)+
  geom_label(data=beam_labels,aes(label_x,label_y,label=beam),size=2.6,fill='white',linewidth=0,label.padding=grid::unit(.08,'lines'))+
  geom_sf(data=chosen_point,aes(shape='Inset waveform'),colour='#A7296B',size=3.2)+
  scale_fill_manual(values=c('Other ground tracks'='#2A7087','Profile beam'='#C65D2B'),name=NULL)+
  scale_shape_manual(values=c('Inset waveform'=8),name=NULL)+
  coord_sf(xlim=zb[c(1,3)],ylim=zb[c(2,4)],datum=NA,expand=FALSE)+
  annotate('segment',x=zb['xmin']+200,xend=zb['xmin']+1200,y=zb['ymin']+200,yend=zb['ymin']+200,linewidth=1.2)+
  annotate('text',x=zb['xmin']+700,y=zb['ymin']+400,label='1 km',size=3.6)+
  annotate('segment',x=zb['xmax']-250,xend=zb['xmax']-250,y=zb['ymax']-900,yend=zb['ymax']-300,arrow=grid::arrow(length=grid::unit(2.4,'mm')))+
  annotate('text',x=zb['xmax']-250,y=zb['ymax']-100,label='N',size=4)+
  labs(title='B  One granule with all eight ground tracks',subtitle='Orbit 02629, part 03; circles are 25 m footprints')+
  theme_void(base_size=11)+theme(legend.position='bottom',plot.title=element_text(face='bold',size=12),legend.text=element_text(size=9))
setorder(wave,height)
polygon <- rbind(data.table(amplitude=0,height=min(wave$height)),wave[,.(amplitude,height)],data.table(amplitude=0,height=max(wave$height)))
wave_inset <- ggplot(wave,aes(amplitude,height))+geom_polygon(data=polygon,fill='#A6CFBB',colour=NA)+
  geom_path(colour='#285F67',linewidth=.65)+geom_hline(yintercept=0,colour='#4F625B',linetype=2,linewidth=.35)+
  geom_hline(yintercept=chosen$rh98,colour='#A7296B',linetype=2,linewidth=.4)+
  coord_cartesian(xlim=c(0,100),ylim=c(-10,60),expand=FALSE)+
  scale_x_continuous(breaks=c(0,50,100))+scale_y_continuous(breaks=c(0,20,40,60))+
  labs(title='Observed L1B waveform',subtitle=sprintf('RH98 %.1f m',chosen$rh98),x='Peak return (%)',y='Height (m)')+
  theme_minimal(base_size=8.5)+theme(panel.grid.minor=element_blank(),plot.background=element_rect(fill='white',colour='#91A69C'),
      panel.grid.major=element_line(colour='#EDF0EC',linewidth=.25),plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=8),plot.margin=margin(4,5,4,5))
profile <- profile+geom_point(data=chosen,aes(x=distance,y=l2a_ground+rh98),shape=8,colour='#A7296B',size=2.2)+
  labs(title=sprintf('C  PAVD profile with a real waveform inset (%d footprints)',nrow(best)))
# Reserve a margin for the inset, so it never obscures measured profiles.
profile_base <- profile+theme(plot.margin=margin(5,180,5,5))
profile <- profile_base+inset_element(wave_inset,left=.80,bottom=.16,right=.99,top=.96,align_to='full')
combined <- (map+zoom+plot_layout(widths=c(1.35,1)))/profile/(vert+box)+plot_layout(heights=c(1.1,1.05,1))
ggsave('readme/fig-forest-characterization.png',combined,width=12.4,height=11.3,dpi=400,bg='white')
ggsave('readme/fig-forest-characterization.pdf',combined,width=12.4,height=11.3,bg='white')
for (name in c('map','zoom','profile','profile_base','wave_inset','vert','box')) {
  width <- if(name %in% c('profile','profile_base')) 12.4 else if(name=='zoom')6 else if(name=='wave_inset')3 else 7
  height <- if(name %in% c('profile','profile_base'))4.2 else if(name=='zoom')6 else if(name=='wave_inset')3 else 5
  ggsave(paste0('readme/data/profile_',name,'.png'),get(name),width=width,height=height,dpi=350,bg='white')
}
cat(sprintf('Verified %d ground tracks and an L1B inset with %d waveform bins.\n',uniqueN(single$beam),nrow(wave)))
