######################################################################
#
#  PNW-TREE vital rates -- FIA data processing
#
#    Rob Smith, phytomosaic@gmail.com, 29 Mar 2020, ver. 21 May 2020
#
##      GNU General Public License, Version 3.0    ###################


### preamble
rm(list=ls())
# devtools::install_github('phytomosaic/ecole')
require(ecole)
require(data.table)
require(raster)
require(FNN) # to assign fire values from nearest neighbors, when NA
setwd('./data_raw/fia_data/')

### timing
t_start <- Sys.time()

### function to read FIADB file
`f` <- function(x, keepcols=NULL, ...) {
  cat('File size (MB):',round(file.info(x)$size/1024^2,1),'--> ')
  cat('Time elapsed:',
      system.time({
        out <- data.table::fread(x, header = T, dec = '.',
                                 fill = T, data.table = T,
                                 select = keepcols,
                                 integer64='character', ...)
      })[[3]], '\n')
  return(out)
}

### read GRM files
d <- rbindlist(lapply(list.files('.',pattern='_GRM_COMPONENT.csv'), f,
                      keepcols=c('TRE_CN','PREV_TRE_CN','PLT_CN',
                                 'DIA_BEGIN','DIA_MIDPT','DIA_END',
                                 'ANN_DIA_GROWTH','STEM_COMPONENT')))
setnames(d, names(d), tolower(names(d))) # cleanup column names

### read TREE files
tr <- rbindlist(lapply(list.files('.', pattern='_TREE.csv'), f,
                       keepcols=c('CN','PLT_CN','INVYR','SPCD',
                                  'DIA','AZIMUTH','DIST','SUBP')))
setnames(tr, names(tr), tolower(names(tr))) # cleanup column names
setnames(tr, 'cn', 'tre_cn')                # rename
# assign subplot identifier
tr$subpid <-paste(formatC(tr$plt_cn,wid=2,format='d',flag='0'),
                  formatC(tr$subp,wid=2,format='d',flag='0'), sep='_')
# calculate stem coordinates, and identify seedlings and edge trees
`stem_coord` <- function(a, d, x0=0, y0=0) {
  data.frame(x = x0 + d * cos(a / 180 * pi),
             y = y0 + d * sin(a / 180 * pi))
}
tr <- cbind(tr, stem_coord(-tr$azimuth, tr$dist,0,0)) ; rm(stem_coord)
is_seed <- (tr$dia < 5
            & !is.na(tr$dia)
            & !is.na(tr$x)
            & !(abs(tr$dist) > 6.8))
tr$x[is_seed] <- tr$x[is_seed] + 12 # offset seedlings by 12 ft
tr$is_seed <- is_seed  ;  rm(is_seed)
`get_dist` <- function(x, y) { sqrt((x - 0)^2 + (y - 0)^2) }
tr[, dist      := get_dist(x,y)] # recalc dist AFTER offset seedlings
tr[, is_macro  := any(abs(dist) > 24), by=plt_cn] # 58.9-ft macroplot?
tr[, edge_dist := ifelse(is_macro, 58.9, 24.0)] # assign size
tr[, is_edge := (abs(edge_dist - dist) < 6)] # arbitrary 6-ft buffer
tr[, c('is_macro','edge_dist') := NULL] # squash unneeded columns

### read PLOT files
p <- rbindlist(lapply(list.files('.', pattern='_PLOT.csv')[-1], f,
                      keepcols=c('CN','STATECD','MEASYEAR',
                                 'LAT','LON','ELEV','ECOSUBCD',
                                 'REMPER')))
setnames(p, names(p), tolower(names(p)))    # cleanup column names
setnames(p, 'cn', 'plt_cn')                 # rename
p[, elev := round(elev * 0.3047851,0)]      # ft to m
p$ecosubcd <- ecole::clean_text(p$ecosubcd) # trim leading whitespace

### read REF_SPECIES file
ref <- read.csv('./REF_SPECIES.csv', stringsAsFactors=F)
ref <- ref[,colnames(ref) %in% c('SPCD','COMMON_NAME','GENUS',
                                 'SPECIES', 'VARIETY','SUBSPECIES')]
names( ref ) <- tolower(names( ref ))     # names to lowercase
ref$spp <- paste0(ref$genus,'_',ref$species,'_',
                  ref$variety,'_',ref$subspecies)
ref$spp <- ecole::clean_text(ref$spp, lower=TRUE)
ref     <- as.data.table(ref[,c('spcd','spp')])

### reference table for ecological subsections (to split east-west)
e <- read.csv('./REF_ECOSUBSECTION.csv', stringsAsFactors=F)
names(e) <- tolower(names(e))
e <- as.data.table(e)
e[, province_name := NULL] # drop one column

### merging
d <- tr[d, on = c(tre_cn = 'tre_cn')]    # get tree attributes at *y2*
d <- p[d, on = c(plt_cn = 'plt_cn')]     # get plot info
d <- ref[d, on = c(spcd = 'spcd')]       # get species names
d <- e[d, on = c(ecosubcd = 'ecosubcd')] # merge on ecosubcd
d[, grep('^i\\.', names(d), value=T) := NULL] # omit dupe columns
setnames(d, 'measyear', 'measyr_2')      # rename
d$measyr_1 <- round(d$measyr_2 - d$remper,0)  # original measyear?
rm(tr, p, ref, f, e, get_dist)

### keep only plots from northern CA, all OR, all WA
d <- d[lat > 38.72,]

### reclassify sides of Cascades, based on elevation (east-west/lo-hi)
d$ewh <- d$ew <- ifelse(d$assigned %in% c('east','eastmont'),
                        'east', 'west')
d$ewh[d$ew == 'west' & d$elev > 1000] <- 'westhi' # 1100 m ~3300 ft
d$ewh[d$ew == 'east' & d$elev > 1500] <- 'easthi' # 1500 m ~4900 ft

### remove some invalid 'species'
d <- d[!c(spp %in% c('tree_broadleaf','tree_evergreen'))]

### remove any infrequent species (< 49 individuals)
rev(sort(table(d$spp)))
d <- d[, frq := .N, by = spp][(frq >= 49)] # remove infrequent spp
d[, frq := NULL ]  # omit frequency column

### change a few names to be compatible with TRY-DB
d$spp[d$spp == 'chrysolepis_chrysophylla_chrysophylla'] <-
  'chrysolepis_chrysophylla'
d$spp[d$spp == 'populus_balsamifera_trichocarpa'] <-
  'populus_trichocarpa'
d$spp[d$spp == 'abies_shastensis'] <-
  'abies_magnifica_var_shastensis'

### temporarily fill NA=0 for new recruits (for next step)
d$dia_begin[is.na(d$dia_begin) & d$stem_component == 'INGROWTH'] <- 0

### *annualize* Y1/Y2 diameters using 'ann_dia_growth'
d$dia_y1   <- d$dia_begin
d$dia_y2   <- d$dia_begin + d$ann_dia_growth

### remove invalid/removed/logged trees
table(d$stem_component)
is_invalid <- d$stem_component %in% c('REMEASURED DEAD','NOT USED',
                                      'UNKNOWN')
was_harvested  <- d$stem_component %in% c('CUT DEAD', 'CUT1', 'CUT2')
sum(c(is_invalid | was_harvested)) / NROW(d) # nearly 20% are invalid
d <- d[!c(is_invalid | was_harvested),]      # CAUTION!

### assign tree fate, recode by recruit/mortality
is_newrecruit  <- d$stem_component == 'INGROWTH'
is_naturaldead <- d$stem_component %in% c('MORTALITY1','MORTALITY2')
is_survivor    <- d$stem_component == 'SURVIVOR'
d$dia_y1[is_newrecruit]  <- NA # if new recruit then expect no Y1
# if dead then expect no Y2 !CAUTION! dead trees still have 'size':
d$dia_y2[is_naturaldead] <- NA
rm(is_invalid,is_naturaldead,is_newrecruit,is_survivor,was_harvested)
# CAUTION: remove if both Y1/Y2 == NA (and therefore untrackable)
sum(is.na(d$dia_y2) & is.na(d$dia_y1)) / NROW(d) # ~16% invalid
d <- d[!(is.na(d$dia_y2) & is.na(d$dia_y1)),]

### tally recruits/mortality/survivors
### identify survivors (if has measurement at both events)
d$surv <- ifelse(c(!is.na(d$dia_y2) & !is.na(d$dia_y1)),1,0)
sum(d$surv)                             # tally survival
sum(is.na(d$dia_y2))                    # tally mortality
sum(is.na(d$dia_y1) & !is.na(d$dia_y2)) # tally new recruits

# remove *negative* diameter growth (only a few hundred cases)
isneg <- c(d$dia_y1 > d$dia_y2)
isneg <- isneg & !is.na(isneg)
sum(isneg) / NROW(d) # 0.2% of cases
d <- d[!isneg,]
rm(isneg)

### get climate info:
a <- d[!duplicated(d$plt_cn),c('plt_cn','lat','lon','elev','measyr_1')]
dim(a) # 13716 unique plots
### --> --> --> go to ClimateNA --> --> --> --> --> -->
prep_climatena(ID1=a$plt_cn, lat=a$lat, long=a$lon, el=a$elev,
               file = 'C:/Users/Rob/Desktop/to_climatena.csv')
### --> --> --> return from ClimateNA --> --> --> --> --> -->
file.copy(from = 'C:/Users/Rob/Desktop/from_climatena.csv',
          to = './from_climatena.csv',
          overwrite = TRUE)   # copy for posterity
### --> --> --> return from ClimateNA --> --> --> --> --> -->
aa <- read.csv('./from_climatena.csv')
aa[aa==(-9999)] <- NA   # replace NA value with NA
all(unique(a$plt_cn) %in% unique(as.character(aa$ID1)))
a <- cbind(a, aa[match(a$plt_cn, as.character(aa$ID1)),])
names(a) <- tolower(names(a))
j <- names(a)[-c(grep(paste( # trim some not-useful columns
  c('dd_0','dd5','dd_18','dd18','bffp','effp','mar','rad',
    'id1','id2','elevation','latitude','longitude'),
  collapse='|'), names(a), value=F))]
a <- a[ , ..j]  ;  rm(j)
d <- cbind(d,a[match(d$plt_cn, a$plt_cn),-c(1:5)]) # match clim data
rm(aa,a)

### squash unneeded columns
d[,c('assigned','province','section','section_name','dia','dia_begin',
     'dia_midpt','dia_end','ann_dia_growth','stem_component',
     'azimuth', 'invyr', 'prev_tre_cn' ) := NULL]

### sort remaining columns
names(d)
j <- c('plt_cn','tre_cn','statecd','lat','lon','elev','subp','subpid',
       'ecosubcd','eco_subsection_name','ew','ewh', 'spcd','spp',
       'dist','x','y','is_seed','is_edge', # 'azimuth',
       'remper','measyr_1','measyr_2','dia_y1','dia_y2','surv',
       'tmax_wt', 'tmax_sp', 'tmax_sm', 'tmax_at', 'tmin_wt',
       'tmin_sp', 'tmin_sm', 'tmin_at', 'tave_wt', 'tave_sp',
       'tave_sm', 'tave_at', 'ppt_wt', 'ppt_sp',  'ppt_sm',  'ppt_at',
       'nffd_wt', 'nffd_sp', 'nffd_sm', 'nffd_at', 'pas_wt','pas_sp',
       'pas_sm',  'pas_at', 'eref_wt','eref_sp','eref_sm','eref_at',
       'cmd_wt',  'cmd_sp', 'cmd_sm', 'cmd_at',
       'rh_wt', 'rh_sp', 'rh_sm', 'rh_at')
d <- d[ , ..j]

### climate PCA
j <- grep('tmax|tmin|tave|ppt|nffd|pas|eref|cmd|rh', names(d))
j <- j[j!=62]
anyNA(d[,..j]) # expect FALSE
(pc <- prcomp(d[,..j], retx=T, center=T, scale.=T, rank.=2))
summary(pc)    # first 2 explain 82.1% variance = 50.6 + 31.4
V <- t(apply(pc$rotation,1,function(a,b){a*b},pc$sdev))[,1:2] # corr
V[,1]    <-  (-1) * V[,1]     # reverse sign of PC1
pc$x[,1] <-  (-1) * pc$x[,1]  # reverse sign of PC1
# write.csv(V, file='./V_climate_pca.csv')
colnames(pc$x) <- tolower(colnames(pc$x))
d <- cbind(d, pc$x[,1:2], stringsAsFactors=F)
### plot the PCs
`h` <- function (x) {
  r <- c("#5E4FA2", "#4F61AA","#4173B3", "#3386BC", "#4198B6",
         "#51ABAE", "#62BEA6", "#77C8A4", "#8ED1A4", "#A4DAA4",
         "#B8E2A1", "#CBEA9D", "#DEF199", "#EAF69F", "#F2FAAC",
         "#FAFDB8", "#FEFAB6", "#FEF0A5", "#FEE695", "#FDD985",
         "#FDC978", "#FDB96A", "#FCA75E", "#F99254", "#F67D4A",
         "#F26943", "#E85A47", "#DE4B4B", "#D33C4E", "#C1284A",
         "#AF1446", "#9E0142")
  image(1:ncol(x), 1:nrow(x), t(x), col=r, axes=F, xlab='', ylab='')
  axis(3, at=1:ncol(x),lab=c('PC1','PC2'),las=1,tick=F,cex.axis=0.6)
  axis(2, at=1:nrow(x), lab=rownames(x), las=1, cex.axis=0.6)
}
png('./../../fig/fig_01_climate_PCA.png',9,3,'in',bg='transparent',
    res=1000)
set_par(3) ; par(font=2, tcl=-0.2, mgp=c(1.6,0.4,0))
h(V)  # loadings
abline(h=seq(1,NROW(V)+4,by=4)-0.5) ; box()
i <- !duplicated(d$plt_cn) # one point per site
plo(d$lon[i],d$lat[i],cex=0.5,col=colvec(d$pc1[i]),asp=1.6,ylab='',
    xlab='')
add_text(0.03,0.95,'PC1 (Thermal)', cex=1.2)
plo(d$lon[i],d$lat[i],cex=0.5,col=colvec(d$pc2[i]),asp=1.6,ylab='',
    xlab='')
add_text(0.03,0.95,'PC2 (Aridity)', cex=1.2)
dev.off()
rm(pc, V, h)
# PC1 (50.6% varexpl) = thermal: less montane, hotter temps, little
#                                   precip as snow, greater evapn
# PC2 (31.4% varexpl) = aridity: larger moisture deficit, lower RH,
#                                   less precip

### calc ANNUAL MEAN values for some climate variables
j     <- grep('nffd', names(d))
d$nffd <- rowSums(d[,..j], na.rm=T) # n frost-free days
j     <- grep('ppt', names(d))
d$maplog <- log10(1+rowSums(d[,..j], na.rm=T)) # mean ann precip
j     <- grep('tave', names(d))
d$mat <- rowMeans(d[,..j], na.rm=T) # mean seasonal temp
j     <- grep('cmd', names(d))
d$cmd <- rowMeans(d[,..j], na.rm=T) # mean climatic moisture deficit
d$cmdlog  <- log10(d$cmd + 1)       # log(x+1) of CMD
d$td  <- d[, tave_sm - tave_wt]     # continentality
rm(j)


###################################################################
####    BEGIN fire probabilities    #########################
###################################################################
#  Preceding steps (in 2019) to calculate burn probabilities:
#    0. Get burn probability data as .GDB from:
#    Short KC, Finney MA, Scott JH, Gilbertson-Day JW, Grenfell IC.
#         2016. Spatial dataset of probabilistic wildfire risk
#         components for the conterminous United States. 1st Edition.
#         Fort Collins, CO: Forest Service Research Data Archive.
#         https://doi.org/10.2737/RDS-2016-0034
#    1. use ARCMAP to export tiffs from the .gdb
#    2. calc PRODUCT of burn probability and conditional flame length,
#         where probs = prob of burn of that intensity
#    3. SUM all 'impactful' burn probabilities where FL > 1.2 m
#    4. SMOOTH that raster as moving average in radius = 2000 m = 2 km
#    5. extract values for each plot location (code below).
#
### get 13716 unique plot locations, and get the fire raster:
u <- unique(d[,c('plt_cn','lat','lon','cmdlog')], by=c('plt_cn'))
r <- raster('../gis/fire/FL_smooth.tif') # burn prob raster
trg_prj <- projection(r)                 # Albers equal-area proj
`make_xy` <- function(xy, crs, ...) {    # reproject points
  xy <- as.data.frame(xy)
  xy <- xy[!(is.na(xy$lon) | is.na(xy$lat)),]  # rm NAs
  coordinates(xy) <- ~lon+lat
  proj4string(xy) <- CRS('+init=epsg:4269')    # NAD83 dec deg (FIA)
  return(spTransform(xy, crs))
}
xy <- make_xy(u, crs=trg_prj) # reproject Albers equal-area proj

### extract values at query points, per chunk
`extract_raster` <- function(r, xy, breaks=10, ...){
  if(class(r)!='RasterLayer') {
    stop('r must be `RasterLayer` object')
  }
  if(!inherits(xy,'SpatialPoints')) {
    stop('xy must be `SpatialPoints` object')
  }
  # divide extent into equal area chunks for speed
  `make_chunk` <- function(e, breaks, ...){
    m  <- matrix(NA, nrow=breaks, ncol=4,
                 dimnames=list(NULL,
                               c('xmin','xmax','ymin','ymax')))
    rng <- range(c(e[1],e[2]))  # x dimension
    b   <- seq(rng[1], rng[2], length=breaks+1)
    for(i in 1:(breaks)){
      m[i,1:2] <- c(b[i], b[i+1])
    }
    rng <- range(c(e[3],e[4]))  # y dimension
    b   <- seq(rng[1], rng[2], length=breaks+1)
    for(i in 1:(breaks)){
      m[i,3:4] <- c(b[i], b[i+1])
    }
    m <- as.matrix(merge(m[,1:2],m[,3:4]))
    m
  }
  ck     <- make_chunk(e=extent(r), breaks=breaks)
  nchunk <- dim(ck)[[1]]
  npts   <- dim(xy@coords)[[1]]
  p      <- matrix(NA, nrow=npts, ncol=nchunk)
  pb     <- pbCreate(nchunk, progress='text', style=3, ...)
  for(i in 1:nchunk){
    pbStep(pb, i)
    p[,i] <- raster::extract(crop(r, extent(ck[i,])), xy,
                             progress='text')
  }
  # collect chunks into one vector 'tmp'
  nonNA <- !is.na(p)
  tmp   <- vector('numeric', npts)
  for(i in 1:npts){
    if( all(!nonNA[i,]) ){
      tmp[i] <- NA   # some points were not queryable
    } else {
      tmp[i] <- p[i, min(which(nonNA[i,]))]
    }
  }
  # assemble matrix
  p <- matrix(c(1:npts, tmp), nrow=npts, ncol=2,
              dimnames=list(NULL, c('pt', 'val')))
  pbClose(pb)
  return(p)
}
v      <- extract_raster(r, xy, breaks=10)   # fire prob values
u$fire <- v[,2]                              # assign values
# d$fire <- u$fire[match(d$plt_cn, u$plt_cn)]  # match
rm(r,v,xy,make_xy,extract_raster,trg_prj)

### assign nearest-neighbor value if fire is NA (mostly at boundaries)
i <- is.na(u$fire)
x <- FNN::get.knn(u[,.(lon,lat)], k=50)$nn.index  # index 50 neighbors
x <- matrix(u$fire[x], nrow=NROW(x), ncol=NCOL(x))[i,]
j <- apply(x, 1, function(x) which.min(is.na(x))) # col index nearest
u$fire[i] <- x[cbind(1:NROW(x), j)]    # assign nearest non-NA
u$firelog <- u$fire^(1/3)              # cube-root transformation
d <- cbind(d, u[match(d$plt_cn, u$plt_cn), .(fire, firelog)]) # match

### plot the environmental gradients
png('./../../fig/fig_02_fire_and_cmd.png', 6.5, 3.0, 'in',
    bg='transparent', res=1000)
set_par_mercury(3, mgp=c(1.8,0.2,0), oma=c(0,1,0,0))
x1 <- expression(log[10] ~ CMD ~ (mm ~ y^{-1}))
x2 <- expression(Fire ~ probability^{1/3})
cc <- colvec(u[,firelog], alpha=0.7)
plot(u[,.(cmdlog,firelog)], col=cc, cex=0.4, xlab=x1, ylab=x2)
plot(u[,.(lon, lat)], col=colvec(u[,cmdlog], alpha=0.7), cex=0.3,
     xlab='', ylab='', asp=1.2)
add_text(0.7, 0.2, 'CMD')
plot(u[,.(lon, lat)], col=cc, cex=0.3, xlab='', ylab='', asp=1.2)
add_text(0.7, 0.2, 'Fire\nprobs')
dev.off()
rm(i,j,x,u,cc) # cleanup
###################################################################
####    END fire probabilities    ###########################
###################################################################


###################################################################
####    BEGIN nci coefficients    ###################
###################################################################
`calc_nci` <- function(d, dist_cutoff=NULL, size_cutoff=NULL) {
  message('Expected time = ',
          round( length(unique(d$plt_cn)) * 0.04 / 60, 1),
          ' minutes')
  t_start <- Sys.time()
  ### check: use STARTING dia_y1 if available, otherwise dia_y2
  isna           <- is.na(d$dia_y1)
  d$dia_y1[isna] <- d$dia_y2[isna]
  d$DBH          <- d$dia_y1 ^ 2   # new column DBH squared
  ### subset based on distance or size cutoffs
  if (is.numeric(dist_cutoff)) { # keep in buffer = 18.0 ft?
    d <- d[d$dist <= dist_cutoff,]
  }
  if (is.numeric(size_cutoff)) { # keep non-seedlings = 5.0 in?
    d <- d[d$DBH > size_cutoff,]
  }
  ### begin loop for each subplot ! ! ! TIMEWARN ! ! !
  i_nci <- sapply(sort(unique(d$subpid)), function(k) {
    x <- d[d$subpid==k,]
    nr <- NROW(x)
    if(nr == 1) {
      return(0) # zero if only 1 individual in the subplot
    } else if (nr > 1) {
      m       <- as.matrix(dist(x[,c('x','y')],upper=T)) ^ 2
      diag(m) <- NA # ignore self
      DEN     <- rowSums(m, na.rm=T) # sums of neighbor dists
      m       <- matrix(x$DBH, NROW(x), NROW(x), byrow=T)
      diag(m) <- NA # ignore self
      NUM     <- rowSums(m, na.rm=T) # sums of neighbor DBH
      return(NUM/DEN) # NCI for each focal tree
    }
  })
  d$i_nci <- unname(unlist(i_nci)) # vector of INDIVIDUAL NCIs
  d$i_nci[is.infinite(d$i_nci) | is.nan(d$i_nci)] <- NA
  print(Sys.time() - t_start)
  return(d$i_nci) # return INDIVIDUAL nci (aggregate by SPP/PLOT afterwards)
}
d$i_nci    <- calc_nci(d, dist_cutoff=NULL, size_cutoff=NULL)
d$i_ncilog <- log10(1 + d$i_nci) # log10 transform individual values

### plot the nci coefficients vs gradients
png('./../../fig/fig_03a_competitioncoeff.png', 6.5, 4.0, 'in',
    bg='transparent', res=700)
ecole::set_par_mercury(6, mgp=c(1.8,0.2,0), oma=c(0,1,0,0))
x1 <- expression(Fire ~ probability^{1/3})
x2 <- expression(log[10] ~ CMD ~ (mm ~ y^{-1}))
x3 <- expression(log[10] ~ MAP ~ (mm ~ y^{-1}))
yl <- 'Individual NCI'
i  <- d$i_nci < 4.0 # remove extreme values just for plotting
u  <- d[i,]
ecole::set_par_mercury(6)
plot(u$firelog, u$i_nci, col='#00000020', xlab=x1, ylab=yl)
plot(u$cmdlog,  u$i_nci, col='#00000020', xlab=x2, ylab=yl)
plot(u$maplog,  u$i_nci, col='#00000020', xlab=x3, ylab=yl)
plot(u$fire,    u$i_nci, col='#00000020', xlab='Fire prob', ylab=yl)
plot(u$cmd,     u$i_nci, col='#00000020', xlab='CMD', ylab=yl)
plot(10^u$maplog,u$i_nci,col='#00000020', xlab='MAP', ylab=yl)
rm(yl,u,calc_nci) # cleanup
dev.off()

png('./../../fig/fig_03b_logcompetitioncoeff.png', 6.5, 4.0, 'in',
    bg='transparent', res=700)
ecole::set_par_mercury(6, mgp=c(1.8,0.2,0), oma=c(0,1,0,0))
yl <- 'Indiv NCI (log10)'
i  <- d$i_ncilog < 0.7 # remove extreme values just for plotting
u  <- d[i,]
ecole::set_par_mercury(6)
plot(u$firelog, u$i_ncilog, col='#00000020', xlab=x1, ylab=yl)
plot(u$cmdlog,  u$i_ncilog, col='#00000020', xlab=x2, ylab=yl)
plot(u$maplog,  u$i_ncilog, col='#00000020', xlab=x3, ylab=yl)
plot(u$fire,    u$i_ncilog, col='#00000020', xlab='Fire prob',ylab=yl)
plot(u$cmd,     u$i_ncilog, col='#00000020', xlab='CMD', ylab=yl)
plot(10^u$maplog,u$i_ncilog,col='#00000020', xlab='MAP', ylab=yl)
rm(x1,x2,x3,yl,u) # cleanup
dev.off()
###################################################################
####    END nci coefficients    #####################
###################################################################


###################################################################
####     save     #################################################
###################################################################
pnw_tree <- d
save(pnw_tree, file='../../data/pnw_tree.rda')
###################################################################
###################################################################
###################################################################

### timing
t_end  <- Sys.time()
t_diff <- round(difftime(t_end, t_start, units="mins"),1)
cat(paste0('Data processing completed at:  ', t_end),
    '\n   after time elapsed of: ', (t_diff), 'minutes')

####    END    ####