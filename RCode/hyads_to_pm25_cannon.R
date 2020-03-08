rm( list = ls())

source( '/n/home03/lhenneman/repos/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

#======================================================================#
## Load meteorology as list of months
#======================================================================#
#define the layer names, do the actual downloading
Sys.setenv(TZ='UTC')
layer.names <- c( "air.2m.mon.mean.nc",
                  "apcp.mon.mean.nc",
                  "rhum.2m.mon.mean.nc",
                  "vwnd.10m.mon.mean.nc",
                  "uwnd.10m.mon.mean.nc")
names( layer.names) <- c( "temp", "apcp", "rhum", "vwnd", "uwnd")

# do the data downloading
# set destination parameter to where you want the data downloaded,
# for example, destination = '~/Desktop'
list.met <- lapply( layer.names,
                    downloader.fn, 
                    destination = '/n/zigler_lab/lhenneman/HyADS_to_pm25/inputdata/met',
                    dataset = 'NARR')

# take over US
mets2005   <- usa.functioner( 2005, list.met, dataset = 'NARR', return.usa.sub = F)
mets2006   <- usa.functioner( 2006, list.met, dataset = 'NARR', return.usa.sub = F)
mets2005.m <- usa.functioner( 2005, list.met, dataset = 'NARR', avg.period = 'month', return.usa.sub = F)
mets2006.m <- usa.functioner( 2006, list.met, dataset = 'NARR', avg.period = 'month', return.usa.sub = F)

# combine monthly rasters into single list
mets.m.all <- append( append( mets2005.m, mets2006.m))

#======================================================================#
## Load ddm as month
#======================================================================#
ddm2005.m <- ddm_to_zip( 
  ddm_coal_file = '/n/zigler_lab/lhenneman/HyADS_to_pm25/inputdata/CMAQ_DDM/COAL_impacts_2005_update.csv',
  Year = 2005, avg.period = 'month')
ddm2006.m <- ddm_to_zip( 
  ddm_coal_file = '/n/zigler_lab/lhenneman/HyADS_to_pm25/inputdata/CMAQ_DDM/COAL_impacts_2006_update.csv',
  Year = 2006, avg.period = 'month')

names( ddm2005.m) <- names( mets2005.m)
names( ddm2006.m) <- names( mets2006.m)

# combine into single list
ddm.m.all <- stack( ddm2005.m, ddm2006.m)

#======================================================================#
## Load ddm as annual
#======================================================================#
ddm2005 <- ddm_to_zip(
  ddm_coal_file = '/n/zigler_lab/lhenneman/HyADS_to_pm25/inputdata/CMAQ_DDM/COAL_impacts_2005_update.csv',
  Year = 2005)
ddm2006 <- ddm_to_zip( 
  ddm_coal_file = '/n/zigler_lab/lhenneman/HyADS_to_pm25/inputdata/CMAQ_DDM/COAL_impacts_2006_update.csv',
  Year = 2006)
names( ddm2005) <- 'cmaq.ddm'
names( ddm2006) <- 'cmaq.ddm'

#======================================================================#
## Load monthly hyads
#======================================================================#
# read monthly grid files
hyads.dir <- '/n/zigler_lab/lhenneman/diseperseR/main/output/exp'

hyads2005.m.l <- lapply( file.path( hyads.dir, paste0( 'grids_exposures_total_2005_',
                                                        formatC( 1:12, flag = '0', width = 2), '.fst')),
                          read.fst, as.data.table = T)
hyads2006.m.l <- lapply( file.path( hyads.dir, paste0( 'grids_exposures_total_2006_',
                                                        formatC( 1:12, flag = '0', width = 2), '.fst')),
                          read.fst, as.data.table = T)

# create lists from monthly grid objects
names( hyads2005.m.l) <- names( mets2005.m)
names( hyads2006.m.l) <- names( mets2006.m)

# create lists of monthly rasters
HyADSrasterizer <- function( X){
  r <- rasterFromXYZ( X[, .( x, y, hyads)], crs = p4s)
  r[is.na( r)] <- 0
  return( r)
}

hyads2005.m <- lapply( hyads2005.m.l, HyADSrasterizer)
hyads2006.m <- lapply( hyads2006.m.l, HyADSrasterizer)

# combine into single list
hyads.m.all <- stack( stack( hyads2005.m), stack( hyads2006.m))


#======================================================================#
## Load anuual hyads
#======================================================================#
hyads2005.dt <- na.omit( read.fst( file.path( hyads.dir, paste0( 'grids_exposures_total_2005.fst')), as.data.table = T))
hyads2006.dt <- na.omit( read.fst( file.path( hyads.dir, paste0( 'grids_exposures_total_2006.fst')), as.data.table = T))
hyads2005 <- rasterFromXYZ( hyads2005.dt[, .( x, y, hyads)], crs = p4s)
hyads2006 <- rasterFromXYZ( hyads2006.dt[, .( x, y, hyads)], crs = p4s)

## ========================================================= ##
##                Read in emissions data
## ========================================================= ##
d_cems_cmaq.f <- "/n/zigler_lab/lhenneman/HyADS_to_pm25/inputdata/NEI/2005_cemsum.txt"
d_nonegu.f <- "/n/zigler_lab/lhenneman/HyADS_to_pm25/inputdata/NEI/ptinv_ptnonipm_xportfrac_cap2005v2_2005cs_orl_06jan2011_v4_orl_COAL.txt"

d_cmaq <- fread( d_cems_cmaq.f)
d_nonegu <- fread( d_nonegu.f, skip = "06029", header = F)[,1:63]
d_nonegu.names <- unlist( fread( d_nonegu.f, skip = 'FIPS,PLANTID,', header = F, nrows = 1))
names( d_nonegu) <- d_nonegu.names
d_nonegu.slim <- d_nonegu[ POLCODE == 'SO2', .( XLOC, YLOC, ANN_EMIS)]

## Convert to spatial object, take over CMAQ raster
d_nonegu.sp <- SpatialPointsDataFrame( d_nonegu.slim[, .( XLOC, YLOC)], 
                                       data.frame( d_nonegu.slim[, ANN_EMIS]),
                                       proj4string = CRS( "+proj=longlat +datum=WGS84 +no_defs"))
d_nonegu.sp <- spTransform( d_nonegu.sp, CRS( p4s))
d_nonegu.r <- rasterize( d_nonegu.sp, ddm.m.all)$d_nonegu.slim...ANN_EMIS.
d_nonegu.r[is.na(d_nonegu.r[])] <- 0


#======================================================================#
# stack up and project annual data
#======================================================================#
# get usa mask for masking
# download USA polygon from rnaturalearth
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
us_states <- st_transform( USAboundaries::us_states(), p4s)
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]

dats2005.a <- project_and_stack( ddm2005, hyads2005,  
                                 mets2005, d_nonegu.r, mask.use = mask.usa)
dats2006.a <- project_and_stack( ddm2006, hyads2006,  
                                 mets2006, d_nonegu.r, mask.use = mask.usa)

#======================================================================#
## Combine into raster stack, train model
#======================================================================#
cov.names = c( "temp", "rhum", "vwnd", "uwnd", "wspd")

# predict each month in 2006 using model trained in 2005
preds.mon.hyads06w05 <- mapply( month.trainer, names( mets2005.m), names( mets2006.m),
                                MoreArgs = list( name.x = 'hyads', y.m = hyads.m.all,
                                                 ddm.m = ddm.m.all, mets.m = mets.m.all,
                                                 idwe.m = hyads.m.all, emiss.m = d_nonegu.r, 
                                                 .mask.use = mask.usa, cov.names = cov.names))

# predict annual 2006 using model trained in 2005
preds.ann.hyads06w05 <- lm.hyads.ddm.holdout( dat.stack = dats2005.a, dat.stack.pred = dats2006.a, 
                                              name.idwe = 'idwe', x.name = 'hyads',
                                              ho.frac = 0, covars.names = cov.names, return.mods = T)

#======================================================================#
## Save data
#======================================================================#
# annual stacks, 
# monthly stacks
# annual model
# monthly models
save( dats2005.a, dats2006.a,
      hyads.m.all, ddm.m.all, mets.m.all,
      d_nonegu.r,
      preds.mon.hyads06w05, #preds.mon.hyads05w06,
      preds.ann.hyads06w05, #preds.ann.hyads05w06,
      file = '/n/zigler_lab/lhenneman/HyADS_to_pm25/rdata/hyads_to_cmaq_models.RData')


# do correlation comparisons on quintiles
# scale all 3 on their Z score scale
#======================================================================#
## Annual plots
#======================================================================#
ggplot.a.raster( preds.ann.hyads06w05$Y.ho.hat.raster$y.hat.lm.cv,
                 preds.ann.hyads06w05$Y.ho.hat.raster$y.hat.gam.cv,
                 preds.ann.idwe06w05$Y.ho.hat.raster$y.hat.lm.cv,
                 preds.ann.idwe06w05$Y.ho.hat.raster$y.hat.gam.cv,
                 ncol. = 2, facet.names = c( 'lm - hyads', 'gam - hyads',
                                             'lm - idwe',  'gam - idwe'),
                 mask.raster = mask.usa)

