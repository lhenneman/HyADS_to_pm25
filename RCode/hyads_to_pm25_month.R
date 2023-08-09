rm( list = ls())

source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')

#coordinate reference system projection string for spatial data
p4s_cmaq  <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"
p4s_hyads <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

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
                    downloader.fn, #destination = '~/Desktop'
                    dataset = 'NARR')

# take over US
mets2005   <- usa.functioner( 2005, list.met, dataset = 'NARR', return.usa.sub = F)
mets2006   <- usa.functioner( 2006, list.met, dataset = 'NARR', return.usa.sub = F)
mets2011   <- usa.functioner( 2011, list.met, dataset = 'NARR', return.usa.sub = F)
mets2005.m <- usa.functioner( 2005, list.met, dataset = 'NARR', avg.period = 'month', return.usa.sub = F)
mets2006.m <- usa.functioner( 2006, list.met, dataset = 'NARR', avg.period = 'month', return.usa.sub = F)
mets2011.m <- usa.functioner( 2011, list.met, dataset = 'NARR', avg.period = 'month', return.usa.sub = F)

# combine monthly rasters into single list
mets.m.all <- append( append( mets2005.m, mets2006.m), mets2011.m)

#======================================================================#
## Load ddm as month
#======================================================================#
ddm2005.m <- ddm_to_zip( ddm_coal_file = '~//Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL_impacts_2005_update.csv',
                         Year = 2005, avg.period = 'month')
ddm2006.m <- ddm_to_zip( ddm_coal_file = '~//Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL_impacts_2006_update.csv',
                         Year = 2006, avg.period = 'month')

names( ddm2005.m) <- names( mets2005.m)
names( ddm2006.m) <- names( mets2006.m)

# combine into single list
ddm.m.all <- stack( ddm2005.m, ddm2006.m)

#======================================================================#
## Load ddm as annual
#======================================================================#
ddm2005 <- ddm_to_zip( ddm_coal_file = '~//Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL_impacts_2005_update.csv',
                       Year = 2005)
ddm2006 <- ddm_to_zip( ddm_coal_file = '~//Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL_impacts_2006_update.csv',
                       Year = 2006)
names( ddm2005) <- 'cmaq.ddm'
names( ddm2006) <- 'cmaq.ddm'

#======================================================================#
## Load monthly hyads
#======================================================================#
# read monthly grid files
hyads_exp_loc <- '~/Dropbox/Harvard/ARP/HyADS/hyads_longterm/exp/grids'

# list annual and monthly files
grid.files.m <- list.files( hyads_exp_loc,
                            pattern = 'grids_exposures_total_\\d{4}_\\d{2}\\.fst',
                            full.names = TRUE)
grid.files.yr <- list.files( hyads_exp_loc,
                             pattern = 'grids_exposures_total_\\d{4}\\.fst',
                             full.names = TRUE)

# read select files
grid.dat.yr <- lapply( grid.files.yr,
                       function( f){
                         year.f <- gsub( '^.*_|\\.fst', '', f)
                         if( !( year.f %in% c( '2005', '2006', '2011')))
                           return( data.table( ))
                         
                         in.f <- read.fst( f, as.data.table = T)
                         in.f[, `:=` (year = year.f,
                                      year.E = NULL,
                                      year.D = NULL)]
                       }) %>% rbindlist
grid.dat.m <- lapply( grid.files.m,
                       function( f){
                         year.m <- gsub( '^.*total_|\\.fst', '', f)
                         year.f <- gsub( '_.*', '', year.m)
                         if( !( year.f %in% c( '2005', '2006', '2011')))
                           return( data.table( ))
                         
                         in.f <- read.fst( f, as.data.table = T)
                       }) %>% rbindlist

hyads2005.m.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/gridexposures/HyADS_grid_month_nopbl2005.csv', drop = 'V1')
hyads2006.m.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/gridexposures/HyADS_grid_month_nopbl2006.csv', drop = 'V1')
hyads2011.m.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/gridexposures/HyADS_grid_month_nopbl2011.csv', drop = 'V1')

# create lists from monthly grid objects
hyads.m.l <- split( grid.dat.m, by = 'yearmonth')
hyads.y.l <- split( grid.dat.yr, by = 'year')
names( hyads.m.l) <- names( mets.m.all)

# create lists of monthly rasters
HyADSrasterizer <- function( X){
  r <- rasterFromXYZ( X[, .( x, y, hyads)], crs = p4s_hyads)
  r[is.na( r)] <- 0
  return( r)
}

# combine into single list
hyads.m.all <- lapply( hyads.m.l, HyADSrasterizer) %>% stack()
hyads.y.all <- lapply( hyads.y.l, HyADSrasterizer) %>% stack()

#======================================================================#
## Load anuual hyads
#======================================================================#
hyads2005 <- hyads.y.all$X2005
hyads2006 <- hyads.y.all$X2006
hyads2011 <- hyads.y.all$X2011

names( hyads2005) <- 'hyads'
names( hyads2006) <- 'hyads'
names( hyads2011) <- 'hyads'

## ========================================================= ##
##                Read in emissions data
## ========================================================= ##
d_cems_cmaq.f <- "~/Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL IMPACTS/INVENTORY/CEM/2005_cemsum.txt"
d_nonegu.f <- "~/Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL IMPACTS/INVENTORY/NONEGU COAL/ptinv_ptnonipm_xportfrac_cap2005v2_2005cs_orl_06jan2011_v4_orl_COAL.txt"

d_cmaq <- fread( d_cems_cmaq.f)
d_nonegu <- fread( d_nonegu.f, skip = "06029", header = F)[,1:63]
d_nonegu.names <- unlist( fread( d_nonegu.f, skip = 'FIPS,PLANTID,', header = F, nrows = 1))
names( d_nonegu) <- d_nonegu.names
d_nonegu.slim <- d_nonegu[ POLCODE == 'SO2', .( XLOC, YLOC, ANN_EMIS)]

## Convert to spatial object, take over CMAQ raster
d_nonegu.sp <- SpatialPointsDataFrame( d_nonegu.slim[, .( XLOC, YLOC)], 
                                       data.frame( d_nonegu.slim[, ANN_EMIS]),
                                       proj4string = CRS( "+proj=longlat +datum=WGS84 +no_defs"))
d_nonegu.sp <- spTransform( d_nonegu.sp, CRS( p4s_cmaq))
d_nonegu.r <- rasterize( d_nonegu.sp, ddm.m.all)$d_nonegu.slim...ANN_EMIS.
d_nonegu.r[is.na(d_nonegu.r[])] <- 0


## ========================================================= ##
##                Source inverse distance weighted raster
## ========================================================= ##
idwe.m.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/ampd_dists_sox_weighted.csv', drop = 'V1')
idwe.m.l <- split( idwe.m.dt, by = 'yearmon')

# create lists of monthly rasters
IDWErasterizer <- function( X){
  r <- rasterFromXYZ( X[, .( x, y, tot.sum)], crs = p4s_cmaq)
  r[is.na( r)] <- 0
  return( r)
}

idwe.m <- stack( lapply( idwe.m.l, IDWErasterizer))
names( idwe.m) <- names( hyads.m.all)

## ========================================================= ##
##                SOx inverse distance by year
## ========================================================= ##
idwe2005.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/ampd_dists_sox_weighted_2005_total.csv', drop = 'V1')
idwe2006.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/ampd_dists_sox_weighted_2006_total.csv', drop = 'V1')
idwe2011.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/ampd_dists_sox_weighted_2011_total.csv', drop = 'V1')
idwe2005 <- rasterFromXYZ( idwe2005.dt, crs = p4s_cmaq)
idwe2006 <- rasterFromXYZ( idwe2006.dt, crs = p4s_cmaq)
idwe2011 <- rasterFromXYZ( idwe2011.dt, crs = p4s_cmaq)
names( idwe2005) <- 'idwe'
names( idwe2006) <- 'idwe'
names( idwe2011) <- 'idwe'

summary(( hyads2006 - hyads2005) / hyads2005)
summary(( ddm2006 - ddm2005) / ddm2005)
summary(( idwe2006 - idwe2005) / idwe2005)

## ========================================================= ##
##                Plots
## ========================================================= ##
# get usa mask for masking
# download USA polygon from rnaturalearth
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
us_states <- st_transform( USAboundaries::us_states(), p4s_hyads)
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]

plot( ( hyads2006 - hyads2005) / hyads2005)
plot(mask.usa, add = T)
plot( (( ddm2006 - ddm2005) / ddm2005))
plot( (( idwe2006 - idwe2005) / idwe2005))


plot( hyads.m.all$X2005.07.01)
plot(mask.usa, add = T)

plot( idwe.m$X2005.07.01)
plot(mask.usa, add = T)

plot( ddm.m.all$X2005.06.01)
plot(mask.usa, add = T)

# plot( data.table( values( project_and_stack( hyads.m.all$X2005.08.01, ddm.m.all$X2005.08.01))))
# plot( data.table( values( project_and_stack( hyads.m.all$X2005.12.01, ddm.m.all$X2005.12.01))))
# plot( data.table( values( project_and_stack( idwe.m$X2005.12.01, ddm.m.all$X2005.12.01))))
#======================================================================#
# stack up and project annual data
#======================================================================#
dats2005.a <- project_and_stack( hyads2005, ddm2005, idwe2005, 
                                 mets2005, d_nonegu.r, mask.use = mask.usa)
dats2006.a <- project_and_stack( hyads2006, ddm2006, idwe2006, 
                                 mets2006, d_nonegu.r, mask.use = mask.usa)
dats2011.a <- project_and_stack( hyads2011, ddm2006, idwe2011, 
                                 mets2011, d_nonegu.r, mask.use = mask.usa)
dats2011.a$cmaq.ddm <- NA

summary( dats2006.a - dats2005.a)
summary( dats2011.a - dats2005.a)

cor( values( dats2005.a), use = 'complete.obs')
cor( values( dats2006.a), use = 'complete.obs')

dats2005.v <- data.table( values( dats2005.a))
dats2006.v <- data.table( values( dats2006.a))
plot( dats2005.v[, .(cmaq.ddm, hyads, idwe)])
plot( dats2006.v[, .(cmaq.ddm, hyads, idwe)])
plot( dats2005.a$cmaq.ddm < 1.2 & dats2005.a$hyads > 1.5e8)
plot( dats2006.a$cmaq.ddm < 1.2 & dats2006.a$hyads > 1.5e8)
d2005.red <- which( dats2005.v$cmaq.ddm < 1.2 & dats2005.v$hyads > 1.5e8)
plot( dats2005.v[d2005.red,.(cmaq.ddm, hyads, idwe)], col = 'red')

cor( dats2005.v[!d2005.red], use = 'complete.obs')

plot( dats2005.v[!d2005.red,.(cmaq.ddm, hyads, idwe)])

#======================================================================#
## Combine into raster stack, train model
#======================================================================#
cov.names = c( "temp", "rhum", "vwnd", "uwnd", "wspd")

# predict each month in 2006 using model trained in 2005
preds.mon.hyads06w05 <- mapply( month.trainer, names( mets2005.m), names( mets2006.m),
                                MoreArgs = list( name.x = 'hyads', y.m = hyads.m.all,
                                                 ddm.m = ddm.m.all, mets.m = mets.m.all,
                                                 idwe.m = idwe.m, emiss.m = d_nonegu.r, 
                                                 .mask.use = mask.usa, cov.names = cov.names))
preds.mon.idwe06w05  <- mapply( month.trainer, names( mets2005.m), names( mets2006.m),
                                MoreArgs = list( name.x = 'idwe', y.m = idwe.m,
                                                 ddm.m = ddm.m.all, mets.m = mets.m.all,
                                                 idwe.m = idwe.m, emiss.m = d_nonegu.r, 
                                                 .mask.use = mask.usa, cov.names = cov.names))
# predict each month in 2006 using model trained in 2005
# preds.mon.hyads05w06 <- mapply( month.trainer, names( mets2006.m), names( mets2005.m),
#                                 MoreArgs = list( name.x = 'hyads', y.m = hyads.m.all,
#                                                  ddm.m = ddm.m.all, mets.m = mets.m.all,
#                                                  idwe.m = idwe.m, emiss.m = d_nonegu.r, 
#                                                  .mask.use = mask.usa, cov.names = cov.names))
# preds.mon.idwe05w06  <- mapply( month.trainer, names( mets2006.m), names( mets2005.m),
#                                 MoreArgs = list( name.x = 'idwe', y.m = idwe.m,
#                                                  ddm.m = ddm.m.all, mets.m = mets.m.all,
#                                                  idwe.m = idwe.m, emiss.m = d_nonegu.r, 
#                                                  .mask.use = mask.usa, cov.names = cov.names))

# predict annual 2006 using model trained in 2005
preds.ann.hyads06w05 <- lm.hyads.ddm.holdout( dat.stack = dats2005.a, dat.stack.pred = dats2006.a, name.idwe = 'idwe', x.name = 'hyads',
                                              ho.frac = 0, covars.names = cov.names, return.mods = T)
preds.ann.idwe06w05  <- lm.hyads.ddm.holdout( dat.stack = dats2005.a, dat.stack.pred = dats2006.a, name.idwe = 'idwe', x.name = 'idwe',
                                              ho.frac = 0, covars.names = cov.names, return.mods = T)

# predict annual 2006 using model trained in 2005
# preds.ann.hyads05w06 <- lm.hyads.ddm.holdout( dat.stack = dats2006.a, dat.stack.pred = dats2005.a, name.idwe = 'idwe', x.name = 'hyads',
#                                               ho.frac = 0, covars.names = cov.names, return.mods = T)
# preds.ann.idwe05w06  <- lm.hyads.ddm.holdout( dat.stack = dats2006.a, dat.stack.pred = dats2005.a, name.idwe = 'idwe', x.name = 'idwe',
#                                               ho.frac = 0, covars.names = cov.names, return.mods = T)

# predict annual 2006 using model trained in 2005 - include inverse distance
# preds.ann.hyads06w05.i <- lm.hyads.ddm.holdout( dat.stack = dats2005.a, dat.stack.pred = dats2006.a, 
#                                               ho.frac = 0, covars.names = c( cov.names, 'idwe'), return.mods = T)

#======================================================================#
## Save data
#======================================================================#
# annual stacks, 
# monthly stacks
# annual model
# monthly models
save( dats2005.a, dats2006.a, dats2011.a,
      hyads.m.all, ddm.m.all, mets.m.all,
      idwe.m, d_nonegu.r,
      preds.mon.hyads06w05, #preds.mon.hyads05w06,
      preds.mon.idwe06w05,  #preds.mon.idwe05w06,
      preds.ann.hyads06w05, #preds.ann.hyads05w06,
      preds.ann.idwe06w05,  #preds.ann.idwe05w06,
      file = '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/hyads_to_cmaq_models3.RData')


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

# preds.ann.hyads06w05$metrics
# preds.ann.idwe06w05$metrics

#======================================================================#
## Extract data, summarize, and plot
#======================================================================#
## things we should show by month
# r (or R^2)
# spatial map of error by month
# each month's holdout?
# plot contributions of inputs


gg_out <- ggplot.a.raster( subset( ddm.m.all, 'X2005.07.01'),
                           preds.mon.hyads05w06 ['Y.ho.hat.raster','X2006.07.01'][[1]]$y.hat.gam.cv,
                           preds.mon.idwe05w06 ['Y.ho.hat.raster','X2006.07.01'][[1]]$y.hat.gam.cv,
                           mask.raster = mask.usa, facet.names = c( 'CMAQ', 'HyADS', 'IDWE'),
                           bounds = c( 0,8), ncol. = 1)

ggsave( '~/Dropbox/Harvard/Meetings_and_People/CMAS_2019/HyADS_pred_model_July.png', gg_out,
        height = 8, width = 3.5, scale = .7)

# plots of monthly predictions
ggplot.a.raster( preds.mon.hyads05w06['Y.ho.hat.raster','X2006.01.01'][[1]]$y.hat.gam.cv,
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.02.01'][[1]]$y.hat.gam.cv,
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.03.01'][[1]]$y.hat.gam.cv,
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.04.01'][[1]]$y.hat.gam.cv,
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.05.01'][[1]]$y.hat.gam.cv,
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.06.01'][[1]]$y.hat.gam.cv,
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.07.01'][[1]]$y.hat.gam.cv,
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.08.01'][[1]]$y.hat.gam.cv,
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.09.01'][[1]]$y.hat.gam.cv,
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.10.01'][[1]]$y.hat.gam.cv,
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.11.01'][[1]]$y.hat.gam.cv,
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.12.01'][[1]]$y.hat.gam.cv,
                 bounds = c( 0,8), ncol. = 3, facet.names = month.name,
                 mask.raster = mask.usa)

ggplot.a.raster( preds.mon.idwe['Y.ho.hat.raster','X2005.01.01'][[1]]$y.hat.lm.cv,
                 preds.mon.idwe['Y.ho.hat.raster','X2005.02.01'][[1]]$y.hat.lm.cv,
                 preds.mon.idwe['Y.ho.hat.raster','X2005.03.01'][[1]]$y.hat.lm.cv,
                 preds.mon.idwe['Y.ho.hat.raster','X2005.04.01'][[1]]$y.hat.lm.cv,
                 preds.mon.idwe['Y.ho.hat.raster','X2005.05.01'][[1]]$y.hat.lm.cv,
                 preds.mon.idwe['Y.ho.hat.raster','X2005.06.01'][[1]]$y.hat.lm.cv,
                 preds.mon.idwe['Y.ho.hat.raster','X2005.07.01'][[1]]$y.hat.lm.cv,
                 preds.mon.idwe['Y.ho.hat.raster','X2005.08.01'][[1]]$y.hat.lm.cv,
                 preds.mon.idwe['Y.ho.hat.raster','X2005.09.01'][[1]]$y.hat.lm.cv,
                 preds.mon.idwe['Y.ho.hat.raster','X2005.10.01'][[1]]$y.hat.lm.cv,
                 preds.mon.idwe['Y.ho.hat.raster','X2005.11.01'][[1]]$y.hat.lm.cv,
                 preds.mon.idwe['Y.ho.hat.raster','X2005.12.01'][[1]]$y.hat.lm.cv,
                 bounds = c( 0,6), ncol. = 3, facet.names = month.name,
                 mask.raster = mask.usa)

# plots of monthly error
ggplot.a.raster( preds.mon.hyads05w06['Y.ho.hat.raster','X2006.01.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.01.01'),
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.02.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.02.01'),
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.03.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.03.01'),
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.04.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.04.01'),
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.05.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.05.01'),
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.06.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.06.01'),
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.07.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.07.01'),
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.08.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.08.01'),
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.09.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.09.01'),
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.10.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.10.01'),
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.11.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.11.01'),
                 preds.mon.hyads05w06['Y.ho.hat.raster','X2006.12.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.12.01'),
                 bounds = c( -2,2), ncol. = 3, facet.names = month.name,
                 mask.raster = mask.usa)

ggplot.a.raster( preds.mon.idwe05w06['Y.ho.hat.raster','X2006.01.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.01.01'),
                 preds.mon.idwe05w06['Y.ho.hat.raster','X2006.02.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.02.01'),
                 preds.mon.idwe05w06['Y.ho.hat.raster','X2006.03.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.03.01'),
                 preds.mon.idwe05w06['Y.ho.hat.raster','X2006.04.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.04.01'),
                 preds.mon.idwe05w06['Y.ho.hat.raster','X2006.05.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.05.01'),
                 preds.mon.idwe05w06['Y.ho.hat.raster','X2006.06.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.06.01'),
                 preds.mon.idwe05w06['Y.ho.hat.raster','X2006.07.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.07.01'),
                 preds.mon.idwe05w06['Y.ho.hat.raster','X2006.08.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.08.01'),
                 preds.mon.idwe05w06['Y.ho.hat.raster','X2006.09.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.09.01'),
                 preds.mon.idwe05w06['Y.ho.hat.raster','X2006.10.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.10.01'),
                 preds.mon.idwe05w06['Y.ho.hat.raster','X2006.11.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.11.01'),
                 preds.mon.idwe05w06['Y.ho.hat.raster','X2006.12.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.12.01'),
                 bounds = c( -2,2), ncol. = 3, facet.names = month.name,
                 mask.raster = mask.usa)

ggplot.a.raster( preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.01.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.01.01'),
                 preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.02.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.02.01'),
                 preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.03.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.03.01'),
                 preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.04.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.04.01'),
                 preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.05.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.05.01'),
                 preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.06.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.06.01'),
                 preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.07.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.07.01'),
                 preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.08.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.08.01'),
                 preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.09.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.09.01'),
                 preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.10.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.10.01'),
                 preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.11.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.11.01'),
                 preds.mon.hyads05w06['Y.ho.hat.bias.raster','X2006.12.01'][[1]]$y.hat.gam.cv / subset( ddm.m.all, 'X2005.12.01'),
                 bounds = c( -1,1), ncol. = 3, facet.names = month.name,
                 mask.raster = mask.usa)

#======================================================================#
## Check out the covariates
#======================================================================#
# plots of monthly covariate contributions
# average of each covariate over the year
# spatial plots of just HyADS/spline covariates over the year
covs.all <- names( preds.ann.hyads06w05$Y.ho.terms.gam.raster)
covs.all.idwe <- names( preds.ann.idwe06w05$Y.ho.terms.gam.raster)
covs.hyads <- covs.all[grep( 'hyads', covs.all)]
covs.hyads.s <- c( covs.hyads, "s.x.y.")
covs.tot.sum  <- covs.all.idwe[grep( 'tot.sum', covs.all.idwe)]
covs.tot.sum.s  <- c( covs.tot.sum, "s.x.y.")
covs.idwe  <- gsub( 'tot.sum', 'idwe', covs.tot.sum)
covs.idwe.s  <- gsub( 'tot.sum', 'idwe', covs.tot.sum.s)

# plot hyads related covariates for different years
ggplot.a.raster( unstack( stack( subset( preds.ann.hyads06w05$Y.ho.terms.gam.raster, covs.hyads.s),
                                 subset( preds.ann.hyads05w06$Y.ho.terms.gam.raster, covs.hyads.s))),
                 bounds = c( -4,4), ncol. = 7, mask.raster = mask.usa,
                 facet.names = paste( covs.hyads.s, rep( c( '06w05', '05w06'), each = 7)))
ggplot.a.raster( unstack( stack( subset( preds.ann.idwe06w05$Y.ho.terms.gam.raster, covs.tot.sum.s),
                                 subset( preds.ann.idwe05w06$Y.ho.terms.gam.raster, covs.tot.sum.s))),
                 bounds = c( -4,4), ncol. = 7, mask.raster = mask.usa,
                 facet.names = paste( covs.tot.sum.s, rep( c( '06w05', '05w06'), each = 7)))

# sum all hyads/tot.sum contributions
hyads_gamters <- stack( sum( subset( preds.ann.hyads06w05$Y.ho.terms.gam.raster, covs.hyads)),
                        subset( preds.ann.hyads06w05$Y.ho.terms.gam.raster, 's.x.y.'),
                        sum( subset( preds.ann.hyads05w06$Y.ho.terms.gam.raster, covs.hyads)),
                        subset( preds.ann.hyads05w06$Y.ho.terms.gam.raster, 's.x.y.'),
                        sum( subset( preds.ann.idwe06w05$Y.ho.terms.gam.raster, covs.tot.sum)),
                        subset( preds.ann.idwe06w05$Y.ho.terms.gam.raster, 's.x.y.'),
                        sum( subset( preds.ann.idwe05w06$Y.ho.terms.gam.raster, covs.tot.sum)),
                        subset( preds.ann.idwe05w06$Y.ho.terms.gam.raster, 's.x.y.'))
ggplot.a.raster( unstack( hyads_gamters),
                 bounds = c( -4,4), ncol. = 2, mask.raster = mask.usa,
                 facet.names = c( paste( c( 'hyads', 'hyads s.x.y.'), 
                                         rep( c( '06w05', '05w06'), each = 2)),
                                  paste( c( 'idwe', 'idwe s.x.y.'), 
                                         rep( c( '06w05', '05w06'), each = 2))))

# plot hyads related covariates for different months
gamters.mon06w05 <- stack( lapply( colnames( preds.mon.hyads06w05), 
                                   function( x) {
                                     subset( preds.mon.hyads06w05['Y.ho.terms.gam.raster', x][[1]], covs.hyads)
                                   }))
gamters.mon06w05.i <- stack( lapply( colnames( preds.mon.idwe06w05), 
                                     function( x) {
                                       subset( preds.mon.idwe06w05['Y.ho.terms.gam.raster', x][[1]], covs.idwe)
                                     }))
gamters.mon05w06 <- stack( lapply( colnames( preds.mon.hyads05w06), 
                                   function( x) {
                                     subset( preds.mon.hyads05w06['Y.ho.terms.gam.raster', x][[1]], covs.hyads)
                                   }))
names.gamters <- paste( covs.hyads, rep( month.abb, each = 7))
names.gamters.i <- paste( covs.idwe, rep( month.abb, each = 7))
ggplot.a.raster( unstack( gamters.mon06w05),
                 bounds = c( -4,4), ncol. = 7, mask.raster = mask.usa,
                 facet.names = names.gamters)
ggplot.a.raster( unstack( gamters.mon06w05.i),
                 bounds = c( -4,4), ncol. = 7, mask.raster = mask.usa,
                 facet.names = names.gamters.i)
ggplot.a.raster( unstack( gamters.mon05w06),
                 bounds = c( -4,4), ncol. = 7, mask.raster = mask.usa,
                 facet.names = names.gamters)

# sum all hyads/tot.sum contributions
gamters.mon06w05.hyadssum <- stack( lapply( colnames( preds.mon.hyads06w05), 
                                            function( x) {
                                              stack( sum( subset( preds.mon.hyads06w05['Y.ho.terms.gam.raster', x][[1]], covs.hyads)),
                                                     subset( preds.mon.hyads06w05['Y.ho.terms.gam.raster', x][[1]], 's.x.y.'))
                                            }))
gamters.mon06w05.idwesum <- stack( lapply( colnames( preds.mon.idwe06w05), 
                                           function( x) {
                                             stack( sum( subset( preds.mon.idwe06w05['Y.ho.terms.gam.raster', x][[1]], covs.idwe)),
                                                    subset( preds.mon.idwe06w05['Y.ho.terms.gam.raster', x][[1]], 's.x.y.'))
                                           }))
names.gamters.hy <- paste( c( 'hyads', 's.x.y.'), rep( month.abb, each = 2))
names.gamters.is <- paste( c( 'idwe', 's.x.y.'), rep( month.abb, each = 2))
ggplot.a.raster( unstack( gamters.mon06w05.hyadssum),
                 bounds = c( -4,4), ncol. = 4, mask.raster = mask.usa,
                 facet.names = names.gamters.hy)
ggplot.a.raster( unstack( gamters.mon06w05.idwesum),
                 bounds = c( -4,4), ncol. = 4, mask.raster = mask.usa,
                 facet.names = names.gamters.is)


ggplot.a.raster( preds.ann.hyads06w05$Y.ho.terms.raster,
                 bounds = c( -4,4), ncol. = 5, mask.raster = mask.usa,
                 facet.names = names( preds.ann.hyads06w05$Y.ho.terms.raster))
ggplot.a.raster( preds.ann.hyads06w05$Y.ho.terms.gam.raster,
                 bounds = c( -4,4), ncol. = 5, mask.raster = mask.usa,
                 facet.names = names( preds.ann.hyads06w05$Y.ho.terms.gam.raster))
ggplot.a.raster( preds.mon.hyads06w05['Y.ho.terms.raster','X2005.01.01'][[1]],
                 bounds = c( -4,4), ncol. = 4, mask.raster = mask.usa,
                 facet.names = names( preds.mon.hyads06w05['Y.ho.terms.raster','X2005.01.01'][[1]]))
ggplot.a.raster( preds.mon.hyads06w05['Y.ho.terms.raster','X2005.07.01'][[1]],
                 bounds = c( -4,4), ncol. = 4, mask.raster = mask.usa,
                 facet.names = names( preds.mon.hyads06w05['Y.ho.terms.raster','X2005.07.01'][[1]]))
ggplot.a.raster( preds.mon.hyads06w05['Y.ho.terms.gam.raster','X2005.01.01'][[1]],
                 bounds = c( -4,4), ncol. = 4, mask.raster = mask.usa,
                 facet.names = names( preds.mon.hyads06w05['Y.ho.terms.gam.raster','X2005.01.01'][[1]]))
ggplot.a.raster( preds.mon.hyads05w06['Y.ho.terms.gam.raster','X2006.07.01'][[1]],
                 bounds = c( -4,4), ncol. = 4, mask.raster = mask.usa,
                 facet.names = names( preds.mon.hyads05w06['Y.ho.terms.gam.raster','X2006.07.01'][[1]]))
ggplot.a.raster( preds.mon.idwe05w06['Y.ho.terms.gam.raster','X2006.07.01'][[1]],
                 bounds = c( -4,4), ncol. = 4, mask.raster = mask.usa,
                 facet.names = names( preds.mon.idwe05w06['Y.ho.terms.gam.raster','X2006.07.01'][[1]]))


#======================================================================#
## Plot the metrics
#======================================================================#
## extract evaluation statistics
## IDWE gets big change from bivariate spline, HyADS does not  
preds.metrics.hyads <- preds.mon.hyads06w05[ 'metrics',]
preds.metrics.idwe  <- preds.mon.idwe06w05[ 'metrics',]

metrics <- data.table( month = c( as.Date( gsub( '\\.', '-', gsub( 'X', '', names( preds.metrics.hyads)))),
                                  as.Date( gsub( '\\.', '-', gsub( 'X', '', names( preds.metrics.idwe))))),
                       model = c( rep( 'hyads', length( names( preds.metrics.hyads))),
                                  rep( 'idwe', length( names( preds.metrics.idwe)))),
                       class = c( rep( 'gam', 2 * length( names( preds.metrics.hyads))),
                                  rep( 'lm', 2 * length( names( preds.metrics.idwe)))),
                       'R^2' = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$R^2),
                                  sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$R^2),
                                  sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$R^2),
                                  sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$R^2)),
                       NMB = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$NMB),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$NMB),
                                sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$NMB),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$NMB)),
                       NME = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$NME),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$NME),
                                sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$NME),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$NME)),
                       RMSE = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$RMSE),
                                 sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$RMSE),
                                 sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$RMSE),
                                 sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$RMSE)))
metrics.m <- melt( metrics, id.vars = c( 'model', 'month', 'class'), variable.name = 'metric')

ggplot( data = metrics.m,
        aes( x = month, y = value, lty = class, color = model)) + 
  geom_line() + geom_point() + 
  facet_wrap( . ~ metric, scales = 'free_y', ncol = 1, labeller = label_parsed) + 
  expand_limits( y = 0)

# metrics - adj. Z score, no model
metrics.Z.only <- data.table( month = c( as.Date( gsub( '\\.', '-', gsub( 'X', '', names( preds.metrics.hyads)))),
                                         as.Date( gsub( '\\.', '-', gsub( 'X', '', names( preds.metrics.idwe))))),
                              model = c( rep( 'hyads', length( names( preds.metrics.hyads))),
                                         rep( 'idwe', length( names( preds.metrics.idwe)))),
                              'R^2' = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'adj.Z.only']$R^2),
                                         sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'adj.Z.only']$R^2)),
                              NMB = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'adj.Z.only']$NMB),
                                       sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'adj.Z.only']$NMB)),
                              NME = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'adj.Z.only']$NME),
                                       sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'adj.Z.only']$NME)),
                              RMSE = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'adj.Z.only']$RMSE),
                                        sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'adj.Z.only']$RMSE)))
metrics.Z.only.m <- melt( metrics.Z.only, id.vars = c( 'model', 'month'), variable.name = 'metric')

ggplot( data = metrics.Z.only.m,
        aes( x = month, y = value, group = model, color = model)) + 
  geom_line() + geom_point() + 
  facet_wrap( . ~ metric, scales = 'free_y', ncol = 1, labeller = label_parsed) + 
  expand_limits( y = 0)

# extract linear model coefficients

#annual comparisons
#5 day avg time
#Check w/ sunni on month/annual etc
#======================================================================#
## Plot changes in evaluation in different areas
#======================================================================#

cors.keep.month.hyads.u05w06 <- rbindlist( preds.mon.hyads05w06['evals.q',], idcol = 'month')[, y := '05w06']
cors.keep.month.hyads.u06w05 <- rbindlist( preds.mon.hyads06w05['evals.q',], idcol = 'month')[, y := '06w05']
cors.keep.month.idwe.u05w06  <- rbindlist( preds.mon.idwe05w06['evals.q',], idcol = 'month')[, y := '05w06']
cors.keep.month.idwe.u06w05  <- rbindlist( preds.mon.idwe06w05['evals.q',], idcol = 'month')[, y := '06w05']
cors.keep.month <- rbind( cors.keep.month.hyads.u05w06, cors.keep.month.hyads.u06w05,
                          cors.keep.month.idwe.u05w06, cors.keep.month.idwe.u06w05)
cors.keep.m <- melt( cors.keep.month, id.vars = c( 'mod.name', 's', 'month', 'y'))
cors.keep.m[, month := month( as.Date( gsub( 'X', '', month), format = '%Y.%m.%d'))]
ggplot( data = cors.keep.m,
        aes( x = s, y = value, color = mod.name, lty = y)) + 
  geom_hline( yintercept = 0) +
  facet_grid( variable ~ month, scales = 'free_y') +  geom_line() 


# plot annual evaluation across s
cors.keep.u06w05 <- rbind( preds.ann.hyads06w05$evals.q, preds.ann.idwe06w05$evals.q)[, y := '06w05']
cors.keep.u05w06 <- rbind( preds.ann.hyads05w06$evals.q, preds.ann.idwe05w06$evals.q)[, y := '05w06']
cors.keep.u <- rbind( cors.keep.u06w05, cors.keep.u05w06)
cors.keep.m <- melt( cors.keep.u, id.vars = c( 'mod.name', 's', 'y'))
ggplot( data = cors.keep.m,
        aes( x = s, y = value, color = mod.name, lty = y)) + 
  geom_hline( yintercept = 0) +
  facet_wrap( . ~ variable, scales = 'free_y') +  geom_line() 


# need somehow to evaluate near vs far sources
# approximate this as high/low tot.sum
# says more about how emissions near sources are handled than
# anything else
# check out wind speed argument --- very key
#  IDWE does better in years with slow windspeed?
# plot cmaq range at each s
# do MSE?
cors.keep <- data.table()
for (y in 2005:2006){
  vals <- values( get( paste0( 'dats', y, '.a')))
  for ( s in seq( 0.01, 1, .01)){
    q <- quantile( vals[,'tot.sum'], s, na.rm = T)
    cors <- cor( vals[vals[,'cmaq.ddm'] < q,], use = 'complete.obs', method = 'spearman') 
    cors.keep <- rbind( cors.keep,
                        data.table( s = s, hyads = cors['cmaq.ddm', 'hyads'],
                                    idwe = cors['cmaq.ddm', 'tot.sum'], year = y))
  }
}
cors.keep.m <- melt( cors.keep, id.vars = c( 's', 'year'))
ggplot( data = cors.keep.m,
        aes( x = s, y = value, color = variable, group = variable)) + 
  geom_line() + facet_wrap( year ~ ., ncol = 2)

cors.keep[which.min( abs( hyads - idwe))]
