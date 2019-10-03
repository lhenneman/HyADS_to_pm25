source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')

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
                    downloader.fn, #destination = '~/Desktop'
                    dataset = 'NARR')

# take over US
mets2005.m <- usa.functioner( 2005, list.met, dataset = 'NARR', avg.period = 'month', return.usa.sub = F)
mets2006.m <- usa.functioner( 2006, list.met, dataset = 'NARR', avg.period = 'month', return.usa.sub = F)

# combine into single list
mets.m.all <- append( mets2005.m, mets2006.m)

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
## Load monthly hyads
#======================================================================#
# read monthly grid files
hyads2005.m.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2005grid/HyADS_grid_month_2005.csv', drop = 'V1')
hyads2006.m.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2006grid/HyADS_grid_month_2006.csv', drop = 'V1')

# create lists from monthly grid objects
hyads2005.m.l <- split( hyads2005.m.dt, by = 'yearmonth')
hyads2006.m.l <- split( hyads2005.m.dt, by = 'yearmonth')
names( hyads2005.m.l) <- names( mets2005.m)
names( hyads2006.m.l) <- names( mets2006.m)

# create lists of monthly rasters
hyads2005.m <- lapply( hyads2005.m.l, 
                       function( X){
                         r <- rasterFromXYZ( X[, .( x, y, hyads)], crs = p4s)
                         r[is.na( r)] <- 0
                         return( r)})
hyads2006.m <- lapply( hyads2006.m.l, 
                       function( X){
                         r <- rasterFromXYZ( X[, .( x, y, hyads)], crs = p4s)
                         r[is.na( r)] <- 0
                         return( r)})

# combine into single list
hyads.m.all <- stack( stack( hyads2005.m), stack( hyads2006.m))

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
d_nonegu.sp <- spTransform( d_nonegu.sp, CRS( p4s))
d_nonegu.r <- rasterize( d_nonegu.sp, ddm.m.all)$d_nonegu.slim...ANN_EMIS.
d_nonegu.r[is.na(d_nonegu.r[])] <- 0

#======================================================================#
## Combine into raster stack, train model
#======================================================================#
month.trainer <- function( name.m = names( mets2005.m)[1], 
                           name.p = names( mets2006.m)[1], 
                           ddm.m = ddm.m.all, 
                           hyads.m = hyads.m.all, 
                           mets.m = mets.m.all,
                           emiss.m = d_nonegu.r,
                           cov.names = c( "temp", "rhum", "vwnd", "uwnd", "wspd", names( d_nonegu.r))){
  # create training dataset
  ddm.use <- ddm.m[[name.m]]
  hyads.use <- projectRaster( hyads.m[[name.m]], ddm.use)
  mets.use <- projectRaster( mets.m[[name.m]], ddm.use)
  emiss.use <- projectRaster( emiss.m, ddm.use)
  
  # create prediction dataset
  ddm.use.p <- ddm.m[[name.p]]
  hyads.use.p <- projectRaster( hyads.m[[name.p]], ddm.use.p)
  mets.use.p <- projectRaster( mets.m[[name.p]], ddm.use.p)
  
  # fix names
  names( ddm.use)     <- 'cmaq.ddm'
  names( ddm.use.p)   <- 'cmaq.ddm'
  names( hyads.use)   <- 'hyads'
  names( hyads.use.p) <- 'hyads'
  
  # combine each dataset as stacks
  dat.s <- stack( ddm.use,   hyads.use,   mets.use,   emiss.use)
  dat.p <- stack( ddm.use.p, hyads.use.p, mets.use.p, emiss.use)
  
  # do the modeling
  pred <- lm.hyads.ddm.holdout( dat.stack = dat.s, dat.stack.pred = dat.p,
                                ho.frac = 0, covars.names = cov.names, return.mods = T)
  
  return( pred)
}


# predict each month in 2006 using model trained in 2005
preds.mon <- mapply( month.trainer, names( mets2005.m), names( mets2006.m))

#======================================================================#
## Extract data, summarize, and plot
#======================================================================#
## things we should show by month
# r (or R^2)
# spatial map of error by month
# each month's holdout?
# plot contributions of inputs

# get usa mask for masking
mask.usa <- usa.functioner( list.met = list.met, return.usa.mask = T)
mask.usa <- spTransform( mask.usa, CRSobj = crs( p4s))

plot( preds.mon['Y.ho.hat.raster','X2005.01.01'][[1]]$y.hat.lm.cv)

# plots of monthly predictions
ggplot.a.raster( preds.mon['Y.ho.hat.raster','X2005.01.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.raster','X2005.02.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.raster','X2005.03.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.raster','X2005.04.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.raster','X2005.05.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.raster','X2005.06.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.raster','X2005.07.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.raster','X2005.08.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.raster','X2005.09.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.raster','X2005.10.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.raster','X2005.11.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.raster','X2005.12.01'][[1]]$y.hat.lm.cv,
                 bounds = c( 0,6), ncol. = 3, facet.names = month.name,
                 mask.raster = mask.usa)

# plots of monthly error
ggplot.a.raster( preds.mon['Y.ho.hat.bias.raster','X2005.01.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.bias.raster','X2005.02.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.bias.raster','X2005.03.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.bias.raster','X2005.04.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.bias.raster','X2005.05.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.bias.raster','X2005.06.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.bias.raster','X2005.07.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.bias.raster','X2005.08.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.bias.raster','X2005.09.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.bias.raster','X2005.10.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.bias.raster','X2005.11.01'][[1]]$y.hat.lm.cv,
                 preds.mon['Y.ho.hat.bias.raster','X2005.12.01'][[1]]$y.hat.lm.cv,
                 bounds = c( -2,2), ncol. = 3, facet.names = month.name,
                 mask.raster = mask.usa)

# plots of monthly covariate contributions
ggplot.a.raster( preds.mon['Y.ho.terms.raster','X2005.01.01'][[1]],
                 bounds = c( -4,4), ncol. = 4, mask.raster = mask.usa,
                 facet.names = names( preds.mon['Y.ho.terms.raster','X2005.01.01'][[1]]))
ggplot.a.raster( preds.mon['Y.ho.terms.raster','X2005.07.01'][[1]],
                 bounds = c( -4,4), ncol. = 4, mask.raster = mask.usa,
                 facet.names = names( preds.mon['Y.ho.terms.raster','X2005.07.01'][[1]]))


## extract evaluation statistics
preds.metrics <- preds.mon[ 'metrics',]

metrics <- data.table( month = as.Date( gsub( '\\.', '-', gsub( 'X', '', names( preds.metrics)))),
                       'R^2' = sapply( preds.metrics, function( dt) dt[ mod.name == 'lm.cv']$R^2),
                       NMB = sapply( preds.metrics, function( dt) dt[ mod.name == 'lm.cv']$NMB),
                       NME = sapply( preds.metrics, function( dt) dt[ mod.name == 'lm.cv']$NME),
                       RMSE = sapply( preds.metrics, function( dt) dt[ mod.name == 'lm.cv']$RMSE))
metrics.m <- melt( metrics, id.vars = 'month', variable.name = 'metric')

ggplot( data = metrics.m,
        aes( x = month, y = value, group = metric)) + 
  geom_line() + geom_point() + 
  facet_wrap( . ~ metric, scales = 'free_y', ncol = 1, labeller = label_parsed) + 
  expand_limits( y = 0)

# extract linear model coefficients


