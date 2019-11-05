rm( list = ls())

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
hyads2005.m.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2005grid/HyADS_grid_month_2005.csv', drop = 'V1')
hyads2006.m.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2006grid/HyADS_grid_month_2006.csv', drop = 'V1')
hyads2011.m.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2011grid/HyADS_grid_month_2011.csv', drop = 'V1')

# create lists from monthly grid objects
hyads2005.m.l <- split( hyads2005.m.dt, by = 'yearmonth')
hyads2006.m.l <- split( hyads2005.m.dt, by = 'yearmonth')
hyads2011.m.l <- split( hyads2011.m.dt, by = 'yearmonth')
names( hyads2005.m.l) <- names( mets2005.m)
names( hyads2006.m.l) <- names( mets2006.m)
names( hyads2011.m.l) <- names( mets2011.m)

# create lists of monthly rasters
HyADSrasterizer <- function( X){
  r <- rasterFromXYZ( X[, .( x, y, hyads)], crs = p4s)
  r[is.na( r)] <- 0
  return( r)
}

hyads2005.m <- lapply( hyads2005.m.l, HyADSrasterizer)
hyads2006.m <- lapply( hyads2006.m.l, HyADSrasterizer)
hyads2011.m <- lapply( hyads2011.m.l, HyADSrasterizer)

# combine into single list
hyads.m.all <- stack( stack( hyads2005.m), stack( hyads2006.m), stack( hyads2011.m))


#======================================================================#
## Load anuual hyads
#======================================================================#
hyads2005.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2005grid/HyADS_grid_2005.csv', drop = 'V1')
hyads2006.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2006grid/HyADS_grid_2006.csv', drop = 'V1')
hyads2011.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2011grid/HyADS_grid_2011.csv', drop = 'V1')
hyads2005 <- rasterFromXYZ( hyads2005.dt[, .( x, y, hyads)], crs = p4s)
hyads2006 <- rasterFromXYZ( hyads2006.dt[, .( x, y, hyads)], crs = p4s)
hyads2011 <- rasterFromXYZ( hyads2011.dt[, .( x, y, hyads)], crs = p4s)


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


## ========================================================= ##
##                Source inverse distance weighted raster
## ========================================================= ##
idwe.m.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/ampd_dists_sox_weighted.csv', drop = 'V1')
idwe.m.l <- split( idwe.m.dt, by = 'yearmon')

# create lists of monthly rasters
IDWErasterizer <- function( X){
  r <- rasterFromXYZ( X[, .( x, y, tot.sum)], crs = p4s)
  r[is.na( r)] <- 0
  return( r)
}

idwe.m <- stack( lapply( idwe.m.l, IDWErasterizer))
names( idwe.m) <- names( hyads.m.all)

## ========================================================= ##
##                SOx inverse distance by year
## ========================================================= ##
idwe.m.dt[, dates := as.Date( paste0( gsub( '_.*', '', yearmon), '_',
                                      formatC( as.numeric( gsub( '.*_', '', yearmon)), width = 2, flag = '0'),
                                      '_01'), 
                              format = '%Y_%m_%d')]
idwe.m.dt[, days := days_in_month( dates)]
idwe.dt <- idwe.m.dt[, .( tot.sum = sum( tot.sum * days) / 365), by = .( x, y, year( dates))]
idwe.a.l <- split( idwe.dt, by = 'year')

idwe.a <- lapply( idwe.a.l, IDWErasterizer)
names( idwe.a) <- unique( idwe.dt$year)

## ========================================================= ##
##                Plots
## ========================================================= ##
# get usa mask for masking
mask.usa <- usa.functioner( list.met = list.met, return.usa.mask = T)
mask.usa <- spTransform( mask.usa, CRSobj = crs( p4s))

plot( hyads.m.all$X2005.07.01)
plot(mask.usa, add = T)

plot( idwe.m$X2005.07.01)
plot(mask.usa, add = T)

plot( ddm.m.all$X2005.07.01)
plot(mask.usa, add = T)
#======================================================================#
# stack up and project annual data
#======================================================================#
dats2005.a <- project_and_stack( ddm2005, hyads2005, idwe.a$`2005`, 
                                 mets2005, d_nonegu.r, mask.use = mask.usa)
dats2006.a <- project_and_stack( ddm2006, hyads2006, idwe.a$`2006`, 
                                 mets2006, d_nonegu.r, mask.use = mask.usa)


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
preds.mon.hyads05w06 <- mapply( month.trainer, names( mets2006.m), names( mets2005.m),
                                MoreArgs = list( name.x = 'hyads', y.m = hyads.m.all,
                                                 ddm.m = ddm.m.all, mets.m = mets.m.all,
                                                 idwe.m = idwe.m, emiss.m = d_nonegu.r, 
                                                 .mask.use = mask.usa, cov.names = cov.names))
preds.mon.idwe05w06  <- mapply( month.trainer, names( mets2006.m), names( mets2005.m),
                                MoreArgs = list( name.x = 'idwe', y.m = idwe.m,
                                                 ddm.m = ddm.m.all, mets.m = mets.m.all,
                                                 idwe.m = idwe.m, emiss.m = d_nonegu.r, 
                                                 .mask.use = mask.usa, cov.names = cov.names))

# predict annual 2006 using model trained in 2005
preds.ann.hyads06w05 <- lm.hyads.ddm.holdout( dat.stack = dats2005.a, dat.stack.pred = dats2006.a, 
                                              ho.frac = 0, covars.names = cov.names, return.mods = T)
preds.ann.idwe06w05  <- lm.hyads.ddm.holdout( dat.stack = dats2005.a, dat.stack.pred = dats2006.a, x.name = 'tot.sum',
                                              ho.frac = 0, covars.names = cov.names, return.mods = T)

# predict annual 2006 using model trained in 2005
preds.ann.hyads05w06 <- lm.hyads.ddm.holdout( dat.stack = dats2006.a, dat.stack.pred = dats2005.a,
                                              ho.frac = 0, covars.names = cov.names, return.mods = T)
preds.ann.idwe05w06  <- lm.hyads.ddm.holdout( dat.stack = dats2006.a, dat.stack.pred = dats2005.a, x.name = 'tot.sum',
                                              ho.frac = 0, covars.names = cov.names, return.mods = T)

#======================================================================#
## Save data
#======================================================================#
# annual stacks, 
# monthly stacks
# annual model
# monthly models
save.list <- c( dats2005.a, )

# do correlation comparisons on quintiles
# scale all 3 on their Z score scale
#======================================================================#
## Extract data, summarize, and plot
#======================================================================#
## things we should show by month
# r (or R^2)
# spatial map of error by month
# each month's holdout?
# plot contributions of inputs


gg_out <- ggplot.a.raster( subset( ddm.m.all, 'X2005.07.01'),
                           preds.mon.hyads05w06 ['Y.ho.hat.raster','X2006.07.01'][[1]]$y.hat.lm.cv,
                           preds.mon.idwe05w06 ['Y.ho.hat.raster','X2006.07.01'][[1]]$y.hat.lm.cv,
                           mask.raster = mask.usa, facet.names = c( 'CMAQ', 'HyADS', 'IDWE'),
                           bounds = c( 0,8), ncol. = 1)

ggsave( '~/Dropbox/Harvard/Meetings_and_People/CMAS_2019/HyADS_pred_model_July.png', gg_out,
        height = 8, width = 3.5, scale = .7)

# plots of monthly predictions
ggplot.a.raster( preds.mon.hyads['Y.ho.hat.raster','X2005.01.01'][[1]]$y.hat.lm.cv,
                 preds.mon.hyads['Y.ho.hat.raster','X2005.02.01'][[1]]$y.hat.lm.cv,
                 preds.mon.hyads['Y.ho.hat.raster','X2005.03.01'][[1]]$y.hat.lm.cv,
                 preds.mon.hyads['Y.ho.hat.raster','X2005.04.01'][[1]]$y.hat.lm.cv,
                 preds.mon.hyads['Y.ho.hat.raster','X2005.05.01'][[1]]$y.hat.lm.cv,
                 preds.mon.hyads['Y.ho.hat.raster','X2005.06.01'][[1]]$y.hat.lm.cv,
                 preds.mon.hyads['Y.ho.hat.raster','X2005.07.01'][[1]]$y.hat.lm.cv,
                 preds.mon.hyads['Y.ho.hat.raster','X2005.08.01'][[1]]$y.hat.lm.cv,
                 preds.mon.hyads['Y.ho.hat.raster','X2005.09.01'][[1]]$y.hat.lm.cv,
                 preds.mon.hyads['Y.ho.hat.raster','X2005.10.01'][[1]]$y.hat.lm.cv,
                 preds.mon.hyads['Y.ho.hat.raster','X2005.11.01'][[1]]$y.hat.lm.cv,
                 preds.mon.hyads['Y.ho.hat.raster','X2005.12.01'][[1]]$y.hat.lm.cv,
                 bounds = c( 0,6), ncol. = 3, facet.names = month.name,
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
ggplot.a.raster( preds.mon.hyads['Y.ho.hat.bias.raster','X2005.01.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.01.01'),
                 preds.mon.hyads['Y.ho.hat.bias.raster','X2005.02.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.02.01'),
                 preds.mon.hyads['Y.ho.hat.bias.raster','X2005.03.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.03.01'),
                 preds.mon.hyads['Y.ho.hat.bias.raster','X2005.04.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.04.01'),
                 preds.mon.hyads['Y.ho.hat.bias.raster','X2005.05.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.05.01'),
                 preds.mon.hyads['Y.ho.hat.bias.raster','X2005.06.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.06.01'),
                 preds.mon.hyads['Y.ho.hat.bias.raster','X2005.07.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.07.01'),
                 preds.mon.hyads['Y.ho.hat.bias.raster','X2005.08.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.08.01'),
                 preds.mon.hyads['Y.ho.hat.bias.raster','X2005.09.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.09.01'),
                 preds.mon.hyads['Y.ho.hat.bias.raster','X2005.10.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.10.01'),
                 preds.mon.hyads['Y.ho.hat.bias.raster','X2005.11.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.11.01'),
                 preds.mon.hyads['Y.ho.hat.bias.raster','X2005.12.01'][[1]]$y.hat.lm.cv / subset( ddm.m.all, 'X2005.12.01'),
                 bounds = c( -1,1), ncol. = 3, facet.names = month.name,
                 mask.raster = mask.usa)

# plots of monthly covariate contributions
ggplot.a.raster( preds.mon['Y.ho.terms.raster','X2005.01.01'][[1]],
                 bounds = c( -4,4), ncol. = 4, mask.raster = mask.usa,
                 facet.names = names( preds.mon['Y.ho.terms.raster','X2005.01.01'][[1]]))
ggplot.a.raster( preds.mon['Y.ho.terms.raster','X2005.07.01'][[1]],
                 bounds = c( -4,4), ncol. = 4, mask.raster = mask.usa,
                 facet.names = names( preds.mon['Y.ho.terms.raster','X2005.07.01'][[1]]))


#======================================================================#
## Plot the metrics
#======================================================================#
## extract evaluation statistics
preds.metrics.hyads <- preds.mon.hyads06w05[ 'metrics',]
preds.metrics.idwe  <- preds.mon.idwe06w05[ 'metrics',]

metrics <- data.table( month = c( as.Date( gsub( '\\.', '-', gsub( 'X', '', names( preds.metrics.hyads)))),
                                  as.Date( gsub( '\\.', '-', gsub( 'X', '', names( preds.metrics.idwe))))),
                       model = c( rep( 'hyads', length( names( preds.metrics.hyads))),
                                  rep( 'idwe', length( names( preds.metrics.idwe)))),
                       'R^2' = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$R^2),
                                  sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$R^2)),
                       NMB = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$NMB),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$NMB)),
                       NME = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$NME),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$NME)),
                       RMSE = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$RMSE),
                                 sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$RMSE)))
metrics.m <- melt( metrics, id.vars = c( 'model', 'month'), variable.name = 'metric')

ggplot( data = metrics.m,
        aes( x = month, y = value, group = model, color = model)) + 
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

cors.keep.month.hyads <- rbindlist( preds.mon.hyads05w06['evals.qq',], idcol = 'month')
cors.keep.month.idwe  <- rbindlist( preds.mon.idwe05w06['evals.qq',], idcol = 'month')
cors.keep.month <- rbind( cors.keep.month.hyads, cors.keep.month.idwe)
cors.keep.m <- melt( cors.keep.month, id.vars = c( 'mod.name', 's', 'month'))
ggplot( data = cors.keep.m,
        aes( x = s, y = value, color = mod.name, group = mod.name)) + 
  geom_hline( yintercept = 0) +
  facet_grid( variable ~ month, scales = 'free') +  geom_line() 


# plot annual evaluation across s
cors.keep.u <- rbind( preds.ann.hyads05w06$evals.q, preds.ann.idwe05w06$evals.q)
cors.keep.m <- melt( cors.keep.u, id.vars = c( 'mod.name', 's'))
ggplot( data = cors.keep.m,
        aes( x = s, y = value, color = mod.name, group = mod.name)) + 
  geom_hline( yintercept = 0) +
  facet_wrap( . ~ variable) +  geom_line() 


