library( data.table)
library( raster)
library( spBayes)
library( hyspdisp)
library( ggplot2)
library( viridis)

#======================================================================#
## define functions
#======================================================================#
`%ni%` <- Negate(`%in%`) 

## ddm
ddm_to_zip <- function( ddm_coal_file,
                        Year,
                        avg.period = 'year'){
  
  p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"
  
  #read data, 
  ddm_coal <- fread(ddm_coal_file)
  
  # melt, extract year, month, day
  ddm_coal.m <- melt( ddm_coal, id.vars = c( 'X', 'Y'),
                      variable.name = 'date.in', value.name = 'coal_pm25')
  ddm_coal.m[, `:=` ( date.in = as.Date( date.in, format = '%m/%d/%y'))]
  ddm_coal.m[, `:=` ( year.in = year( date.in),
                      month.in = month( date.in))]
  
  # remove blow up values (greater than 30)
  ddm_coal.m[ coal_pm25 < 0 | coal_pm25 > 30, coal_pm25 := NA]
  
  #rasterize as brick
  names <- unique( ddm_coal.m$date.in)
  ddm_coal.b <- brick( lapply( names,
                               function( name, dt.m){
                                 r <- rasterFromXYZ(dt.m[ date.in == name,
                                                          .(x = X, y = Y, z = coal_pm25)],
                                                    crs = CRS(p4s))
                                 names( r) <- name
                                 return( r)
                               }, ddm_coal.m))
  
  # fill NA's with linear interpolation across days
  ddm_coal.b <- approxNA( ddm_coal.b, rule=2)
  
  
  
  # take monthly averages
  if( avg.period == 'month'){
    ddm_coal.month <- lapply( 1:12,
                              function( m, ddm_raster.b){
                                M <- formatC( m, width = 2, flag = '0')
                                id <- grep( paste0( '\\.', M, '\\.'), names( ddm_coal.b))
                                
                                ddm_coal.mon <- mean( subset( ddm_raster.b, id))
                                # ddm_coal.mon.z <- raster_to_zip(ddm_raster = ddm_coal.mon,
                                #                                 zcta_shapefile, crosswalk_csv)
                                # ddm_coal.mon.z[, `:=` (month = m, year = year)]
                                return( ddm_coal.mon)
                              }, ddm_coal.b)
    
    ddm_coal.z <- brick( ddm_coal.month)
    names( ddm_coal.z) <- paste( Year, formatC( 1:12, width = 2, flag = '0'), sep = '.')
    
    # take annual averages
  } else if( avg.period == 'year'){
    # take annual average
    ddm_coal.z <- mean( ddm_coal.b)
    names( ddm_coal.z) <- Year
    
  }
  
  return( ddm_coal.z)
}


## functions to get meteorology data
# download the necessary met files, 20th century reanalysis
downloader.fn <- function( filename,
                           dataset = c( '20thC_ReanV2c', 'ncep.reanalysis.derived', 'NARR')){
  if( length( dataset) > 1)
    dataset <- dataset[1]
  fileloc <- file.path('~', 'Dropbox', 'Harvard', 'RFMeval_Local', 'Comparisons_Intermodel', 
                       'Global_meteorology', dataset)
  
  
  # create directory to store in
  dir.create( fileloc, 
              recursive = T, 
              showWarnings = F)
  
  # name variable, filenames
  varname_NOAA <- gsub( "\\..*", "", filename)
  file_NOAA <- file.path( fileloc, filename)
  
  # define URL
  if( dataset == '20thC_ReanV2c'){
    # https://www.esrl.noaa.gov/psd/data/gridded/data.20thC_ReanV2c.monolevel.mm.html
    url_NOAA <- paste0( "ftp://ftp.cdc.noaa.gov/Datasets/20thC_ReanV2c/Monthlies/gaussian/monolevel/", filename)
  } else if( dataset == 'ncep.reanalysis.derived'){
    url_NOAA <- paste0( "ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/", filename)
  }else if( dataset == 'NARR'){
    url_NOAA <- paste0( "ftp://ftp.cdc.noaa.gov/Datasets/NARR/Monthlies/monolevel/", filename)
  }
  
  if( !file.exists( file_NOAA))
    download.file( url = url_NOAA,
                   destfile = file_NOAA)
  
  hpbl_rasterin <- brick( x = file_NOAA, 
                          varname = varname_NOAA)
  
  return( hpbl_rasterin)
  
}

# download the necessary met files, 20th century reanalysis
extract_year.fn <- function( raster.in = list.met[[1]],
                             year.in = 2005,
                             dataset = c( '20thC_ReanV2c', 'ncep.reanalysis.derived', 'NARR')){
  
  # default to 20th cent reanalysis
  if( length( dataset) > 1){
    dataset <- dataset[1]
    print( paste( 'No dataset specified, defaulting to', dataset))
  }
  
  # name months 1:12 for extracting from raster
  names.months <- paste0( year.in, '-',
                          formatC( 1:12, width = 2, flag = '0'), '-',
                          '01')
  
  # extract monthly dates using function from hyspdisp
  raster.sub <- subset_nc_date(  hpbl_brick = raster.in,
                                 vardate = names.months)
  
  # take annual mean
  raster.sub.mean <- stackApply( raster.sub, indices = rep( 1, 12), fun = mean)
  
  #NARR dataset requires rotating
  if( dataset != 'NARR')
    raster.sub.mean <- rotate( raster.sub.mean)
  
  return( raster.sub.mean)
}

# trim data over US, combine into data.table, create spatial object
usa.functioner <- function( year.in = 2005,
                            list.met,
                            dataset = c( '20thC_ReanV2c', 'ncep.reanalysis.derived', 'NARR'),
                            return.usa.mask = F,
                            return.usa.sub = F){
  
  # extract year
  mets <- brick( lapply( list.met,
                         extract_year.fn,
                         year.in = year.in,
                         dataset = dataset))
  crs.usa <- crs( "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000")
  
  # convert temp to celcius
  mets$temp <- mets$temp - 273.15
  
  # calculate windspeed
  # calculate meteorology wind angle (0 is north wind)
  # http://weatherclasses.com/uploads/3/6/2/3/36231461/computing_wind_direction_and_speed_from_u_and_v.pdf
  if( 'uwnd' %in% names( list.met) & 'vwnd' %in% names( list.met)){
    mets$wspd <- sqrt( mets$uwnd ^ 2 + mets$vwnd ^ 2)
    mets$phi <- atan2( mets$uwnd, mets$vwnd) * 180 / pi + 180
  }
  
  # download USA polygon from rnaturalearth
  usa <- rnaturalearth::ne_countries(scale = 110, type = "countries", country = "United States of America", 
                                     geounit = NULL, sovereignty = NULL,
                                     returnclass = c("sp"))
  usa.sub <- disaggregate(usa)[6,]
  usa.sub <- spTransform(usa.sub, CRSobj = proj4string( mets))
  
  if( return.usa.mask){
    usa.sub.p <- spTransform(usa.sub, CRSobj = crs.usa)
    usa.sub.sf <- data.table( st_as_sf( usa.sub.p))
    return( usa.sub.sf)
  }
  
  if( return.usa.sub){
    # crop to USA
    mets.usa <- crop( mask(mets, usa.sub), usa.sub)
    
    # reproject
    # mets.usa <- projectRaster( mets.usa, crs = crs.usa)
    
    # convert rasters to sf - differences
    # mets.usa.sp <- rasterToPolygons( mets.usa)
    
    # convert rasters to sf - annual
    # mets.usa.sf <- data.table( st_as_sf( mets.usa.sp))[, year := year.in]
    
    # merge with coordinates - annual
    # coords <- st_coordinates( st_centroid( mets.usa.sf$geometry))
    # mets.usa.sf <- cbind( mets.usa.sf, coords)
    
    # return sf object
    return( mets.usa)
  }
  
  else{
    mets <- projectRaster( mets, crs = crs.usa)
    return( mets)
  }
}

# define the model function
splM.hyads.ddm <- function( dummy.n = 1,
                            coords = as.matrix( hyads2005.dt[,.( x, y)]),
                            Y = ddm2005.dt$X2005,
                            X = data.table( intercept = 1, hyads = hyads2005.dt$X2005),
                            seed.n = NULL,
                            ...){
  
  set.seed( seed.n)
  
  quants <- function(x){
    quantile(x, prob=c(0.5, 0.025, 0.975))
  }
  
  # holdout fraction
  ho.frac <- .1
  
  # define number of samples
  n.samples <- 5000
  
  # define priors
  starting <- list("tau.sq"=1, "sigma.sq"=1, "phi"=6)
  tuning <- list("tau.sq"=0.01, "sigma.sq"=0.01, "phi"=0.1)
  priors <- list("beta.Flat", "tau.sq.IG"=c(2, 1),
                 "sigma.sq.IG"=c(2, 1), "phi.Unif"=c(3, 30))
  
  # define holdout parameters
  ho <- sample( 1:length( Y), ceiling( ho.frac * length( Y)))
  
  #convert X & Y to matrices
  X.m <- as.matrix( X)
  Y.m <- as.matrix( Y)
  
  # define inputs
  coords.ho <- coords[ho,]
  coords.tr <- coords[-ho,]
  Y.ho <- Y.m[ho]
  Y.tr <- Y.m[-ho]
  X.ho <- X.m[ho,]
  X.tr <- X.m[-ho,]
  
  # define burn in
  burn.in <- floor(0.75*n.samples)
  
  # train the model
  m.i <- spLM( Y.tr ~ X.tr - 1, coords = coords.tr, 
               modified.pp = TRUE, ...,
               starting = starting, tuning = tuning, priors = priors, 
               cov.model = "exponential",
               n.samples = n.samples, n.report = 2500)
  
  # recover estimates of beta and theta
  rt.rec <- system.time( {
    m.i.rec <- spRecover(m.i, start=burn.in, thin=5, n.report=100)
  })
  
  # return simulated y.hat's
  rt.y.hat <- system.time( {
    m.i.pred <- spPredict( m.i, start=burn.in, thin=2, pred.covars = X.ho, 
                           pred.coords=coords.ho, verbose=FALSE)
  })
  
  # find estimates of beta, theta, w.hat, and y.hat
  beta.hat <- round(summary(m.i.rec$p.beta.recover.samples)$quantiles[c(3,1,5)],6)
  theta.hat <- round( summary( window( m.i$p.theta.samples, 
                                       start = burn.in))$quantiles[, c( 3, 1, 5)], 2)
  w.hat <- apply(m.i.rec$p.w.recover.samples, 1, median)
  y.hat <- data.table( t( apply(m.i.pred$p.y.predictive.samples, 1, quants)))
  
  # set up evaluation data.table
  Y.ho.tr.dt <- data.table( cbind( coords.ho, y.hat[,`50%`], Y.ho))
  setnames( Y.ho.tr.dt, 'V3', 'Y.hat')
  
  # rasterize output for plots
  w.hat.raster <- rasterFromXYZ( cbind( coords.tr, w.hat))
  y.hat.raster <- rasterFromXYZ( Y.ho.tr.dt[, .( x, y, Y.hat)])
  Y.tr.raster <- rasterFromXYZ( cbind( coords.tr, Y.tr))
  Y.ho.raster <- rasterFromXYZ( Y.ho.tr.dt[, .( x, y, Y.ho)])
  
  # mcmc plots
  # plot( m.i$p.theta.samples)
  # plot( m.i.rec$p.beta.recover.samples)
  
  # spatial areas of y.hat - Y.ho
  par(mfrow=c(1,3))
  plot( w.hat.raster, main = "w.hat (spatial adjustment term)")
  points( m.i$knot.coords, cex=1)
  plot( y.hat.raster - Y.ho.raster, main = "Y.hat - Y.ho")
  plot( (y.hat.raster - Y.ho.raster) / Y.ho.raster, main = "(Y.hat - Y.ho) / Y.ho")
  par(mfrow=c(1,1))
  
  # check out simpler models - set up inputs
  names.covars <- names( X)
  XY.tr <- data.table( Y = Y.tr, X.tr)
  XY.ho <- data.table( Y = Y.ho, X.ho)
  setnames( XY.tr, names(XY.tr)[names(XY.tr) %ni% 'Y'], names.covars)
  setnames( XY.ho, names(XY.tr)[names(XY.tr) %ni% 'Y'], names.covars)
  
  # check out simpler models - define them
  form <- as.formula( paste( 'Y ~ -1 +', paste( names.covars, collapse = '+')))
  m.lm <- lm( form, data = XY.tr)
  m.mean <- mean( Y.tr / XY.tr$hyads)
  
  # check out simpler models - get predicted Y.ho
  y.hat.lm <- predict( m.lm, newdata = XY.ho)
  y.hat.mean <- XY.ho$hyads * m.mean
  
  # calculate evaluation metrics
  eval.fn <- function( Yhat, Yact, mod.name){
    num.diff <- sum( Yhat - Yact)
    abs.diff <- sum( abs( Yhat - Yact))
    denom <- sum( Yact)
    metrics <- data.table( mod.name = mod.name,
                           NMB = num.diff / denom,
                           NME = abs.diff / denom,
                           MB   = num.diff / length( ho),
                           RMSE = sqrt( sum( ( Yhat - Yact) ^ 2) / length( Y.ho)),
                           R = cor( Yhat, Yact))
    return( metrics)
  }
  metrics.out <- rbind( eval.fn( Y.ho.tr.dt$Y.hat, Y.ho, 'spLM'),
                        eval.fn( y.hat.lm, Y.ho, 'lm'),
                        eval.fn( y.hat.mean, Y.ho, 'mean'))
  
  out <- list( metrics = metrics.out,
               runtimes = data.table( ncells.tr = length( Y.tr),
                                      train = round(m.i$run.time[3]/60,3),
                                      recover = round(rt.rec[3]/60,3),
                                      est.y.hat = round(rt.y.hat[3]/60,3)))
  
  return( out)
}

# define the linear model holdout function
lm.hyads.ddm.holdout <- function( seed.n = NULL,
                                  dat.stack = dats2005.s.small,
                                  dat.stack.pred = NULL,
                                  y.name = 'cmaq.ddm',
                                  x.name = 'hyads',
                                  covars.names = NULL, #c( 'temp', 'apcp'),
                                  ho.frac = .1,
                                  ...){
  if( is.null( covars.names))
    covars.names <- names( dat.stack)[names( dat.stack) %ni% c( x.name, y.name)]
  
  set.seed( seed.n)
  
  # define holdout parameters
  N <- ncell( dat.stack)
  ho <- sample( 1:N, ceiling( ho.frac * N))
  
  # extract coordinates
  dat.coords <- coordinates( dat.stack)
  dat.coords.ho <- data.table( dat.coords[ho,])
  dat.coords.tr <- data.table( dat.coords[-ho,])
  
  # define inputs
  dat.stack.ho <- data.table( values( dat.stack)[ ho,])
  dat.stack.tr <- data.table( values( dat.stack)[-ho,])
  
  # special case for ho is zero
  if( length( ho) == 0){
    # extract coordinates
    dat.coords.ho <- data.table( dat.coords)
    dat.coords.tr <- data.table( dat.coords)
    
    # define inputs
    dat.stack.ho <- data.table( values( dat.stack))
    dat.stack.tr <- data.table( values( dat.stack))
  } 
  if( !is.null( dat.stack.pred)){
    # extract coordinates
    dat.coords.ho <- data.table( coordinates( dat.stack.pred))
    
    # define inputs
    dat.stack.ho <- data.table( values( dat.stack.pred))
  } 
  
  # check out linear regression models - define them
  form.cv <-  as.formula( paste( y.name, '~', paste( c( x.name, covars.names), collapse = '+')))
  form.ncv <- as.formula( paste( y.name, '~', x.name))
  lm.cv <-  lm( form.cv,  data = dat.stack.tr)
  lm.ncv <- lm( form.ncv, data = dat.stack.tr)
  
  # check out linear regression models - define them
  mean.y.over.x <- mean( unlist( dat.stack.tr[,..y.name]) / unlist( dat.stack.tr[,..x.name]), na.rm = T)
  mean.y <- mean( unlist( dat.stack.tr[,..y.name]), na.rm = T)
  mean.x <- mean( unlist( dat.stack.tr[,..x.name]), na.rm = T)
  sd.y <- sd( unlist( dat.stack.tr[,..y.name]), na.rm = T)
  sd.x <- sd( unlist( dat.stack.tr[,..x.name]), na.rm = T)
  
  # check out simpler models - get predicted Y.ho
  y.ho <- unlist( dat.stack.ho[,..y.name])
  y.hat.lm.cv  <- predict( lm.cv,  newdata = dat.stack.ho, se.fit = T)
  y.hat.lm.ncv <- predict( lm.ncv, newdata = dat.stack.ho, se.fit = T)
  y.hat.mean <- unlist( dat.stack.ho[,..x.name] * mean.y.over.x)
  y.hat.Z <-    unlist( (dat.stack.ho[,..x.name] - mean.x) / sd.x * sd.y + mean.y)
  
  # set up evaluation data.table
  Y.ho.hat <- data.table( dat.coords.ho, y.ho, y.hat.lm.cv = y.hat.lm.cv$fit, 
                          y.hat.lm.ncv = y.hat.lm.ncv$fit, y.hat.mean, y.hat.Z)
  Y.ho.hat.bias <- data.table( dat.coords.ho, y.ho, 
                               y.hat.lm.cv = y.hat.lm.cv$fit - y.ho, 
                               y.hat.lm.ncv = y.hat.lm.ncv$fit - y.ho, 
                               y.hat.mean = y.hat.mean - y.ho, 
                               y.hat.Z = y.hat.Z - y.ho)
  Y.ho.hat.se <- data.table( dat.coords.ho, y.hat.lm.cv = y.hat.lm.cv$se.fit, 
                             y.hat.lm.ncv = y.hat.lm.ncv$se.fit)
  
  # rasterize output for plots
  crs.in <- crs( dat.stack)
  Y.ho.hat.raster <- projectRaster( rasterFromXYZ( Y.ho.hat, crs = crs.in), dat.stack)
  Y.ho.hat.se.raster <- projectRaster( rasterFromXYZ( Y.ho.hat.se, crs = crs.in), dat.stack)
  Y.ho.hat.bias.raster <- projectRaster( rasterFromXYZ( Y.ho.hat.bias, crs = crs.in), dat.stack)
  
  # calculate evaluation metrics
  eval.fn <- function( Yhat, Yact, mod.name){
    num.diff <- sum( Yhat - Yact, na.rm = T)
    abs.diff <- sum( abs( Yhat - Yact), na.rm = T)
    denom <- sum( Yact, na.rm = T)
    metrics <- data.table( mod.name = mod.name,
                           NMB = num.diff / denom,
                           NME = abs.diff / denom,
                           MB   = num.diff / length( Yhat),
                           RMSE = sqrt( sum( ( Yhat - Yact) ^ 2, na.rm = T) / length( Yhat)),
                           R = cor( Yhat, Yact, use = 'complete.obs'))
    return( metrics)
  }
  metrics.out <- rbind( eval.fn( Y.ho.hat$y.hat.lm.cv, y.ho, 'lm.cv'),
                        eval.fn( Y.ho.hat$y.hat.lm.ncv, y.ho, 'lm.ncv'),
                        eval.fn( Y.ho.hat$y.hat.mean, y.ho, 'adj.mean'),
                        eval.fn( Y.ho.hat$y.hat.Z, y.ho, 'adj.Z'))
  
  out <- list( metrics = metrics.out, 
               Y.ho.hat.raster = Y.ho.hat.raster, 
               Y.ho.hat.se.raster = Y.ho.hat.se.raster,
               Y.ho.hat.bias.raster = Y.ho.hat.bias.raster)
  
  return( out)
}

# extract a mean from a list of lists of rasters
mean.lol <- function( lol, layer1, layer2, plot.out = TRUE){
  first.lol <- lapply( lol, '[[', layer1)
  second.lol <- lapply( first.lol, '[[', layer2)
  mean.lol <- mean( stack( second.lol), na.rm = T)
  
  if( plot.out)
    plot( mean.lol, main = layer2)
  return( mean.lol)
}

# ggplot a raster
ggplot.a.raster <- function( ..., bounds = NULL, facet.names = NULL){
  
  in.x <- list( ...)
  if( !is.null( facet.names))
    names( in.x) <- facet.names
  
  dat.dt <- rbindlist( lapply( names( in.x), function( x.name, x.list) {
    x <- x.list[[x.name]]
    r_points <- rasterToPoints( x)
    r_dt <- data.table( r_points)[, name.in := x.name]
    setnames( r_dt, names( x), 'z')
    return( r_dt)
  }, in.x))
  
  ggplot( dat.dt) +
    geom_tile( aes( x = x, y = y, fill = z)) + 
    scale_fill_viridis( name = NULL, limits = bounds, oob = scales::squish) + 
    facet_wrap( . ~ name.in) + 
    # expand_limits( fill = 0) +
    theme_bw() + 
    theme( axis.text = element_blank(),
           axis.title = element_blank(),
           axis.ticks = element_blank(),
           panel.grid = element_blank(),
           strip.background = element_blank())
  
}





#======================================================================#
## import data
#======================================================================#
ddm2005 <- ddm_to_zip( ddm_coal_file = '~//Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL_impacts_2005_update.csv',
                       Year = 2005)
ddm2006 <- ddm_to_zip( ddm_coal_file = '~//Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL_impacts_2006_update.csv',
                       Year = 2006)
names( ddm2005) <- 'cmaq.ddm'
names( ddm2006) <- 'cmaq.ddm'

p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"
hyads2005.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2005grid/HyADS_grid_2005.csv', drop = 'V1')
hyads2006.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2006grid/HyADS_grid_2006.csv', drop = 'V1')
hyads2005 <- rasterFromXYZ( hyads2005.dt[, .( x, y, hyads)], crs = p4s)
hyads2006 <- rasterFromXYZ( hyads2006.dt[, .( x, y, hyads)], crs = p4s)

#======================================================================#
## download met data
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
list.met <- lapply( layer.names,
                    downloader.fn,
                    dataset = 'NARR')

# take over US
mets2005 <- usa.functioner( 2005, list.met, dataset = 'NARR', return.usa.sub = F)
mets2006 <- usa.functioner( 2006, list.met, dataset = 'NARR', return.usa.sub = F)

#======================================================================#
## reproject, align on same grid
#======================================================================#
hyads2005.p <- projectRaster( hyads2005, ddm2005)
mets2005.p <- projectRaster( mets2005, ddm2005)
hyads2006.p <- projectRaster( hyads2006, ddm2006)
mets2006.p <- projectRaster( mets2006, ddm2006)

# fill NA's with linear interpolation across days
hyads2005.pa <- hyads2005.p #focal( hyads2005.p, w=matrix(1,1,1), fun = mean, na.rm = T, NAonly = T)
hyads2006.pa <- hyads2006.p #focal( hyads2005.p, w=matrix(1,1,1), fun = mean, na.rm = T, NAonly = T)
names( hyads2005.pa) <- names( hyads2005.p)
names( hyads2006.pa) <- names( hyads2006.p)

#======================================================================#
## stack up!
#======================================================================#
dats2005.s <- stack( ddm2005, hyads2005.pa, mets2005.p)
dats2006.s <- stack( ddm2006, hyads2006.pa, mets2006.p)

# crop to US
usa <- rnaturalearth::ne_countries(scale = 110, type = "countries", country = "United States of America", 
                                   geounit = NULL, sovereignty = NULL,
                                   returnclass = c("sp"))
usa.sub <- disaggregate(usa)[6,]
usa.sub <- spTransform(usa.sub, CRSobj = proj4string( dats2005.s))

dats2005.s <- crop( mask(dats2005.s, usa.sub), usa.sub)
dats2006.s <- crop( mask(dats2006.s, usa.sub), usa.sub)
#======================================================================#
## pick out section for analysis
#======================================================================#
square <- extent( 0e6, 1e6, -5e5, 5e5)
square.sm <- extent( 0.25e6, 0.75e6, -2.5e5, 2.5e5)
square.lg <- extent( -.5e6, 1.5e6, -10e5, 10e5)

#======================================================================#
## crop to test area
#======================================================================#
dats2005.s.square    <- crop( dats2005.s, square)
dats2005.s.small     <- crop( dats2005.s, square.sm)
dats2005.s.large     <- crop( dats2005.s, square.lg)

#======================================================================#
## run linear model holdouts
#======================================================================#
# single.crossval.sm <- lapply( 1:50,
#                               lm.hyads.ddm.holdout,
#                               dat.stack = dats2005.s.small)
# single.crossval.la <- lapply( 1:50,
#                               lm.hyads.ddm.holdout,
#                               dat.stack = dats2005.s.large)
single.crossval <- lapply( 1:500,
                           lm.hyads.ddm.holdout,
                           dat.stack = dats2005.s,
                           covars.names = )

pred2006 <- lm.hyads.ddm.holdout( dat.stack = dats2005.s,
                                  ho.frac = 0,
                                  dat.stack.pred = dats2006.s)

save.image( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_PM25/RData/cmaqddm_to_pm25_crossval.RData')


mean.lol( single.crossval.sm, 'Y.ho.hat.bias.raster', 'y.hat.lm.cv')
mean.lol( single.crossval.sm, 'Y.ho.hat.bias.raster', 'y.hat.lm.ncv')
mean.lol( single.crossval.la, 'Y.ho.hat.bias.raster', 'y.hat.lm.cv')
mean.lol( single.crossval.la, 'Y.ho.hat.bias.raster', 'y.hat.lm.ncv')
mean.lol( single.crossval, 'Y.ho.hat.bias.raster', 'y.hat.lm.cv')
mean.lol( single.crossval, 'Y.ho.hat.bias.raster', 'y.hat.lm.ncv')


single.crossval2 <- lm.hyads.ddm.holdout( dat.stack = dats2005.s,
                                          ho.frac = 0)
plot( single.crossval2$Y.ho.hat.bias.raster$y.hat.lm.cv)
plot( single.crossval2$Y.ho.hat.bias.raster$y.hat.Z)
plot( single.crossval2$Y.ho.hat.bias.raster$y.hat.mean)

cuts=0:6
par(mfrow=c(1,2))
plot( single.crossval2$Y.ho.hat.raster$y.hat.lm.cv, main = 'Y.hat.2006')
plot( dats2005.s$cmaq.ddm, main = 'Y.2006')
par(mfrow=c(1,1))


# to do - 
# 1. limit model to over land
# 2. Build distributions of metrics for different models

# to save - input data and cross validation
#======================================================================#
## extract distributions of metrics
#======================================================================#

metrics.out <- rbindlist( lapply( single.crossval, '[[', 'metrics'))
metrics.out.m <- melt( metrics.out, id.vars = 'mod.name')
metrics.means <- metrics.out.m[, .( med = median( value)), by = .( mod.name, variable)]

ggplot( data = metrics.out.m[ mod.name %ni% 'adj.mean'],
        aes( color = mod.name)) +
  facet_wrap( . ~ variable, ncol = 1, scales = 'free') +
  geom_vline( data = metrics.means[ mod.name %ni% 'adj.mean'], 
              aes( xintercept = med, color = mod.name),
              lty = 2) +
  geom_density( aes( x = value)) + 
  scale_color_discrete( name = 'Model Name') 
#scale_x_continuous( limits = c( .23, .37))

# compare 2005 and 2006

ggplot.a.raster( pred2006$Y.ho.hat.raster$y.hat.lm.cv,
                 dats2006.s$cmaq.ddm, bounds = c( 0,5),
                 facet.names = c( 'Yhat', 'Hybrid CMAQ-DDM'))


#======================================================================#
#======================================================================#
## old analysis - compare spLM and lm
#======================================================================#
#======================================================================#
us_buonds <- sf::as_Spatial( USAboundaries::us_states( ))
us_buonds <- spTransform( us_buonds, CRS( p4s))
par(mfrow=c(1,2))
plot( ddm2005, main = 'DDM, 2005')
plot( square, add = T)
plot( square.sm, add = T)
plot( square.lg, add = T)
plot( us_buonds, add = T, border = 'grey50')
plot( hyads2005.p, main = 'HyADS, 2005')
plot( square, add = T)
plot( square.sm, add = T)
plot( square.lg, add = T)
plot( us_buonds, add = T, border = 'grey50')
par(mfrow=c(1,1))

#======================================================================#
## crop to test area
#======================================================================#
hyads2005.pc    <- crop( hyads2005.p, square)
hyads2005.pc.sm <- crop( hyads2005.p, square.sm)
hyads2005.pc.lg <- crop( hyads2005.p, square.lg)

ddm2005.pc      <- crop( ddm2005, square)
ddm2005.pc.sm   <- crop( ddm2005, square.sm)
ddm2005.pc.lg   <- crop( ddm2005, square.lg)

mets2005.pc     <- crop( mets2005.p, square)
mets2005.pc.sm  <- crop( mets2005.p, square.sm)
mets2005.pc.lg  <- crop( mets2005.p, square.lg)

#======================================================================#
## extract coordinates, values
#======================================================================#
hyads2005.dt    <- data.table( rasterToPoints( hyads2005.pc))
hyads2005.sm.dt <- data.table( rasterToPoints( hyads2005.pc.sm))
hyads2005.lg.dt <- data.table( rasterToPoints( hyads2005.pc.lg))

ddm2005.dt      <- data.table( rasterToPoints( ddm2005.pc))
ddm2005.sm.dt   <- data.table( rasterToPoints( ddm2005.pc.sm))
ddm2005.lg.dt   <- data.table( rasterToPoints( ddm2005.pc.lg))

mets2005.dt     <- data.table( rasterToPoints( mets2005.pc))
mets2005.sm.dt  <- data.table( rasterToPoints( mets2005.pc.sm))
mets2005.lg.dt  <- data.table( rasterToPoints( mets2005.pc.lg))

#======================================================================#
## set up Y & X covariates
#======================================================================#
set.inputs <- function( hyads.dt, ddm.dt, mets.dt = NULL){
  list( Y = ddm.dt$X2005,
        X = data.table( intercept = 1, 
                        hyads = hyads.dt$X2005, 
                        mets.dt),
        coords = as.matrix( hyads.dt[,.( x, y)]))
}

in.cv <- set.inputs( hyads2005.dt, ddm2005.dt, mets2005.dt)
in.sm.cv <- set.inputs( hyads2005.sm.dt, ddm2005.sm.dt, mets2005.sm.dt)
in.lg.cv <- set.inputs( hyads2005.lg.dt, ddm2005.lg.dt, mets2005.lg.dt)
in.no <- set.inputs( hyads2005.dt, ddm2005.dt)


###################################################
### benchmark time tests for different sizes
###################################################
# time goes about with the cube of ncells
nin <- 1
# normal size square
# with covariates
results_sq.cv <- splM.hyads.ddm( coords = in.cv$coords, Y = in.cv$Y, X = in.cv$X,
                                 seed.n = 1)
# without covariates
results.sq.on <- splM.hyads.ddm( coords = in.no$coords, Y = in.no$Y, X = in.no$X,
                                 seed.n = 1)

# with covariates & knots
results_sq.cv.5k <- splM.hyads.ddm( coords = in.cv$coords, Y = in.cv$Y, X = in.cv$X,
                                    knots = c( 5, 5, 0), seed.n = 1)
results_sq.cv.10k <- splM.hyads.ddm( coords = in.cv$coords, Y = in.cv$Y, X = in.cv$X,
                                     knots = c( 10, 10, 0), seed.n = 1)

# with small, large & knots
results_sm.cv <- splM.hyads.ddm( coords = in.sm.cv$coords, Y = in.sm.cv$Y, X = in.sm.cv$X,
                                 seed.n = 1)
results_sm.cv.5k <- splM.hyads.ddm( coords = in.sm.cv$coords, Y = in.sm.cv$Y, X = in.sm.cv$X,
                                    knots = c( 5, 5, 0), seed.n = 1)
results_lg.cv.5k <- splM.hyads.ddm( coords = in.lg.cv$coords, Y = in.lg.cv$Y, X = in.lg.cv$X,
                                    knots = c( 5, 5, 0), seed.n = 1)
results_lg.cv.10k <- splM.hyads.ddm( coords = in.lg.cv$coords, Y = in.lg.cv$Y, X = in.lg.cv$X,
                                     knots = c( 10, 10, 0), seed.n = 1)
results_lg.cv.20k <- splM.hyads.ddm( coords = in.lg.cv$coords, Y = in.lg.cv$Y, X = in.lg.cv$X,
                                     knots = c( 20, 20, 0), seed.n = 1)

###################################################
### put together benchmarks
###################################################
results.list <- list( results_sq.cv, results.sq.on, results_sq.cv.5k, results_sq.cv.10k,
                      results_sm.cv, results_sm.cv.5k, results_lg.cv.5k, results_lg.cv.10k, results_lg.cv.20k)

results.l.dt <- rbindlist( lapply( seq_along(results.list), 
                                   function( n, l, cv, kn){
                                     lm.n <- copy( l[[n]]$metrics)
                                     lt.n <- copy( l[[n]]$runtimes)
                                     lm.n[, `:=` (ncells.tr = lt.n$ncells.tr,
                                                  cv = cv[n])]
                                     lm.n[ mod.name == 'spLM', `:=` ( kn = kn[n],
                                                                      time = sum( lt.n$train,
                                                                                  lt.n$recover,
                                                                                  lt.n$est.y.hat))]
                                     # lm.n[ mod.name == 'lm', `:=` ( cv = cv[n])]
                                     return( lm.n)                          
                                   }, results.list, 
                                   cv = c( 'yes', 'no', rep( 'yes', 7)),
                                   kn = c( 0, 0, 5, 10, 0, 5, 5, 10, 20)))

ggplot( results.l.dt) + 
  geom_point( aes( color = kn, y = time, x = ncells.tr))
ggplot( results.l.dt) + 
  geom_point( aes( color = ncells.tr, y = RMSE, x = kn))
ggplot( results.l.dt) + 
  geom_point( aes( color = mod.name, y = NME, x = ncells.tr), alpha = 0.5)
ggplot( results.l.dt) + 
  geom_point( aes( color = mod.name, y = R, x = ncells.tr), alpha = 0.5)
ggplot( results.l.dt) + 
  geom_point( aes( color = mod.name, y = MB, x = ncells.tr), alpha = 0.5)
ggplot( results.l.dt) + 
  geom_point( aes( color = cv, y = MB, x = mod.name), alpha = 0.5)

###################################################
### old tests
###################################################
set.seed( 1)
results.sq1a <- splM.hyads.ddm( dummy.n,
                                knots = c( 5, 5),
                                coords = as.matrix( hyads2005.dt[,.( x, y)]),
                                Y = ddm2005.dt$X2005,
                                X = cbind( 1, hyads2005.dt$X2005, mets2005.dt$temp, mets2005.dt$apcp))
set.seed( 1)
results.sq2 <- splM.hyads.ddm( dummy.n,
                               knots = c( 20, 20),
                               coords = as.matrix( hyads2005.dt[,.( x, y)]),
                               Y = ddm2005.dt$X2005,
                               X = cbind( 1, hyads2005.dt$X2005, mets2005.dt$temp, mets2005.dt$apcp))
# benchmark different sizes
results.sm <- splM.hyads.ddm( dummy.n,
                              coords = as.matrix( hyads2005.sm.dt[,.( x, y)]),
                              Y = ddm2005.sm.dt$X2005,
                              X = cbind( 1, hyads2005.sm.dt$X2005))
# knots 10 => 15 went from time 3 => 15
results.lg <- splM.hyads.ddm( dummy.n,
                              knots = c( 20, 20, 0),
                              coords = as.matrix( hyads2005.lg.dt[,.( x, y)]),
                              Y = ddm2005.lg.dt$X2005,
                              X = cbind( 1, hyads2005.lg.dt$X2005, mets2005.lg.dt$temp, mets2005.lg.dt$apcp))

# set the seed
set.seed( 12345)
holdouts.eval <- rbindlist( lapply( 1, splM.hyads.ddm))
# thoughts
# 1. log first?
# 2. what to do with the uncertainties in y.hat?
# 3. what about uncertainties in HyADS and CMAQ-DDM?
# 
# To-dp
# 1. benchmark times for different sizes
# 2. can it be run with multiple X variables?

###################################################
### code chunk number 3: spRecover1
###################################################
m.i.rec <- spRecover(m.i, start=burn.in, thin=5, n.report=100)
plot( m.i.rec$p.beta.recover.samples)
beta.hat <- round(summary(m.i.rec$p.beta.recover.samples)$quantiles[c(3,1,5)],6)

###################################################
### code chunk number 4: wHatSurf1
###################################################
# The w(si)’s provide local adjustment (with structured dependence) to the mean
# and capturing the effect of unmeasured or unobserved regressors with spatial
# pattern.
w.hat <- apply(m.i.rec$p.w.recover.samples, 1, median)

w.hat.raster <- rasterFromXYZ( cbind( coords, w.hat))
plot( w.hat.raster)

theta.raster <- rasterFromXYZ( cbind( coords, m.i.rec$p.theta.recover.samples))

m.i.pred <- spPredict(m.i, start=burn.in, thin=2, pred.covars=X.ho, 
                      pred.coords=coords, verbose=FALSE)
y.hat <- apply(m.i.pred$p.y.predictive.samples, 1, quants)
par(mar=c(5,5,5,5))
plot(Y, y.hat[1,], pch=19, cex=0.5, xlab="Observed y", ylab="Predicted y", 
     ylim=range(y.hat), xlim=range(y.hat), cex.lab=2, cex.axis=2)
arrows(Y, y.hat[1,], Y, y.hat[2,], angle=90, length=0.05)
arrows(Y, y.hat[1,], Y, y.hat[3,], angle=90, length=0.05)
lines(-20:20,-20:20, col="blue")

###################################################
### code chunk number 5: runTime1
###################################################
run.time <- round(m.i$run.time[3]/60,3)

###################################################
### code chunk number 7: ppModI
###################################################
m.i <- spLM( Y ~ X - 1, coords=coords, knots=c(5, 5, 0), starting=starting,
             tuning=tuning, priors=priors, cov.model="exponential",
             modified.pp=FALSE, n.samples=n.samples, n.report=2500)


###################################################
### code chunk number 8: ppModII-IV
###################################################
m.ii <- spLM(Y~X-1, coords=coords, knots=c(5, 5, 0), starting=starting,
             tuning=tuning, priors=priors, cov.model="exponential",
             modified.pp=TRUE, n.samples=n.samples, verbose=FALSE)

m.iii <- spLM(Y~X-1, coords=coords, knots=c(10, 10, 0), starting=starting,
              tuning=tuning, priors=priors, cov.model="exponential",
              modified.pp=FALSE, n.samples=n.samples, verbose=FALSE)

m.iv <- spLM(Y~X-1, coords=coords, knots=c(10, 10, 0), starting=starting,
             tuning=tuning, priors=priors, cov.model="exponential",
             modified.pp=TRUE, n.samples=n.samples, verbose=FALSE)

quants <- function(x){
  quantile(x, prob=c(0.5, 0.025, 0.975))
}

ci.print <- function(x, digits=2){
  apply(round(apply(x, 2, quants), digits), 2, function(x){paste(x[1]," (",x[2],", ",x[3],")",  sep="")})
}

m.i <- spRecover(m.i, start=burn.in, thin=5, verbose=FALSE)
m.ii <- spRecover(m.ii, start=burn.in, thin=5, verbose=FALSE)
m.iii <- spRecover(m.iii, start=burn.in, thin=5, verbose=FALSE)
m.iv <- spRecover(m.iv, start=burn.in, thin=5, verbose=FALSE)

sub <- burn.in:n.samples

B <- as.matrix(c(1,5))
sigma.sq <- 2
tau.sq <- 1
phi <- 6

ests <- cbind(c(B[1], B[2], sigma.sq, tau.sq, phi), 
              c(ci.print(m.i$p.beta.samples[sub,]), ci.print(m.i$p.theta.samples[sub,])),
              c(ci.print(m.ii$p.beta.samples[sub,]), ci.print(m.ii$p.theta.samples[sub,])),
              c(ci.print(m.iii$p.beta.samples[sub,]), ci.print(m.iii$p.theta.samples[sub,])),
              c(ci.print(m.iv$p.beta.samples[sub,]), ci.print(m.iv$p.theta.samples[sub,])))

run.times <- c("", round(c(m.i$run.time[3]/60, m.ii$run.time[3]/60, m.iii$run.time[3]/60,  m.iv$run.time[3]/60), 2))

rel.run.times <- c("", round(c(m.i$run.time[3]/m.iv$run.time[3], m.ii$run.time[3]/m.iv$run.time[3], m.iii$run.time[3]/m.iv$run.time[3],  m.iv$run.time[3]/m.iv$run.time[3]), 2))

tab <- rbind(ests, run.times, rel.run.times)

rownames(tab) <- c("$\\beta_0$","$\\beta_1$","$\\sigma^2$", "$\\tau^2$", "$\\phi$", "Time", "Rel. time")

colnames(tab) <- c("true", "i", "ii", "iii", "iv")


###################################################
### code chunk number 9: tab1
###################################################
print(xtable(tab, 
             caption="Candidate predictive process models' parameter estimates, run-time (wall time) in minutes, and run-time relative to model $iv$. Parameter posterior summary 50 (2.5, 97.5) percentiles.", label="tab1", align="cccccc"),
      table.placement="!ht",  caption.placement="bottom", sanitize.text.function=function(x){x})



crs.usa <- crs( "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000")

# download USA polygon from rnaturalearth
usa <- rnaturalearth::ne_countries(scale = 110, type = "countries", country = "United States of America", 
                                   geounit = NULL, sovereignty = NULL,
                                   returnclass = c("sp"))
usa.sub <- disaggregate(usa)[6,]
usa.sub.t <- spTransform(usa.sub, CRSobj = crs.usa)
ext <- extent(usa.sub)

r <- raster( resolution = 1, ext = extent(usa), crs = crs( usa))
crop( r, usa.sub)
crop( r, usa.sub.t)
ext.proj <- projectExtent( usa.sub, p4s)
ext.proj.t <- projectExtent( usa.sub.t, crs( usa))
crop( r, ext.proj.t)

projectExtent( raster( resolution = 1, ext, crs = crs( usa)), crs = crs( p4s))

r <- raster(xmn=-110, xmx=-90, ymn=40, ymx=60, ncols=40, nrows=40)

