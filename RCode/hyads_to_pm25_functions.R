library( data.table)
library( raster)
library( spBayes)
library( hyspdisp)
library( ggplot2)
library( viridis)
library( lubridate)

`%ni%` <- Negate(`%in%`) 

#======================================================================#
## ddm to grid raster
#======================================================================#
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
  names.ddm <- unique( ddm_coal.m$date.in)
  ddm_coal.b <- brick( lapply( names.ddm,
                               function( name, dt.m){
                                 r <- rasterFromXYZ(dt.m[ date.in == name,
                                                          .(x = X, y = Y, z = coal_pm25)],
                                                    crs = CRS(p4s))
                                 names( r) <- name
                                 return( r)
                               }, ddm_coal.m))
  
  # fill NA's with linear interpolation across days
  ddm_coal.b <- approxNA( ddm_coal.b, rule=2)
  names( ddm_coal.b) <- names.ddm
  
  # take monthly averages
  if( avg.period == 'month'){
    ddm_coal.month <- lapply( 1:12,
                              function( m, ddm_raster.b){
                                names.dates <- as.Date( gsub( '\\.', '-', gsub( 'X', '', names( ddm_raster.b))))
                                id <- which( month( names.dates) == m)
                                
                                ddm_coal.mon <- mean( subset( ddm_raster.b, id))
                                
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


#======================================================================#
## functions to get meteorology data
# download the necessary met files, 20th century reanalysis
#======================================================================#
downloader.fn <- function( filename,
                           destination = file.path('~', 'Dropbox', 'Harvard', 'RFMeval_Local', 
                                                   'Comparisons_Intermodel', 'Global_meteorology'),
                           dataset = c( '20thC_ReanV2c', 'ncep.reanalysis.derived', 'NARR')){
  if( length( dataset) > 1)
    dataset <- dataset[1]
  fileloc <- file.path( destination, dataset)
  
  
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

#======================================================================#
## functions to get meteorology data
# extract the year of interest, average by year or return months
#======================================================================#
extract_year.fn <- function( raster.in = list.met[[1]],
                             year.in = 2005,
                             avg.period = 'year',
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
  raster.sub <- brick( subset_nc_date(  hpbl_brick = raster.in,
                                        vardate = names.months))
  
  #NARR dataset requires rotating
  if( dataset != 'NARR')
    raster.sub <- rotate( raster.sub)
  
  # take annual mean
  if( avg.period == 'year'){
    out <- stackApply( raster.sub, indices = rep( 1, 12), fun = mean)
  } else
    out <- raster.sub
  
  return( out)
}

#======================================================================#
## functions to get meteorology data
# trim data over US, create raster object
#======================================================================#
usa.functioner <- function( year.in = 2005,
                            list.met,
                            dataset = c( '20thC_ReanV2c', 'ncep.reanalysis.derived', 'NARR'),
                            avg.period = 'year',
                            return.usa.mask = F,
                            return.usa.sub = T){
  
  # extract year
  if( avg.period == 'year'){
    mets <- brick( lapply( list.met,
                           extract_year.fn,
                           year.in = year.in,
                           avg.period = avg.period,
                           dataset = dataset))
  } else
    mets <- lapply( list.met,
                    extract_year.fn,
                    year.in = year.in,
                    avg.period = avg.period,
                    dataset = dataset)
  
  crs.str <- projection( list.met[[1]])
  crs.usa <- crs( crs.str)
  
  # convert temp to celcius
  mets$temp <- mets$temp - 273.15
  
  # calculate windspeed
  # calculate meteorology wind angle (0 is north wind)
  # http://weatherclasses.com/uploads/3/6/2/3/36231461/computing_wind_direction_and_speed_from_u_and_v.pdf
  if( 'uwnd' %in% names( list.met) & 'vwnd' %in% names( list.met)){
    mets$wspd <- sqrt( mets$uwnd ^ 2 + mets$vwnd ^ 2)
    mets$phi <- atan2( mets$uwnd, mets$vwnd) * 180 / pi + 180
    names( mets$wspd) <- names(mets$vwnd)
    names( mets$phi) <- names(mets$vwnd)
  }
  
  # download USA polygon from rnaturalearth
  usa <- rnaturalearth::ne_countries(scale = 110, type = "countries", country = "United States of America", 
                                     geounit = NULL, sovereignty = NULL,
                                     returnclass = c("sp"))
  usa.sub <- disaggregate(usa)[6,]
  usa.sub <- spTransform(usa.sub, CRSobj = crs.str) #proj4string( ))
  
  if( return.usa.mask){
    usa.sub.p <- spTransform(usa.sub, CRSobj = crs.usa)
    # usa.sub.sf <- sf::st_as_sf( usa.sub.p)
    return( usa.sub.p)
  }
  
  if( return.usa.sub){
    # crop to USA
    if( avg.period == 'year'){
      mets.out <- crop( mask(mets, usa.sub), usa.sub)
    } else{ 
      mets.mask <- lapply( mets, mask, usa.sub)
      mets.usa <- lapply( mets.mask, crop, usa.sub)
      
      # reorganize - list of months
      names.months <- names( mets.usa[[1]])
      mets.out <- lapply( names.months, 
                          function( name.ext, X){
                            brick( lapply( X, '[[', name.ext))
                          }, mets.usa)
      names( mets.out) <- names.months
    }
    
    
    return( mets.out)
  } else{
    if( avg.period == 'year'){
      mets.out <- mets
    } else{ 
      
      # reorganize - list of months
      names.months <- names( mets[[1]])
      mets.out <- lapply( names.months, 
                          function( name.ext, X){
                            brick( lapply( X, '[[', name.ext))
                          }, mets)
      names( mets.out) <- names.months
    }
    return( mets.out)
  }
}

#======================================================================#
# define the spBayes model function
#======================================================================#
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

#======================================================================#
# define the linear model holdout function
#======================================================================#
lm.hyads.ddm.holdout <- function( seed.n = NULL,
                                  dat.stack = dats2005.s.small,
                                  dat.stack.pred = NULL,
                                  y.name = 'cmaq.ddm',
                                  x.name = 'hyads',
                                  name.idwe = 'tot.sum',
                                  covars.names = NULL, #c( 'temp', 'apcp'),
                                  ho.frac = .1,
                                  return.mods = F,
                                  ...){
  set.seed( seed.n)
  
  # if no covar names provided, use all covariates that aren't x or y
  if( is.null( covars.names))
    covars.names <- names( dat.stack)[names( dat.stack) %ni% c( x.name, y.name)]
  
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
  form.cv <-  as.formula( paste( y.name, '~ (', paste( c( x.name, covars.names), 
                                                       collapse = '+'), ') ^2'))
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
  
  # calculate covariate contributions
  Y.ho.terms <- predict( lm.cv,  newdata = dat.stack.ho, se.fit = T, type = 'terms')
  Y.ho.terms.dt <- data.table( dat.coords.ho, Y.ho.terms$fit)
  
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
  Y.ho.terms.raster <- projectRaster( rasterFromXYZ( Y.ho.terms.dt, crs = crs.in), dat.stack)
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
                           R = cor( Yhat, Yact, use = 'complete.obs'),
                           R.s = cor( Yhat, Yact, use = 'complete.obs', method = 'spearman'))
    return( metrics)
  }
  metrics.out <- rbind( eval.fn( Y.ho.hat$y.hat.lm.cv, y.ho, 'lm.cv'),
                        eval.fn( Y.ho.hat$y.hat.lm.ncv, y.ho, 'lm.ncv'),
                        eval.fn( Y.ho.hat$y.hat.mean, y.ho, 'adj.mean'),
                        eval.fn( Y.ho.hat$y.hat.Z, y.ho, 'adj.Z'))
  
  # calculate correlations along quantiles of idwe
  evals.q <- data.table()
  if( ho.frac == 0){
    vals.idwe <- unlist( data.table( values( dat.stack))[, ..name.idwe])
    vals.eval <- data.table( values( Y.ho.hat.raster))
    for ( s in seq( 0.01, 1, .01)){
      q <- quantile( vals.idwe, s, na.rm = T)
      vals.use <- vals.eval[vals.idwe < q,]
      evals <- eval.fn( vals.use$y.hat.lm.cv, vals.use$y.ho, x.name)[, s := s]
      evals.q <- rbind( evals.q, evals)
    }
  }
  
  # listify the models
  if( return.mods)
    out <- list( metrics = metrics.out, 
                 evals.q = evals.q,
                 model.cv = lm.cv,
                 model.ncv = lm.ncv,
                 Y.ho.hat.raster = Y.ho.hat.raster, 
                 Y.ho.hat.se.raster = Y.ho.hat.se.raster,
                 Y.ho.hat.bias.raster = Y.ho.hat.bias.raster,
                 Y.ho.terms.raster = Y.ho.terms.raster)
  if( !return.mods)
    out <- list( metrics = metrics.out, 
                 evals.q = evals.q,
                 Y.ho.hat.raster = Y.ho.hat.raster, 
                 Y.ho.hat.se.raster = Y.ho.hat.se.raster,
                 Y.ho.hat.bias.raster = Y.ho.hat.bias.raster,
                 Y.ho.terms.raster = Y.ho.terms.raster)
  
  return( out)
}

#======================================================================#
# extract a mean from a list of lists of rasters
#======================================================================#
mean.lol <- function( lol, layer1, layer2, plot.out = TRUE){
  first.lol <- lapply( lol, '[[', layer1)
  second.lol <- lapply( first.lol, '[[', layer2)
  mean.lol <- mean( stack( second.lol), na.rm = T)
  
  if( plot.out)
    plot( mean.lol, main = layer2)
  return( mean.lol)
}

#======================================================================#
# ggplot a raster
#======================================================================#
ggplot.a.raster <- function( ..., bounds = NULL, facet.names = NULL, 
                             mask.raster = NULL, nrow. = NULL, ncol. = NULL){
  
  in.x <- list( ...)
  
  if( length( in.x) == 1)
    in.x <- in.x[[1]]
  
  if( is.null( names( in.x)) & !is.null( facet.names))
    names( in.x) <- facet.names
  
  if( !is.null( mask.raster) & is.list( in.x))
    in.x.crop <- lapply( in.x, function( X, mask.raster.){
      X.mask <- mask( X, mask.raster.)
      X.crop <- crop( X.mask, mask.raster.)
      return( X.crop)
    }, mask.raster)
  if( !is.null( mask.raster) & !is.list( in.x)){
    in.x.mask <- mask( in.x, mask.raster)
    in.x.crop <- crop( in.x.mask, mask.raster)
  }
  
  dat.dt <- rbindlist( lapply( names( in.x.crop), function( x.name, x.list) {
    x <- x.list[[x.name]]
    r_points <- rasterToPoints( x)
    r_dt <- data.table( r_points)[, name.in := x.name]
    setnames( r_dt, names( x), 'z')
    return( r_dt)
  }, in.x.crop))
  
  if( !is.null( facet.names))
    dat.dt[, name.in := factor( name.in, levels = facet.names)]
  
  ggplot( dat.dt) +
    geom_tile( aes( x = x, y = y, fill = z)) + 
    scale_fill_viridis( name = NULL, limits = bounds, oob = scales::squish) + 
    facet_wrap( . ~ name.in, nrow = nrow., ncol = ncol.) + 
    # expand_limits( fill = 0) +
    theme_bw() + 
    theme( axis.text = element_blank(),
           axis.title = element_blank(),
           axis.ticks = element_blank(),
           legend.position = 'bottom',
           panel.grid = element_blank(),
           strip.background = element_blank())
}

#======================================================================#
# do the predictions for many months
#======================================================================#
month.trainer <- function( name.m = names( mets2005.m)[1], 
                           name.p = names( mets2006.m)[1], 
                           name.x,
                           y.m, 
                           ddm.m = ddm.m.all, 
                           mets.m = mets.m.all,
                           emiss.m = d_nonegu.r,
                           idwe.m = idwe.m,
                           .mask.use = NULL,
                           cov.names = c( "temp", "rhum", "vwnd", "uwnd", "wspd")){ #, names( d_nonegu.r))){
  
  # create training dataset
  ddm.use <- ddm.m[[name.m]]
  hyads.use <- y.m[[name.m]]
  mets.use <- mets.m[[name.m]]
  emiss.use <- emiss.m
  idwe.use <- idwe.m[[name.m]]
  
  # create prediction dataset
  ddm.use.p <- ddm.m[[name.p]]
  hyads.use.p <- y.m[[name.p]]
  mets.use.p <- mets.m[[name.p]]
  idwe.use.p <- idwe.m[[name.p]]
  
  # fix names
  names( ddm.use)     <- 'cmaq.ddm'
  names( ddm.use.p)   <- 'cmaq.ddm'
  names( hyads.use)   <- name.x
  names( hyads.use.p) <- name.x
  names( idwe.use)   <- 'tot.sum.idwe'
  names( idwe.use.p) <- 'tot.sum.idwe'
  
  # combine each dataset as stacks
  dat.s <- project_and_stack( ddm.use,   hyads.use,   idwe.use,   
                              mets.use, emiss.use, mask.use = .mask.use)
  dat.p <- project_and_stack( ddm.use.p, hyads.use.p, idwe.use.p, 
                              mets.use.p, emiss.use, mask.use = .mask.use)
  
  # do the modeling
  pred <- lm.hyads.ddm.holdout( dat.stack = dat.s, dat.stack.pred = dat.p, x.name = name.x,
                                ho.frac = 0, covars.names = cov.names, 
                                name.idwe = 'tot.sum.idwe', return.mods = T)
  
  return( pred)
}


#======================================================================#
# project and stack in one command
#======================================================================#
project_and_stack <- function( ..., mask.use = NULL){
  list.r <- list( ...)
  for( r in 2: length( list.r)){
    list.r[[r]] <- projectRaster( list.r[[r]], list.r[[1]])
  }
  
  # mask over usa
  if( !is.null( mask.use)){
    list.out <- lapply( list.r, mask, mask.use)
  } else
    list.out <- list.r
  
  return( stack( list.out))
}
